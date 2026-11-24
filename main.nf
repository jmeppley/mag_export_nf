#!/usr/bin/env nextflow
/*
 * mag_export_nf
 * For each configured table: pull MAG fasta names from a column, copy those
 * fastas from params.mag_dir into the table's dest folder, and (optionally)
 * predict genes with pyrodigal and annotate them with eggNOG-mapper.
 *
 * Each MAG is annotated once even if several tables list it; results are
 * copied into every annotating table's dest.
 */

nextflow.enable.dsl = 2

/* ---------- processes ---------- */

process EXTRACT_MAG_NAMES {
    tag { meta.id }
    label 'small'
    conda "${moduleDir}/envs/python.yml"
    publishDir path: { meta.dest }, mode: 'copy'

    input:
    tuple val(meta), path(table)

    output:
    tuple val(meta), path('mags.txt'), path('missing_mags.txt')

    script:
    def header = meta.header ? '--header' : ''
    """
    extract_mag_names.py ${table} \\
        --column '${meta.column}' ${header} \\
        --sep '${meta.sep}' \\
        --mag-dir '${meta.mag_dir}' \\
        --mag-ext '${params.mag_ext}' \\
        --on-missing ${params.on_missing}
    """
}

process COPY_MAGS {
    tag { meta.id }
    label 'small'
    // mode 'copy' resolves the staged symlinks, so dest gets real files
    publishDir path: { "${meta.dest}/mags" }, mode: 'copy', saveAs: { it.tokenize('/')[-1] }

    input:
    tuple val(meta), path(fastas, stageAs: 'mags/*')

    output:
    tuple val(meta), path('mags/*', includeInputs: true)

    script:
    """
    echo "Copying \$(ls mags | wc -l) MAGs for ${meta.id}"
    """
}

process DOWNLOAD_EGGNOG_DB {
    label 'small'
    conda "${moduleDir}/envs/eggnog.yml"
    // skipped automatically once the files exist in params.eggnog_db
    storeDir { params.eggnog_db }

    output:
    path 'eggnog.db'
    path 'eggnog.taxa.db'
    path 'eggnog.taxa.db.traverse.pkl'
    path 'eggnog_proteins.dmnd'

    script:
    """
    download_eggnog_data.py -y --data_dir .
    """

    stub:
    """
    touch eggnog.db eggnog.taxa.db eggnog.taxa.db.traverse.pkl eggnog_proteins.dmnd
    """
}

process PYRODIGAL {
    tag { mag_id }
    label 'small'
    conda "${moduleDir}/envs/pyrodigal.yml"
    // CPUs come in via an env var so changing them keeps -resume cache hits
    beforeScript { "export TASK_CPUS=${task.cpus}" }

    input:
    tuple val(mag), path(fasta)

    output:
    tuple val(mag), path("${mag_id}.faa"), path("${mag_id}.ffn"), path("${mag_id}.gff")

    script:
    mag_id = fasta.baseName
    """
    run_pyrodigal() {
        pyrodigal -i ${fasta} -p \$1 -j \$TASK_CPUS \\
            -a ${mag_id}.faa -d ${mag_id}.ffn -f gff -o ${mag_id}.gff
    }
    # single mode needs >= 20 kb of sequence to train; fall back to meta mode
    if [ "${params.pyrodigal_mode}" = "single" ]; then
        run_pyrodigal single || { echo "single mode failed for ${mag_id}; using meta" >&2; run_pyrodigal meta; }
    else
        run_pyrodigal ${params.pyrodigal_mode}
    fi
    """

    stub:
    mag_id = fasta.baseName
    """
    printf ">${mag_id}_1\\nMKV\\n" > ${mag_id}.faa
    printf ">${mag_id}_1\\nATGAAAGTT\\n" > ${mag_id}.ffn
    touch ${mag_id}.gff
    """
}

process EGGNOG_MAPPER {
    tag { "batch ${batch_id} (${faas instanceof List ? faas.size() : 1} MAGs)" }
    label 'big'
    conda "${moduleDir}/envs/eggnog.yml"
    // CPUs come in via an env var so changing them keeps -resume cache hits
    beforeScript { "export TASK_CPUS=${task.cpus}" }

    input:
    tuple val(batch_id), path(faas, stageAs: 'faa/*')
    path db

    output:
    path 'split/*.emapper.annotations'

    script:
    """
    # merge the batch, prefixing protein IDs with the MAG id so they stay unique
    for f in faa/*.faa; do
        m=\$(basename "\$f" .faa)
        awk -v m="\$m" '/^>/{sub(/^>/, ">" m "|")} {print}' "\$f"
    done > batch.faa

    mkdir -p tmp
    emapper.py -i batch.faa --itype proteins -m diamond \\
        --data_dir ${db} --cpu \$TASK_CPUS --temp_dir tmp \\
        -o batch --output_dir . --override ${params.emapper_args}

    split_emapper.py batch.emapper.annotations --outdir split \\
        \$(for f in faa/*.faa; do basename "\$f" .faa; done)
    """

    stub:
    """
    mkdir -p split
    for f in faa/*.faa; do touch split/\$(basename "\$f" .faa).emapper.annotations; done
    """
}

process ANNOTATION_TABLE {
    tag { mag_id }
    label 'small'
    conda "${moduleDir}/envs/python.yml"

    input:
    tuple val(mag), path(gff), path(annotations)

    output:
    tuple val(mag), path("${mag_id}.annotations.tsv")

    script:
    mag_id = gff.baseName
    """
    annotation_table.py ${gff} ${annotations} -o ${mag_id}.annotations.tsv
    """
}

process COPY_ANNOTATIONS {
    tag { meta.id }
    label 'small'
    publishDir path: { "${meta.dest}/annotation" }, mode: 'copy', saveAs: { fn ->
        def name = fn.tokenize('/')[-1]
        def id = name.replaceAll(/\.(faa|ffn|gff|emapper\.annotations|annotations\.tsv)$/, '')
        "${id}/${name}"
    }

    input:
    tuple val(meta), path(files, stageAs: 'annotation/*')

    output:
    tuple val(meta), path('annotation/*', includeInputs: true)

    script:
    """
    echo "Copying annotations for ${meta.id}"
    """
}

process COMBINE_TABLES {
    tag { meta.id }
    label 'small'
    publishDir path: { meta.dest }, mode: 'copy'

    input:
    tuple val(meta), val(mags), path(tables, stageAs: 'tables/*')

    output:
    tuple val(meta), path('annotations.tsv')

    script:
    // mags[i] is the fasta name (with extension) for tables[i]; sort by name
    def tlist = tables instanceof List ? tables : [tables]
    def pairs = [mags, tlist].transpose().sort { it[0] }
    def rows = pairs.collect { m, t -> "add_rows '${m}' '${t}'" }.join('\n    ')
    """
    # one table for all MAGs in this dest, with a leading "mag" column
    set -o pipefail
    add_rows() { tail -n +2 "\$2" | awk -v m="\$1" '{print m "\\t" \$0}'; }
    {
    head -n 1 '${pairs[0][1]}' | awk '{print "mag\\t" \$0}'
    ${rows}
    } > annotations.tsv
    """
}

/* ---------- helpers ---------- */

def parseTables(tables, magDir) {
    if (!tables) error "params.tables is empty; give a -params-file (see params.example.yaml)"
    def ids = [] as Set
    tables.withIndex().collect { t, i ->
        ['table', 'column', 'dest'].each { k ->
            if (t[k] == null) error "tables[${i}] is missing required key '${k}'"
        }
        def tableFile = file(t.table, checkIfExists: true)
        def id = tableFile.baseName
        if (id in ids) id = "${id}_${i + 1}"
        ids << id
        def meta = [
            id      : id,
            column  : t.column.toString(),
            header  : (t.header ?: false) as boolean,
            sep     : t.sep ?: '\\t',
            dest    : file(t.dest).toString(),
            annotate: (t.annotate ?: false) as boolean,
            mag_dir : magDir.toString(),
        ]
        [meta, tableFile]
    }
}

/* ---------- workflow ---------- */

workflow {
    if (!params.mag_dir) error "Set params.mag_dir (folder of MAG fasta files)"
    if (!(params.on_missing in ['warn', 'fail'])) error "params.on_missing must be 'warn' or 'fail'"
    def magDir = file(params.mag_dir, checkIfExists: true)
    def tables = parseTables(params.tables, magDir)

    // ---- extract names, copy MAGs ----
    names_ch = EXTRACT_MAG_NAMES(channel.fromList(tables))

    names_ch.subscribe { meta, found, missing ->
        def n = missing.readLines().findAll { it }.size()
        if (n) log.warn "${meta.id}: ${n} MAG(s) not found in ${params.mag_dir}; listed in ${meta.dest}/missing_mags.txt"
    }

    // [meta, [mag names]]
    members_ch = names_ch.map { meta, found, missing -> [meta, found.readLines().findAll { it }] }

    COPY_MAGS(
        members_ch
            .filter { meta, mags -> mags }
            .map { meta, mags -> [meta, mags.collect { magDir.resolve(it) }] }
    )

    // ---- annotation ----
    if (tables.any { meta, t -> meta.annotate }) {
        if (!params.eggnog_db) error "Set params.eggnog_db (eggNOG-mapper data dir); add --download_eggnog_db true to fetch it (~50 GB)"

        if (params.download_eggnog_db) {
            db_ch = DOWNLOAD_EGGNOG_DB()[0].map { it.parent }
        } else {
            def db = file(params.eggnog_db)
            def absent = ['eggnog.db', 'eggnog_proteins.dmnd'].findAll { !db.resolve(it).exists() }
            if (absent)
                error "Missing ${absent.join(', ')} in ${db}; add --download_eggnog_db true to fetch the DB (~50 GB)"
            db_ch = channel.value(db)
        }

        // [mag name, meta] for every annotating table
        annot_members = members_ch
            .filter { meta, mags -> meta.annotate }
            .flatMap { meta, mags -> mags.collect { [it, meta] } }

        // Each MAG once, assigned to a fixed eggNOG batch by sorted name, so a
        // batch holds the same MAGs on every run and -resume can reuse it.
        // [mag name, batch number, batch size]
        def batchSize = params.eggnog_batch_size as int
        batches_ch = annot_members
            .map { mag, meta -> mag }
            .unique()
            .toSortedList()
            .flatMap { mags ->
                mags.collate(batchSize).withIndex().collectMany { grp, i ->
                    grp.collect { [it, i + 1, grp.size()] }
                }
            }

        genes_ch = PYRODIGAL(batches_ch.map { mag, b, n -> [mag, magDir.resolve(mag)] })

        // a batch starts as soon as all of its MAGs have genes
        emapper_ch = EGGNOG_MAPPER(
            genes_ch
                .join(batches_ch)
                .map { mag, faa, ffn, gff, b, n -> [groupKey(b, n), faa] }
                .groupTuple()
                .map { b, faas -> [b.getGroupTarget(), faas.sort { it.name }] },
            db_ch
        )
            .flatten()
            .map { f -> [f.name - ~/\.emapper\.annotations$/, f] }

        // [mag name, faa, ffn, gff, annotations]
        annots_ch = genes_ch
            .map { mag, faa, ffn, gff -> [faa.baseName, mag, faa, ffn, gff] }
            .join(emapper_ch)
            .map { id, mag, faa, ffn, gff, ann -> [mag, faa, ffn, gff, ann] }

        // per-gene summary table: [mag name, tsv]
        table_ch = ANNOTATION_TABLE(annots_ch.map { mag, faa, ffn, gff, ann -> [mag, gff, ann] })

        per_table_ch = annot_members
            .combine(annots_ch.join(table_ch), by: 0)
            .map { mag, meta, faa, ffn, gff, ann, tsv -> [meta, [mag, faa, ffn, gff, ann, tsv]] }
            .groupTuple()
            // fixed order, so the downstream tasks' cache keys don't depend on finish order
            .map { meta, recs -> [meta, recs.sort { it[0] }] }

        COPY_ANNOTATIONS(per_table_ch.map { meta, recs -> [meta, recs.collect { it[1..5] }.flatten()] })

        // dest/annotations.tsv: all of this table's MAGs in one file, "mag" = fasta name
        COMBINE_TABLES(per_table_ch.map { meta, recs -> [meta, recs.collect { it[0] }, recs.collect { it[5] }] })
    }
}
