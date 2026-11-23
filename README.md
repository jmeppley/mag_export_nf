# mag_export_nf

A Nextflow workflow that, for each of one or more tables:

1. pulls MAG fasta file names from one column (de-duplicated, first-seen order),
2. copies those fastas from a MAG directory into the table's destination folder, and
3. optionally predicts genes (pyrodigal) and annotates them (eggNOG-mapper).

A MAG listed in several tables is annotated only once; its results are copied into
every table that asked for annotation.

## Setup (once per server)

```bash
git clone <this repo> ~/work/mag_export_nf    # or copy the folder
cd ~/work/mag_export_nf
mamba env create -f environment.yml           # nextflow + java
mamba activate nextflow                       # sets JAVA_HOME for nextflow
```

The tool environments (python, pyrodigal, eggnog-mapper) in `envs/` are built
automatically by Nextflow the first time you run with `-profile mamba` (or
`-profile conda`). They are cached in `~/.cache/mag_export_nf/conda`.

### eggNOG database (only needed if any table has `annotate: true`)

The DB is about 50 GB. Either point `eggnog_db` at an existing eggNOG-mapper data dir
(it must contain `eggnog.db` and `eggnog_proteins.dmnd`), or let the workflow
download it on the first run:

```bash
nextflow run main.nf -profile mamba -params-file params.yaml \
    --eggnog_db /path/to/eggnog_db --download_eggnog_db true
```

The download step is skipped on later runs once the files exist in `eggnog_db`.

## Configure

Copy `params.example.yaml` to `params.yaml` and edit it:

```yaml
mag_dir: /mnt/delong/scratch2/for_John
eggnog_db: /path/to/eggnog_db
on_missing: warn          # warn | fail
tables:
  - table: 795_decreasing_prev.tsv
    header: false
    column: 1             # 1-based index when there is no header
    dest: results/decreasing_795
    annotate: true
  - table: combined_maaslin3_results.marine_MAGs.taxon_annotated.2.tsv
    header: true
    column: feature       # column name when there is a header
    dest: results/maaslin3_all
    annotate: false
```

Per-table keys: `table`, `column`, `dest` (required), plus `header` (default false),
`annotate` (default false) and `sep` (default tab). Relative paths are resolved from
the folder you launch in.

Global params (defaults in `nextflow.config`):

| param | default | meaning |
|---|---|---|
| `mag_dir` | — | folder of MAG fasta files (required) |
| `mag_ext` | `''` | suffix to add to names that lack it, e.g. `.fa` |
| `on_missing` | `warn` | `warn`: log and continue; `fail`: stop the run |
| `eggnog_db` | — | eggNOG-mapper data dir |
| `download_eggnog_db` | `false` | download the DB into `eggnog_db` if absent |
| `pyrodigal_mode` | `single` | `single` (falls back to `meta` for small MAGs) or `meta` |
| `eggnog_batch_size` | `50` | MAGs per eggNOG-mapper run |
| `emapper_cpus` | `16` | CPUs per eggNOG-mapper run |
| `emapper_args` | `''` | extra `emapper.py` options |

## Run

```bash
mamba activate nextflow
nextflow run ~/work/mag_export_nf/main.nf -profile mamba -params-file params.yaml -resume
```

Use `screen`/`tmux` for long runs. Nextflow runs tasks in parallel on the local
machine and may use all its CPUs. To share the server, cap it with a small config:

```bash
echo 'executor.cpus = 32' > local.config
nextflow run ~/work/mag_export_nf/main.nf -profile mamba -params-file params.yaml -c local.config -resume
```

## Outputs (per table `dest`)

```
dest/
├── mags.txt                 # MAG names found in mag_dir
├── missing_mags.txt         # names in the table but not in mag_dir
├── mags/                    # copied fasta files
└── annotation/              # only if annotate: true
    └── <mag>/
        ├── <mag>.faa        # proteins (pyrodigal)
        ├── <mag>.ffn        # gene nucleotide sequences
        ├── <mag>.gff        # gene coordinates
        └── <mag>.emapper.annotations
```

Intermediate files are in `work/`; delete it after a successful run.

## Test

```bash
# copy step only, synthetic data
nextflow run main.nf -profile test,mamba
# full wiring including annotation, without running the tools
nextflow run main.nf -params-file tests/params.test.yaml -stub \
    --eggnog_db tests/stub_db --download_eggnog_db true --emapper_cpus 2
```
