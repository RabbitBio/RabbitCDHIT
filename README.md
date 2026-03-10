# RabbitCD-HIT

RabbitCD-HIT is a scalable greedy incremental clustering tool designed for massive protein sequences. It supports running on distributed clusters to efficiently process massive protein sequence datasets. RabbitCD-HIT combines inter-node MPI parallelism with intra-node multithreading to achieve superior speed and scalability while maintaining full compatibility with the CD-HIT output format.



## 1) Build

```bash
git clone --recursive https://github.com/RabbitBio/RabbitCDHIT.git
cd RabbitCDHIT
make
```

## 2) Workflow (must follow this order)

1. Run `cdhit-preprocess` to generate preprocessed files and `info.json`.
2. Run `cdhit-mpi` for clustering using the generated files.

## 3) Preprocess stage

```bash
./cdhit-preprocess -i DB.fa -N 1 -NT 128 -T 64 -ST 64 -nT 32 -tmp /absolute/path/to/tmp_runs -pre_out /absolute/path/to/preprocess_output
```

### Important preprocess arguments

- `-i`: input FASTA file (supports `.gz`)
- `-tmp`: temp run-file directory (recommend local fast disk)
- `-pre_out`: preprocess output directory (`info.json` and `_proc*.fa` will be written here)
- `-ST`: physical cores per socket (optional, recommended for multi-machine runs)
- `-nT`: physical cores per NUMA node (optional, recommended for multi-machine runs)

## 4) Clustering stage

Use `-i` in `cdhit-mpi` to point to the same directory as preprocess `-pre_out`.

```bash
mpirun -np <NP_FROM_JSON> ./cdhit-mpi -i /absolute/path/to/preprocess_output -o DB_output -T <T_FROM_JSON>
```

Work-stealing example:

```bash
mpirun -np <NP_FROM_JSON> ./cdhit-mpi -i /absolute/path/to/preprocess_output -o DB_output -T <T_FROM_JSON> -stealing 1
```

## 5) Critical requirement for `cdhit-mpi`

When running `./cdhit-mpi`, these two values MUST match `info.json`:

- `mpirun -np` **must equal** `info.total_mpi_num`
- `cdhit-mpi -T` **must equal** `info.threads_per_node`

If they do not match, the run may fail or show unstable performance/results.

## 6) Read `-np` and `-T` from `info.json`

Read `info.total_mpi_num` and `info.threads_per_node` from:

- `/absolute/path/to/preprocess_output/info.json`

Then run:

```bash
mpirun -np <value_from_info.total_mpi_num> ./cdhit-mpi -i /absolute/path/to/preprocess_output -o DB_output -T <value_from_info.threads_per_node>
```

## 7) Output files

- `*.txt`: representative sequence names and sequences
- `*.clstr`: cluster membership information
