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
2. Run `cdhit-mpi` for clustering and **you must use the output of Step 1 as input**.

## 3) Preprocess stage

```bash
./cdhit-preprocess -i DB.fa -N 1 -NT 128 -T 64 -tmp /absolute/path/to/tmp_runs -pre_out /absolute/path/to/preprocess_output
```

### Important preprocess arguments

- `-i`: input FASTA file (supports `.gz`)
- `-N`: number of nodes for the MPI workflow. 
- `-NT`: total threads per node for the MPI workflow.
- `-T`: threads used by `cdhit-preprocess` itself.
- `-tmp`: temp run-file directory.
- `-pre_out`: preprocess output directory.

## 4) Clustering stage
Read `info.total_mpi_num` and `info.threads_per_node` from either:

- `/absolute/path/to/preprocess_output/info.json`
- the terminal output printed by the Preprocess stage

Then run:

```bash
mpirun -np <value_from_info.total_mpi_num> ./cdhit-mpi -i /absolute/path/to/preprocess_output -o DB_output -T <value_from_info.threads_per_rank>
```

When running `./cdhit-mpi`, these two values MUST match `info.json`:

- `mpirun -np` **must equal** `info.total_mpi_num`
- `cdhit-mpi -T` **must equal** `info.threads_per_rank`

If they do not match, the run may fail.

Work-stealing example:

```bash
mpirun -np <value_from_info.total_mpi_num> ./cdhit-mpi -i /absolute/path/to/preprocess_output -o DB_output -T <value_from_info.threads_per_rank> -stealing 1
```

## 5) Output files

- `*.txt`: representative sequence names and sequences
- `*.clstr`: cluster membership information
