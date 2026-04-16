# RabbitCD-HIT

RabbitCD-HIT is a scalable greedy incremental clustering tool designed for massive protein sequences. It supports running on distributed clusters to efficiently process massive protein sequence datasets. RabbitCD-HIT combines inter-node MPI parallelism with intra-node multithreading to achieve superior speed and scalability while maintaining full compatibility with the CD-HIT output format.

## 1) Build

```bash
git clone --recursive https://github.com/RabbitBio/RabbitCDHIT.git
cd RabbitCDHIT
make
```

By default, the project uses `mpiicpx` as the MPI C++ compiler.  
If `mpiicpx` is not available on your machine, edit `CC` in `Makefile` and switch it to `mpicxx`, then rebuild.  
When using `mpicxx`, also change `-Ofast` to `-O2` in `Makefile` to avoid potential correctness issues.

If TBB (Threading Building Blocks) is not found during build, set `TBB_ROOT` in `Makefile` to point to your TBB installation:

```makefile
TBB_ROOT = /path/to/your/tbb
```

If your CPU supports AVX512, compile with:

```bash
make AVX512=yes
```

---

## 2) Single-node quick start

For single-node runs, use the provided `run_single.sh` script to complete both stages in one command:

```bash
bash run_single.sh -i DB.fa -o result
```

The script automatically detects the machine's CPU topology (`nproc`, `lscpu`) and sets `-NT` (total cores) and `-T` (cores per NUMA node) accordingly; both can also be specified manually. It also reads `info.json` produced by the preprocess stage and passes the correct `-np` and `-T` to `cdhit-mpi`.

### Options


| Option | Default      | Description                       |
| ------ | ------------ | --------------------------------- |
| `-i`   | *(required)* | Input FASTA file (supports `.gz`) |
| `-o`   | *(required)* | Output file prefix                |
| `-c`   | `0.9`        | Sequence identity threshold       |


Any additional `cdhit-mpi` options (`-g`, `-fo`, `-stealing`) can be appended directly:

```bash
bash run_single.sh -i DB.fa -o result -c 0.9 -g 1 -fo 1 -stealing 1
```

---

## 3) Manual workflow

> **When running on a cluster (multi-node), this manual workflow must be used.**  
> The two stages below **must be run in order**.

### 3.1) Preprocess stage

```bash
./cdhit-preprocess -i DB.fa -N 1 -NT 128 -T 64 -tmp /absolute/path/to/tmp_runs -pre_out /absolute/path/to/preprocess_output
```


| Argument   | Description                                                           |
| ---------- | --------------------------------------------------------------------- |
| `-i`       | Input FASTA file (supports `.gz`)                                     |
| `-N`       | Number of nodes                                                       |
| `-NT`      | Total threads per node for the MPI workflow                           |
| `-T`       | Threads used by `cdhit-preprocess` itself (one NUMA node recommended) |
| `-tmp`     | Temp run-file directory                                               |
| `-pre_out` | Preprocess output directory                                           |


### 3.2) Clustering stage

Read `total_mpi_num` and `threads_per_rank` from the generated `info.json`, then run:

```bash
mpirun -np <total_mpi_num> ./cdhit-mpi -i /absolute/path/to/preprocess_output -o DB_output -T <threads_per_rank>
```

> `**mpirun -np` and `cdhit-mpi -T` MUST match `info.json`, otherwise the run may fail.**

Work-stealing example:

```bash
mpirun -np <total_mpi_num> ./cdhit-mpi -i /absolute/path/to/preprocess_output -o DB_output -T <threads_per_rank> -stealing 1
```

---

## 4) Output files

- `*.txt`: representative sequence names and sequences
- `*.clstr`: cluster membership information

