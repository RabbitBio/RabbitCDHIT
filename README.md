# RabbitCD-HIT

RabbitCD-HIT is a scalable greedy incremental clustering tool designed for massive protein sequences. It supports running on distributed clusters to efficiently process massive protein sequence datasets. RabbitCD-HIT combines inter-node MPI parallelism with intra-node multithreading to achieve superior speed and scalability while maintaining full compatibility with the CD-HIT output format.

## (1) Dependencies

RabbitCD-HIT requires an MPI implementation for distributed parallel execution.

### Required dependency

- **MPI**: required for compiling and running `cdhit-mpi`.
  - Intel MPI is recommended when using `mpiicpx`.
  - OpenMPI or MPICH can also be used through `mpicxx`.

### Optional dependency

- **TBB (Threading Building Blocks)**: optional.
  - TBB can be enabled to accelerate some intra-node operations, such as parallel sorting.
  - By default, RabbitCD-HIT is built **without TBB**.
  - To enable TBB support, both `tbb=yes` and `TBB_ROOT` must be specified during compilation.

Please make sure that an MPI implementation is installed and available in your build environment before compiling RabbitCD-HIT.

---

## (2) Build

Clone the repository and build RabbitCD-HIT with the default configuration:

```bash
git clone --recursive https://github.com/RabbitBio/RabbitCDHIT.git
cd RabbitCDHIT
make
```

By default, the project uses `mpiicpx` as the MPI C++ compiler.

If `mpiicpx` is not available on your machine, you can build with `mpicxx`:

```bash
make CC=mpicxx
```

When using `mpicxx`, it is recommended to use `-O2` instead of aggressive optimization flags such as `-Ofast`, which may cause potential correctness issues on some platforms.

---

### Build with TBB

TBB support is disabled by default. To enable TBB, the TBB installation path must be specified explicitly through `TBB_ROOT`:

```bash
make tbb=yes TBB_ROOT=/path/to/your/tbb
```

For example, when using Intel oneAPI TBB:

```bash
make tbb=yes TBB_ROOT=/opt/intel/oneapi/tbb/latest
```

Using only `make tbb=yes` without specifying `TBB_ROOT` is not supported.

---

### Build with AVX512

If your CPU supports AVX512, compile with:

```bash
make AVX512=yes
```

AVX512 support can also be combined with other build options:

```bash
make AVX512=yes tbb=yes TBB_ROOT=/opt/intel/oneapi/tbb/latest
```

or:

```bash
make CC=mpicxx AVX512=yes
```

---

### Clean build files

```bash
make clean
```

## (3) Single-node quick start

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


Any additional `cdhit-mpi` options can be appended directly:

```bash
bash run_single.sh -i DB.fa -o result -c 0.9 -g 1 -fo 1 -load_all 1
```

| Option      | Default | Description |
| ----------- | ------- | ----------- |
| `-g`        | `0`     | `1`: accurate mode — assign to the most similar cluster. `0`: fast mode — assign to the first qualifying cluster. |
| `-fo`       | `0`     | `1`: enforce strict input file order when comparing sequences of equal length. |
| `-load_all` | `0`     | `0` (default): double-buffered streaming, minimises peak memory. `1`: load all sequences into RAM. |

---

## (4) Manual workflow

> **When running on a cluster (multi-node), this manual workflow must be used.**  
> The two stages below **must be run in order**.

### 4.1 Preprocess stage

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


### 4.2 Clustering stage

Read `total_mpi_num` and `threads_per_rank` from the generated `info.json`, then run:

```bash
mpirun -np <total_mpi_num> ./cdhit-mpi -i /absolute/path/to/preprocess_output -o DB_output -T <threads_per_rank>
```

> `**mpirun -np` and `cdhit-mpi -T` MUST match `info.json`, otherwise the run may fail.**

---

## (5) Output files

- `*.txt`: representative sequence names and sequences
- `*.clstr`: cluster membership information

