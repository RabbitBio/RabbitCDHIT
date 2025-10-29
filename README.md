# `RabbitCD-HIT`
RabbitCD-HIT is a scalable greedy incremental clustering tool designed for massive protein sequences.
It supports running on distributed clusters to efficiently process massive protein sequence datasets.
RabbitCD-HIT combines inter-node MPI parallelism with intra-node multithreading to achieve superior speed and scalability while maintaining full compatibility with the CD-HIT output format.

#### Compile and install
```bash
git clone --recursive https://github.com/RabbitBio/RabbitCDHIT.git
cd RabbitTClust
make
```
## Usage
```bash
Usage: ./cdhit-preprocess [Options]
#cdhit-preprocess, external sorting for protein sequences
Options

   -i   input filename in fasta format, required, can be in .gz format
   -o   output filename, required
   -N   number of nodes, required
   -NT  number of threads per node, required
   -tmp output directory for sequences after external sorting, default tmp
   -T   number of threads, default 1; with 0, all CPUs will be used
   -l   length of throw_away_sequences, default 10
   -h   print this help


Usage: ./cdhit-mpi [Options]
#chit-mpi, distributed clustering
Options

   -c   sequence identity threshold, default 0.9
        this is the default cd-hit's "global sequence identity" calculated as:
        number of identical amino acids or bases in alignment
        divided by the full length of the shorter sequence
   -G   use global sequence identity, default 1
        if set to 0, then use local sequence identity, calculated as :
        number of identical amino acids or bases in alignment
        divided by the length of the alignment
        NOTE!!! don't use -G 0 unless you use alignment coverage controls
        see options -aL, -AL, -aS, -AS
   -b   band_width of alignment, default 20
   -T   number of threads per mpi process, required, must be consistent with the JSON file.
   -n   word_length, default 5, see user's guide for choosing it
   -t   tolerance for redundance, default 2
   -s   length difference cutoff, default 0.0
        if set to 0.9, the shorter sequences need to be
        at least 90% length of the representative of the cluster
   -S   length difference cutoff in amino acid, default 999999
        if set to 60, the length difference between the shorter sequences
        and the representative of the cluster can not be bigger than 60
   -aL  alignment coverage for the longer sequence, default 0.0
        if set to 0.9, the alignment must covers 90% of the sequence
   -AL  alignment coverage control for the longer sequence, default 99999999
        if set to 60, and the length of the sequence is 400,
        then the alignment must be >= 340 (400-60) residues
   -aS  alignment coverage for the shorter sequence, default 0.0
        if set to 0.9, the alignment must covers 90% of the sequence
   -AS  alignment coverage control for the shorter sequence, default 99999999
        if set to 60, and the length of the sequence is 400,
        then the alignment must be >= 340 (400-60) residues
   -A   minimal alignment coverage control for the both sequences, default 0
        alignment must cover >= this value for both sequences
   -uL  maximum unmatched percentage for the longer sequence, default 1.0
        if set to 0.1, the unmatched region (excluding leading and tailing gaps)
        must not be more than 10% of the sequence
   -uS  maximum unmatched percentage for the shorter sequence, default 1.0
        if set to 0.1, the unmatched region (excluding leading and tailing gaps)
        must not be more than 10% of the sequence
   -U   maximum unmatched length, default 99999999
        if set to 10, the unmatched region (excluding leading and tailing gaps)
        must not be more than 10 bases
   -g   1 or 0, default 0
        by cd-hit's default algorithm, a sequence is clustered to the first
        cluster that meet the threshold (fast cluster). If set to 1, the program
        will cluster it into the most similar cluster that meet the threshold
        (accurate but slow mode)
        but either 1 or 0 won't change the representatives of final clusters
   -h   print this help
```
## Example:
```bash
# The clustering task runs on a single node with 128 threads.
# The output files from external sorting are stored in the tmp directory.  
# The external sorting process uses 64 threads.  
 ./cdhit-preprocess -i DB.fa  -T 64 -N 1 -NT 128 -tmp tmp
```
After preprocessing, a JSON file will be generated in the tmp directory:
```json
 {
    "files": {
        "out_prefix": "_proc",
        "output_dir": "huge_tmp/"
    },
    "info": {
        "chunk_bytes": 35900000,
        "chunk_size": 100000,
        "chunks_num": 18,
        "first_chunk_size": 2000,
        "len_n50": 455,
        "max_idf": 14,
        "max_len": 35808,
        "min_len": 11,
        "num_procs": 3,
        "threads_per_node": 32,
        "total_chunk": 17,
        "total_desc": 20839413,
        "total_letter": 600968830,
        "total_mpi_num": 4,
        "total_num": 1670150
    }
}
```

Here, total_mpi_num and threads_per_node correspond to the number of MPI processes and threads required for clustering,
which are equivalent to the parameters -np and -T, respectively.
The -tmp option should be set to the same directory specified above.
```bash
mpirun -np 4 ./cdhit-mpi -o output -T 32 -M 0 -tmp huge_tmp
```

