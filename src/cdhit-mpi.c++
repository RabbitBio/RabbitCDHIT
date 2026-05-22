

#include <mpi.h>

#include <fstream>
#include <numeric>

#include "cdhit-common.h"

Options options;
SequenceDB seq_db;

////////////////////////////////////  MAIN /////////////////////////////////////
int main(int argc, char* argv[]) {
    sleep(0);
    string db_in;
    string db_out;
    vector<SequenceMeta> meta_table;
    std::vector<std::string> run_files;
    vector<pair<int, int>>& all_chunks = seq_db.all_chunks;
    vector<pair<int, int>>& my_chunks = seq_db.my_chunks;
    vector<int>& chunks_id = seq_db.chunks_id;
    auto start = std::chrono::high_resolution_clock::now();
    bool master = true;
    bool worker = false;

    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    if (provided < MPI_THREAD_MULTIPLE) {
        cerr << "ERROR: Need MPI_THREAD_MULTIPLE" << endl;
        return 1;
    }
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm worker_comm;
    if (rank != 0) {
        MPI_Comm_split(MPI_COMM_WORLD, 1, rank, &worker_comm);
    } else {
        MPI_Comm_split(MPI_COMM_WORLD, 0, rank, &worker_comm);
    }
    int worker_rank, worker_size;
    MPI_Comm_rank(worker_comm, &worker_rank);
    MPI_Comm_size(worker_comm, &worker_size);
    if (rank != 0) {
        worker = true;
        master = false;
    }
    int total_seqs = 0;
    // ***********************************    parse command line and open file
    if (argc < 5) print_usage_mpi(argv[0]);
    if (options.SetOptions(argc, argv,rank) == 0) print_usage_mpi(argv[0]);
    if (options.input.empty()) bomb_error("no input directory for cdhit-mpi (-i)");
    if (options.output.empty()) bomb_error("no output prefix for cdhit-mpi (-o)");
    // cdhit-mpi: use -i as preprocess input directory.
    if (!options.input.empty()) {
        options.preprocess_dir = options.input;
    }
    options.Validate();

    // The purpose is merely to facilitate git.
    // const char* RUN_DIR = "output";
    // std::string tmp_prefix = std::string(RUN_DIR);
    // if (tmp_prefix.back() != '/') {
    //     tmp_prefix += '/';
    // }
    if (!options.output.empty() && options.output[0] == '/') {
        options.output = options.output.substr(1);
    }
    // options.output = tmp_prefix + options.output;
    // ------------------------------------------------------------

    if (options.output.size() == 0) bomb_error("no output file");
    // Modify by MGL: remove the rank directory
    // Place both the representative sequence file and the specific cluster file in the root directory of output.
    // if (master) {
    //     mkdir(options.output.c_str(), 0755);
    // }
    db_in = options.input;
    db_out = options.output;

    InitNAA(MAX_UAA);
    options.NAAN = NAAN_array[options.NAA];
    seq_db.NAAN = NAAN_array[options.NAA];

    string preprocess_output_dir = options.preprocess_dir;
    if (!preprocess_output_dir.empty() && preprocess_output_dir.back() != '/' &&
        preprocess_output_dir.back() != '\\') {
        preprocess_output_dir += '/';
    }
    seq_db.ReadJsonInfo("info.json", preprocess_output_dir, options, master);
    if (size != seq_db.total_mpi_num) 
    {
        #ifdef DEBUG
        cerr<<"size "<<size<<endl;
        cerr<<"seq_db.total_mpi_num "<<seq_db.total_mpi_num<<endl;
        #endif
bomb_error("Number of processes does not match");
    }

    if (!master) {
        seq_db.read_sorted_files(preprocess_output_dir, rank, size, false, worker_comm, options);
        MPI_Barrier(worker_comm);
    }

    seq_db.DoClustering_MPI(options, rank, master, worker, worker_rank, db_out.c_str(), worker_comm);
    MPI_Barrier(MPI_COMM_WORLD);
    if (master) {
        
        cout << "Cluster is Finished" << endl;
        printf("writing new database\n");
        cout << "program completed !" << endl << endl;
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        std::cout << "Cluster time:    " << elapsed.count() << " second\n";

    }
    MPI_Barrier(MPI_COMM_WORLD);
    // cerr<<"this"<<endl;
    MPI_Finalize();
    return 0;
}
