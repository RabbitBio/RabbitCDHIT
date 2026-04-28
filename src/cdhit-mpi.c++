

#include <mpi.h>

#include <fstream>
#include <numeric>
#include <sys/stat.h>

#include "cdhit-common.h"

Options options;
SequenceDB seq_db;

// ---------------------------------------------------------------------------
// Auto-detect whether memory_control should be enabled.
//
// Reads /proc/meminfo (Linux) to get MemAvailable, then compares the
// estimated per-worker memory requirement against the per-worker memory
// budget.  Returns true if memory_control should be turned on.
//
//   file_path    : worker's _proc{rank-1}.fa file
//   num_workers  : total number of worker MPI processes (used as divisor for
//                  the node's available memory – conservative on multi-node
//                  runs but always safe)
//   safety_factor: fraction of MemAvailable considered usable (default 0.80)
// ---------------------------------------------------------------------------
static bool auto_detect_memory_control(const std::string &file_path,
                                        int num_workers,
                                        double safety_factor = 0.80)
{
    // 1. Estimate memory requirement from file size.
    //    File ≈ raw sequence bytes; add ~50 % for Sequence structs + identifiers.
    struct stat st;
    if (stat(file_path.c_str(), &st) != 0)
        return false;  // can't stat – skip detection
    long long estimated_bytes = (long long)(st.st_size * 1.5);

    // 2. Query available system memory from /proc/meminfo (Linux).
    long long avail_bytes = 0;
#ifdef __linux__
    FILE *mf = fopen("/proc/meminfo", "r");
    if (mf) {
        char key[64];
        long long val;
        while (fscanf(mf, "%63s %lld %*s\n", key, &val) >= 2) {
            if (strcmp(key, "MemAvailable:") == 0) {
                avail_bytes = val * 1024LL;   // kB → bytes
                break;
            }
        }
        fclose(mf);
    }
#endif
    if (avail_bytes == 0)
        return false;  // can't determine available memory – skip detection

    // 3. Each worker gets an equal share of the node's available memory.
    long long per_worker_avail =
        (long long)(avail_bytes * safety_factor) / max(num_workers, 1);

    return estimated_bytes > per_worker_avail;
}

////////////////////////////////////  MAIN /////////////////////////////////////
int main(int argc, char* argv[]) {
    // sleep(0);
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
        // Auto-detect whether memory is tight enough to warrant memory_control.
        // Always runs unless the user explicitly passed -memory_control <0|1>.
        // If memory is tight AND -stealing is active the program aborts, because
        // the two modes are hard mutual exclusive.
        if (!options.memory_control_explicit) {
            std::string proc_file = preprocess_output_dir + "_proc"
                                    + std::to_string(rank - 1) + ".fa";
            // Use workers_per_node from info.json so the divisor matches the
            // number of processes actually sharing this node's RAM.  Fall back
            // to total workers if the field is absent (older info.json).
            int num_workers = (seq_db.workers_per_node > 0)
                                  ? seq_db.workers_per_node
                                  : (size - 1);
            if (auto_detect_memory_control(proc_file, num_workers)) {
                if (options.stealing) {
                    fprintf(stderr,
                            "[rank %d] ERROR: memory is insufficient for the current "
                            "workload but -stealing is active. "
                            "-memory_control and -stealing are mutually exclusive. "
                            "Either disable -stealing or pass -memory_control 0 to "
                            "skip this check.\n", rank);
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
                options.memory_control = true;
                if (rank == 1)
                    fprintf(stderr,
                            "[INFO] memory_control auto-enabled: estimated per-worker "
                            "requirement exceeds available memory budget "
                            "(nodes=%d, workers_per_node=%d). "
                            "Use -memory_control 0 to suppress this check.\n",
                            seq_db.num_nodes, num_workers);
            }
        }
        // Hard mutual-exclusion guard: Validate() catches the explicit case;
        // the auto-detect block above catches the implicit case.
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
