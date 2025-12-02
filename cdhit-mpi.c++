// =============================================================================
// CD-HIT
// http://cd-hit.org/
// http://bioinformatics.burnham-inst.org/cd-hi
//
// program written by
//                    Weizhong Li
//                    UCSD, San Diego Supercomputer Center
//                    La Jolla, CA, 92093
//                    Email liwz@sdsc.edu
//                 at
//                    Adam Godzik's lab
//                    The Burnham Institute
//                    La Jolla, CA, 92037
//                    Email adam@burnham-inst.org
//
// Modified by:
//                    Limin Fu
//                    Center for Research in Biological Systems (CRBS), UCSD
//                    La Jolla, CA, 92093
//                    Email: l2fu@ucsd.edu, fu@daovm.net
// =============================================================================

#include "cdhit-common.h"
#include <numeric>
#include <mpi.h>
#include <fstream>

Options options;
SequenceDB seq_db;



////////////////////////////////////  MAIN /////////////////////////////////////
int main(int argc, char *argv[])
{
	string db_in;
	string db_out;
	vector<SequenceMeta> meta_table;
	std::vector<std::string> run_files;
	vector<pair<int, int>>& all_chunks = seq_db.all_chunks;
	vector<pair<int, int>>& my_chunks = seq_db.my_chunks;
	vector<int>& chunks_id = seq_db.chunks_id;
	int total_chunk=seq_db.total_chunk;
	float begin_time = current_time();
	float end_time;
	bool master = true;
	bool worker = false;
	// sleep(0);
	//初始化MPI
	MPI_Init(&argc, &argv);
	int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm worker_comm;
    if (rank != 0) {
        MPI_Comm_split(MPI_COMM_WORLD, 1, rank, &worker_comm);  // 工作进程（非主进程）创建一个新的group
    } else {
        MPI_Comm_split(MPI_COMM_WORLD, 0, rank, &worker_comm);  // 主进程不属于工作进程组
    }
	int worker_rank, worker_size;
    MPI_Comm_rank(worker_comm, &worker_rank);
	MPI_Comm_size(worker_comm, &worker_size);
	if(rank != 0){
		worker = true;
		master = false;
	}
	int total_seqs = 0;
	// ***********************************    parse command line and open file
	if (argc < 5) print_usage_mpi(argv[0]);
	if (options.SetOptions( argc, argv ) == 0) print_usage_mpi(argv[0]);
	options.Validate();
	if (options.output.size()  == 0) bomb_error("no output file");
	
	db_in = options.input;
	db_out = options.output;
	
	InitNAA( MAX_UAA );
	options.NAAN = NAAN_array[options.NAA];
	seq_db.NAAN = NAAN_array[options.NAA];
	//printf( "%i  %i  %i\n", sizeof(NVector<IndexCount>), seq_db.NAAN, sizeof(NVector<IndexCount>) * seq_db.NAAN );
	// MPI_Barrier(MPI_COMM_WORLD);
	string temp_dir = options.tmp_dir;
	if (!temp_dir.empty() && temp_dir.back() != '/' && temp_dir.back() != '\\')
	{
		temp_dir += '/';
	}
	seq_db.ReadJsonInfo("info.json", temp_dir, options, master);
	if (size != seq_db.total_mpi_num)
		bomb_error("Number of processes does not match");
	if (options.threads != seq_db.Production_threads)
		bomb_error("Number of threads does not match");
	if (!master)
	{
		seq_db.read_sorted_files(temp_dir,rank, size, false,worker_comm,options);
		MPI_Barrier(worker_comm);
		// exit(0);
	}

	// if (rank == 0) {
		


	
	seq_db.DoClustering_MPI(options, rank, master, worker, worker_rank,db_out.c_str(),worker_comm);
	MPI_Barrier(MPI_COMM_WORLD);
	if (master) {
		cout << "Cluster is Finished" << endl;
		// seq_db.checkRepSeqs();
		printf("writing new database\n");
		// seq_db.WriteClustersSort(db_in.c_str(), db_out.c_str(), options);
		// seq_db.WriteClusterDetail(options);
		cout << "program completed !" << endl << endl;
	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	return 0;

	

}  