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
	int worker_rank = -1;
	
	//初始化MPI
	MPI_Init(&argc, &argv);
	int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
	if(rank != 0){
		worker_rank = rank - 1;
		worker = true;
		master = false;
	}
	int total_seqs = 0;
	// ***********************************    parse command line and open file
	if (argc < 5) print_usage(argv[0]);
	if (options.SetOptions( argc, argv ) == 0) print_usage(argv[0]);
	options.Validate();
	if (options.output.size()  == 0) bomb_error("no output file");
	
	db_in = options.input;
	db_out = options.output;
	
	InitNAA( MAX_UAA );
	options.NAAN = NAAN_array[options.NAA];
	seq_db.NAAN = NAAN_array[options.NAA];
	//printf( "%i  %i  %i\n", sizeof(NVector<IndexCount>), seq_db.NAAN, sizeof(NVector<IndexCount>) * seq_db.NAAN );

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
		seq_db.read_sorted_files(temp_dir,rank, size, false);
	}

	// if (rank == 0) {
		


	// 	//外部排序
	// 	 auto start = std::chrono::high_resolution_clock::now();
	// 	seq_db.GenerateSorted_Parallel(db_in.c_str(), 500 * 1024 * 1024, run_files,options); 
		
	// 	seq_db.MergeSortedRuns_KWay(run_files, "output/",size-1);
	// 	auto end = std::chrono::high_resolution_clock::now();
    //     std::chrono::duration<double> elapsed = end - start;
    //     std::cout << "外部排序耗时:    " << elapsed.count() << " 秒\n";
	// }

  
	// else {
		
	// 	seq_db.read_sorted_files(rank,size);
		
    // }
	
	seq_db.DoClustering_MPI(options, rank, master, worker, worker_rank,db_out.c_str());
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