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
	int chunks_size=seq_db.chunks_size;
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

	db_in = options.input;
	db_out = options.output;
	
	InitNAA( MAX_UAA );
	options.NAAN = NAAN_array[options.NAA];
	seq_db.NAAN = NAAN_array[options.NAA];
	//printf( "%i  %i  %i\n", sizeof(NVector<IndexCount>), seq_db.NAAN, sizeof(NVector<IndexCount>) * seq_db.NAAN );
	
	if (rank == 0) {
		// 元数据排序
		// seq_db.Read( db_in.c_str(), options ,meta_table);
        // seq_db.SortDivideMetaTable(meta_table, options );
		// seq_db.GenerateSortedRuns(db_in.c_str(), meta_table, 10 * 1024 * 1024, run_files);  // 50MB per chunk
		// seq_db.MergeRuns_Sequential(run_files, "final_sorted_1.fa");

		//外部排序
		seq_db.GenerateSorted_Parallel(db_in.c_str(), 500 * 1024 * 1024, run_files,options); 
		seq_db.MergeSortedRuns_KWay(run_files, "output/",size-1,50000);
		sleep(10);
	
		MPI_Barrier(MPI_COMM_WORLD);

	}


	else {
		

		// sleep(10);
		seq_db.read_sorted_files(rank,size);
		MPI_Barrier(MPI_COMM_WORLD);
		

    }
	MPI_Barrier(MPI_COMM_WORLD);
	seq_db.DoClustering_MPI(options, rank, master, worker, worker_rank);
	MPI_Barrier(MPI_COMM_WORLD);

	

}
