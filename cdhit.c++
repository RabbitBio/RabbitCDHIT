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
	float begin_time = current_time();
	float end_time;
	//初始化MPI
	MPI_Init(&argc, &argv);
	int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
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
	vector<Sequence*> all_sequences;
	vector<vector<Sequence*>> chunks;
	//printf( "%i  %i  %i\n", sizeof(NVector<IndexCount>), seq_db.NAAN, sizeof(NVector<IndexCount>) * seq_db.NAAN );
	
	if (rank == 0) {
		// 元数据排序
		// seq_db.Read( db_in.c_str(), options ,meta_table);
        // seq_db.SortDivideMetaTable(meta_table, options );
		// seq_db.GenerateSortedRuns(db_in.c_str(), meta_table, 10 * 1024 * 1024, run_files);  // 50MB per chunk
		// seq_db.MergeRuns_Sequential(run_files, "final_sorted_1.fa");

		//外部排序
		seq_db.GenerateSorted_Parallel(db_in.c_str(), 500 * 1024 * 1024, run_files,options); 
		seq_db.MergeSortedRuns_KWay(run_files, "output/",4,10000);
		std::cout << "Rank " << rank << " entered Barrier" << std::endl;
		MPI_Barrier(MPI_COMM_WORLD);
		std::cout << "Rank " << rank << " exiting normally." << std::endl;
		MPI_Finalize();
		return 0;
	}


	else {
		std::cout << "Rank " << rank << " entered Barrier" << std::endl;
		MPI_Barrier(MPI_COMM_WORLD);  // 等所有进程都输出完
		std::cout << "Rank " << rank << " exiting normally." << std::endl;
		MPI_Finalize();
		exit(0);

    }
	


	

// 	MPI_Barrier(MPI_COMM_WORLD);  // 等所有进程同步

// // 让每个进程都输出一个信息以确认它们活着
// printf("Rank %d finished preprocessing\n", rank);

// MPI_Finalize();
// return 0;
	// seq_db.Read( db_in.c_str(), options );
	// cout << "total seq: " << seq_db.sequences.size() << endl;

	// seq_db.SortDivide( options );

	// seq_db.DoClustering( options );

	// printf( "writing new database\n" );
	// seq_db.WriteClusters( db_in.c_str(), db_out.c_str(), options );

	// // write a backup clstr file in case next step crashes
	// seq_db.WriteExtra1D( options );
	// cout << "program completed !" << endl << endl;
	// end_time = current_time();
	// printf( "Total CPU time %.2f\n", end_time - begin_time );
	// return 0;
} // END int main
