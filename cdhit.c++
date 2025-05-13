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
#include <mpi.h>
#include <iomanip>
#include <cmath>

Options options;
SequenceDB seq_db;


int digitCount(int num) {
    if (num == 0) return 1;
    return static_cast<int>(log10(abs(num))) + 1;
}

////////////////////////////////////  MAIN /////////////////////////////////////
int main(int argc, char* argv[])
{
	// 	modify by mgl version 1:
	// 	First, we initially set each chunk to have the same number of sequences 
	//		rather than roughly the same number of bytes

	int my_rank, rank_size;
	int worker_size;
	int worker_rank = -1;
	bool master = true;
	bool worker = false;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &rank_size);

	if(my_rank != 0){
		worker_rank = my_rank - 1;
		worker = true;
		master = false;
	}
	worker_size = rank_size - 1;

	string db_in;
	string db_out;

	if (argc < 5) print_usage(argv[0]);
	if (options.SetOptions(argc, argv) == 0)
		print_usage(argv[0]);
	options.Validate();

	db_in = options.input;
	db_out = options.output;
	InitNAA(MAX_UAA);
	options.NAAN = NAAN_array[options.NAA];
	seq_db.NAAN = NAAN_array[options.NAA];

	seq_db.Read(db_in.c_str(), options);
	size_t num_seqs = seq_db.sequences.size();
	if(master)
		cout << "total seq: " << num_seqs << endl;

	seq_db.SortDivide(options);

	MPI_Barrier(MPI_COMM_WORLD);
	// To solve for chunk_size, because in the previous CDHIT each table should not contain more than 10000 m values.
	// In the following code, 'node_chunks' indicates how many chunks should be allocated to a node,
	//		and 'chunk_size' indicates the size of a chunk.
	int node_chunks = 0, chunk_size = 0;
	node_chunks = num_seqs / (worker_size * 20000) + 1;
	if (num_seqs % (worker_size * node_chunks + 1))
		chunk_size = num_seqs / (worker_size * node_chunks + 1) + 1;
	else chunk_size = num_seqs / (worker_size * node_chunks + 1);

	vector<pair<int, int>>& chunks = seq_db.chunks;
	vector<int>& chunks_id = seq_db.chunks_id;

	if (master) {
		int num_chunks = node_chunks * worker_size + 1;
		cout << "\t>>> Chunk size: " << chunk_size << endl;
		cout << "\t>>> Total chunks: " << num_chunks << endl;
		cout << "\t>>> Number chunks in each nodes" << node_chunks << endl;
		int temp_target = 0;
		// TODO: Here perhaps OpenMP could be used for multithreading?
		chunks_id.resize(num_chunks, -1);
		chunks.resize(num_chunks);
		for (int i = 0;i < num_chunks;i++) {
			chunks_id[i] = i;
			chunks[i] = make_pair(i * chunk_size, min((i + 1) * chunk_size, num_seqs));
			if (i == 0) continue;
			int target = temp_target + 1;
			// First send the chunk_id of this chunk, stored in the vector chunk_id
			MPI_Send(&i, 1, MPI_INT, target, 0, MPI_COMM_WORLD);
			// Then send the begin and the end of this chunk,stored in vector chunk of that node
			MPI_Send(&chunks[i].first, 1, MPI_INT, target, 0, MPI_COMM_WORLD);
			MPI_Send(&chunks[i].second, 1, MPI_INT, target, 0, MPI_COMM_WORLD);
			// Update the temp_target
			temp_target = (temp_target + 1) % (rank_size - 1);
		}
		printf("In master node, we built [ %d ] chunks and allocated [ %d ] chunks to each worker node.\n",
			num_chunks, node_chunks);
		MPI_Barrier(MPI_COMM_WORLD);
		cerr << ">>> Splitting of the chunk has been completed!" << endl;
	}

	if (worker) {
		int source = 0;
		chunks_id.resize(node_chunks, -1);
		chunks.resize(node_chunks);
		for (int i = 0;i < node_chunks;i++) {
			int receive_id = 0, begin = 0, end = 0;
			MPI_Recv(&receive_id, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			chunks_id[i] = receive_id;
			MPI_Recv(&begin, 1, MPI_INT, source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Recv(&end, 1, MPI_INT, source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			chunks[i] = make_pair(begin, end);
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	int num_chunks = node_chunks * worker_size + 1;
	for (int p = 0;p < rank_size;p++) {
		MPI_Barrier(MPI_COMM_WORLD);
		if (worker_rank == p-1 && my_rank != 0) {
			cout << "Worker " << worker_rank << " received chunks: " << endl;
			for (int i = 0;i < node_chunks;i++) {
				cout << "    Chunk id: " << setw(digitCount(num_chunks) + 1) << chunks_id[i]
					<< "  Begin: " << setw(digitCount(num_seqs) + 1) << chunks[i].first
					<< "  End: " << setw(digitCount(num_seqs) + 1) << chunks[i].second
					<< "  Chunk_size: " << setw(digitCount(chunk_size) + 1) << chunks[i].second - chunks[i].first << endl;
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	seq_db.DoClustering_MPI(options, my_rank, master, worker, worker_rank);
	MPI_Barrier(MPI_COMM_WORLD);
	if (master) {
		cout << "Cluster is Finished" << endl;
		printf("writing new database\n");
		seq_db.WriteClusters(db_in.c_str(), db_out.c_str(), options);

		// write a backup clstr file in case next step crashes
		seq_db.WriteExtra1D(options);
		cout << "program completed !" << endl << endl;
	}
	// MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	return 0;
} // END int main
