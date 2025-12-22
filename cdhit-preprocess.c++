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
#include <fstream>
#include <dirent.h>
#include <sys/types.h>
#include <fstream>
#include <set>
#include <string>
#include <map>
#include <cstdio>
#undef min
#include <regex>

Options options;
SequenceDB seq_db;

static bool read_int_file(const std::string& path, int& val) {
    std::ifstream f(path);
    if (!f) return false;
    f >> val;
    return true;
}
////////////////////////////////////  MAIN /////////////////////////////////////
int main(int argc, char *argv[])
{
	string db_in;
	string db_out;
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
	
	sleep(0);
	
	int total_seqs = 0;
	// ***********************************    parse command line and open file
	if (argc < 5) print_usage_preprocess(argv[0]);
	if (options.SetOptions( argc, argv ) == 0) print_usage_preprocess(argv[0]);
	options.Validate();

	db_in = options.input;
	db_out = options.output;
	
	InitNAA( MAX_UAA );
	options.NAAN = NAAN_array[options.NAA];
	seq_db.NAAN = NAAN_array[options.NAA];
	//printf( "%i  %i  %i\n", sizeof(NVector<IndexCount>), seq_db.NAAN, sizeof(NVector<IndexCount>) * seq_db.NAAN );
    if (options.input.size()  == 0) bomb_error("no input file");
    if (options.NodeNum  == 0) bomb_error("no NodeNum");
    if (options.threads_per_node  == 0) bomb_error("no threads_per_node");
	    std::regex cpu_dir("^cpu([0-9]+)$");
    std::map<int, std::set<int>> socket_coreids;  // socket -> {core_id set}

    DIR* d = opendir("/sys/devices/system/cpu");
    if (!d) { perror("opendir"); return 1; }
    while (dirent* e = readdir(d)) {
        std::cmatch m;
        if (!std::regex_match(e->d_name, m, cpu_dir)) continue;
        std::string base = std::string("/sys/devices/system/cpu/") + e->d_name + "/topology";
        int pkg = -1, core = -1;
        if (!read_int_file(base + "/physical_package_id", pkg)) continue;
        if (!read_int_file(base + "/core_id", core)) continue;
        socket_coreids[pkg].insert(core);
    }
    closedir(d);

    // for (auto& [socket, cores] : socket_coreids) {
    //     printf("Socket %d: physical cores = %zu\n", socket, cores.size());
    // }
    int core_size = socket_coreids[0].size();
			// 外部排序
			size_t min_file_size = 512ull * 1024 * 1024;
			// seq_db.GenerateSorted_Parallel(db_in.c_str(), min_file_size , run_files,options);
			auto start = std::chrono::high_resolution_clock::now();
			seq_db.Pipeline_External_Sort(db_in.c_str(), min_file_size, run_files, options,core_size);
			mkdir(options.tmp_dir.c_str(), 0755);
			string temp_dir = options.tmp_dir;
			if (!temp_dir.empty() && temp_dir.back() != '/' && temp_dir.back() != '\\')
			{
				temp_dir += '/'; 
			}
			seq_db.MergeSortedRuns_KWay(run_files, temp_dir);
			auto end = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> elapsed = end - start;
			std::cout << "external sorting cost:    " << elapsed.count() << " second\n";
	
	return 0;

	

}  