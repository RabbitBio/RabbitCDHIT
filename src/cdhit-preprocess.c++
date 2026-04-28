

#include <dirent.h>
#include <sys/types.h>

#include <cstdio>
#include <fstream>
#include <map>
#include <numeric>
#include <set>
#include <string>

#include "cdhit-common.h"
#undef min
#undef max
#include <sys/stat.h>

#include <regex>

Options options;
SequenceDB seq_db;

static bool read_text_file(const string& path, string& out) {
    ifstream ifs(path);
    if (!ifs) return false;
    getline(ifs, out);
    return true;
}
static bool read_int_file(const string& path, int& out) {
    ifstream ifs(path);
    if (!ifs) return false;
    ifs >> out;
    return bool(ifs);
}
static vector<int> parse_cpulist(const string& s) {
    vector<int> cpus;
    size_t i = 0;
    while (i < s.size()) {
        while (i < s.size() && (s[i] == ' ' || s[i] == '\t' || s[i] == '\n' || s[i] == ',')) i++;
        if (i >= s.size()) break;
        int a = 0; while (i < s.size() && isdigit((unsigned char)s[i])) { a = a*10 + (s[i]-'0'); i++; }
        int b = a;
        if (i < s.size() && s[i] == '-') {
            i++; b = 0;
            while (i < s.size() && isdigit((unsigned char)s[i])) { b = b*10 + (s[i]-'0'); i++; }
        }
        if (a <= b) for (int x = a; x <= b; x++) cpus.push_back(x);
        else        for (int x = a; x >= b; x--) cpus.push_back(x);
    }
    sort(cpus.begin(), cpus.end());
    cpus.erase(unique(cpus.begin(), cpus.end()), cpus.end());
    return cpus;
}
////////////////////////////////////  MAIN /////////////////////////////////////
int main(int argc, char* argv[]) {
    // sleep(0);
    string db_in;
    string db_out;
    std::vector<std::string> run_files;
    vector<pair<int, int>>& all_chunks = seq_db.all_chunks;
    vector<pair<int, int>>& my_chunks = seq_db.my_chunks;
    vector<int>& chunks_id = seq_db.chunks_id;
    int total_chunk = seq_db.total_chunk;
    float begin_time = current_time();
    float end_time;
    bool master = true;
    bool worker = false;
    int worker_rank = -1;

    int total_seqs = 0;
    // ***********************************    parse command line and open file
    if (argc < 5) print_usage_preprocess(argv[0]);
    if (options.SetOptions(argc, argv,0) == 0) print_usage_preprocess(argv[0]);
    options.Validate();

    db_in = options.input;
    db_out = options.output;

    InitNAA(MAX_UAA);
    options.NAAN = NAAN_array[options.NAA];
    seq_db.NAAN = NAAN_array[options.NAA];

    if (options.input.size() == 0) bomb_error("no input file");
    if (options.NodeNum == 0) bomb_error("no NodeNum");
    if (options.threads_per_node == 0) bomb_error("no threads_per_node");
    const char *RUN_DIR = "output";
    if (mkdir(RUN_DIR, 0755) != 0 && errno != EEXIST)
    {
        perror("mkdir output");
        return 0;
    }
    int core_size;
    if (options.core_size)
        core_size = options.core_size;
    else
    {
        std::regex cpu_dir("^cpu([0-9]+)$");
        std::map<int, std::set<int>> socket_coreids;

        DIR *d = opendir("/sys/devices/system/cpu");
        if (!d)
        {
            perror("opendir");
            return 1;
        }
        while (dirent *e = readdir(d))
        {
            std::cmatch m;
            if (!std::regex_match(e->d_name, m, cpu_dir))
                continue;
            std::string base = std::string("/sys/devices/system/cpu/") + e->d_name + "/topology";
            int pkg = -1, core = -1;
            if (!read_int_file(base + "/physical_package_id", pkg))
                continue;
            if (!read_int_file(base + "/core_id", core))
                continue;
            socket_coreids[pkg].insert(core);
        }
        closedir(d);
        core_size = socket_coreids[0].size();
    }
    int numa_size;
    if (options.numa_size)
        numa_size = options.numa_size;
    else
    {
        int node = 0;

        string cpulist;
        if (!read_text_file("/sys/devices/system/node/node" + to_string(node) + "/cpulist", cpulist))
        {
            cerr << "cannot read cpulist\n";
            return 1;
        }
        auto cpus = parse_cpulist(cpulist);

        // 统计 unique physical cores：用 (package_id, core_id) 去重最稳
        set<pair<int, int>> unique_cores;

        for (int cpu : cpus)
        {
            string base = "/sys/devices/system/cpu/cpu" + to_string(cpu) + "/topology";
            int pkg = -1, core = -1;
            if (!read_int_file(base + "/physical_package_id", pkg))
                continue;
            if (!read_int_file(base + "/core_id", core))
                continue;
            unique_cores.insert({pkg, core});
        }

        // cout << "NUMA node " << node << " logical CPUs  = " << cpus.size() << "\n";
        cout << " physical cores = " << unique_cores.size() << "\n";
        numa_size = unique_cores.size();
    }
    #ifdef DEBUG
    cerr<<"core_size "<<core_size<<" numa_size "<<numa_size<<endl;
    #endif
    size_t min_file_size = 512ull * 1024 * 1024;

    auto start = std::chrono::high_resolution_clock::now();
    seq_db.Pipeline_External_Sort(db_in.c_str(), min_file_size, run_files, options, core_size,numa_size);

    mkdir(options.preprocess_dir.c_str(), 0755);
    string preprocess_output_dir = options.preprocess_dir;
    if (!preprocess_output_dir.empty() && preprocess_output_dir.back() != '/' &&
        preprocess_output_dir.back() != '\\') {
        preprocess_output_dir += '/';
    }
    cout << "preprocess_output_dir: " << preprocess_output_dir << endl;
    seq_db.MergeSortedRuns_KWay(run_files, preprocess_output_dir, options.threads);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "external sorting cost:    " << elapsed.count() << " second\n";

    return 0;
}