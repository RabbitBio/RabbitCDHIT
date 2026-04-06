#include <sys/stat.h>
#include <unistd.h>

#include <cerrno>
#include <iostream>
#include <string>

#include "cdhit-common.h"

Options options;
SequenceDB seq_db;

int main(int argc, char *argv[]) {
    sleep(0);

    if (argc < 5) print_usage_mpi(argv[0]);
    if (options.SetOptions(argc, argv, 0) == 0) print_usage_mpi(argv[0]);
    if (options.input.empty()) bomb_error("no input directory for master-kernel-test (-i)");
    if (options.output.empty()) bomb_error("no output prefix for master-kernel-test (-o)");

    // Reuse cdhit-mpi convention: -i is preprocess output directory.
    options.preprocess_dir = options.input;
    options.Validate();

    const char *RUN_DIR = "output";
    if (mkdir(RUN_DIR, 0755) != 0 && errno != EEXIST) {
        perror("mkdir output");
        return 1;
    }

    std::string tmp_prefix = std::string(RUN_DIR);
    if (tmp_prefix.back() != '/') tmp_prefix += '/';
    if (!options.output.empty() && options.output[0] == '/') options.output = options.output.substr(1);
    options.output = tmp_prefix + options.output;

    InitNAA(MAX_UAA);
    options.NAAN = NAAN_array[options.NAA];
    seq_db.NAAN = NAAN_array[options.NAA];

    std::string preprocess_output_dir = options.preprocess_dir;
    if (!preprocess_output_dir.empty() && preprocess_output_dir.back() != '/' &&
        preprocess_output_dir.back() != '\\') {
        preprocess_output_dir += '/';
    }

    seq_db.ReadJsonInfo("info.json", preprocess_output_dir, options, true);
    seq_db.DoClustering_MasterKernel_Test(options, options.output.c_str());

    std::cout << "Master kernel test completed." << std::endl;
    return 0;
}
