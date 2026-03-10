// =============================================================================
// CD-HI/CD-HIT
//
// Cluster Database at High Identity
//
// CD-HIT clusters protein sequences at high identity threshold.
// This program can remove the high sequence redundance efficiently.
//
// program written by
//                    Weizhong Li
//                    UCSD, San Diego Supercomputer Center
//                    La Jolla, CA, 92093
//                    Email liwz@sdsc.edu
//
//                 at
//                    Adam Godzik's lab
//                    The Burnham Institute
//                    La Jolla, CA, 92037
//                    Email adam@burnham-inst.org
//
// modified by:
//                    Limin Fu
//                    Center for Research in Biological Systems (CRBS), UCSD
//                    La Jolla, CA, 92093
//                    Email: l2fu@ucsd.edu, fu@daovm.net
// ==========================================================================

#include "cdhit-common.h"

#include <assert.h>
#include <immintrin.h>
#include <limits.h>
#include <mpi.h>
#include <stdint.h>
#include <sys/time.h>
#include <unistd.h>

#include <valarray>
#ifndef NO_OPENMP

#include <omp.h>
// #include<mpi.h>

#define WITH_OPENMP " (+OpenMP)"
#define WITH_OPEMMP_AND_MPI " (+OpenMP and MPI)"
#else

#define WITH_OPENMP ""

#define omp_set_num_threads(T) (T = T)
#define omp_get_thread_num() 0

#endif
#define MAX_AA_2 529
// class function definition
const char aa[] = {"ARNDCQEGHILKMFPSTWYVBZX"};
//{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,2,6,20};
int aa2idx[] = {0, 2, 4, 3, 6, 13, 7, 8, 9, 20, 11, 10, 12, 2, 20, 14, 5, 1, 15, 16, 20, 19, 17, 20, 18, 6};
// idx for  A  B  C  D  E  F  G  H  I  J  K  L  M  N  O  P
//          Q  R  S  T  U  V  W  X  Y  Z
// so  aa2idx[ X - 'A'] => idx_of_X, eg aa2idx['A' - 'A'] => 0,
// and aa2idx['M'-'A'] => 12
std::atomic<bool> progress_running{true};
double get_time() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double) tv.tv_sec + (double) tv.tv_usec / 1000000;
}
alignas(256) static const uint8_t REVERSE_MASK_LUT[256] = {
    0x00, 0x80, 0x40, 0xC0, 0x20, 0xA0, 0x60, 0xE0, 0x10, 0x90, 0x50, 0xD0, 0x30, 0xB0, 0x70, 0xF0, 0x08, 0x88, 0x48,
    0xC8, 0x28, 0xA8, 0x68, 0xE8, 0x18, 0x98, 0x58, 0xD8, 0x38, 0xB8, 0x78, 0xF8, 0x04, 0x84, 0x44, 0xC4, 0x24, 0xA4,
    0x64, 0xE4, 0x14, 0x94, 0x54, 0xD4, 0x34, 0xB4, 0x74, 0xF4, 0x0C, 0x8C, 0x4C, 0xCC, 0x2C, 0xAC, 0x6C, 0xEC, 0x1C,
    0x9C, 0x5C, 0xDC, 0x3C, 0xBC, 0x7C, 0xFC, 0x02, 0x82, 0x42, 0xC2, 0x22, 0xA2, 0x62, 0xE2, 0x12, 0x92, 0x52, 0xD2,
    0x32, 0xB2, 0x72, 0xF2, 0x0A, 0x8A, 0x4A, 0xCA, 0x2A, 0xAA, 0x6A, 0xEA, 0x1A, 0x9A, 0x5A, 0xDA, 0x3A, 0xBA, 0x7A,
    0xFA, 0x06, 0x86, 0x46, 0xC6, 0x26, 0xA6, 0x66, 0xE6, 0x16, 0x96, 0x56, 0xD6, 0x36, 0xB6, 0x76, 0xF6, 0x0E, 0x8E,
    0x4E, 0xCE, 0x2E, 0xAE, 0x6E, 0xEE, 0x1E, 0x9E, 0x5E, 0xDE, 0x3E, 0xBE, 0x7E, 0xFE, 0x01, 0x81, 0x41, 0xC1, 0x21,
    0xA1, 0x61, 0xE1, 0x11, 0x91, 0x51, 0xD1, 0x31, 0xB1, 0x71, 0xF1, 0x09, 0x89, 0x49, 0xC9, 0x29, 0xA9, 0x69, 0xE9,
    0x19, 0x99, 0x59, 0xD9, 0x39, 0xB9, 0x79, 0xF9, 0x05, 0x85, 0x45, 0xC5, 0x25, 0xA5, 0x65, 0xE5, 0x15, 0x95, 0x55,
    0xD5, 0x35, 0xB5, 0x75, 0xF5, 0x0D, 0x8D, 0x4D, 0xCD, 0x2D, 0xAD, 0x6D, 0xED, 0x1D, 0x9D, 0x5D, 0xDD, 0x3D, 0xBD,
    0x7D, 0xFD, 0x03, 0x83, 0x43, 0xC3, 0x23, 0xA3, 0x63, 0xE3, 0x13, 0x93, 0x53, 0xD3, 0x33, 0xB3, 0x73, 0xF3, 0x0B,
    0x8B, 0x4B, 0xCB, 0x2B, 0xAB, 0x6B, 0xEB, 0x1B, 0x9B, 0x5B, 0xDB, 0x3B, 0xBB, 0x7B, 0xFB, 0x07, 0x87, 0x47, 0xC7,
    0x27, 0xA7, 0x67, 0xE7, 0x17, 0x97, 0x57, 0xD7, 0x37, 0xB7, 0x77, 0xF7, 0x0F, 0x8F, 0x4F, 0xCF, 0x2F, 0xAF, 0x6F,
    0xEF, 0x1F, 0x9F, 0x5F, 0xDF, 0x3F, 0xBF, 0x7F, 0xFF};

alignas(32) static const int32_t LEFT_SHIFT_TABLE[8][8] = {
    {0, 1, 2, 3, 4, 5, 6, 7}, {7, 0, 1, 2, 3, 4, 5, 6}, {7, 7, 0, 1, 2, 3, 4, 5}, {7, 7, 7, 0, 1, 2, 3, 4},
    {7, 7, 7, 7, 0, 1, 2, 3}, {7, 7, 7, 7, 7, 0, 1, 2}, {7, 7, 7, 7, 7, 7, 0, 1}, {7, 7, 7, 7, 7, 7, 7, 0}};

alignas(32) static const int32_t RIGHT_SHIFT_TABLE[8][8] = {
    {7, 6, 5, 4, 3, 2, 1, 0}, {6, 5, 4, 3, 2, 1, 0, 7}, {5, 4, 3, 2, 1, 0, 7, 7}, {4, 3, 2, 1, 0, 7, 7, 7},
    {3, 2, 1, 0, 7, 7, 7, 7}, {2, 1, 0, 7, 7, 7, 7, 7}, {1, 0, 7, 7, 7, 7, 7, 7}, {0, 7, 7, 7, 7, 7, 7, 7}};
alignas(32) static const int32_t SHUFFLE_TABLE[8][8] = {
    {0, 1, 2, 3, 4, 5, 6, 7}, {1, 2, 3, 4, 5, 6, 7, 0}, {2, 3, 4, 5, 6, 7, 0, 1}, {3, 4, 5, 6, 7, 0, 1, 2},
    {4, 5, 6, 7, 0, 1, 2, 3}, {5, 6, 7, 0, 1, 2, 3, 4}, {6, 7, 0, 1, 2, 3, 4, 5}, {7, 0, 1, 2, 3, 4, 5, 6}};
#ifndef NO_AVX512
void _mm256_load_char_array_forward(__m256i &vec_index, const char *arr, __m256i &seq_vals, __mmask8 mask) {
    __mmask8 copy_mask = mask;
    uint32_t lowbit = copy_mask & (-copy_mask);
    int pos = _mm_popcnt_u32(lowbit - 1);
    __m256i shuffle_perm = _mm256_load_si256((__m256i *) SHUFFLE_TABLE[pos]);
    vec_index = seq_vals = _mm256_permutevar8x32_epi32(vec_index, shuffle_perm);
    int base_index = _mm256_extract_epi32(vec_index, 0);
    copy_mask >>= pos;
    __mmask16 mask_16 = copy_mask;
    __m128i bytes = _mm_maskz_loadu_epi8(mask_16, arr + base_index);
    seq_vals = _mm256_cvtepu8_epi32(bytes);
    __m256i perm = _mm256_load_si256((__m256i *) LEFT_SHIFT_TABLE[pos]);
    seq_vals = _mm256_permutevar8x32_epi32(seq_vals, perm);
    seq_vals = _mm256_maskz_mov_epi32(mask, seq_vals);
}

void _mm256_load_char_array_backward(__m256i &vec_index, const char *arr, __m256i &seq_vals, __mmask8 mask) {
    __mmask8 reverse_mask = REVERSE_MASK_LUT[mask];
    __mmask8 copy_mask = reverse_mask;
    uint32_t lowbit = copy_mask & (-copy_mask);
    int pos = _mm_popcnt_u32(lowbit - 1);
    __m256i shuffle_perm = _mm256_load_si256((__m256i *) SHUFFLE_TABLE[7 - pos]);
    vec_index = seq_vals = _mm256_permutevar8x32_epi32(vec_index, shuffle_perm);
    int base_index = _mm256_extract_epi32(vec_index, 0);
    copy_mask >>= pos;
    __mmask16 mask_16 = copy_mask;
    __m128i bytes = _mm_maskz_loadu_epi8(mask_16, arr + base_index);
    seq_vals = _mm256_cvtepu8_epi32(bytes);
    __m256i perm = _mm256_load_si256((__m256i *) RIGHT_SHIFT_TABLE[pos]);
    seq_vals = _mm256_permutevar8x32_epi32(seq_vals, perm);
    seq_vals = _mm256_maskz_mov_epi32(mask, seq_vals);
}
#endif
int BLOSUM62[] = {
    4,                                                                                     // A
    -1, 5,                                                                                 // R
    -2, 0,  6,                                                                             // N
    -2, -2, 1,  6,                                                                         // D
    0,  -3, -3, -3, 9,                                                                     // C
    -1, 1,  0,  0,  -3, 5,                                                                 // Q
    -1, 0,  0,  2,  -4, 2,  5,                                                             // E
    0,  -2, 0,  -1, -3, -2, -2, 6,                                                         // G
    -2, 0,  1,  -1, -3, 0,  0,  -2, 8,                                                     // H
    -1, -3, -3, -3, -1, -3, -3, -4, -3, 4,                                                 // I
    -1, -2, -3, -4, -1, -2, -3, -4, -3, 2,  4,                                             // L
    -1, 2,  0,  -1, -3, 1,  1,  -2, -1, -3, -2, 5,                                         // K
    -1, -1, -2, -3, -1, 0,  -2, -3, -2, 1,  2,  -1, 5,                                     // M
    -2, -3, -3, -3, -2, -3, -3, -3, -1, 0,  0,  -3, 0,  6,                                 // F
    -1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7,                             // P
    1,  -1, 1,  0,  -1, 0,  0,  0,  -1, -2, -2, 0,  -1, -2, -1, 4,                         // S
    0,  -1, 0,  -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1,  5,                     // T
    -3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1,  -4, -3, -2, 11,                // W
    -2, -2, -2, -3, -2, -1, -2, -3, 2,  -1, -1, -2, -1, 3,  -3, -2, -2, 2,  7,             // Y
    0,  -3, -3, -3, -1, -2, -2, -3, -3, 3,  1,  -2, 1,  -1, -2, -2, 0,  -3, -1, 4,         // V
    -2, -1, 3,  4,  -3, 0,  1,  -1, 0,  -3, -4, 0,  -3, -3, -2, 0,  -1, -4, -3, -3, 4,     // B
    -1, 0,  0,  1,  -3, 3,  4,  -2, 0,  -3, -3, 1,  -1, -3, -1, 0,  -1, -3, -2, -2, 1,  4, // Z
    0,  -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2, 0,  0,  -2, -1, -1, -1, -1,
    -1 // X
    // A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X
    // 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19  2  6 20
};

int na2idx[] = {0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4};
// idx for  A  B  C  D  E  F  G  H  I  J  K  L  M  N  O  P
//          Q  R  S  T  U  V  W  X  Y  Z
// so aa2idx[ X - 'A'] => idx_of_X, eg aa2idx['A' - 'A'] => 0,
// and aa2idx['M'-'A'] => 4
int BLOSUM62_na[] = {
    2,                     // A
    -2, 2,                 // C
    -2, -2, 2,             // G
    -2, -2, -2, 2,         // T
    -2, -2, -2, 1,  2,     // U
    -2, -2, -2, -2, -2, 1, // N
    0,  0,  0,  0,  0,  0,
    1 // X
    // A  C  G  T  U  N  X
    // 0  1  2  3  3  4  5
};

void setaa_to_na();

struct Slot {
    // 头信息：固定 7 个 long（你原代码用 malloc 了，其实可直接数组）
    long info[7];
    // 变长数据缓冲
    long *cluster_id = nullptr;
    size_t cluster_n = 0; // count (elements)
    long *seqs_suffix = nullptr;
    size_t suffix_n = 0;
    long long *prefix = nullptr;
    size_t prefix_n = 0; // 注意 long long
    long *indexCount = nullptr;
    size_t indexCount_n = 0;

    // 非阻塞请求（把一次“一个块”的所有 Ibcast 请求都挂在这）
    std::vector<MPI_Request> reqs;

    // 析构/释放
    void release() {
        free(cluster_id);
        cluster_id = nullptr;
        cluster_n = 0;
        free(seqs_suffix);
        seqs_suffix = nullptr;
        suffix_n = 0;
        free(prefix);
        prefix = nullptr;
        prefix_n = 0;
        free(indexCount);
        indexCount = nullptr;
        indexCount_n = 0;
        reqs.clear();
    }
};

struct TempFile {
    FILE *file;
    char buf[512];

    TempFile(const char *dir = NULL) {
        int len = dir ? strlen(dir) : 0;
        assert(len < 400);
        buf[0] = 0;
        if (len) {
            strcat(buf, dir);
            if (buf[len - 1] != '/' && buf[len - 1] != '\\') {
                buf[len] = '/';
                len += 1;
            }
        }
        strcat(buf, "cdhit.temp.");
        len += 11;
        sprintf(buf + len, "%p", this);
        file = fopen(buf, "w+");
    }
    ~TempFile() {
        if (file) {
            fclose(file);
            remove(buf);
        }
    }
};

struct TempFiles {
    NVector<TempFile *> files;

    ~TempFiles() { Clear(); }

    void Clear() {
        int i;
#pragma omp critical
        {
            //  cerr << "files.size     " << files.size << endl;
            for (i = 0; i < files.size; i++) {
                if (files[i]) {
                    cerr << "Deleting temporary file: " << files[i]->buf << endl; // 输出临时文件的路径
                    delete files[i];
                }
            }
        }
    }
};

const char *temp_dir = "";
TempFiles temp_files;

static std::string MakeAbsolutePathFromCwd(const std::string &path, const char *default_name) {
    std::string target = path.empty() ? std::string(default_name) : path;
    if (!target.empty() && target[0] == '/') return target;

    char cwd_buf[PATH_MAX];
    if (getcwd(cwd_buf, sizeof(cwd_buf)) == nullptr) {
        return target;
    }
    std::string cwd(cwd_buf);
    if (!cwd.empty() && cwd.back() != '/') cwd.push_back('/');
    return cwd + target;
}

FILE *OpenTempFile(const char *dir = NULL) {
    TempFile *file = new TempFile(dir);
#pragma omp critical
    {
        // cerr<<"caonima"<<endl;
        temp_files.files.Append(file);
    }
    return file->file;
}
static void CleanUpTempFiles() {
    if (temp_files.files.Size()) printf("Clean up temporary files ...\n");
    temp_files.Clear();
}

int NAA1;
int NAA2;
int NAA3;
int NAA4;
int NAA5;
int NAA6;
int NAA7;
int NAA8;
int NAA9;
int NAA10;
int NAA11;
int NAA12;
int NAAN_array[13] = {1};

void InitNAA(int max) {
    NAA1 = NAAN_array[1] = max;
    NAA2 = NAAN_array[2] = NAA1 * NAA1;
    NAA3 = NAAN_array[3] = NAA1 * NAA2;
    NAA4 = NAAN_array[4] = NAA2 * NAA2;
    NAA5 = NAAN_array[5] = NAA2 * NAA3;
    NAA6 = NAAN_array[6] = NAA3 * NAA3;
    NAA7 = NAAN_array[7] = NAA3 * NAA4;
    NAA8 = NAAN_array[8] = NAA4 * NAA4;
    NAA9 = NAAN_array[9] = NAA4 * NAA5;
    NAA10 = NAAN_array[10] = NAA5 * NAA5;
    NAA11 = NAAN_array[11] = NAA5 * NAA6;
    NAA12 = NAAN_array[12] = NAA6 * NAA6;
}

extern Options options;
ScoreMatrix mat;
Vector<int> Comp_AAN_idx;

void make_comp_iseq(int len, char *iseq_comp, char *iseq) {
    int i, c[6] = {3, 2, 1, 0, 4, 5};
    for (i = 0; i < len; i++) iseq_comp[i] = c[(int) iseq[len - i - 1]];
} // make_comp_iseq

bool Options::SetOptionCommon(const char *flag, const char *value) {
    int intval = atoi(value);
    if (strcmp(flag, "-i") == 0)
        input = value;
    else if (strcmp(flag, "-j") == 0)
        input_pe = value;
    else if (strcmp(flag, "-o") == 0)
        output = value;
    else if (strcmp(flag, "-op") == 0)
        output_pe = value;
    else if (strcmp(flag, "-pre_out") == 0)
        preprocess_dir = value;
    else if (strcmp(flag, "-stealing") == 0)
        stealing = true;
    else if (strcmp(flag, "-M") == 0)
        max_memory = atoll(value) * 1000000;
    else if (strcmp(flag, "-l") == 0)
        min_length = intval;
    else if (strcmp(flag, "-c") == 0)
        cluster_thd = atof(value), useIdentity = true;
    else if (strcmp(flag, "-D") == 0)
        distance_thd = atof(value), useDistance = true;
    else if (strcmp(flag, "-b") == 0)
        band_width = intval;
    else if (strcmp(flag, "-n") == 0)
        NAA = intval;
    else if (strcmp(flag, "-d") == 0)
        des_len = intval;
    else if (strcmp(flag, "-s") == 0)
        diff_cutoff = atof(value);
    else if (strcmp(flag, "-S") == 0)
        diff_cutoff_aa = intval;
    else if (strcmp(flag, "-B") == 0)
        store_disk = intval;
    else if (strcmp(flag, "-P") == 0)
        PE_mode = intval;
    else if (strcmp(flag, "-cx") == 0)
        trim_len = intval;
    else if (strcmp(flag, "-cy") == 0)
        trim_len_R2 = intval;
    else if (strcmp(flag, "-ap") == 0)
        align_pos = intval;
    else if (strcmp(flag, "-sc") == 0)
        sort_output = intval;
    else if (strcmp(flag, "-sf") == 0)
        sort_outputf = intval;
    else if (strcmp(flag, "-p") == 0)
        print = intval;
    else if (strcmp(flag, "-g") == 0)
        cluster_best = intval;
    else if (strcmp(flag, "-G") == 0)
        global_identity = intval;
    else if (strcmp(flag, "-aL") == 0)
        long_coverage = atof(value);
    else if (strcmp(flag, "-AL") == 0)
        long_control = intval;
    else if (strcmp(flag, "-aS") == 0)
        short_coverage = atof(value);
    else if (strcmp(flag, "-AS") == 0)
        short_control = intval;
    else if (strcmp(flag, "-A") == 0)
        min_control = intval;
    else if (strcmp(flag, "-uL") == 0)
        long_unmatch_per = atof(value);
    else if (strcmp(flag, "-uS") == 0)
        short_unmatch_per = atof(value);
    else if (strcmp(flag, "-U") == 0)
        unmatch_len = intval;
    else if (strcmp(flag, "-tmp") == 0)
        temp_dir = value;
    else if (strcmp(flag, "-N") == 0)
        NodeNum = atof(value);
    else if (strcmp(flag, "-NT") == 0)
        threads_per_node = atof(value);
    else if (strcmp(flag, "-ST") == 0)
        core_size = atof(value);
    else if (strcmp(flag, "-nT") == 0)
        numa_size = atof(value);
    else if (strcmp(flag, "-bak") == 0)
        backupFile = intval;
    else if (strcmp(flag, "-ready") == 0)
        ready = true;
    else if (strcmp(flag, "-T") == 0) {
#ifndef NO_OPENMP
        int cpu = omp_get_num_procs();
        threads = intval;
        if (threads > cpu) {
            threads = cpu;
            printf("Warning: total number of CPUs in the system is %i\n", cpu);
        } else if (threads < 0) {
            threads += cpu;
            if (threads < 0) threads = 0;
        }
        if (threads == 0) {
            threads = cpu;
            printf("total number of CPUs in the system is %i\n", cpu);
        }
        if (threads != intval) printf("Actual number of CPUs to be used: %i\n\n", threads);
#else
        printf("Option -T is ignored: multi-threading with OpenMP is NOT enabled!\n");
#endif
    } else
        return false;
    return true;
}
bool Options::SetOption(const char *flag, const char *value) {
    if (is454) {
        if (strcmp(flag, "-s") == 0)
            return false;
        else if (strcmp(flag, "-S") == 0)
            return false;
        else if (strcmp(flag, "-G") == 0)
            return false;
        else if (strcmp(flag, "-A") == 0)
            return false;
        else if (strcmp(flag, "-r") == 0)
            return false;
        else if (strcmp(flag, "-D") == 0) {
            max_indel = atoi(value);
            return true;
        }
    }
    if (SetOptionCommon(flag, value)) return true;
    if (strcmp(flag, "-t") == 0)
        tolerance = atoi(value);
    else if (strcmp(flag, "-F") == 0)
        frag_size = atoi(value);
    else if (has2D && SetOption2D(flag, value))
        return true;
    else if (isEST && SetOptionEST(flag, value))
        return true;
    else
        return false;
    return true;
}
bool Options::SetOption2D(const char *flag, const char *value) {
    if (SetOptionCommon(flag, value)) return true;
    if (strcmp(flag, "-i2") == 0)
        input2 = value;
    else if (strcmp(flag, "-j2") == 0)
        input2_pe = value;
    else if (strcmp(flag, "-s2") == 0)
        diff_cutoff2 = atof(value);
    else if (strcmp(flag, "-S2") == 0)
        diff_cutoff_aa2 = atoi(value);
    else
        return false;
    return true;
}
bool Options::SetOptionEST(const char *flag, const char *value) {
    NAA_top_limit = 12;
    if (SetOptionCommon(flag, value)) return true;
    if (strcmp(flag, "-r") == 0)
        option_r = atoi(value);
    else if (strcmp(flag, "-gap") == 0)
        mat.gap = MAX_SEQ * atoi(value);
    else if (strcmp(flag, "-gap-ext") == 0)
        mat.ext_gap = MAX_SEQ * atoi(value);
    else if (strcmp(flag, "-match") == 0)
        mat.set_match(atoi(value));
    else if (strcmp(flag, "-mismatch") == 0)
        mat.set_mismatch(atoi(value));
    else if (strcmp(flag, "-mask") == 0) {
        string letters = value;
        int i, n = letters.size();
        for (i = 0; i < n; i++) {
            char ch = toupper(letters[i]);
            if (ch < 'A' || ch > 'Z') continue;
            na2idx[ch - 'A'] = 5;
        }
        setaa_to_na();
    } else
        return false;
    return true;
}
bool Options::SetOptions(int argc, char *argv[], int rank, bool twod, bool est) {
    int i, n;
    char date[100];
    strcpy(date, __DATE__);
    n = strlen(date);
    for (i = 1; i < n; i++)
        if (date[i - 1] == ' ' && date[i] == ' ') date[i] = '0';

    if (rank == 0) {
        printf("================================================================\n");
        printf("Program: RabbitCD-HIT," WITH_OPEMMP_AND_MPI ", %s, " __TIME__ "\n", date);

        printf("Command:");
        n = 9;
        for (i = 0; i < argc; i++) {
            n += (int)strlen(argv[i]) + 1;
            if (n >= 64) {
                printf("\n         %s", argv[i]);
                n = (int)strlen(argv[i]) + 9;
            } else {
                printf(" %s", argv[i]);
            }
        }
        printf("\n\n");

        time_t tm = time(NULL);
        printf("Started: %s", ctime(&tm));
        printf("================================================================\n");
        printf("                            Output                              \n");
        printf("----------------------------------------------------------------\n");
        fflush(stdout);
    }

    has2D = twod;
    isEST = est;

    for (i = 1; i + 1 < argc; i += 2)
        if (SetOption(argv[i], argv[i + 1]) == 0) return false;
    if (i < argc) return false;

    atexit(CleanUpTempFiles);
    return true;
}

void Options::Validate() {
    if (useIdentity and useDistance) bomb_error("can not use both identity cutoff and distance cutoff");
    if (useDistance) {
        if ((distance_thd > 1.0) || (distance_thd < 0.0)) bomb_error("invalid distance threshold");
    } else if (isEST) {
        if ((cluster_thd > 1.0) || (cluster_thd < 0.8)) bomb_error("invalid clstr threshold, should >=0.8");
    } else {
        if ((cluster_thd > 1.0) || (cluster_thd < 0.4)) bomb_error("invalid clstr");
    }

    // if (input.size()  == 0) bomb_error("no input file");
    // if (output.size() == 0) bomb_error("no output file");
    if (PE_mode) {
        if (input_pe.size() == 0) bomb_error("no input file for R2 sequences in PE mode");
        if (output_pe.size() == 0) bomb_error("no output file for R2 sequences in PE mode");
    }
    if (isEST && (align_pos == 1)) option_r = 0;

    if (band_width < 1) bomb_error("invalid band width");
    if (NAA < 2 || NAA > NAA_top_limit) bomb_error("invalid word length");
    if (des_len < 0) bomb_error("too short description, not enough to identify sequences");
    if (not isEST && (tolerance < 0 || tolerance > 5)) bomb_error("invalid tolerance");
    if ((diff_cutoff < 0) || (diff_cutoff > 1)) bomb_error("invalid value for -s");
    if (diff_cutoff_aa < 0) bomb_error("invalid value for -S");
    if (has2D) {
        if ((diff_cutoff2 < 0) || (diff_cutoff2 > 1)) bomb_error("invalid value for -s2");
        if (diff_cutoff_aa2 < 0) bomb_error("invalid value for -S2");
        if (PE_mode) {
            if (input2_pe.size() == 0) bomb_error("no input file for R2 sequences for 2nd db in PE mode");
        }
    }
    if (global_identity == 0) print = 1;
    if (short_coverage < long_coverage) short_coverage = long_coverage;
    if (short_control > long_control) short_control = long_control;
    if ((global_identity == 0) && (short_coverage == 0.0) && (min_control == 0))
        bomb_error("You are using local identity, but no -aS -aL -A option");
    if (frag_size < 0) bomb_error("invalid fragment size");
    preprocess_dir = MakeAbsolutePathFromCwd(preprocess_dir, "preprocess_output");

#if 0
	if( useDistance ){
		/* when required_aan becomes zero */
		if( distance_thd * NAA >= 1 )
			bomb_warning( "word length is too long for the distance cutoff" );
	}else{
		/* when required_aan becomes zero */
		if( cluster_thd <= 1.0 - 1.0 / NAA )
			bomb_warning( "word length is too long for the identity cutoff" );
	}
#endif

    const char *message = "Your word length is %i, using %i may be faster!\n";
    if (not isEST && tolerance) {
        int i, clstr_idx = (int) (cluster_thd * 100) - naa_stat_start_percent;
        int tcutoff = naa_stat[tolerance - 1][clstr_idx][5 - NAA];

        if (tcutoff < 5)
            bomb_error(
                "Too low cluster threshold for the word length.\n"
                "Increase the threshold or the tolerance, or decrease the word length.");
        for (i = 5; i > NAA; i--) {
            if (naa_stat[tolerance - 1][clstr_idx][5 - i] > 10) {
                printf(message, NAA, i);
                break;
            }
        }
    } else if (isEST) {
        if (cluster_thd > 0.9 && NAA < 8)
            printf(message, NAA, 8);
        else if (cluster_thd > 0.87 && NAA < 5)
            printf(message, NAA, 5);
        else if (cluster_thd > 0.80 && NAA < 4)
            printf(message, NAA, 4);
        else if (cluster_thd > 0.75 && NAA < 3)
            printf(message, NAA, 3);
    } else {
        if (cluster_thd > 0.85 && NAA < 5)
            printf(message, NAA, 5);
        else if (cluster_thd > 0.80 && NAA < 4)
            printf(message, NAA, 4);
        else if (cluster_thd > 0.75 && NAA < 3)
            printf(message, NAA, 3);
    }

    if ((min_length + 1) < NAA) bomb_error("Too short -l, redefine it");
}
void Options::Print() {
    printf("isEST = %i\n", isEST);
    printf("has2D = %i\n", has2D);
    printf("NAA = %i\n", NAA);
    printf("NAA_top_limit = %i\n", NAA_top_limit);
    printf("min_length = %i\n", min_length);
    printf("cluster_best = %i\n", cluster_best);
    printf("global_identity = %i\n", global_identity);
    printf("cluster_thd = %g\n", cluster_thd);
    printf("diff_cutoff = %g\n", diff_cutoff);
    printf("diff_cutoff_aa = %i\n", diff_cutoff_aa);
    printf("tolerance = %i\n", tolerance);
    printf("long_coverage = %g\n", long_coverage);
    printf("long_control = %i\n", long_control);
    printf("short_coverage = %g\n", short_coverage);
    printf("short_control = %i\n", short_control);
    printf("frag_size = %i\n", frag_size);
    printf("option_r = %i\n", option_r);
    printf("print = %i\n", print);
}

void bomb_error(const char *message) {
    fprintf(stderr, "\nFatal Error:\n%s\nProgram halted !!\n\n", message);
    temp_files.Clear();
    exit(1);
} // END void bomb_error

void bomb_error(const char *message, const char *message2) {
    fprintf(stderr, "\nFatal Error:\n%s %s\nProgram halted !!\n\n", message, message2);
    temp_files.Clear();
    exit(1);
} // END void bomb_error

void bomb_warning(const char *message) {
    fprintf(stderr, "\nWarning:\n%s\nNot fatal, but may affect results !!\n\n", message);
} // END void bomb_warning

void bomb_warning(const char *message, const char *message2) {
    fprintf(stderr, "\nWarning:\n%s %s\nNot fatal, but may affect results !!\n\n", message, message2);
} // END void bomb_warning

void format_seq(char *seq) {
    int i, j;
    char c1;
    int len = strlen(seq);

    for (i = 0, j = 0; i < len; i++) {
        c1 = toupper(seq[i]);
        if (isalpha(c1)) seq[j++] = c1;
    }
    seq[j] = 0;
} // END void format_seq

void strrev(char *p) {
    char *q = p;
    while (q && *q) ++q;
    for (--q; p < q; ++p, --q) *p = *p ^ *q, *q = *p ^ *q, *p = *p ^ *q;
}

////For smiple len1 <= len2, len2 is for existing representative
////walk along all diag path of two sequences,
////find the diags with most aap
////return top n diags
////added on 2006 11 13
////band 0                      XXXXXXXXXXXXXXXXXX               seq2, rep seq
////                            XXXXXXXXXXXXXXX                  seq1
////band 1                      XXXXXXXXXXXXXXXXXX               seq2, rep seq
////                             XXXXXXXXXXXXXXX                 seq1
////extreme right (+)           XXXXXXXXXXXXXXXXXX               seq2, rep seq
////    band = len2-1                            XXXXXXXXXXXXXXX seq1
////band-1                      XXXXXXXXXXXXXXXXXX               seq2, rep seq
////                           XXXXXXXXXXXXXXX                   seq1
////extreme left (-)            XXXXXXXXXXXXXXXXXX               seq2, rep seq
////              XXXXXXXXXXXXXXX   band = -(len1-1)             seq1
////index of diag_score = band+len1-1;
int diag_test_aapn(int NAA1, char iseq2[], int len1, int len2, WorkingBuffer &buffer, int &best_sum, int band_width,
                   int &band_left, int &band_center, int &band_right, int required_aa1) {
    int i, i1, j, k;
    int *pp;
    int nall = len1 + len2 - 1; // 总“对角线”数
    Vector<int> &taap = buffer.taap;
    Vector<INTs> &aap_begin = buffer.aap_begin;
    Vector<INTs> &aap_list = buffer.aap_list;
    Vector<int> &diag_score = buffer.diag_score;   // 记命中次数
    Vector<int> &diag_score2 = buffer.diag_score2; // 记加权命中次数

    if (nall > MAX_DIAG) bomb_error("in diag_test_aapn, MAX_DIAG reached");
    for (pp = &diag_score[0], i = nall; i; i--, pp++) *pp = 0;
    for (pp = &diag_score2[0], i = nall; i; i--, pp++) *pp = 0;

    int c22, cpx;
    INTs *bip;
    int len11 = len1 - 1;
    int len22 = len2 - 1;
    i1 = len11;
    for (i = 0; i < len22; i++, i1++) {
        c22 = iseq2[i] * NAA1 + iseq2[i + 1];
        cpx = 1 + (iseq2[i] != iseq2[i + 1]);
        if ((j = taap[c22]) == 0) continue;
        int m = aap_begin[c22];
        for (int k = 0; k < j; k++) {
            diag_score[i1 - aap_list[m + k]]++;
            diag_score2[i1 - aap_list[m + k]] += cpx;
        }
    }

    // find the best band range
    //   int band_b = required_aa1;
    int band_b = required_aa1 - 1 >= 0 ? required_aa1 - 1 : 0; // on dec 21 2001
    int band_e = nall - band_b;

    int band_m = (band_b + band_width - 1 < band_e) ? band_b + band_width - 1 : band_e;
    int best_score = 0;
    int best_score2 = 0;
    int max_diag = 0;
    int max_diag2 = 0;
    int imax_diag = 0;
    for (i = band_b; i <= band_m; i++) {
        best_score += diag_score[i];
        best_score2 += diag_score2[i];
        if (diag_score2[i] > max_diag2) {
            max_diag2 = diag_score2[i];
            max_diag = diag_score[i];
            imax_diag = i;
        }
    }
    int from = band_b;
    int end = band_m;
    int score = best_score;
    int score2 = best_score2;
    for (k = from, j = band_m + 1; j < band_e; j++, k++) {
        score -= diag_score[k];
        score += diag_score[j];
        score2 -= diag_score2[k];
        score2 += diag_score2[j];
        if (score2 > best_score2) {
            from = k + 1;
            end = j;
            best_score = score;
            best_score2 = score2;
            if (diag_score2[j] > max_diag2) {
                max_diag2 = diag_score2[j];
                max_diag = diag_score[j];
                imax_diag = j;
            }
        }
    }
    int mlen = imax_diag;
    if (imax_diag > len1) mlen = nall - imax_diag;
    int emax = int((1.0 - options.cluster_thd) * mlen) + 1;
    for (j = from; j < imax_diag; j++) { // if aap pairs fail to open gap
        if ((imax_diag - j) > emax || diag_score[j] < 1) {
            best_score -= diag_score[j];
            from++;
        } else
            break;
    }
    for (j = end; j > imax_diag; j--) { // if aap pairs fail to open gap
        if ((j - imax_diag) > emax || diag_score[j] < 1) {
            best_score -= diag_score[j];
            end--;
        } else
            break;
    }

    //  delete [] diag_score;
    band_left = from - len1 + 1;
    band_right = end - len1 + 1;
    band_center = imax_diag - len1 + 1;
    best_sum = best_score;
    return OK_FUNC;
}
// END diag_test_aapn

int diag_test_aapn_est(int NAA1, char iseq2[], int len1, int len2, WorkingBuffer &buffer, int &best_sum, int band_width,
                       int &band_left, int &band_center, int &band_right, int required_aa1) {
    int i, i1, j, k;
    int *pp, *pp2;
    int nall = len1 + len2 - 1;
    int NAA2 = NAA1 * NAA1;
    int NAA3 = NAA2 * NAA1;
    Vector<int> &taap = buffer.taap;
    Vector<INTs> &aap_begin = buffer.aap_begin;
    Vector<INTs> &aap_list = buffer.aap_list;
    Vector<int> &diag_score = buffer.diag_score;
    Vector<int> &diag_score2 = buffer.diag_score2;

    if (nall > MAX_DIAG) bomb_error("in diag_test_aapn_est, MAX_DIAG reached");
    pp = &diag_score[0];
    pp2 = &diag_score2[0];
    for (i = nall; i; i--, pp++, pp2++) *pp = *pp2 = 0;

    INTs *bip;
    int c22, cpx;
    int len22 = len2 - 3;
    i1 = len1 - 1;
    for (i = 0; i < len22; i++, i1++, iseq2++) {
        unsigned char c0 = iseq2[0];
        unsigned char c1 = iseq2[1];
        unsigned char c2 = iseq2[2];
        unsigned char c3 = iseq2[3];
        if ((c0 >= 4) || (c1 >= 4) || (c2 >= 4) || (c3 >= 4)) continue; // skip N

        c22 = c0 * NAA3 + c1 * NAA2 + c2 * NAA1 + c3;
        if ((j = taap[c22]) == 0) continue;
        cpx = 1 + (c0 != c1) + (c1 != c2) + (c2 != c3);
        bip = &aap_list[aap_begin[c22]]; //    bi = aap_begin[c22];
        for (; j; j--, bip++) {
            diag_score[i1 - *bip]++;
            diag_score2[i1 - *bip] += cpx;
        }
    }
#if 0
	int mmax = 0;
	int immax = 0;
	for(i=0; i<=nall; i++){
		if( i%len2 ==0 or i == nall ) printf( "\n" );
		printf( "%3i ", diag_score[i] );
		if( diag_score[i] > mmax ){
			mmax = diag_score[i];
			immax = i;
		}
	}
#endif

    // find the best band range
    //   int band_b = required_aa1;
    int band_b = required_aa1 - 1 >= 0 ? required_aa1 - 1 : 0; // on dec 21 2001
    int band_e = nall - band_b;

    if (options.is454) {
        band_b = len1 - band_width;
        band_e = len1 + band_width;
        if (band_b < 0) band_b = 0;
        if (band_e > nall) band_e = nall;
    }

    int band_m = (band_b + band_width - 1 < band_e) ? band_b + band_width - 1 : band_e;
    int best_score = 0;
    int best_score2 = 0;
    int max_diag = 0;
    int max_diag2 = 0;
    int imax_diag = 0;
    for (i = band_b; i <= band_m; i++) {
        best_score += diag_score[i];
        best_score2 += diag_score2[i];
        if (diag_score2[i] > max_diag2) {
            max_diag2 = diag_score2[i];
            max_diag = diag_score[i];
            imax_diag = i;
        }
    }
    int from = band_b;
    int end = band_m;
    int score = best_score;
    int score2 = best_score2;

    for (k = from, j = band_m + 1; j < band_e; j++, k++) {
        score -= diag_score[k];
        score += diag_score[j];
        score2 -= diag_score2[k];
        score2 += diag_score2[j];
        if (score2 > best_score2) {
            from = k + 1;
            end = j;
            best_score = score;
            best_score2 = score2;
            if (diag_score2[j] > max_diag2) {
                max_diag2 = diag_score2[j];
                max_diag = diag_score[j];
                imax_diag = j;
            }
        }
    }
#if 0
	printf( "%i %i\n", required_aa1, from );
	printf( "max=%3i  imax=%3i; band:  %3i  %3i  %i\n", max_diag, imax_diag, band_b, band_e, band_m );
	printf( "best: %i\n", best_score );
	printf( "from: %i, end: %i,  best: %i\n", from, end, best_score );
#endif
    int mlen = imax_diag;
    if (imax_diag > len1) mlen = nall - imax_diag;
    int emax = int((1.0 - options.cluster_thd) * mlen) + 1;
    for (j = from; j < imax_diag; j++) { // if aap pairs fail to open gap
        if ((imax_diag - j) > emax || diag_score[j] < 1) {
            best_score -= diag_score[j];
            from++;
        } else
            break;
    }
    for (j = end; j > imax_diag; j--) { // if aap pairs fail to open gap
        if ((j - imax_diag) > emax || diag_score[j] < 1) {
            best_score -= diag_score[j];
            end--;
        } else
            break;
    }

    band_left = from - len1 + 1;
    band_right = end - len1 + 1;
    band_center = imax_diag - len1 + 1;
    best_sum = best_score;
    if (options.is454) {
        if (band_left > 0) best_sum = 0;
        if (band_right < 0) best_sum = 0;
    }
#if 0
	printf( "%3i:  best: %i,  %i  %i  %i\n", required_aa1, best_score, band_left, band_right, band_width );
	printf( "max=%3i  imax=%3i; band:  %3i  %3i  %i\n", mmax, immax, band_b, band_e, band_m );
#endif
    return OK_FUNC;
}
// END diag_test_aapn_est

/*
local alignment of two sequence within a diag band
for band 0 means direction (0,0) -> (1,1)
         1 means direction (0,1) -> (1,2)
        -1 means direction (1,0) -> (2,1)
added on 2006 11 13
band 0                      XXXXXXXXXXXXXXXXXX               seq2, rep seq
                            XXXXXXXXXXXXXXX                  seq1
band 1                      XXXXXXXXXXXXXXXXXX               seq2, rep seq
                             XXXXXXXXXXXXXXX                 seq1
extreme right (+)           XXXXXXXXXXXXXXXXXX               seq2, rep seq
    band = len2-1                            XXXXXXXXXXXXXXX seq1
band-1                      XXXXXXXXXXXXXXXXXX               seq2, rep seq
                           XXXXXXXXXXXXXXX                   seq1
extreme left (-)            XXXXXXXXXXXXXXXXXX               seq2, rep seq
              XXXXXXXXXXXXXXX   band = -(len1-1)             seq1
iseq len are integer sequence and its length,
mat is matrix, return ALN_PAIR class

       band:  -101   seq2 len2 = 17
                \\\1234567890123456
              0  \xxxxxxxxxxxxxxxxx
              1   xxxxxxxxxxxxxxxxx\ most right band = len2-1
              2   xxxxxxxxxxxxxxxxx
    seq1      3   xxxxxxxxxxxxxxxxx
    len1 = 11 4   xxxxxxxxxxxxxxxxx
              5   xxxxxxxxxxxxxxxxx
              6   xxxxxxxxxxxxxxxxx
              7   xxxxxxxxxxxxxxxxx
              8   xxxxxxxxxxxxxxxxx
              9   xxxxxxxxxxxxxxxxx
              0   xxxxxxxxxxxxxxxxx
                  \
                   most left band = -(len1-1)

*/
#ifndef NO_AVX512
int rotation_band_align_AVX512(char iseq1[], char iseq2[], int len1, int len2, ScoreMatrix &mat, int &best_score,
                               int &iden_no, int &alnln, float &dist, int *alninfo, int band_left, int band_center,
                               int band_right, WorkingBuffer &buffer) {
    int i, j, k, j1;
    int jj, kk;
    int x, y;
    int iden_no1;
    int64_t best_score1;
    iden_no = 0;

    if ((band_right >= len2) || (band_left <= -len1) || (band_left > band_right)) return FAILED_FUNC;

    int band_width = band_right - band_left + 1;
    int band_width1 = 17;

    MatrixInt64 &score_mat = buffer.score_mat;
    MatrixInt &back_mat = buffer.back_mat;

    int L = band_left, R = band_right;
    int kmin = (R < 0) ? -R : (L > 0) ? L : 0;
    int kmax = (R < len2 - len1) ? (R + 2 * len1) : (L > len2 - len1) ? (2 * len2 - L) : (len1 + len2);

    int lenY = kmax - kmin + 1;

    if (score_mat.size() <= lenY) {
        VectorInt row(band_width1, 0);
        VectorInt64 row2(band_width1, 0);
        while (score_mat.size() <= lenY) {
            score_mat.Append(row2);
            back_mat.Append(row);
        }
    }
    for (int i = 0; i <= lenY; i++) {
        if (score_mat[i].size < band_width1) score_mat[i].Resize(band_width1);
        if (back_mat[i].size < band_width1) back_mat[i].Resize(band_width1);
    }

    best_score = 0;

    // 初始化边界
    if (L < 0) {
        int T = (R < 0) ? R : 0;
        for (int X = L; X <= T; ++X) {
            int I = -X, J = 0;
            if (I < 0 || I > len1) continue;
            int n = (I + J) - kmin;
            int m = (X - L + 1) >> 1;
            score_mat[n][m] = (int64_t) mat.ext_gap * I;
            back_mat[n][m] = DP_BACK_TOP;
        }
        back_mat[(-T) - kmin][(T - L + 1) / 2] = DP_BACK_NONE;
    }

    if (R >= 0) {
        int T = (L > 0) ? L : 0;
        for (int X = T; X <= R; ++X) {
            int I = 0, J = X;
            if (J < 0 || J > len2) continue;
            int n = (I + J) - kmin;
            int m = (X - L + 1) >> 1;
            score_mat[n][m] = (int64_t) mat.ext_gap * J;
            back_mat[n][m] = DP_BACK_LEFT;
        }
        back_mat[T - kmin][(T - L + 1) / 2] = DP_BACK_NONE;
    }

    int gap_open[2] = {mat.gap, mat.ext_gap};
    int max_diag = band_center - band_left;
    int extra_score[4] = {4, 3, 2, 1};

    // ============ AVX-512 向量化变量定义 ============
    const int SIMD_WIDTH = 8;

    __m512i vec_0 = _mm512_setzero_si512();
    __m512i vec_1 = _mm512_set1_epi64(1);
    __m512i vec_3 = _mm512_set1_epi64(3);
    __m512i vec_4 = _mm512_set1_epi64(4);
    __m256i vec_DP_BACK_LEFT = _mm256_set1_epi32(DP_BACK_LEFT);
    __m256i vec_DP_BACK_TOP = _mm256_set1_epi32(DP_BACK_TOP);
    __m256i vec_DP_BACK_LEFT_TOP = _mm256_set1_epi32(DP_BACK_LEFT_TOP);
    __m512i vec_offset = _mm512_setr_epi64(0, 2, 4, 6, 8, 10, 12, 14);

    __m512i vec_len1 = _mm512_set1_epi64(len1);
    __m512i vec_len2 = _mm512_set1_epi64(len2);
    __m512i vec_lenY = _mm512_set1_epi64(lenY);
    __m512i vec_R = _mm512_set1_epi64(R);

    __m512i vec_band_center = _mm512_set1_epi64(band_center);

    __m512i vec_ext_gap = _mm512_set1_epi64(mat.ext_gap);
    __m512i vec_gap_open = _mm512_set1_epi64(mat.gap);

    size_t mat_size = MAX_AA_2;
    char *base = (char *) mat.flat_matrix;
    for (size_t offset = 0; offset < mat_size; offset += 64) {
        _mm_prefetch(base + offset, _MM_HINT_T0);
    }

    // ============ 主循环 ============
    for (int y = kmin + 1; y <= kmax; y++) {
        int offset = (abs(y + L)) % 2;
        int x_start = L + offset;

        __m512i vec_y = _mm512_set1_epi64(y);
        __m512i vec_y_minus_kmin = _mm512_set1_epi64(y - kmin);

        int num_elements = 16;
        int vec_iterations = num_elements / SIMD_WIDTH;

        int index_y = y - kmin;
        int index_y1 = y - 1 - kmin;
        int index_y2 = y - 2 - kmin;

        // size_t mat_size = MAX_AA_2;
        // char* base = (char*)mat.flat_matrix;
        // for (size_t offset = 0; offset < mat_size; offset += 64) {
        //     _mm_prefetch(base + offset, _MM_HINT_T0);
        // }

        if (y + 1 <= kmax) {
            _mm_prefetch((const char *) &score_mat[y + 1 - kmin][0], _MM_HINT_T0);
            _mm_prefetch((const char *) &back_mat[y + 1 - kmin][0], _MM_HINT_T0);
        }

        // ============ 向量化主循环 ============
        for (int vec_idx = 0; vec_idx < vec_iterations; vec_idx++) {
            int x_base = x_start + vec_idx * SIMD_WIDTH * 2;
            int index_x = (x_base - L + 1) >> 1;
            int index_x_L = (x_base - 1 - L + 1) >> 1;
            int index_x_R = (x_base + 1 - L + 1) >> 1;
            // cout<<"x_base:"<<x_base<<endl;

            __m512i vec_x = _mm512_add_epi64(_mm512_set1_epi64(x_base), vec_offset);

            __m512i vec_i = _mm512_sub_epi64(vec_y, vec_x);
            vec_i = _mm512_srai_epi64(vec_i, 1);
            __m512i vec_j = _mm512_add_epi64(vec_x, vec_y);
            vec_j = _mm512_srai_epi64(vec_j, 1);

            __mmask8 vec_i_gt_0 = _mm512_cmpgt_epi64_mask(vec_i, vec_0);
            __mmask8 vec_i_le_len1 = _mm512_cmple_epi64_mask(vec_i, vec_len1);
            __mmask8 vec_j_gt_0 = _mm512_cmpgt_epi64_mask(vec_j, vec_0);
            __mmask8 vec_j_le_len2 = _mm512_cmple_epi64_mask(vec_j, vec_len2);
            __mmask8 mask = vec_i_gt_0 & vec_i_le_len1 & vec_j_gt_0 & vec_j_le_len2;
            int count = _mm_popcnt_u32((unsigned int) mask);
            if (count == 0) continue;

#ifdef GATHER
            // printf(" (0x%02X)\n", (unsigned int)mask);
            __m512i vec_i_minus_1 = _mm512_sub_epi64(vec_i, vec_1);

            __m512i vec_j_minus_1 = _mm512_sub_epi64(vec_j, vec_1);

            vec_i_minus_1 = _mm512_mask_blend_epi64(mask, vec_0, vec_i_minus_1);
            vec_j_minus_1 = _mm512_mask_blend_epi64(mask, vec_0, vec_j_minus_1);

            __m256i vec_i_index = _mm512_cvtepi64_epi32(vec_i_minus_1);
            __m256i vec_j_index = _mm512_cvtepi64_epi32(vec_j_minus_1);

            __m256i seq1_i32 =
                _mm256_mmask_i32gather_epi32(_mm256_setzero_si256(), mask, vec_i_index, (const int *) iseq1, 1);

            __m256i seq2_i32 =
                _mm256_mmask_i32gather_epi32(_mm256_setzero_si256(), mask, vec_j_index, (const int *) iseq2, 1);
            __m256i seq1_vals = _mm256_and_si256(seq1_i32, _mm256_set1_epi32(0xFF));
            __m256i seq2_vals = _mm256_and_si256(seq2_i32, _mm256_set1_epi32(0xFF));
#else

            __m512i vec_i64_index = _mm512_sub_epi64(vec_i, vec_1);
            __m512i vec_j64_index = _mm512_sub_epi64(vec_j, vec_1);
            __m256i vec_i_index = _mm512_cvtepi64_epi32(vec_i64_index);
            __m256i vec_j_index = _mm512_cvtepi64_epi32(vec_j64_index);

            __m256i seq1_vals, seq2_vals;
            _mm256_load_char_array_backward(vec_i_index, iseq1, seq1_vals, mask);
            _mm256_load_char_array_forward(vec_j_index, iseq2, seq2_vals, mask);
#endif

            __m256i matrix_indices = _mm256_mullo_epi32(seq1_vals, _mm256_set1_epi32(MAX_AA));
            matrix_indices = _mm256_add_epi32(matrix_indices, seq2_vals);

            __m512i vec_sij =
                _mm512_mask_i32gather_epi64(_mm512_setzero_si512(), mask, matrix_indices, mat.flat_matrix, 8);

            // _mm512_store_epi64(i_arr, vec_i);
            // _mm512_store_epi64(j_arr, vec_j);
            // // double t1 =get_time();
            // for(int index=0; index<SIMD_WIDTH; index++){
            // 	int bit = (mask >> index) & 1;
            // 	if(!bit) sij_arr[index] = 0;
            // 	else sij_arr[index] = mat.matrix[iseq1[i_arr[index]-1]][iseq2[j_arr[index]-1]];
            // }
            // __m512i vec_sij = _mm512_load_si512(sij_arr);

            // double t2 = get_time();
            // t += t2 -t1;

            // int extra = extra_score[abs(x-band_center)&3];
            __m512i vec_extra = _mm512_sub_epi64(vec_x, vec_band_center);
            vec_extra = _mm512_abs_epi64(vec_extra);
            vec_extra = _mm512_and_si512(vec_extra, vec_3);
            vec_extra = _mm512_sub_epi64(vec_4, vec_extra);

            __mmask8 mask_sij = _mm512_cmpgt_epi64_mask(vec_sij, vec_0);
            mask_sij = mask_sij & mask;
            vec_extra = _mm512_mask_blend_epi64(mask_sij, vec_0, vec_extra);
            vec_sij = _mm512_add_epi64(vec_sij, vec_extra);

            __m512i vec_best_score1 = _mm512_loadu_si512(&score_mat[index_y][index_x]);
            __m256i vec_back = _mm256_loadu_si256((__m256i *) &back_mat[index_y][index_x]);

            // left-top
            if (y - 2 >= kmin) {
                __m512i vec_score_y2 = _mm512_loadu_si512(&score_mat[index_y2][index_x]);
                vec_best_score1 = _mm512_add_epi64(vec_score_y2, vec_sij);
                vec_back = _mm256_set1_epi32(DP_BACK_LEFT_TOP);
            }

            __mmask8 mask_gap = _mm512_cmpeq_epi64_mask(vec_i, vec_len1) | _mm512_cmpeq_epi64_mask(vec_j, vec_len2);
            __m512i vec_gap0 = _mm512_mask_blend_epi64(mask_gap, vec_gap_open, vec_ext_gap);
            __m512i vec_gap;
            __m512i vec_score;

            // left
            if (y - 1 >= kmin) {
                vec_gap = vec_gap0;
                __m512i vec_score_y1;
                __m256i vec_back_y1;

                vec_score_y1 = _mm512_loadu_si512(&score_mat[index_y1][index_x_L]);
                vec_back_y1 = _mm256_loadu_si256((__m256i *) &back_mat[index_y1][index_x_L]);

                __mmask8 gap_flag = _mm256_cmpeq_epi32_mask(vec_back_y1, vec_DP_BACK_LEFT);

                vec_gap = _mm512_mask_blend_epi64(gap_flag, vec_gap0, vec_ext_gap);

                vec_score = _mm512_add_epi64(vec_score_y1, vec_gap);

                __mmask8 modify_flag =
                    _mm512_cmpgt_epi64_mask(vec_score, vec_best_score1) & (0xFE | ((offset) | (vec_idx)));
                vec_best_score1 = _mm512_mask_blend_epi64(modify_flag, vec_best_score1, vec_score);
                vec_back = _mm256_mask_blend_epi32(modify_flag, vec_back, vec_DP_BACK_LEFT);
            }

            // top
            if (y - 1 >= kmin) {
                vec_gap = vec_gap0;
                __m512i vec_score_y1;
                __m256i vec_back_y1;
                vec_score_y1 = _mm512_loadu_si512(&score_mat[index_y1][index_x_R]);
                vec_back_y1 = _mm256_loadu_si256((__m256i *) &back_mat[index_y1][index_x_R]);

                __mmask8 gap_flag = _mm256_cmpeq_epi32_mask(vec_back_y1, vec_DP_BACK_TOP);
                vec_gap = _mm512_mask_blend_epi64(gap_flag, vec_gap0, vec_ext_gap);

                vec_score = _mm512_add_epi64(vec_score_y1, vec_gap);
                __mmask8 modify_flag = _mm512_cmpgt_epi64_mask(vec_score, vec_best_score1) &
                                       _mm512_cmple_epi64_mask(_mm512_add_epi64(vec_x, vec_1), vec_R);
                vec_best_score1 = _mm512_mask_blend_epi64(modify_flag, vec_best_score1, vec_score);
                vec_back = _mm256_mask_blend_epi32(modify_flag, vec_back, vec_DP_BACK_TOP);
            }
            _mm512_mask_storeu_epi64(&score_mat[index_y][index_x], mask, vec_best_score1);
            _mm256_mask_storeu_epi32((__m256i *) &back_mat[index_y][index_x], mask, vec_back);
        }
        // if (y >= kmin + 3) exit(0);
        // for(int x=L+(abs(y+L))%2;x<=R;x+=2){
        // 	// avx<<score_mat[y-kmin][(x-L+1)>>1]<<" ";
        // 	cout << score_mat[y-kmin][(x-L+1)>>1]<<" ";
        // }
        // // avx<<endl;
        // cout << endl;
        // if(y > kmin+20) exit(0);
    }
    // cout<<t<<endl;

    // ============ 回溯部分（保持原样）============
    x = (R < len2 - len1) ? R : (L > len2 - len1) ? L : len2 - len1;
    y = kmax;
    i = (-x + y) >> 1;
    j = (x + y) >> 1;
    // printf("\n(%d,%d)\n", i, j);
    best_score = score_mat[y - kmin][(x - L + 1) >> 1];
    best_score1 = score_mat[y - kmin][(x - L + 1) >> 1];

    int back = back_mat[y - kmin][(x - L + 1) >> 1];
    int last = back;

    int count = 0, count2 = 0, count3 = 0;
    int match, begin1, begin2, end1, end2;
    int gbegin1 = 0, gbegin2 = 0, gend1 = 0, gend2 = 0;
    int64_t score, smin = best_score1, smax = best_score1 - 1;
    int posmin, posmax, pos = 0;
    int bl, dlen = 0, dcount = 0;
    posmin = posmax = 0;
    begin1 = begin2 = end1 = end2 = 0;

    int masked = 0;
    int indels = 0;
    int max_indels = 0;

    while (back != DP_BACK_NONE) {
        switch (back) {
            case DP_BACK_TOP:
                bl = (last != back) & (j != 1) & (j != len2);
                dlen += bl;
                dcount += bl;
                score = score_mat[y - kmin][(x - L + 1) >> 1];
                if (score < smin) {
                    count2 = 0;
                    smin = score;
                    posmin = pos - 1;
                    begin1 = i;
                    begin2 = j;
                }
                i -= 1;
                x = x + 1;
                y = y - 1;
                break;
            case DP_BACK_LEFT:
                bl = (last != back) & (i != 1) & (i != len1);
                dlen += bl;
                dcount += bl;
                score = score_mat[y - kmin][(x - L + 1) >> 1];
                if (score < smin) {
                    count2 = 0;
                    smin = score;
                    posmin = pos - 1;
                    begin1 = i;
                    begin2 = j;
                }
                j -= 1;
                x = x - 1;
                y = y - 1;
                break;
            case DP_BACK_LEFT_TOP:
                if (alninfo && true) {
                    if (i == 1 || j == 1) {
                        gbegin1 = i - 1;
                        gbegin2 = j - 1;
                    } else if (i == len1 || j == len2) {
                        gend1 = i - 1;
                        gend2 = j - 1;
                    }
                }
                score = score_mat[y - kmin][(x - L + 1) >> 1];
                i -= 1;
                j -= 1;
                y -= 2;
                match = iseq1[i] == iseq2[j];
                if (score > smax) {
                    count = 0;
                    smax = score;
                    posmax = pos;
                    end1 = i;
                    end2 = j;
                }
                if (false && (iseq1[i] > 4 || iseq2[j] > 4)) {
                    masked += 1;
                } else {
                    dlen += 1;
                    dcount += !match;
                    count += match;
                    count2 += match;
                    count3 += match;
                }
                if (score < smin) {
                    int mm = match == 0;
                    count2 = 0;
                    smin = score;
                    posmin = pos - mm;
                    begin1 = i + mm;
                    begin2 = j + mm;
                }
                break;
            default:
                printf("%i\n", back);
                break;
        }
        pos += 1;
        last = back;
        back = back_mat[y - kmin][(x - L + 1) >> 1];
    }

    iden_no = true ? count3 : count - count2;
    alnln = posmin - posmax + 1 - masked;
    dist = dcount / (float) dlen;

    int umtail1 = len1 - 1 - end1;
    int umtail2 = len2 - 1 - end2;
    int umhead = begin1 < begin2 ? begin1 : begin2;
    int umtail = umtail1 < umtail2 ? umtail1 : umtail2;
    int umlen = umhead + umtail;

    if (umlen > 99999999) return FAILED_FUNC;
    if (umlen > len1 * 1.0) return FAILED_FUNC;
    if (umlen > len2 * 1.0) return FAILED_FUNC;

    if (alninfo) {
        alninfo[0] = begin1;
        alninfo[1] = end1;
        alninfo[2] = begin2;
        alninfo[3] = end2;
        alninfo[4] = masked;
        if (true) {
            alninfo[0] = gbegin1;
            alninfo[1] = gend1;
            alninfo[2] = gbegin2;
            alninfo[3] = gend2;
        }
    }
    return OK_FUNC;
}
#endif
int local_band_align(char iseq1[], char iseq2[], int len1, int len2, ScoreMatrix &mat, int &best_score, int &iden_no,
                     int &alnln, float &dist, int *alninfo, int band_left, int band_center, int band_right,
                     WorkingBuffer &buffer) {
    int i, j, k, j1;
    int jj, kk;
    int iden_no1;
    int64_t best_score1;
    iden_no = 0;

    if ((band_right >= len2) || (band_left <= -len1) || (band_left > band_right)) return FAILED_FUNC;

    // allocate mem for score_mat[len1][len2] etc
    int band_width = band_right - band_left + 1;
    int band_width1 = band_width + 1;

    // score_mat, back_mat [i][j]: i index of seqi (0 to len(seqi)-1), j index of band (0 to band_width-1)
    MatrixInt64 &score_mat = buffer.score_mat;
    MatrixInt &back_mat = buffer.back_mat;

    // printf( "%i  %i\n", band_right, band_left );

    if (score_mat.size() <= len1) {
        VectorInt row(band_width1, 0);
        VectorInt64 row2(band_width1, 0);
        while (score_mat.size() <= len1) {
            score_mat.Append(row2);
            back_mat.Append(row);
        }
    }
    for (i = 0; i <= len1; i++) {
        if (score_mat[i].Size() < band_width1) score_mat[i].Resize(band_width1);
        if (back_mat[i].Size() < band_width1) back_mat[i].Resize(band_width1);
    }

    best_score = 0;
    /* seq1 is query, seq2 is rep
                  seq2    len2 = 17       seq2    len2 = 17    seq2    len2 = 17
                  01234567890123456       01234567890123456    01234567890123456
       0          xxxxxxxxxxxxxxxxx \\\\\\XXXxxxxxxxxxxxxxx    xXXXXXXXxxxxxxxxx
       1     \\\\\Xxxxxxxxxxxxxxxxx  \\\\\Xxx\xxxxxxxxxxxxx    xx\xxxxx\xxxxxxxx
       2      \\\\X\xxxxxxxxxxxxxxx   \\\\Xxxx\xxxxxxxxxxxx    xxx\xxxxx\xxxxxxx
  seq1 3       \\\Xx\xxxxxxxxxxxxxx    \\\Xxxxx\xxxxxxxxxxx    xxxx\xxxxx\xxxxxx
  len1 4        \\Xxx\xxxxxxxxxxxxx     \\Xxxxxx\xxxxxxxxxx    xxxxx\xxxxx\xxxxx
  = 11 5         \Xxxx\xxxxxxxxxxxx      \Xxxxxxx\xxxxxxxxx    xxxxxx\xxxxx\xxxx
       6          Xxxxx\xxxxxxxxxxx       Xxxxxxxx\xxxxxxxx    xxxxxxx\xxxxx\xxx
       7          x\xxxx\xxxxxxxxxx       x\xxxxxxx\xxxxxxx    xxxxxxxx\xxxxx\xx
       8          xx\xxxx\xxxxxxxxx       xx\xxxxxxx\xxxxxx    xxxxxxxxx\xxxxx\x
       9          xxx\xxxx\xxxxxxxx       xxx\xxxxxxx\xxxxx    xxxxxxxxxx\xxxxx\
       0          xxxx\xxxx\xxxxxxx       xxxx\xxxxxxx\xxxx    xxxxxxxxxxx\xxxxx
                  band_left < 0 (-6)      band_left < 0 (-6)   band_left >=0 (1)
                  band_right < 0 (-1)     band_right >=0 (2)   band_right >=0(7)
                  band_width 6            band_width 9         band_width 7
       init score_mat, and iden_mat (place with upper 'X')
     */

    if (band_left < 0) { // set score to left border of the matrix within band
        int tband = (band_right < 0) ? band_right : 0;
        // for (k=band_left; k<tband; k++)
        for (k = band_left; k <= tband; k++) { // fixed on 2006 11 14
            i = -k;
            j1 = k - band_left;
            // penalty for leading gap opening = penalty for gap extension
            // each of the left side query hunging residues give ext_gap (-1)
            score_mat[i][j1] = (int64_t) mat.ext_gap * (int64_t) i;
            back_mat[i][j1] = DP_BACK_TOP;
        }
        back_mat[-tband][tband - band_left] = DP_BACK_NONE;
    }

    if (band_right >= 0) { // set score to top border of the matrix within band
        int tband = (band_left > 0) ? band_left : 0;
        for (j = tband; j <= band_right; j++) {
            j1 = j - band_left;
            score_mat[0][j1] = (int64_t) mat.ext_gap * (int64_t) j;
            back_mat[0][j1] = DP_BACK_LEFT;
        }
        back_mat[0][tband - band_left] = DP_BACK_NONE;
    }

    int gap_open[2] = {mat.gap, mat.ext_gap};
    int max_diag = band_center - band_left;
    int extra_score[4] = {4, 3, 2, 1};
    for (i = 1; i <= len1; i++) {
        int J0 = 1 - band_left - i;
        int J1 = len2 - band_left - i;
        if (J0 < 0) J0 = 0;
        if (J1 >= band_width) J1 = band_width;
        for (j1 = J0; j1 <= J1; j1++) {
            j = j1 + i + band_left;

            int ci = iseq1[i - 1];
            int cj = iseq2[j - 1];
            int sij = mat.matrix[ci][cj];
            // int iden_ij = (ci == cj);
            int s1, k0, back;

            /* extra score according to the distance to the best diagonal */
            int extra = extra_score[abs(j1 - max_diag) & 3]; // max distance 3
            sij += extra * (sij > 0);
            back = DP_BACK_LEFT_TOP;
            best_score1 = score_mat[i - 1][j1] + sij;
            int gap0 = gap_open[(i == len1) | (j == len2)];
            int gap = 0;
            int64_t score;

            if (j1 > 0) {
                gap = gap0;
                if (back_mat[i][j1 - 1] == DP_BACK_LEFT) gap = mat.ext_gap;
                if ((score = score_mat[i][j1 - 1] + gap) > best_score1) {
                    back = DP_BACK_LEFT;
                    best_score1 = score;
                }
            }
            if (j1 + 1 < band_width) {
                gap = gap0;
                if (back_mat[i - 1][j1 + 1] == DP_BACK_TOP) gap = mat.ext_gap;
                if ((score = score_mat[i - 1][j1 + 1] + gap) > best_score1) {
                    back = DP_BACK_TOP;
                    best_score1 = score;
                }
            }
            score_mat[i][j1] = best_score1;
            back_mat[i][j1] = back;
            // printf( "%2i(%2i) ", best_score1, iden_no1 );
        }
        // printf( "\n" );
    }

    i = j = 0;
    if (len2 - band_left < len1) {
        i = len2 - band_left;
        j = len2;
    } else if (len1 + band_right < len2) {
        i = len1;
        j = len1 + band_right;
    } else {
        i = len1;
        j = len2;
    }
    j1 = j - i - band_left;
    best_score = score_mat[i][j1];
    best_score1 = score_mat[i][j1];
#if 1
    const char *letters = "acgtnx";
    const char *letters2 = "ACGTNX";
#else
    const char *letters = "arndcqeghilkmfpstwyvbzx";
    const char *letters2 = "ARNDCQEGHILKMFPSTWYVBZX";
#endif
    int back = back_mat[i][j1];
    int last = back;
    int count = 0, count2 = 0, count3 = 0;
    int match, begin1, begin2, end1, end2;
    int gbegin1 = 0, gbegin2 = 0, gend1 = 0, gend2 = 0;
    int64_t score, smin = best_score1, smax = best_score1 - 1;
    int posmin, posmax, pos = 0;
    int bl, dlen = 0, dcount = 0;
    posmin = posmax = 0;
    begin1 = begin2 = end1 = end2 = 0;

#ifdef PRINT
#define PRINT
    printf("%i %i\n", best_score, score_mat[i][j1]);
    printf("%i %i %i\n", band_left, band_center, band_right);
    printf("%i %i %i %i\n", i, j, j1, len2);
#endif
#ifdef MAKEALIGN
#define MAKEALIGN
    char AA[MAX_SEQ], BB[MAX_SEQ];
    int NN = 0;
    int IA, IB;
    for (IA = len1; IA > i; IA--) {
        AA[NN] = letters[iseq1[IA - 1]];
        BB[NN++] = '-';
    }
    for (IB = len2; IB > j; IB--) {
        AA[NN] = '-';
        BB[NN++] = letters[iseq2[IB - 1]];
    }
#endif

    int masked = 0;
    int indels = 0;
    int max_indels = 0;
    while (back != DP_BACK_NONE) {
        switch (back) {
            case DP_BACK_TOP:
#ifdef PRINT
                printf("%5i: %c %c %9i\n", pos, letters[iseq1[i - 1]], '|', score_mat[i][j1]);
#endif
#ifdef MAKEALIGN
                AA[NN] = letters[iseq1[i - 1]];
                BB[NN++] = '-';
#endif
                bl = (last != back) & (j != 1) & (j != len2);
                dlen += bl;
                dcount += bl;
                score = score_mat[i][j1];
                if (score < smin) {
                    count2 = 0;
                    smin = score;
                    posmin = pos - 1;
                    begin1 = i;
                    begin2 = j;
                }
                i -= 1;
                j1 += 1;
                break;
            case DP_BACK_LEFT:
#ifdef PRINT
                printf("%5i: %c %c %9i\n", pos, '|', letters[iseq2[j - 1]], score_mat[i][j1]);
#endif
#ifdef MAKEALIGN
                AA[NN] = '-';
                BB[NN++] = letters[iseq2[j - 1]];
#endif
                bl = (last != back) & (i != 1) & (i != len1);
                dlen += bl;
                dcount += bl;
                score = score_mat[i][j1];
                if (score < smin) {
                    count2 = 0;
                    smin = score;
                    posmin = pos - 1;
                    begin1 = i;
                    begin2 = j;
                }
                j1 -= 1;
                j -= 1;
                break;
            case DP_BACK_LEFT_TOP:
#ifdef PRINT
                if (iseq1[i - 1] == iseq2[j - 1]) {
                    printf("%5i: %c %c %9i\n", pos, letters2[iseq1[i - 1]], letters2[iseq2[j - 1]], score_mat[i][j1]);
                } else {
                    printf("%5i: %c %c %9i\n", pos, letters[iseq1[i - 1]], letters[iseq2[j - 1]], score_mat[i][j1]);
                }
#endif
#ifdef MAKEALIGN
                if (iseq1[i - 1] == iseq2[j - 1]) {
                    AA[NN] = letters2[iseq1[i - 1]];
                    BB[NN++] = letters2[iseq2[j - 1]];
                } else {
                    AA[NN] = letters[iseq1[i - 1]];
                    BB[NN++] = letters[iseq2[j - 1]];
                }
#endif
                if (alninfo && options.global_identity) {
                    if (i == 1 || j == 1) {
                        gbegin1 = i - 1;
                        gbegin2 = j - 1;
                    } else if (i == len1 || j == len2) {
                        gend1 = i - 1;
                        gend2 = j - 1;
                    }
                }
                score = score_mat[i][j1];
                i -= 1;
                j -= 1;
                match = iseq1[i] == iseq2[j];
                if (score > smax) {
                    count = 0;
                    smax = score;
                    posmax = pos;
                    end1 = i;
                    end2 = j;
                }
                if (options.isEST && (iseq1[i] > 4 || iseq2[j] > 4)) {
                    masked += 1;
                } else {
                    dlen += 1;
                    dcount += !match;
                    count += match;
                    count2 += match;
                    count3 += match;
                }
                if (score < smin) {
                    int mm = match == 0;
                    count2 = 0;
                    smin = score;
                    posmin = pos - mm;
                    begin1 = i + mm;
                    begin2 = j + mm;
                }
                break;
            default:
                printf("%i\n", back);
                break;
        }
        if (options.is454) {
            if (back == DP_BACK_LEFT_TOP) {
                if (indels > max_indels) max_indels = indels;
                indels = 0;
            } else {
                if (last == DP_BACK_LEFT_TOP) {
                    indels = 1;
                } else if (indels) {
                    indels += 1;
                }
            }
        }
        pos += 1;
        last = back;
        back = back_mat[i][j1];
    }
    if (options.is454 and max_indels > options.max_indel) return FAILED_FUNC;
    iden_no = options.global_identity ? count3 : count - count2;
    alnln = posmin - posmax + 1 - masked;
    dist = dcount / (float) dlen;

    // dist = - 0.75 * log( 1.0 - dist * 4.0 / 3.0 );
    int umtail1 = len1 - 1 - end1;
    int umtail2 = len2 - 1 - end2;
    int umhead = begin1 < begin2 ? begin1 : begin2;
    int umtail = umtail1 < umtail2 ? umtail1 : umtail2;
    int umlen = umhead + umtail;
    if (umlen > options.unmatch_len) return FAILED_FUNC;
    if (umlen > len1 * options.short_unmatch_per) return FAILED_FUNC;
    if (umlen > len2 * options.long_unmatch_per) return FAILED_FUNC;
    if (alninfo) {
        alninfo[0] = begin1;
        alninfo[1] = end1;
        alninfo[2] = begin2;
        alninfo[3] = end2;
        alninfo[4] = masked;
        if (options.global_identity) {
            alninfo[0] = gbegin1;
            alninfo[1] = gend1;
            alninfo[2] = gbegin2;
            alninfo[3] = gend2;
        }
    }
#ifdef PRINT
    printf("%6i %6i:  %4i %4i %4i %4i\n", alnln, iden_no, begin1, end1, begin2, end2);
    printf("%6i %6i:  %4i %4i\n", posmin, posmax, posmin - posmax, count - count2);
    printf("smin = %9i, smax = %9i\n", smin, smax);
    printf("dlen = %5i, dcount = %5i, dist = %.3f\n", dlen, dcount, dcount / (float) dlen);
#endif
#ifdef MAKEALIGN
    float identity = iden_no / (float) (options.global_identity ? (len1 - masked) : alnln);
    if (identity < options.cluster_thd) return OK_FUNC;
    while (i--) {
        AA[NN] = letters[iseq1[i - 1]];
        BB[NN++] = '-';
    }
    while (j--) {
        AA[NN] = '-';
        BB[NN++] = letters[iseq2[j - 1]];
    }
    AA[NN] = '\0';
    BB[NN] = '\0';
    for (i = 0; i < NN / 2; i++) {
        char aa = AA[i], bb = BB[i];
        AA[i] = AA[NN - i - 1];
        BB[i] = BB[NN - i - 1];
        AA[NN - i - 1] = aa;
        BB[NN - i - 1] = bb;
    }
    static int fcount = 0;
    fcount += 1;
    FILE *fout = fopen("alignments.txt", "a");
    if (fout == NULL) {
        if (fcount <= 1) printf("alignment files open failed\n");
        return OK_FUNC;
    }
    fprintf(fout, "\n\n######################################################\n");
    fprintf(fout, "# length X = %i\n", len2);
    fprintf(fout, "# length Y = %i\n", len1);
    fprintf(fout, "# best align X: %i-%i\n", begin2 + 1, end2 + 1);
    fprintf(fout, "# best align Y: %i-%i\n", begin1 + 1, end1 + 1);
    if (alninfo) {
        fprintf(fout, "# align X: %i-%i\n", alninfo[2] + 1, alninfo[3] + 1);
        fprintf(fout, "# align Y: %i-%i\n", alninfo[0] + 1, alninfo[1] + 1);
    }
    fprintf(fout, "# alignment length: %i\n", alnln);
    fprintf(fout, "# identity count: %i\n", iden_no);
    fprintf(fout, "# identity: %g\n", identity);
    fprintf(fout, "# distance: %g\n", dist);
    if (options.is454) fprintf(fout, "# max indel: %i\n", max_indels);
#if 0
	fprintf( fout, "%i %s\n", seq1->index, AA );
	fprintf( fout, "%i %s\n", seq2->index, BB );
#else
    bool printaa = true;
    IB = IA = 0;
    fprintf(fout, "\n\nX ");
    while (IA < NN) {
        if (printaa) {
            fprintf(fout, "%c", BB[IB]);
            IB += 1;
            if (IB % 75 == 0 or IB == NN) printaa = false, fprintf(fout, "\nY ");
        } else {
            fprintf(fout, "%c", AA[IA]);
            IA += 1;
            if (IA % 75 == 0) printaa = true, fprintf(fout, "\n\nX ");
        }
    }
#endif
    fclose(fout);
#endif

    return OK_FUNC;
} // END int local_band_align

void setaa_to_na() {
    int i;
    for (i = 0; i < 26; i++) aa2idx[i] = na2idx[i];
} // END void setaa_to_na

/////////////////
ScoreMatrix::ScoreMatrix() {
    size_t req_size = MAX_AA * MAX_AA * sizeof(int64_t);
    size_t real_size = (req_size + 63) / 64 * 64;

    flat_matrix = (int64_t *) aligned_alloc(64, real_size);
    if (!flat_matrix) {
        fprintf(stderr, "Failed to allocate flat_matrix\n");
        exit(1);
    }
    init();
}

void ScoreMatrix::init() {
    set_gap(-11, -1);
    set_matrix(BLOSUM62);
}

void ScoreMatrix::set_gap(int gap1, int ext_gap1) {
    int i;
    gap = MAX_SEQ * gap1;
    ext_gap = MAX_SEQ * ext_gap1;
}

void ScoreMatrix::set_matrix(int *mat1) {
    int i, j, k;
    k = 0;
    for (i = 0; i < MAX_AA; i++)
        for (j = 0; j <= i; j++) matrix[j][i] = matrix[i][j] = MAX_SEQ * mat1[k++];
    update_flat_matrix();
}

void ScoreMatrix::update_flat_matrix() {
    for (int i = 0; i < MAX_AA; i++) {
        for (int j = 0; j < MAX_AA; j++) {
            flat_matrix[i * MAX_AA + j] = (int64_t) matrix[i][j];
        }
    }
}

void ScoreMatrix::set_to_na() {
    set_gap(-6, -1);
    set_matrix(BLOSUM62_na);
}
// Only for est
void ScoreMatrix::set_match(int score) {
    int i;
    for (i = 0; i < 5; i++) matrix[i][i] = MAX_SEQ * score;
    // matrix[3][4] = matrix[4][3] = MAX_SEQ * score;
}
// Only for est
void ScoreMatrix::set_mismatch(int score) {
    int i, j;
    for (i = 0; i < MAX_AA; i++)
        for (j = 0; j < i; j++) matrix[j][i] = matrix[i][j] = MAX_SEQ * score;
    matrix[3][4] = matrix[4][3] = MAX_SEQ;
}

WordTable::WordTable(int naa, int naan) {
    NAA = 0;
    NAAN = 0;
    is_aa = 1;
    size = 0;
    frag_count = 0;
    Init(naa, naan);
}

void WordTable::SetDNA() { is_aa = 0; }

void WordTable::Init(int naa, int naan) {
    NAA = naa;
    NAAN = naan;
    indexCounts.resize(NAAN);
}

void WordTable::Clear() {
    int i;
#if 0
	int n1 = 0, n2 = 0, n3 = 0, ns = 0;
	for(i=0; i<NAAN; i++){
		NVector<IndexCount> & ics = indexCounts[i];
		for(int j=0; j<ics.size; j++){
			IndexCount ic = ics[j];
			n1 += ic.count == 1;
			n2 += ic.count == 2;
			n3 += ic.count == 3;
			ns += ic.count >= 4;
		}
	}
	printf( "%9i %9i %9i %9i\n", n1, n2, n3, ns );
#endif
    size = 0;
    frag_count = 0;
    sequences.clear();
    for (i = 0; i < NAAN; i++) indexCounts[i].size = 0; // Clear();
}

int WordTable::AddWordCounts(NVector<IndexCount> &counts, Sequence *seq, bool skipN) {
    int aan_no = counts.Size();
    int i, j, k, idx = sequences.size();
    for (i = 0; i < aan_no; i++) {
        IndexCount ic = counts[i];
        if ((k = ic.count)) {
            j = ic.index;
            if (skipN && j < 0) continue; // for those has 'N'
            NVector<IndexCount> &row = indexCounts[j];
            ic.index = idx;
            row.Append(ic);
            size += 1;
        }
    }
    sequences.Append(seq);
    return OK_FUNC;
}
int WordTable::AddWordCountsFrag(NVector<IndexCount> &counts, int frag, int frag_size, int repfrag) { return 0; }
// 建立索引表
int WordTable::AddWordCounts(int aan_no, Vector<int> &word_encodes, Vector<INTs> &word_encodes_no, int idx,
                             bool skipN) {
    int i, j, k;
    // printf( "seq %6i: ", idx );
    for (i = 0; i < aan_no; i++) {
        if ((k = word_encodes_no[i])) {
            // assert(k<1000);
            j = word_encodes[i];
            if (skipN && j < 0) continue; // for those has 'N'
            NVector<IndexCount> &row = indexCounts[j];
            row.Append(IndexCount(idx, k));
            size += 1;
            // if( k >1 ) printf( " %3i", k );
        }
    }
    // printf( "\n" );
    return OK_FUNC;
}

int WordTable::AddWordCountsFrag(int aan_no, Vector<int> &word_encodes, Vector<INTs> &word_encodes_no, int frag,
                                 int frag_size) {
    int i, j, k, i1, k1, fra;

    for (i = 0; i < frag; i++) {
        k = (i + 1) * frag_size < aan_no ? (i + 1) * frag_size - 1 : aan_no - 1;
        // quick_sort(&word_encodes[0], i*frag_size, k);
        std::sort(word_encodes.begin() + i * frag_size, word_encodes.begin() + k + 1);
    }
    for (j = aan_no - 1; j; j--) {
        if (word_encodes[j] == word_encodes[j - 1]) {
            word_encodes_no[j - 1] += word_encodes_no[j];
            word_encodes_no[j] = 0;
        }
    }
    // END check_word_encodes

    for (i = 0; i < aan_no; i += frag_size) {
        k = frag_size < (aan_no - i) ? frag_size : (aan_no - i);
        fra = i / frag_size;
        // AddWordCounts(k, word_encodes+i, word_encodes_no+i, NR90f_no+fra);
        for (i1 = i; i1 < i + k; i1++) {
            if ((k1 = word_encodes_no[i1])) {
                j = word_encodes[i1];
                NVector<IndexCount> &row = indexCounts[j];
                row.Append(IndexCount(frag_count + fra, k1));
                size += 1;
            }
        }
    }
    frag_count += frag;

    return 0;
}

void WordTable::PrintAll() {
    int i, j, k;
    int cols = 0;
    long long total_words = 0;
    k = 0;
    for (i = 0; i < NAAN; i++) {
        int size = indexCounts[i].Size();
        if (size == 0) continue;
        cols++;
        cout << k << "\t" << i << "\tsize:" << size << "\t";
        for (j = 0; j < size; j++) {
            cout << indexCounts[i][j].index << "," << indexCounts[i][j].count << " ";
            total_words += indexCounts[i][j].count;
        }
        cout << endl;
        k++;
    }

    cout << "total cols: " << cols << " total words: " << total_words << endl;
}

/* Quick Sort.
 * Adam Drozdek: Data Structures and Algorithms in C++, 2nd Edition.
 */
void PartialQuickSort(IndexCount *data, int first, int last, int partial) {
    int lower = first + 1, upper = last;
    IndexCount pivot;
    IndexCount val;
    if (first >= last) return;
    val = data[first];
    data[first] = data[(first + last) / 2];
    data[(first + last) / 2] = val;
    pivot = data[first];

    while (lower <= upper) {
        while (lower <= last && data[lower].count < pivot.count) lower++;
        while (pivot.count < data[upper].count) upper--;
        if (lower < upper) {
            val = data[lower];
            data[lower] = data[upper];
            data[upper] = val;
            upper--;
        }
        lower++;
    }
    val = data[first];
    data[first] = data[upper];
    data[upper] = val;
    if (first < upper - 1) PartialQuickSort(data, first, upper - 1, partial);
    if (upper >= partial) return;
    if (upper + 1 < last) PartialQuickSort(data, upper + 1, last, partial);
}
int WordTable::CountWords(int aan_no, Vector<int> &word_encodes, Vector<INTs> &word_encodes_no,
                          NVector<IndexCount> &lookCounts, NVector<uint32_t> &indexMapping, bool est, int min) {
    int S = frag_count ? frag_count : sequences.size();
    int j, k, j0, j1, k1, m;
    int ix1, ix2, ix3, ix4;
    IndexCount tmp;

    IndexCount *ic = lookCounts.items;
    for (j = 0; j < lookCounts.size; j++, ic++) indexMapping[ic->index] = 0;
    lookCounts.size = 0;

    int *we = &word_encodes[0];
    j0 = 0;
    if (est)
        while (*we < 0) j0++, we++; // if met short word has 'N'
    INTs *wen = &word_encodes_no[j0];
    // printf( "\nquery : " );
    for (; j0 < aan_no; j0++, we++, wen++) {
        j = *we;
        j1 = *wen;
        // if( j1 >1 ) printf( " %3i", j1 );
        if (j < 0 || j >= indexCounts.size())continue;
        if (j1 == 0) continue;
        NVector<IndexCount> &one = indexCounts[j];
        k1 = one.Size();
        IndexCount *ic = one.items;
        // cerr<<"lookcounts size  "<<lookCounts.size<<endl;
        int rest = aan_no - j0 + 1;
        for (k = 0; k < k1; k++, ic++) {
            int c = ic->count < j1 ? ic->count : j1;
            uint32_t *idm = indexMapping.items + ic->index;
            // if(my_rank==1){
            // 	cerr<<"lookCounts.size     "<<lookCounts.size<<endl;
            // 	cerr<<"ic->index     "<<ic->index<<endl;
            // 	cerr<<"*idm     "<<*idm<<endl;
            // }

            if (*idm == 0) {
                if (rest < min) continue;
                IndexCount *ic2 = lookCounts.items + lookCounts.size;
                lookCounts.size += 1;
                *idm = lookCounts.size;
                ic2->index = ic->index;
                ic2->count = c;
            } else {
                lookCounts[*idm - 1].count += c;
            }
        }
    }
    // printf( "%6i %6i\n", S, lookCounts.size );
    lookCounts[lookCounts.size].count = 0;
    // printf( "\n\n" );
    return OK_FUNC;
}
int WordTable::CountWords(int aan_no, int qid, Vector<int> &word_encodes, Vector<INTs> &word_encodes_no,
                          NVector<IndexCount> &lookCounts, NVector<uint32_t> &indexMapping, bool est, int min) {
    int S = frag_count ? frag_count : sequences.size();
    int j, k, j0, j1, k1, m;
    int ix1, ix2, ix3, ix4;
    IndexCount tmp;

    IndexCount *ic = lookCounts.items;
    for (j = 0; j < lookCounts.size; j++, ic++) indexMapping[ic->index] = 0;
    lookCounts.size = 0;

    int *we = &word_encodes[0];
    j0 = 0;
    if (est)
        while (*we < 0) j0++, we++; // if met short word has 'N'
    INTs *wen = &word_encodes_no[j0];
    // printf( "\nquery : " );
    for (; j0 < aan_no; j0++, we++, wen++) {
        j = *we;
        j1 = *wen;
        // if( j1 >1 ) printf( " %3i", j1 );
        if (j < 0 || j >= indexCounts.size())continue;
        if (j1 == 0) continue;
        NVector<IndexCount> &one = indexCounts[j];
        k1 = one.Size();
        IndexCount *ic = one.items;
        // cerr<<"lookcounts size  "<<lookCounts.size<<endl;
        int rest = aan_no - j0 + 1;
        for (k = 0; k < k1; k++, ic++) {
            if (ic->index >= qid) break;
            int c = ic->count < j1 ? ic->count : j1;
            uint32_t *idm = indexMapping.items + ic->index;

            // if(my_rank==1){
            // cerr<<"lookCounts.size     "<<lookCounts.size<<endl;
            // cerr<<"ic->index     "<<ic->index<<endl;
            // cerr<<"*idm     "<<*idm<<endl;
            // }

            if (*idm == 0) {
                if (rest < min) continue;
                IndexCount *ic2 = lookCounts.items + lookCounts.size;
                lookCounts.size += 1;
                *idm = lookCounts.size;
                ic2->index = ic->index;
                ic2->count = c;
            } else {
                lookCounts[*idm - 1].count += c;
            }
        }
    }
    // printf( "%6i %6i\n", S, lookCounts.size );
    lookCounts[lookCounts.size].count = 0;
    // printf( "\n\n" );
    return OK_FUNC;
}
Sequence::Sequence() {
    memset(this, 0, sizeof(Sequence));
    distance = 2.0;
}
Sequence::Sequence(const Sequence &other) {
    int i;
    // printf( "new: %p  %p\n", this, & other );
    memcpy(this, &other, sizeof(Sequence));
    distance = 2.0;
    if (other.data && other.master_flag == 1) {
        size = bufsize = other.size;

        data = new char[size + 1];
        true_data = new char[size + 1];
        // printf( "data: %p  %p\n", data, other.data );
        data[size] = 0;
        true_data[size] = 0;
        memcpy(data, other.data, size);
        memcpy(true_data, other.data, size);

        // for (i=0; i<size; i++) data[i] = other.data[i];
    }
    if (other.data && other.master_flag == 0) {
        size = bufsize = other.size;

        data = new char[size + 1];
        data[size] = 0;
        memcpy(data, other.data, size);
    }
    if (other.identifier) {
        int len = strlen(other.identifier);
        identifier = new char[len + 1];
        memcpy(identifier, other.identifier, len);
        identifier[len] = 0;
    }
}

// back to back merge for PE
// R1 -> XXXXXXABC ------------------- NMLYYYYYY <--R2
// >R1           >R2
// XXXXXXABC     YYYYYYLMN =====> Merge into
// >R12
// NMLYYYYYYXXXXXXABC
// Sequence::Sequence(const Sequence &other, const Sequence &other2, int mode) {
//     int i;
//     if (mode != 1) bomb_error("unknown mode");

//     // printf( "new: %p  %p\n", this, & other );
//     memcpy(this, &other, sizeof(Sequence));
//     distance = 2.0;

//     if (other.data && other2.data) {
//         size = bufsize = (other.size + other2.size);
//         size_R2 = other2.size;
//         data = new char[size + 1];
//         // printf( "data: %p  %p\n", data, other.data );
//         data[size] = 0;
//         data[size_R2] = 0;
//         memcpy(data, other2.data, size_R2);                 // copy R2 first
//         strrev(data);                                       // reverse R2 on data
//         memcpy(data + size_R2, other.data, size - size_R2); // copy R1 to end of R2
//         // for (i=0; i<size; i++) data[i] = other.data[i];
//         des_begin2 = other2.des_begin;
//         tot_length2 = other2.tot_length;
//     } else if (other.data || other2.data) {
//         bomb_error("Not both PE sequences have data");
//     }

//     if (other.identifier) { // only use R1
//         int len = strlen(other.identifier);
//         identifier = new char[len + 1];
//         memcpy(identifier, other.identifier, len);
//         identifier[len] = 0;
//     }
// }

Sequence::~Sequence() {
    // printf( "delete: %p\n", this );
    if (data) delete[] data;
    if (identifier) delete[] identifier;
    if (true_data) delete[] true_data;
}

void Sequence::Clear() {
    if (data) delete[] data;
    if (true_data) delete[] true_data;
    bufsize = 0;
    data = nullptr;
    true_data = nullptr;
}
void Sequence::worker_Clear() {
    if (data) delete[] data;
    if (identifier) delete[] identifier;
    if (true_data) delete[] true_data;
    bufsize = 0;
    data = nullptr;
    true_data = nullptr;
    identifier = nullptr;
}
void Sequence::operator=(const char *s) {
    size = 0; // avoid copying;
    Resize(strlen(s));
    strcpy(data, s);
    strcpy(true_data, s);
}
void Sequence::operator+=(const char *s) {
    int i, m = size, n = strlen(s);
    Reserve(m + n);
    memcpy(data + m, s, n);
}
void Sequence::Resize(int n) {
    int i, m = size < n ? size : n;
    size = n;
    if (size != bufsize) {
        char *old = data;
        bufsize = size;
        data = new char[bufsize + 1];
        if (data == NULL) bomb_error("Memory");
        if (old) {
            memcpy(data, old, m);
            delete[] old;
        }
        if (size) data[size] = 0;
    }
}
void Sequence::Reserve(int n) {
    int i, m = size < n ? size : n;
    size = n;
    if (size > bufsize) {
        char *old = data;
        bufsize = size + size / 5 + 1;
        data = new char[bufsize + 1];
        if (data == NULL) bomb_error("Memory");
        if (old) {
            memcpy(data, old, m);
            delete[] old;
        }
    }
    if (size) data[size] = 0;
}
void Sequence::trim(int trim_len) {
    if (trim_len >= size) return;
    size = trim_len;
    if (size) data[size] = 0;
}
void Sequence::ConvertBases() {
    int i;
    for (i = 0; i < size; i++) data[i] = aa2idx[data[i] - 'A'];
}

void Sequence::Swap(Sequence &other) {
    Sequence tmp;
    memcpy(&tmp, this, sizeof(Sequence));
    memcpy(this, &other, sizeof(Sequence));
    memcpy(&other, &tmp, sizeof(Sequence));
    memset(&tmp, 0, sizeof(Sequence));
}
int Sequence::Format() {
    int i, j = 0, m = 0;
    while (size && isspace(data[size - 1])) size--;
    if (size && data[size - 1] == '*') size--;
    if (size) data[size] = 0;
    for (i = 0; i < size; i++) {
        char ch = data[i];
        m += !(isalpha(ch) | isspace(ch));
    }
    if (m) return m;
    for (i = 0; i < size; i++) {
        char ch = data[i];
        if (isalpha(ch)) data[j++] = toupper(ch);
    }
    data[j] = 0;
    size = j;
    return 0;
}

void Sequence::PrintInfo(int id, FILE *fout, const Options &options, char *buf) {
    const char *tag = options.isEST ? "nt" : "aa";
    bool print = options.print != 0;
    bool strand = options.isEST;
    fprintf(fout, "%i\t%i%s, >%s...", id, size, tag, identifier + 1);
    if (identity) {
        int *c = coverage;
        fprintf(fout, " at ");
        if (print) fprintf(fout, "%i:%i:%i:%i/", c[0], c[1], c[2], c[3]);
        if (strand) fprintf(fout, "%c/", (state & IS_MINUS_STRAND) ? '-' : '+');
        fprintf(fout, "%.2f%%", identity * 100);
        if (options.useDistance) fprintf(fout, "/%.2f%%", distance * 100);
        fprintf(fout, "\n");
    } else {
        fprintf(fout, " *\n");
    }
}

static void write_run_fasta(const std::vector<std::pair<std::string, std::string>> &chunk, const std::string &path) {
    FILE *fout = std::fopen(path.c_str(), "wb");
    if (!fout) {
        perror(("fopen" + path).c_str());
        exit(1);
    }
    // Reduce syscall overhead when writing large run files.
    std::setvbuf(fout, nullptr, _IOFBF, 1 << 20);
    for (const auto &p : chunk) {
        std::fputc('>', fout);
        if (!p.first.empty()) std::fwrite(p.first.data(), 1, p.first.size(), fout);
        if (p.first.empty() || p.first.back() != '\n') std::fputc('\n', fout);

        if (!p.second.empty()) std::fwrite(p.second.data(), 1, p.second.size(), fout);
        if (p.second.empty() || p.second.back() != '\n') std::fputc('\n', fout);
    }
    fclose(fout);
}

void SequenceDB::Pipeline_External_Sort(const char *file, size_t chunk_size_bytes, std::vector<std::string> &run_files,
                                        Options &options, size_t core_num, size_t numa_size) {
    int option_l = options.min_length;
    first_chunk_size = 2000;
    total_num = 0;
    long long total_num_divide = 0;
    std::string max_name, min_name;
    total_letter = 0;
    total_desc = 0;
    max_len = 0;
    min_len = (size_t) -1;
    max_idf = 0;
    int T = options.threads;
    struct stat st;
    if (stat(file, &st) == 0) {
        std::cout << "File size: " << st.st_size << " bytes" << std::endl;
        size_t perhaps_file_size = st.st_size / 1024;
        chunk_size_bytes = max(perhaps_file_size, chunk_size_bytes);
        std::cout << "chunk_size_bytes: " << chunk_size_bytes << " bytes" << std::endl;
    } else
        perror("stat");
    std::string run_dir = MakeAbsolutePathFromCwd(temp_dir, "tmp_runs");
    if (mkdir(run_dir.c_str(), 0755) != 0 && errno != EEXIST) {
        perror(("mkdir " + run_dir).c_str());
        return;
    }

    std::cout << "kseq init" << endl;

    gzFile fp = gzopen(file, "rb");
    if (!fp) {
        perror("gzopen");
        return;
    }
    kseq_t *seq = kseq_init(fp);

    static constexpr int MAX_INFLIGHT = 8;
    static constexpr int MAX_WRITE_PAR = 1;

    std::atomic<int> inflight{0};
    std::atomic<int> write_slots{MAX_WRITE_PAR};
    std::atomic<int> run_id{0};
    run_files.clear();
    run_files.shrink_to_fit();
    std::vector<std::pair<std::string, std::string>> current_chunk;
    size_t current_chunk_size = 0;

    using clock = std::chrono::steady_clock;
    auto read_start = clock::now();

    omp_set_dynamic(0);
    omp_set_num_threads(T);
#pragma omp parallel
    {
#pragma omp single nowait
        {
            while (true) {
                while (inflight.load(std::memory_order_relaxed) >= MAX_INFLIGHT) {
#pragma omp taskyield
                }
                int len = kseq_read(seq);
                if (len < 0) break;
                if (len <= option_l) continue;
                std::string data = seq->seq.s ? seq->seq.s : "";
                std::string idf = seq->name.s ? seq->name.s : "";
                if (options.trim_len > 0) {
                    if (options.trim_len >= len) continue;
                    len = options.trim_len;
                    if (len) data[len] = 0;
                }
                total_num++;
                total_letter += len;
                total_desc += seq->name.l;
                if ((size_t) len > max_len) {
                    max_len = (size_t) len;
                    max_name = idf;
                }
                if (seq->name.l > max_idf) max_idf = seq->name.l;
                if ((size_t) len < min_len) {
                    min_len = (size_t) len;
                    min_name = idf;
                }
                current_chunk.emplace_back(std::move(idf), std::move(data));
                current_chunk_size += seq->name.l + (size_t) len;
                if (current_chunk_size >= chunk_size_bytes) {
                    auto read_end = clock::now();
                    auto read_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(read_end - read_start).count();
                    total_num_divide += current_chunk.size();
                    auto *ch = new std::vector<std::pair<std::string, std::string>>(std::move(current_chunk));
                    current_chunk.clear();
                    current_chunk_size = 0;
                    int my_run = run_id.fetch_add(1, std::memory_order_relaxed);
                    inflight.fetch_add(1, std::memory_order_relaxed);
#pragma omp task firstprivate(ch) depend(out : ch[0 : 1]) priority(1)
                    {
                        std::stable_sort(ch->begin(), ch->end(), [](auto const &a, auto const &b) {
                            return a.second.size() > b.second.size();
                        });
                    }
#pragma omp task firstprivate(ch, my_run, read_ns) depend(in : ch[0 : 1])
                    {
                        while (true) {
                            int slots = write_slots.load(std::memory_order_relaxed);
                            if (slots > 0 && write_slots.compare_exchange_weak(slots, slots - 1)) break;
#pragma omp taskyield
                        }
                        auto sort_write_start = clock::now();
                        std::string path = run_dir + "/run_" + std::to_string(my_run) + ".fa";
                        write_run_fasta(*ch, path);
                        auto sort_write_end = clock::now();
                        auto sort_write_ns =
                            std::chrono::duration_cast<std::chrono::nanoseconds>(sort_write_end - sort_write_start)
                                .count();
                        std::cout << read_ns << "\tReading\n" << sort_write_ns << "\tSort+Write\n";

                        write_slots.fetch_add(1, std::memory_order_relaxed);
                        delete ch;
                        inflight.fetch_sub(1, std::memory_order_relaxed);
                    }
                }
            }
            if (current_chunk_size > 0) {
                total_num_divide += current_chunk.size();
                auto *ch = new std::vector<std::pair<std::string, std::string>>(std::move(current_chunk));
                current_chunk.clear();
                current_chunk_size = 0;
                int my_run = run_id.fetch_add(1, std::memory_order_relaxed);
                inflight.fetch_add(1, std::memory_order_relaxed);

#pragma omp task firstprivate(ch) depend(out : ch[0 : 1]) priority(1)
                {
                    std::stable_sort(ch->begin(), ch->end(),
                                     [](auto const &a, auto const &b) { return a.second.size() > b.second.size(); });
                }
#pragma omp task firstprivate(ch, my_run) depend(in : ch[0 : 1])
                {
                    while (true) {
                        int slots = write_slots.load(std::memory_order_relaxed);
                        if (slots > 0 && write_slots.compare_exchange_weak(slots, slots - 1)) break;
#pragma omp taskyield
                    }
                    std::string path = run_dir + "/run_" + std::to_string(my_run) + ".fa";
                    write_run_fasta(*ch, path);
                    write_slots.fetch_add(1, std::memory_order_relaxed);
                    delete ch;
                    inflight.fetch_sub(1, std::memory_order_relaxed);
                }
            }
#pragma omp taskwait
        }
    }
    kseq_destroy(seq);
    gzclose(fp);
    options.max_entries = max_len * MAX_TABLE_SEQ;
    if (max_len >= 65536 && sizeof(INTs) <= 2) bomb_warning("Some seqs longer than 65536, you may define LONG_SEQ");
    if (max_len > MAX_SEQ)
        bomb_warning("Some seqs are too long, please rebuild the program with make parameter MAX_SEQ=...");
    for (int i = 0; i < run_id; i++) {
        std::string path = run_dir + "/run_" + std::to_string(i) + ".fa";
        run_files.push_back(path);
    }
    std::cout << "longest and shortest : " << max_name << " and " << min_name << "\n";
    std::cout << "longest and shortest : " << max_len << " and " << min_len << "\n";
    std::cout << "longest name : " << max_idf << "\n";
    std::cout << "Total letters: " << total_letter << "\n";
    std::cout << "total_num_divede: " << total_num_divide << "\n";
    std::cout << "Total number: " << total_num << "\n";
    if (total_num > 25000000) {
        chunk_size = 500000;
    } else {
        chunk_size = total_num / 50;
    }
    if (chunk_size < 100000) chunk_size = 100000;
    int total_threads = options.NodeNum * options.threads_per_node;
    int Consumption_threads;
    int Production_rate;
    int Consumption_rate;
    int factor = 2;
    // if (chunk_size >= 500000) factor = 4;
    // if (chunk_size <= 100000) factor = 2;
    int numa_num = options.threads_per_node/numa_size;
    int MPI_max = max(8,numa_num);
    cerr<<"MPI_max "<<MPI_max<<endl;
    for (int t = 4; t <= MPI_max; t = t * 2) {
        Production_threads = options.threads_per_node / t;
        Consumption_threads = total_threads - Production_threads;
        Production_rate = chunk_size / Production_threads;
        Consumption_rate = (total_num / 2) / Consumption_threads;
        cerr << "Consumption_rate " << Consumption_rate << endl;
        cerr << "Production_rate " << Production_rate << endl;
        if (Production_rate * factor < Consumption_rate)
            mpi_size = t;
        else
            break;
    }

    Production_threads = options.threads_per_node / mpi_size;
    // if (Production_threads > core_num || options.NodeNum > 4) {
    if (Production_threads > core_num ) {
        Production_threads = core_num;
        mpi_size = options.threads_per_node / Production_threads;
    }

    total_mpi_num = options.NodeNum * mpi_size;
    chunk_bytes = total_letter / total_num * chunk_size;
    std::cout << "chunk_size: " << chunk_size << std::endl;
    std::cout << "threads_per_node : " << Production_threads << std::endl;
    std::cout << "total_mpi_num : " << total_mpi_num << std::endl;
}

void SequenceDB::WriteToJSON(const std::string &file, const std::string &output_dir, const std::string &output_prefix,
                             int num_procs) {
    // 构造 JSON 对象
    using json = nlohmann::json;
    json j;

    j["info"] = {
        {"max_len", max_len},
        {"max_idf", max_idf},
        {"chunks_num", chunks_num},
        {"chunk_size", chunk_size},
        {"total_num", total_num},
        {"min_len", min_len},
        {"num_procs", num_procs},
        {"total_letter", total_letter},
        {"total_desc", total_desc},
        {"total_chunk", total_chunk},
        {"chunk_bytes", chunk_bytes},
        {"first_chunk_size", first_chunk_size},
        {"threads_per_node", Production_threads},
        {"total_mpi_num", total_mpi_num},

    };

    j["files"] = {{"output_dir", output_dir}, {"out_prefix", output_prefix}};

    std::ofstream ofs(output_dir + file, std::ios::binary);
    if (!ofs) return;
    ofs << j.dump(4) << "\n";
}

void SequenceDB::ReadJsonInfo(const std::string &file, const std::string &output_dir, Options &options, bool master) {
    using json = nlohmann::json;
    std::string path = output_dir + file;
    std::ifstream ifs(path);
    if (!ifs.is_open()) {
        throw std::runtime_error("Cannot open JSON file: " + path);
    }
    json j;
    ifs >> j;

    max_len = j["info"]["max_len"].get<size_t>();
    max_idf = j["info"]["max_idf"].get<size_t>();
    chunks_num = j["info"]["chunks_num"].get<int>();
    chunk_size = j["info"]["chunk_size"].get<int>();
    total_num = j["info"]["total_num"].get<long long>();
    min_len = j["info"]["min_len"].get<size_t>();
    chunk_bytes = j["info"]["chunk_bytes"].get<long long>();
    first_chunk_size = j["info"]["first_chunk_size"].get<int>();
    Production_threads = j["info"]["threads_per_node"].get<int>();
    total_mpi_num = j["info"]["total_mpi_num"].get<int>();
    string temp_dir = j["files"]["output_dir"].get<string>();
    if (master) {
        total_letter = j["info"]["total_letter"].get<long long>();
        total_desc = j["info"]["total_desc"].get<long long>();
        total_chunk = j["info"]["total_chunk"].get<int>();
    }
    options.max_entries = max_len * MAX_TABLE_SEQ;
}

void SequenceDB::MergeSortedRuns_KWay(const std::vector<std::string> &run_files, const std::string &output_prefix) {
    if (run_files.empty()) return;
    int num_procs = total_mpi_num - 1;
    std::priority_queue<FastaRecord> pq;
    std::vector<FILE *> fps(run_files.size(), nullptr);

    for (size_t i = 0; i < run_files.size(); ++i) {
        FILE *fp = fopen(run_files[i].c_str(), "r");
        if (!fp) {
            std::cerr << "Failed to open run file: " << run_files[i] << std::endl;
            continue;
        }
        fps[i] = fp;

        char desc[MAX_LINE_SIZE], seq[MAX_LINE_SIZE];
        if (fgets(desc, MAX_LINE_SIZE, fp) && fgets(seq, MAX_LINE_SIZE, fp)) {
            FastaRecord rec;
            rec.desc = desc;
            rec.seq = seq;
            rec.file_id = i;
            pq.push(rec);
        }
    }

    struct ProcFile {
        FILE *fp;
    };

    std::vector<ProcFile> proc_files(num_procs);
    for (int i = 0; i < num_procs; ++i) {
        std::string filename = output_prefix + "_proc" + std::to_string(i) + ".fa";

        proc_files[i].fp = fopen(filename.c_str(), "w");
        if (!proc_files[i].fp) {
            fprintf(stderr, "FATAL: Failed to create output file %s (%s)\n", filename.c_str(), strerror(errno));
            exit(EXIT_FAILURE);
        }
    }

    size_t global_chunk_id = -1;
    size_t current_chunk_size = 0;
    size_t current_chunk_bytes = 0;
    int current_proc = -1;
    auto rotate_chunk = [&] {
        if (current_proc != -1) {
            auto &pf = proc_files[current_proc];
            fflush(pf.fp);
        }

        global_chunk_id++;
        current_proc = (global_chunk_id) % num_procs;
        current_chunk_size = 0;
        current_chunk_bytes = 0;

        auto &pf = proc_files[current_proc];
    };

    rotate_chunk();

    while (!pq.empty()) {
        FastaRecord rec = pq.top();
        pq.pop();

        auto &pf = proc_files[current_proc];
        fwrite(rec.desc.data(), 1, rec.desc.size(), pf.fp);
        fwrite(rec.seq.data(), 1, rec.seq.size(), pf.fp);
        current_chunk_bytes += rec.seq.length() - 1;
        current_chunk_size++;
        int fid = rec.file_id;
        char desc[MAX_LINE_SIZE], seq[MAX_LINE_SIZE];
        if (fgets(desc, MAX_LINE_SIZE, fps[fid]) && fgets(seq, MAX_LINE_SIZE, fps[fid])) {
            FastaRecord new_rec;
            new_rec.desc = desc;
            new_rec.seq = seq;
            new_rec.file_id = fid;
            pq.push(new_rec);
        }

        if ((current_chunk_bytes > chunk_bytes) || (chunks_num == 0 && current_chunk_size >= first_chunk_size) ||
            (current_chunk_size >= chunk_size)) {
            if (current_chunk_bytes > chunk_bytes && chunks_num == 0) first_chunk_size = current_chunk_size;
            rotate_chunk();
            chunks_num++;
        }
    }
    if (current_chunk_bytes > 0) {
        auto &pf = proc_files[current_proc];
        chunks_num++;
        fflush(pf.fp);
    }
    total_chunk = global_chunk_id;

    for (auto &pf : proc_files) fclose(pf.fp);

    for (auto *fp : fps)
        if (fp) fclose(fp);

    for (const auto &fname : run_files) remove(fname.c_str());
    proc_files.clear();
    std::vector<ProcFile>().swap(proc_files);
    fps.clear();
    std::vector<FILE *>().swap(fps);
    pq = std::priority_queue<FastaRecord>();
    WriteToJSON("info.json", output_prefix, "_proc", num_procs);
    std::cout << "Successfully write info.json!\n";
    std::cout << "chunk_num: " << chunks_num << std::endl;
    if (chunks_num < total_mpi_num) std::cout << "Warring:There is a waste of computing resources  " << std::endl;
}

void SequenceDB::read_sorted_files(const std::string &temp_dir, int rank, int rank_size, bool mpi_status,
                                   MPI_Comm worker_comm, Options &options) {
    int file_index = rank;
    std::string file = temp_dir + "_proc" + std::to_string(rank - 1) + ".fa";
    int start_my_id = sequences.size();
    int chunk_id = (rank - 1);
    long long now_bytes = 0;
    long now_num = 0;
    gzFile fp = gzopen(file.c_str(), "r");
    kseq_t *seq = kseq_init(fp);
    int len;
    Sequence one;
    char *id_ptr = new char[max_idf + 1];
    char *data_ptr = new char[max_len + 1];
    while ((len = kseq_read(seq)) >= 0) {
        memcpy(id_ptr, seq->name.s, seq->name.l);
        id_ptr[seq->name.l] = 0;
        one.identifier = id_ptr;

        memcpy(data_ptr, seq->seq.s, seq->seq.l);
        data_ptr[seq->seq.l] = 0;
        one.data = data_ptr;
        one.size = len;
        one.master_flag = 0;
        sequences.Append(new Sequence(one));
        now_bytes += len;
        now_num++;
        if ((now_bytes > chunk_bytes) || (rank == 1 && chunks_id.size() == 0 && sequences.size() >= first_chunk_size) ||
            (now_num >= chunk_size)) {
            my_chunks.push_back(make_pair(start_my_id, sequences.size() - 1));
            cerr << "chunk_id    " << chunk_id << endl;
            chunks_id.push_back(chunk_id);
            chunk_id = (rank_size - 1) + chunk_id;
            start_my_id = sequences.size();
            now_bytes = 0;
            now_num = 0;
        }
    }
    one.identifier = nullptr;
    one.data = nullptr;
    delete[] id_ptr;
    delete[] data_ptr;
    if (now_bytes > 0) {
        my_chunks.push_back(make_pair(start_my_id, sequences.size() - 1));
        chunks_id.push_back(chunk_id);
        cerr << "chunk_id    " << chunk_id << endl;
    }
    kseq_destroy(seq);
    gzclose(fp);

#pragma omp parallel for num_threads(options.threads)
    for (int i = 0; i < sequences.size(); i++) {
        Sequence *seq = sequences[i];
        seq->ConvertBases();
    }

    if (options.stealing) {
        // SUB = chunk_size/10000;
        for (size_t ci = 0; ci < my_chunks.size(); ++ci) {
            int L = my_chunks[ci].first;
            int R = my_chunks[ci].second;

            const int n = R - L + 1;
            int base = n / SUB;
            int rem = n % SUB;

            int cur = L;
            for (int s = 0; s < SUB; ++s) {
                int len = base + (s < rem ? 1 : 0);
                if (len == 0) continue;
                int l = cur;
                int r = cur + len - 1;
                cur = r + 1;

                sub_chunks.emplace_back(l, r);
                sub_chunks_id.push_back(chunks_id[ci] * SUB + s);
            }
        }
        cerr << "SUB  " << SUB << endl;

        tasks_local_.clear();
        tasks_local_.reserve(sub_chunks.size());
        tasks_flag.assign(sub_chunks.size(), 0);
        for (auto &pr : sub_chunks) tasks_local_.push_back(Task{pr.first, pr.second});

        ctrl_[0] = 0;
        ctrl_[1] = (int) tasks_local_.size() - 1;
        ctrl_[2] = 0;
        MPI_Win_create(tasks_local_.empty() ? MPI_BOTTOM : (void *) tasks_local_.data(),
                       (MPI_Aint) tasks_local_.size() * sizeof(Task), sizeof(Task), MPI_INFO_NULL, worker_comm,
                       &win_tasks_);
        // MPI_Win_create(tasks_flag.empty() ? MPI_BOTTOM : (void *) tasks_flag.data(), tasks_flag.size() * sizeof(int),
        //                sizeof(int), MPI_INFO_NULL, worker_comm, &win_tasks_flag_);

        MPI_Win_create((void *)ctrl_, 3 * (MPI_Aint)sizeof(int), sizeof(int), MPI_INFO_NULL, worker_comm, &win_ctrl_);
        ptr_ctrl_ = (volatile int*) ctrl_;
        MPI_Win_lock_all(0, win_tasks_);
        MPI_Win_lock_all(0, win_ctrl_);
        // MPI_Win_lock_all(0, win_tasks_flag_);
        const size_t Nseq = sequences.size();
        meta_.clear();
        meta_.resize(Nseq);
        pool_data_.clear();
        size_t est_data = 0;
        for (size_t k = 0; k < Nseq; ++k) {
            Sequence *s = sequences[k];
            est_data += (size_t) s->size;
        }
        pool_data_.reserve(est_data);
        size_t off_d = 0;
        for (size_t k = 0; k < Nseq; ++k) {
            Sequence *s = sequences[k];
            SeqMeta m{};

            // data
            m.data_off = off_d;
            m.data_len = (int32_t) s->size;
            if (s->size > 0 && s->data) {
                pool_data_.insert(pool_data_.end(), (uint8_t *) s->data, (uint8_t *) s->data + s->size);
                off_d += s->size;
            }
            m.size = s->size;
            m.state = s->state;
            m.cluster_id = s->cluster_id;
            m.identity = s->identity;
            m.distance = s->distance;
            for (int t = 0; t < 4; ++t) m.coverage[t] = s->coverage[t];

            meta_[k] = m;
        }
        MPI_Win_create(pool_data_.empty() ? MPI_BOTTOM : (void *) pool_data_.data(), (MPI_Aint) pool_data_.size(), 1,
                       MPI_INFO_NULL, worker_comm, &win_pool_d_);
        MPI_Win_create(meta_.empty() ? MPI_BOTTOM : (void *) meta_.data(), (MPI_Aint) meta_.size() * sizeof(SeqMeta),
                       sizeof(SeqMeta), MPI_INFO_NULL, worker_comm, &win_meta_);
        MPI_Win_lock_all(0, win_meta_);
        MPI_Win_lock_all(0, win_pool_d_);
    }
}
void SequenceDB::Encodeseqs(Sequence *seq, int NAA, int id, bool est) {
    char *seqi = seq->data;
    int len = seq->size;
    // check_word_encodes
    int aan_no = len - NAA + 1;
    int i, j, i0, i1;
    int skip = 0;
    unsigned char k, k1;
    vector<int> word_encodes(len);
    vector<INTs> word_encodes_no(len);
    for (j = 0; j < aan_no; j++) {
        char *word = seqi + j;
        int encode = 0;
        for (k = 0, k1 = NAA - 1; k < NAA; k++, k1--) encode += word[k] * NAAN_array[k1];
        word_encodes[j] = encode;
    }

    if (est) {
        for (j = 0; j < len; j++) {
            if (seqi[j] >= 4) { // here N is 4
                i0 = (j - NAA + 1 > 0) ? j - NAA + 1 : 0;
                i1 = j < aan_no ? j : aan_no - 1;
                for (i = i0; i <= i1; i++) word_encodes[i] = -1;
            }
        }
        for (j = 0; j < aan_no; j++) skip += (word_encodes[j] == -1);
    }

    std::sort(word_encodes.begin(), word_encodes.begin() + aan_no);
    for (j = 0; j < aan_no; j++) word_encodes_no[j] = 1;
    for (j = aan_no - 1; j; j--) {
        if (word_encodes[j] == word_encodes[j - 1]) {
            word_encodes_no[j - 1] += word_encodes_no[j];
            word_encodes_no[j] = 0;
        }
    }
    for (j = 0; j < aan_no; j++) {
        if (word_encodes_no[j] != 0) {
            total_encodes[id].emplace_back(word_encodes[j]);
            total_encodes_no[id].emplace_back(word_encodes_no[j]);
        }
    }
    total_encodes[id].emplace_back(0);
    total_encodes_no[id].emplace_back(0);
}

void SequenceDB::WriteClustersSort(const char *input, const char *output, const Options &options) {
    ofstream fout(output);
    int i, j, n = rep_seqs.size();
    int count, rest;
    char *buf = new char[MAX_LINE_SIZE + 1];
    vector<uint64_t> sorting(n);
    if (!fout) throw runtime_error("Cannot open output file");
    for (i = 0; i < n; i++) {
        Sequence *seq = sequences[rep_seqs[i]];
        fout << ">" << seq->identifier << "\n";
        // fout << seq->data << "\n";
    }
    fout.close();
    return;
}

void SequenceDB::WriteClusterDetail(const Options &options) {
    string output = options.output + ".clstr";
    string db_clstr_bak = options.output + ".bak.clstr";
    int i, i0, k, N = total_num;
    FILE *fout;
    char *buf = new char[MAX_DES + 1];

    if (options.backupFile) {
        fout = fopen(db_clstr_bak.c_str(), "w+");
        for (i = 0; i < N; i++) {
            Sequence *seq = sequences[i];
            seq->PrintInfo(seq->cluster_id, fout, options, buf);
        }
        fclose(fout);
    }

    cout << "writing clustering information" << endl;
    int rep_size = rep_seqs.size();
    vector<vector<int>> clusters(rep_size);
    for (int i = 0; i < total_num; i++) {
        int id = sequences[i]->cluster_id;
        clusters[id].push_back(i);
    }
    fout = fopen(output.c_str(), "w+");

    if (options.sort_output) {
        int *clstr_size = new int[rep_size];
        int *clstr_idx1 = new int[rep_size];

        for (i = 0; i < rep_size; i++) {
            clstr_size[i] = (int) clusters[i].size();
            clstr_idx1[i] = i;
        }
        quick_sort_idxr(clstr_size, clstr_idx1, 0, rep_size - 1);

        for (i = 0; i < rep_size; i++) {
            i0 = clstr_idx1[i];
            fprintf(fout, ">Cluster %i\n", i);
            for (k = 0; k < (int) clusters[i0].size(); k++)
                sequences[clusters[i0][k]]->PrintInfo(k, fout, options, buf);
        }
    } else {
        for (i = 0; i < rep_size; i++) {
            fprintf(fout, ">Cluster %i\n", i);
            for (k = 0; k < (int) clusters[i].size(); k++) sequences[clusters[i][k]]->PrintInfo(k, fout, options, buf);
        }
    }

    delete[] buf;
}

void SequenceDB::WriteExtra1D(const Options &options) {
    string db_clstr = options.output + ".clstr";
    string db_clstr_bak = options.output + ".bak.clstr";
    int i, i0, k, N = sequences.size();
    vector<long long> sorting(N);
    for (i = 0; i < N; i++) sorting[i] = ((long long) sequences[i]->index << 32) | i;
    std::sort(sorting.begin(), sorting.end());

    FILE *fout;
    char *buf = new char[MAX_DES + 1];

    if (options.backupFile) {
        fout = fopen(db_clstr_bak.c_str(), "w+");
        for (i = 0; i < N; i++) {
            Sequence *seq = sequences[sorting[i] & 0xffffffff];
            seq->PrintInfo(seq->cluster_id, fout, options, buf);
        }
        fclose(fout);
    }

    cout << "writing clustering information" << endl;
    int M = rep_seqs.size();
    Vector<Vector<int>> clusters(M);
    for (i = 0; i < N; i++) {
        int k = sorting[i] & 0xffffffff;
        int id = sequences[k]->cluster_id;
        clusters[id].Append(k);
    }

    fout = fopen(db_clstr.c_str(), "w+");

    if (options.sort_output) {
        int *clstr_size = new int[M];
        int *clstr_idx1 = new int[M];

        for (i = 0; i < M; i++) {
            clstr_size[i] = (int) clusters[i].size();
            clstr_idx1[i] = i;
        }
        quick_sort_idxr(clstr_size, clstr_idx1, 0, M - 1);

        for (i = 0; i < M; i++) {
            i0 = clstr_idx1[i];
            fprintf(fout, ">Cluster %i\n", i);
            for (k = 0; k < (int) clusters[i0].size(); k++)
                sequences[clusters[i0][k]]->PrintInfo(k, fout, options, buf);
        }
    } else {
        for (i = 0; i < M; i++) {
            fprintf(fout, ">Cluster %i\n", i);
            for (k = 0; k < (int) clusters[i].size(); k++) sequences[clusters[i][k]]->PrintInfo(k, fout, options, buf);
        }
    }

    delete[] buf;
}
void SequenceDB::WriteExtra2D(SequenceDB &other, const Options &options) {
    string db_clstr = options.output + ".clstr";
    string db_clstr_bak = options.output + ".bak.clstr";
    int i, k, N = other.sequences.size();
    int N2 = sequences.size();
    vector<long long> sorting(N);
    for (i = 0; i < N; i++) sorting[i] = ((long long) other.sequences[i]->index << 32) | i;
    std::sort(sorting.begin(), sorting.end());

    FILE *fout;
    char *buf = new char[MAX_DES + 1];
    if (options.backupFile) {
        fout = fopen(db_clstr_bak.c_str(), "w+");
        for (i = 0; i < N; i++) {
            Sequence *seq = other.sequences[sorting[i] & 0xffffffff];
            seq->PrintInfo(seq->cluster_id, fout, options, buf);
        }
        for (i = 0; i < N2; i++) {
            Sequence *seq = sequences[i];
            if (seq->state & IS_REDUNDANT) seq->PrintInfo(seq->cluster_id, fout, options, buf);
        }
        fclose(fout);
    }

    cout << "writing clustering information" << endl;
    Vector<Vector<int>> clusters(N);
    for (i = 0; i < N2; i++) {
        int id = sequences[i]->cluster_id;
        if (sequences[i]->state & IS_REDUNDANT) clusters[id].Append(i);
    }

    fout = fopen(db_clstr.c_str(), "w+");
    for (i = 0; i < N; i++) {
        Sequence *seq = other.sequences[i];
        fprintf(fout, ">Cluster %i\n", i);
        seq->PrintInfo(0, fout, options, buf);
        for (k = 0; k < (int) clusters[i].size(); k++) sequences[clusters[i][k]]->PrintInfo(k + 1, fout, options, buf);
    }
    delete[] buf;
}
void WorkingParam::ControlShortCoverage(int len, const Options &options) {
    len_eff = len;
    aln_cover_flag = 0;
    if ((options.short_coverage > 0.0) || (options.min_control > 0)) { // has alignment coverage control
        aln_cover_flag = 1;
        min_aln_lenS = (int) (double(len) * options.short_coverage);
        if (len - options.short_control > min_aln_lenS) min_aln_lenS = len - options.short_control;
        if (options.min_control > min_aln_lenS) min_aln_lenS = options.min_control;
    }
    if (options.global_identity == 0) len_eff = min_aln_lenS; // global_identity==0
}
void WorkingParam::ControlLongCoverage(int len2, const Options &options) {
    if (aln_cover_flag) {
        min_aln_lenL = (int) (double(len2) * options.long_coverage);
        if (len2 - options.long_control > min_aln_lenL) min_aln_lenL = len2 - options.long_control;
        if (options.min_control > min_aln_lenL) min_aln_lenL = options.min_control;
    }
}

// when alignment coverage such as -aL is specified
// if a existing rep is too long, it won't be qulified
int upper_bound_length_rep(int len, double opt_s, int opt_S, double opt_aL, int opt_AL) {
    int len_upper_bound = 99999999;
    double r1 = (opt_s > opt_aL) ? opt_s : opt_aL;
    int a2 = (opt_S < opt_AL) ? opt_S : opt_AL;
    if (r1 > 0.0) len_upper_bound = (int) (((float) len) / r1);
    if ((len + a2) < len_upper_bound) len_upper_bound = len + a2;

    return len_upper_bound;
} // END upper_bound_length_rep
int upper_bound_length_rep(int len, const Options &options) {
    double opt_s = options.diff_cutoff;
    int opt_S = options.diff_cutoff_aa;
    double opt_aL = options.long_coverage;
    int opt_AL = options.long_control;
    return upper_bound_length_rep(len, opt_s, opt_S, opt_aL, opt_AL);
}

void cal_aax_cutoff(double &aa1_cutoff, double &aa2_cutoff, double &aan_cutoff, double cluster_thd, int tolerance,
                    int naa_stat_start_percent, int naa_stat[5][61][4], int NAA) {
    aa1_cutoff = cluster_thd;
    aa2_cutoff = 1 - (1 - cluster_thd) * 2;
    aan_cutoff = 1 - (1 - cluster_thd) * NAA;
    if (tolerance == 0) return;

    int clstr_idx = (int) (cluster_thd * 100) - naa_stat_start_percent;
    if (clstr_idx < 0) clstr_idx = 0;
    double d2 = ((double) (naa_stat[tolerance - 1][clstr_idx][3])) / 100;
    double dn = ((double) (naa_stat[tolerance - 1][clstr_idx][5 - NAA])) / 100;
    aa2_cutoff = d2 > aa2_cutoff ? d2 : aa2_cutoff;
    aan_cutoff = dn > aan_cutoff ? dn : aan_cutoff;
    return;
} // END cal_aax_cutoff

void update_aax_cutoff(double &aa1_cutoff, double &aa2_cutoff, double &aan_cutoff, int tolerance,
                       int naa_stat_start_percent, int naa_stat[5][61][4], int NAA, double cluster_thd) {
    if (cluster_thd > 1.0) cluster_thd = 1.00;

    double aa1_t, aa2_t, aan_t;
    cal_aax_cutoff(aa1_t, aa2_t, aan_t, cluster_thd, tolerance, naa_stat_start_percent, naa_stat, NAA);
    if (aa1_t > aa1_cutoff) aa1_cutoff = aa1_t;
    if (aa2_t > aa2_cutoff) aa2_cutoff = aa2_t;
    if (aan_t > aan_cutoff) aan_cutoff = aan_t;
    return;
} // END update_aax_cutoff
// 聚类参数初始化
void WorkingParam::ComputeRequiredBases(int NAA, int ss, const Options &option) {
    // d: distance, fraction of errors;
    // e: number of errors;
    // g: length of the maximum gap;
    // m: word length;
    // n: sequence length;
    // alignment length = n - g + 1;
    // d = e / (n - g + 1);
    // e >= 1, so that, g <= n + 1 - 1/d
    // word count = (n - g - m + 1) - (e - 1)*m;
    //            = (n - g - m + 1) - (d*(n - g + 1) - 1)*m
    //            = (n - g + 1) - d*m*(n - g + 1)
    //            = (n - g + 1)*(1 - d*m)
    // minimum word count is reached when g == n + 1 - 1/d
    // so, minimum word count = 1/d - m.
    // if g == band_width: word count = (n - band + 1)*(1 - d*m);
    if (options.useDistance) {
        int band = options.band_width + 1;
        int invd = int(1.0 / (options.distance_thd + 1E-9));
        int k = len_eff < invd ? len_eff : invd;
        int ks = len_eff - ss + 1;
        int kn = len_eff - NAA + 1;
        int ks2 = invd - ss;
        int kn2 = invd - NAA;
        int ks3 = int((len_eff - band + 1.0) * (1.0 - options.distance_thd * ss));
        int kn3 = int((len_eff - band + 1.0) * (1.0 - options.distance_thd * NAA));
        // if( ks3 > ks2 ) ks2 = ks3;
        // if( kn3 > kn2 ) kn2 = kn3;
        required_aa1 = required_aas = (ks2 < ks ? ks2 : ks);
        required_aan = kn2 < kn ? kn2 : kn;
        if (required_aa1 <= 0) required_aa1 = required_aas = 1;
        if (required_aan <= 0) required_aan = 1;
        // required_aa1 = required_aas = required_aan = 0;
        return;
    }
    // (N-K)-K*(1-C)*N = C*K*N-(K-1)*N-K = (C*K-K+1)*N-K
    required_aa1 = (len_eff - ss) - int(ss * ceil((1.0 - aa1_cutoff) * len_eff));
    if (required_aa1 < 0) required_aa1 = 0;
    required_aas = required_aa1;
    required_aan = (len_eff - NAA) - int(NAA * ceil((1.0 - aa1_cutoff) * len_eff));
    // printf( "%i %i\n", required_aa1, required_aan );
    if (required_aan < 0) required_aan = 0;

    int aa1_old = int(aa1_cutoff * (double) len_eff) - ss + 1;
    int aas_old = int(aas_cutoff * (double) len_eff);
    int aan_old = int(aan_cutoff * (double) len_eff);

    double thd = option.cluster_thd;
    // double rest = (len_eff - ss) / double(len_eff * ss);
    double rest = (len_eff - NAA) / double(len_eff * NAA);
    double thd0 = 1.0 - rest;
    double fnew = 0;
    double fold = 1;
    if (thd > thd0) {
        fnew = (thd - thd0) / rest;
        fold = 1.0 - fnew;
    }
    // printf( "%g %g %g\n", thd, thd0, fnew );

    required_aa1 = (int) (fnew * required_aa1 + fold * aa1_old);
    required_aas = (int) (fnew * required_aas + fold * aas_old);
    required_aan = (int) (fnew * required_aan + fold * aan_old);
}
// 编码 k-mer NAA-kmer 5
int WorkingBuffer::EncodeWords(Sequence *seq, int NAA, bool est) {
    char *seqi = seq->data;
    int len = seq->size;
    // check_word_encodes
    int aan_no = len - NAA + 1;
    int i, j, i0, i1;
    int skip = 0;
    // unsigned char k, k1;
    // for (j = 0; j < aan_no; j++) {
    //     char *word = seqi + j;
    //     int encode = 0;
    //     for (k = 0, k1 = NAA - 1; k < NAA; k++, k1--) encode += word[k] * NAAN_array[k1];
    //     word_encodes[j] = word_encodes_backup[j] = encode;
    // }
    const int base    = NAAN_array[1];
    const int basePow = NAAN_array[NAA - 1];

    // 1) 初始化第一个窗口编码：O(NAA)
    int enc = 0;
    for (int k = 0; k < NAA; ++k) {
        const int k1 = NAA - 1 - k;
        enc += (int)seqi[k] * NAAN_array[k1];
    }

    // 2) rolling 生成所有窗口编码：O(aan_no)
    // enc(j+1) = (enc(j) - seq[j]*basePow) * base + seq[j+NAA]
    for (int j = 0; j < aan_no; ++j) {
        word_encodes[j] = word_encodes_backup[j] = enc;

        // update to next window
        if (j + NAA < len) {
            enc = (enc - (int)seqi[j] * basePow) * base + (int)seqi[j + NAA];
        }
    }


    if (est) {
        for (j = 0; j < len; j++) {
            if (seqi[j] >= 4) { // here N is 4
                i0 = (j - NAA + 1 > 0) ? j - NAA + 1 : 0;
                i1 = j < aan_no ? j : aan_no - 1;
                for (i = i0; i <= i1; i++) word_encodes[i] = -1;
            }
        }
        for (j = 0; j < aan_no; j++) skip += (word_encodes[j] == -1);
    }
    // assert(aan_no<35808);
    //ips4o::sort(word_encodes.begin(), word_encodes.begin() + aan_no);
    std::sort(word_encodes.begin(), word_encodes.begin() + aan_no);
    for (j = 0; j < aan_no; j++) word_encodes_no[j] = 1;
    for (j = aan_no - 1; j; j--) {
        if (word_encodes[j] == word_encodes[j - 1]) {
            word_encodes_no[j - 1] += word_encodes_no[j];
            word_encodes_no[j] = 0;
        }
    }
    return skip;
    // END check_word_encodes
}

void WorkingBuffer::ComputeAAP(const char *seqi, int size) {
    int len1 = size - 1;
    int sk, j1, mm, c22;
    for (sk = 0; sk < NAA2; sk++) taap[sk] = 0;
    for (j1 = 0; j1 < len1; j1++) {
        c22 = seqi[j1] * NAA1 + seqi[j1 + 1];
        taap[c22]++;
    }
    for (sk = 0, mm = 0; sk < NAA2; sk++) {
        aap_begin[sk] = mm;
        mm += taap[sk];
        taap[sk] = 0;
    }
    for (j1 = 0; j1 < len1; j1++) {
        c22 = seqi[j1] * NAA1 + seqi[j1 + 1];
        aap_list[aap_begin[c22] + taap[c22]++] = j1;
    }
}
void WorkingBuffer::ComputeAAP2(const char *seqi, int size) {
    int len1 = size - 3;
    int sk, j1, mm, c22;
    for (sk = 0; sk < NAA4; sk++) taap[sk] = 0;
    for (j1 = 0; j1 < len1; j1++) {
        if ((seqi[j1] >= 4) || (seqi[j1 + 1] >= 4) || (seqi[j1 + 2] >= 4) || (seqi[j1 + 3] >= 4)) continue; // skip N
        c22 = seqi[j1] * NAA3 + seqi[j1 + 1] * NAA2 + seqi[j1 + 2] * NAA1 + seqi[j1 + 3];
        taap[c22]++;
    }
    for (sk = 0, mm = 0; sk < NAA4; sk++) {
        aap_begin[sk] = mm;
        mm += taap[sk];
        taap[sk] = 0;
    }
    for (j1 = 0; j1 < len1; j1++) {
        if ((seqi[j1] >= 4) || (seqi[j1 + 1] >= 4) || (seqi[j1 + 2] >= 4) || (seqi[j1 + 3] >= 4)) continue; // skip N
        c22 = seqi[j1] * NAA3 + seqi[j1 + 1] * NAA2 + seqi[j1 + 2] * NAA1 + seqi[j1 + 3];
        aap_list[aap_begin[c22] + taap[c22]++] = j1;
    }
}
void SequenceDB::ClusterOne(Sequence *seq, int id, WordTable &table, WorkingParam &param, WorkingBuffer &buffer,
                            const Options &options) {
    if (seq->state & IS_REDUNDANT) return;
    int frag_size = options.frag_size;
    int NAA = options.NAA;
    int len = seq->size;
    int len_bound = upper_bound_length_rep(len, options);
    param.len_upper_bound = len_bound;
    int flag = CheckOne(seq, table, param, buffer, options);

    if (flag == 0) {
        if ((seq->identity > 0) && (options.cluster_best)) {
            // because of the -g option, this seq is similar to seqs in old SEGs
            seq->state |= IS_REDUNDANT;
            seq->Clear();
        } else { // else add to NR90 db
            int aan_no = len - NAA + 1;
            int size = rep_seqs.size();
            rep_seqs.Append(id);
            seq->cluster_id = size;
            seq->identity = 0;
            seq->state |= IS_REP;
            if (frag_size) { /* not used for EST */
                int frg1 = (len - NAA) / frag_size + 1;
                table.AddWordCountsFrag(aan_no, buffer.word_encodes_backup, buffer.word_encodes_no, frg1, frag_size);
            } else {
                table.AddWordCounts(aan_no, buffer.word_encodes, buffer.word_encodes_no, table.sequences.size(),
                                    options.isEST);
            }
            table.sequences.Append(seq);
            if (frag_size) {
                while (table.sequences.size() < table.frag_count) table.sequences.Append(seq);
            }
        }
    }
    if ((id + 1) % 1000 == 0) {
        int size = rep_seqs.size();
        printf(".");
        fflush(stdout);
        if ((id + 1) % 10000 == 0) printf("\r..........%9i  finished  %9i  clusters\n", id + 1, size);
    }
}
void SequenceDB::ClusterOne(Sequence *seq, int id, WordTable &table, WorkingParam &param, WorkingBuffer &buffer,
                            const Options &options, int my_rank) {
    if (seq->state & IS_REDUNDANT) return;
    int frag_size = options.frag_size;
    int NAA = options.NAA;
    int len = seq->size;
    int len_bound = upper_bound_length_rep(len, options);
    param.len_upper_bound = len_bound;
    int flag = CheckOne(seq, table, param, buffer, options);

    if (flag == 0) {
        if ((seq->identity > 0) && (options.cluster_best)) {
            // because of the -g option, this seq is similar to seqs in old SEGs
            seq->state |= IS_REDUNDANT;
            seq->Clear();
        } else { // else add to NR90 db
            int aan_no = len - NAA + 1;
            int size = rep_seqs.size();
            rep_seqs.Append(id);
            seq->cluster_id = size;
            seq->identity = 0;
            seq->state |= IS_REP;
            if (frag_size) { /* not used for EST */
                int frg1 = (len - NAA) / frag_size + 1;
                table.AddWordCountsFrag(aan_no, buffer.word_encodes_backup, buffer.word_encodes_no, frg1, frag_size);
            } else {
                table.AddWordCounts(aan_no, buffer.word_encodes, buffer.word_encodes_no, table.sequences.size(),
                                    options.isEST);
            }
            table.sequences.Append(seq);
            if (frag_size) {
                while (table.sequences.size() < table.frag_count) table.sequences.Append(seq);
            }
        }
    }
    if ((id + 1) % 1000 == 0) {
        int size = rep_seqs.size();
        printf(".");
        fflush(stdout);
        if ((id + 1) % 10000 == 0) printf("\r..........%9i  finished  %9i  clusters\n", id + 1, size);
    }
}
void SequenceDB::ClusterOne(Sequence *seq, int id, WordTable &local_table, WordTable &table, WorkingParam &param,
                            WorkingBuffer &buffer, const Options &options, int my_rank) {
    if (seq->state & IS_REDUNDANT) return;
    int frag_size = options.frag_size;
    int NAA = options.NAA;
    int len = seq->size;
    int len_bound = upper_bound_length_rep(len, options);
    param.len_upper_bound = len_bound;
    int flag = CheckOne(seq, local_table, param, buffer, options);

    if (flag == 0) {
        if ((seq->identity > 0) && (options.cluster_best)) {
            // because of the -g option, this seq is similar to seqs in old SEGs
            seq->state |= IS_REDUNDANT;
            seq->Clear();
        } else { // else add to NR90 db
            int aan_no = len - NAA + 1;
            int size = rep_seqs.size();
            rep_seqs.Append(id);
            seq->cluster_id = size;
            seq->identity = 0;
            seq->state |= IS_REP;
            if (frag_size) { /* not used for EST */
                int frg1 = (len - NAA) / frag_size + 1;
                local_table.AddWordCountsFrag(aan_no, buffer.word_encodes_backup, buffer.word_encodes_no, frg1,
                                              frag_size);
                table.AddWordCountsFrag(aan_no, buffer.word_encodes_backup, buffer.word_encodes_no, frg1, frag_size);
            } else {
                local_table.AddWordCounts(aan_no, buffer.word_encodes, buffer.word_encodes_no,
                                          local_table.sequences.size(), options.isEST);
                table.AddWordCounts(aan_no, buffer.word_encodes, buffer.word_encodes_no, table.sequences.size(),
                                    options.isEST);
            }
            table.sequences.Append(seq);
            local_table.sequences.Append(seq);
            if (frag_size) {
                while (local_table.sequences.size() < local_table.frag_count) local_table.sequences.Append(seq);
                while (table.sequences.size() < table.frag_count) table.sequences.Append(seq);
            }
        }
    }
    if ((id + 1) % 1000 == 0) {
        int size = rep_seqs.size();
        printf(".");
        fflush(stdout);
        if ((id + 1) % 10000 == 0) printf("\r..........%9i  finished  %9i  clusters\n", id + 1, size);
    }
}
void SequenceDB::ClusterOne_single(Sequence *seq, int id, WordTable &word_table, WorkingParam &param,
                                   WorkingBuffer &buffer, const Options &options, int &centers) {
    int len = seq->size;
    int NAA = options.NAA;
    int len_bound = upper_bound_length_rep(len, options);
    param.len_upper_bound = len_bound;
    int flag = CheckOne_single(seq, id, word_table, param, buffer, options);
    if (flag == 1) {
        seq->state |= IS_REDUNDANT;
    } else {
        // buffer.EncodeWords( seq, options.NAA, false );
        int size = rep_seqs.size();
        rep_seqs.Append(seq->index);
        seq->cluster_id = size;
        seq->identity = 0;
        seq->state |= IS_REP;
        seq->table_idx = centers;
        centers++;
        int aan_no = len - NAA + 1;
        for (int j = 0; j < aan_no; ++j) {
            int bucket = buffer.word_encodes[j];
            int count = buffer.word_encodes_no[j];

            if (count > 0) {
                NVector<IndexCount> &row = word_table.indexCounts[bucket];
                row.Append(IndexCount(id, count));
            }
        }
    }
}
#include <assert.h>
size_t SequenceDB::MinimalMemory(int frag_no, int bsize, int T, const Options &options, size_t extra) {
    int N = total_num;
    int F = frag_no < MAX_TABLE_SEQ ? frag_no : MAX_TABLE_SEQ;
    size_t mem_need = 0;
    size_t mem, mega = 1000000;
    int table = T > 1 ? 2 : 1;

    printf("\nApproximated minimal memory consumption:\n");
    mem = N * sizeof(Sequence) + total_desc + N + extra;
    if (options.store_disk == false) mem += total_letter + N;
    printf("%-16s: %zuM\n", "Sequence", mem / mega);
    mem_need += mem;

    mem = bsize;
    printf("%-16s: %i X %zuM = %zuM\n", "Buffer", T, mem / mega, T * mem / mega);
    mem_need += T * mem;

    mem = F * (sizeof(Sequence *) + sizeof(IndexCount)) + NAAN * sizeof(NVector<IndexCount>);
    printf("%-16s: %i X %zuM = %zuM\n", "Table", table, mem / mega, table * mem / mega);
    mem_need += table * mem;

    mem = sequences.capacity() * sizeof(Sequence *) + N * sizeof(int);
    mem += Comp_AAN_idx.size() * sizeof(int);
    printf("%-16s: %zuM\n", "Miscellaneous", mem / mega);
    mem_need += mem;

    printf("%-16s: %zuM\n\n", "Total", mem_need / mega);

    if (options.max_memory and options.max_memory < mem_need + 50 * table) {
        char msg[200];
        sprintf(msg, "not enough memory, please set -M option greater than %zu\n", 50 * table + mem_need / mega);
        bomb_error(msg);
    }
    return mem_need;
}
size_t MemoryLimit(size_t mem_need, const Options &options) {
    size_t mem_limit = (options.max_memory - mem_need) / sizeof(IndexCount);

    // printf( "Table limit with the given memory limit:\n" );
    if (options.max_memory == 0) {
        mem_limit = options.max_entries;
        if (mem_limit > MAX_TABLE_SIZE) mem_limit = MAX_TABLE_SIZE;
    }
    // printf( "Max number of representatives: %zu\n", mem_limit );
    // printf( "Max number of word counting entries: %zu\n\n", mem_limit );
    return mem_limit;
}
void Options::ComputeTableLimits(int min_len, int max_len, int typical_len, size_t mem_need) {
    // liwz Fri Jan 15 15:44:47 PST 2016
    // T=1 scale=1
    // T=2 scale=0.6035
    // T=4 scale=0.375
    // T=8 scale=0.2392
    // T=16 scale=0.1562
    // T=32 scale=0.104
    // T=64 scale=0.0703

    double scale = 0.5 / threads + 0.5 / sqrt(threads);
    max_sequences = (size_t) (scale * MAX_TABLE_SEQ);
    max_entries = (size_t) (scale * (500 * max_len + 500000 * typical_len + 50000000));
    if (max_memory) {
        double frac = max_sequences / (double) max_entries;
        max_entries = (options.max_memory - mem_need) / sizeof(IndexCount);
        max_sequences = (size_t) (max_entries * frac);
        if (max_sequences < MAX_TABLE_SEQ / 100) max_sequences = MAX_TABLE_SEQ / 100;
        if (max_sequences > MAX_TABLE_SEQ) max_sequences = MAX_TABLE_SEQ;
    }
    printf("Table limit with the given memory limit:\n");
    printf("Max number of representatives: %zu\n", max_sequences);
    printf("Max number of word counting entries: %zu\n\n", max_entries);
}

void post_ibcasts_for_this_block(Slot &s, int source, MPI_Comm comm) {
    MPI_Bcast((void *) s.info, 7, MPI_LONG, source, MPI_COMM_WORLD);
    s.cluster_n = (size_t) s.info[1];
    s.suffix_n = (size_t) s.info[1];
    s.cluster_id = (long *) realloc(s.cluster_id, sizeof(long) * s.cluster_n);
    s.seqs_suffix = (long *) realloc(s.seqs_suffix, sizeof(long) * s.suffix_n);
    MPI_Bcast((void *) s.cluster_id, (int) s.cluster_n, MPI_LONG, source, MPI_COMM_WORLD);

    MPI_Bcast((void *) s.seqs_suffix, (int) s.suffix_n, MPI_LONG, source, MPI_COMM_WORLD);
}
void post_ibcasts_for_next_block(Slot &s, int source, MPI_Comm comm) {
    MPI_Request head;
    MPI_Ibcast((void *) s.info, 7, MPI_LONG, source, comm, &head);
    MPI_Wait(&head, MPI_STATUS_IGNORE);

    s.cluster_n = (size_t) s.info[1];
    s.suffix_n = (size_t) s.info[1];
    s.reqs.clear();
    s.cluster_id = (long *) realloc(s.cluster_id, sizeof(long) * s.cluster_n);
    s.seqs_suffix = (long *) realloc(s.seqs_suffix, sizeof(long) * s.suffix_n);

    s.reqs.emplace_back(MPI_REQUEST_NULL);
    MPI_Ibcast((void *) s.cluster_id, (int) s.cluster_n, MPI_LONG, source, comm, &s.reqs.back());
    s.reqs.emplace_back(MPI_REQUEST_NULL);
    MPI_Ibcast((void *) s.seqs_suffix, (int) s.suffix_n, MPI_LONG, source, comm, &s.reqs.back());
    // double t1 = get_time();
    // MPI_Waitall(s.reqs.size(), s.reqs.data(), MPI_STATUSES_IGNORE);
    // double t2 = get_time();
    // cerr << "..................wait  time " << t2 - t1 << endl;
}

void wait_all(Slot &s) {
    if (!s.reqs.empty()) MPI_Waitall((int) s.reqs.size(), s.reqs.data(), MPI_STATUSES_IGNORE);
    s.reqs.clear();
}

void SequenceDB::encode_WordTable(long *&info_buf, int chunk_id, int start, int end, long *&cluster_id_buf,
                                  long *&suffix_buf, long *&indexCount_buf, long long *&prefix_buf,
                                  long long &indexCount_buf_size, long &prefix_size, int send_file_index,
                                  int start_global_id) {
    int T = options.threads;
    int len = end - start;

    info_buf = (long *) malloc(7 * sizeof(long));
    info_buf[0] = chunk_id;
    info_buf[1] = len;

    suffix_buf = (long *) malloc(len * sizeof(long));
    cluster_id_buf = (long *) malloc(len * sizeof(long));
    for (int i = start; i < end; i++) {
        int index = i - start;
        suffix_buf[index] = rep_seqs[i];
        cluster_id_buf[index] = sequences[rep_seqs[i] - start_global_id]->cluster_id;
    }
    info_buf[5] = send_file_index;
    info_buf[6] = start_global_id;

    // 计算每行 size
}
void SequenceDB::prepare_to_decode(WordTable &table, long *&info_buf, long *&cluster_id_buf, long *&suffix_buf,
                                   long *&indexCount_buf, long long *&prefix_buf, long long &indexCount_buf_size) {
    int len = info_buf[1];
    cluster_id_buf = (long *) malloc(len * sizeof(long));
    suffix_buf = (long *) malloc(len * sizeof(long));
}
void SequenceDB::decode_WordTable(WordTable &table, int start, Slot &s, const Options &options)
{
    int T = options.threads;
    int len = s.info[1];
    int start_id = s.info[6];
    // int my_rank;
    // MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    table.sequences.resize(len);
    if (chunks_id[start] == s.info[0])
    {
        int table_start_id = my_chunks[start].first;
#pragma omp parallel for num_threads(T)
        for (int i = 0; i < len; i++)
        {
            int index = i;

            Sequence *seq = rep_sequences[s.seqs_suffix[index] - start_id];
            Sequence *seq1 = sequences[s.seqs_suffix[index] - start_id + table_start_id];
            if (options.stealing)
            {
                auto &m = meta_[s.seqs_suffix[index] - start_id + table_start_id];
                m.cluster_id = s.cluster_id[index];
                m.state = IS_REP;
                m.identity = 0;
            }
            else
            {
                seq1->cluster_id = s.cluster_id[index];
                seq1->state |= IS_REP;
                seq1->identity = 0;
            }
            seq->cluster_id = s.cluster_id[index];
            seq->identity = 0;
            seq->table_idx = i;
            seq->state |= IS_REP;

            seq1->Clear();

            table.sequences[i] = seq;
        }
    }
    else
    {
#pragma omp parallel for num_threads(T)
        for (int i = 0; i < len; i++) {
            int index = i;
            Sequence *seq = rep_sequences[s.seqs_suffix[index] - start_id];
            seq->cluster_id = s.cluster_id[index];
            seq->identity = 0;
            seq->table_idx = i;
            seq->state |= IS_REP;
            table.sequences[i] = seq;
        }
    }
}

char *SequenceDB::str_copy(const char *str) {
    if (!str) return nullptr;
    char *data = new char[strlen(str) + 1];
    strcpy(data, str);
    return data;
}
void mpi_progress_thread() {
    while (progress_running) {
        int dummy;
        MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &dummy, MPI_STATUS_IGNORE);
        std::this_thread::sleep_for(std::chrono::microseconds(50));
    }
}
void SequenceDB::send_cluster(const std::vector<std::vector<std::string>> &clusters_identifier,
                              const std::vector<std::vector<int>> &clusters_size,
                              const std::vector<std::vector<float>> &clusters_identity,
                              const std::vector<std::vector<int>> &clusters_coverage, int *&prefix_seq, int *&flat_size,
                              float *&flat_identity, int *&flat_coverage, char *&flat_identifier, int &C, int &N,
                              int &IDLEN) {
    int *seq_cnt = (int *) malloc(C * sizeof(int));
    prefix_seq = (int *) malloc((C + 1) * sizeof(int));
    int T = options.threads;
#pragma omp parallel for num_threads(T)
    for (int c = 0; c < C; ++c) seq_cnt[c] = clusters_identifier[c].size();

    prefix_seq[0] = 0;
    for (int c = 0; c < C; ++c) prefix_seq[c + 1] = prefix_seq[c] + seq_cnt[c];

    N = prefix_seq[C];

    // 分配定长 flat_identifier
    flat_identifier = (char *) malloc(N * IDLEN);
    flat_size = (int *) malloc(N * sizeof(int));
    flat_identity = (float *) malloc(N * sizeof(float));
    flat_coverage = (int *) malloc(N * 4 * sizeof(int));

#pragma omp parallel for num_threads(T)
    for (int c = 0; c < C; ++c) {
        int base = prefix_seq[c];
        for (int j = 0; j < seq_cnt[c]; ++j) {
            int g = base + j;

            // identifier 拷贝 + 补零
            const std::string &id = clusters_identifier[c][j];
            char *dst = flat_identifier + g * IDLEN;
            memcpy(dst, id.data(), id.size());
            memset(dst + id.size(), 0, IDLEN - id.size());

            flat_size[g] = clusters_size[c][j];
            flat_identity[g] = clusters_identity[c][j];

            for (int k = 0; k < 4; ++k) flat_coverage[g * 4 + k] = clusters_coverage[c][j * 4 + k];
        }
    }

    free(seq_cnt);
    seq_cnt = nullptr;
}
void SequenceDB::ClusterOne_worker(Sequence *seq, int id, WordTable &table, WorkingParam &param, WorkingBuffer &buffer,
                                   const Options &options, PaddedLock *locks, int num_locks, int lock_mask) {
    int len = seq->size;
    int NAA = options.NAA;
    buffer.EncodeWords(seq, options.NAA, false);
    int aan_no = len - NAA + 1;
    for (int j = 0; j < aan_no; ++j) {
        int bucket = buffer.word_encodes[j];
        int count = buffer.word_encodes_no[j];
        if (count > 0) {
            NVector<IndexCount> &row = table.indexCounts[bucket];
            int lock_idx = bucket & lock_mask;
            omp_set_lock(&locks[lock_idx].lock);
            row.Append(IndexCount(id, count));
            omp_unset_lock(&locks[lock_idx].lock);
        }
    }
}
void SequenceDB::ClusterOne_master(Sequence *seq, int id, std::vector<std::vector<std::pair<int, int>>> &word_table,
                                   WorkingParam &param, WorkingBuffer &buffer, const Options &options,
                                   PaddedLock *locks, int num_locks, int lock_mask) {
    int len = seq->size;
    int NAA = options.NAA;
    buffer.EncodeWords(seq, options.NAA, false);
    int aan_no = len - NAA + 1;
    for (int j = 0; j < aan_no; ++j) {
        int bucket = buffer.word_encodes[j];
        int count = buffer.word_encodes_no[j];
        if (count > 0) {
            int lock_idx = bucket & lock_mask;
            omp_set_lock(&locks[lock_idx].lock);
            word_table[bucket].emplace_back(id, count);
            omp_unset_lock(&locks[lock_idx].lock);
        }
    }
}
int SequenceDB::CheckOne_single(Sequence *seq, int qid, WordTable &word_table, WorkingParam &param, WorkingBuffer &buf,
                                const Options &options) {
    int len = seq->size;
    param.len_upper_bound = upper_bound_length_rep(len, options);
    // if( options.isEST ) return CheckOneEST( seq, word_table, param, buf, options );
    return CheckOneAA_single(seq, qid, word_table, param, buf, options);
}
int SequenceDB::CheckOneAA_single(Sequence *seq, int qid, WordTable &word_table, WorkingParam &param,
                                  WorkingBuffer &buf,
                                  const Options &options) { // Todo  uint32_t
    NVector<IndexCount> &lookCounts = buf.lookCounts;

    NVector<uint32_t> &indexMapping = buf.indexMapping;
    Vector<INTs> &word_encodes_no = buf.word_encodes_no;
    Vector<INTs> &aap_list = buf.aap_list;
    Vector<INTs> &aap_begin = buf.aap_begin;
    Vector<int> &word_encodes = buf.word_encodes;
    Vector<int> &taap = buf.taap;
    double aa1_cutoff = param.aa1_cutoff;
    double aa2_cutoff = param.aas_cutoff;
    double aan_cutoff = param.aan_cutoff;

    char *seqi = seq->data;
    int j, k, j1, len = seq->size;
    int flag = 0;
    int frag_size = options.frag_size;
    int &aln_cover_flag = param.aln_cover_flag;
    int &required_aa1 = param.required_aa1;
    int &required_aa2 = param.required_aas;
    int &required_aan = param.required_aan;
    int &min_aln_lenS = param.min_aln_lenS;
    int &min_aln_lenL = param.min_aln_lenL;

    int NAA = options.NAA;
    int len_eff = len;
    param.ControlShortCoverage(len_eff, options);
    param.ComputeRequiredBases(options.NAA, 2, options);

    buf.EncodeWords(seq, options.NAA, false);

    // if minimal alignment length > len, return
    // I can not return earlier, because I need to calc the word_encodes etc
    if (options.min_control > len) return 0; // return flag=0

    // lookup_aan
    int aan_no = len - options.NAA + 1;
    // int M = frag_size ? table.frag_count : S;
    word_table.CountWords(aan_no, qid, word_encodes, word_encodes_no, lookCounts, indexMapping, false, required_aan);

    // contained_in_old_lib()
    int len_upper_bound = param.len_upper_bound;
    int len_lower_bound = param.len_lower_bound;
    int band_left, band_right, best_score, band_width1, best_sum, len2, alnln, len_eff1;
    int tiden_no, band_center;
    float tiden_pc, distance = 0;
    int talign_info[5];
    int best1, sum;
    INTs *lookptr;
    char *seqj;
    int frg2 = frag_size ? (len - NAA + options.band_width) / frag_size + 1 + 1 : 0;
    int lens;
    int has_aa2 = 0;
    IndexCount *ic = lookCounts.items;
    ic = lookCounts.items;
    for (; ic->count; ic++) {
        if (!frag_size) {
            // indexMapping[ ic->index ] = 0;
            if (ic->count < required_aan) continue;
        }

        // Sequence *rep = table.sequences[ ic->index ];
        Sequence *rep = sequences[ic->index];

        len2 = rep->size;
        if (len2 > len_upper_bound) continue;
        if (options.has2D && len2 < len_lower_bound) continue;

        param.ControlLongCoverage(len2, options);

        if (has_aa2 == 0) { // calculate AAP array
            buf.ComputeAAP(seqi, seq->size);
            has_aa2 = 1;
        }
        seqj = rep->data; // NR_seq[NR90_idx[j]];

        band_width1 = (options.band_width < len + len2 - 2) ? options.band_width : len + len2 - 2;
        diag_test_aapn(NAA1, seqj, len, len2, buf, best_sum, band_width1, band_left, band_center, band_right,
                       required_aa1);
        if (best_sum < required_aa2) continue;

        int rc = FAILED_FUNC;
#ifndef NO_AVX512

        if (options.print || aln_cover_flag) // return overlap region
            rc = rotation_band_align_AVX512(seqi, seqj, len, len2, mat, best_score, tiden_no, alnln, distance,
                                            talign_info, band_left, band_center, band_right, buf);
        else
            rc = rotation_band_align_AVX512(seqi, seqj, len, len2, mat, best_score, tiden_no, alnln, distance,
                                            talign_info, band_left, band_center, band_right, buf);
#else
        // auto t0 = std::chrono::high_resolution_clock::now();
        if (options.print || aln_cover_flag) // return overlap region
            rc = local_band_align(seqi, seqj, len, len2, mat, best_score, tiden_no, alnln, distance, talign_info,
                                  band_left, band_center, band_right, buf);
        else
            rc = local_band_align(seqi, seqj, len, len2, mat, best_score, tiden_no, alnln, distance, talign_info,
                                  band_left, band_center, band_right, buf);
#endif
        // auto t1 = std::chrono::high_resolution_clock::now();
        // param.local_align_total_time_ns += std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count();
        // param.local_align_call_count++;
        if (rc == FAILED_FUNC) continue;
        if (tiden_no < required_aa1) continue;
        lens = len;
        if (options.has2D && len > len2) lens = len2;
        len_eff1 = (options.global_identity == 0) ? alnln : (lens - talign_info[4]);
        tiden_pc = tiden_no / (float) len_eff1;
        if (options.useDistance) {
            if (distance > options.distance_thd) continue;
            if (distance >= seq->distance) continue; // existing distance
        } else {
            if (tiden_pc < options.cluster_thd) continue;
            if (tiden_pc <= seq->identity) continue; // existing iden_no
        }
        if (aln_cover_flag) {
            if (talign_info[3] - talign_info[2] + 1 < min_aln_lenL) continue;
            if (talign_info[1] - talign_info[0] + 1 < min_aln_lenS) continue;
        }
        flag = 1;
        break;
    }
    return flag;
}
int SequenceDB::CheckOne_master(Sequence *seq, int qid, WordTable &table, WorkingParam &param, WorkingBuffer &buf,
                                const Options &options) {
    int len = seq->size;
    param.len_upper_bound = upper_bound_length_rep(len, options);
    // if( options.isEST ) return CheckOneEST( seq, word_table, param, buf, options );
    return CheckOneAA_master(seq, qid, table, param, buf, options);
}
int SequenceDB::CheckOneAA_master(Sequence *seq, int qid, WordTable &table, WorkingParam &param, WorkingBuffer &buf,
                                  const Options &options) { // Todo  uint32_t
    NVector<IndexCount> &lookCounts = buf.lookCounts;

    NVector<uint32_t> &indexMapping = buf.indexMapping;
    Vector<INTs> &word_encodes_no = buf.word_encodes_no;
    Vector<INTs> &aap_list = buf.aap_list;
    Vector<INTs> &aap_begin = buf.aap_begin;
    Vector<int> &word_encodes = buf.word_encodes;
    Vector<int> &taap = buf.taap;
    double aa1_cutoff = param.aa1_cutoff;
    double aa2_cutoff = param.aas_cutoff;
    double aan_cutoff = param.aan_cutoff;

    char *seqi = seq->data;
    int j, k, j1, len = seq->size;
    int flag = 0;
    int frag_size = options.frag_size;
    int &aln_cover_flag = param.aln_cover_flag;
    int &required_aa1 = param.required_aa1;
    int &required_aa2 = param.required_aas;
    int &required_aan = param.required_aan;
    int &min_aln_lenS = param.min_aln_lenS;
    int &min_aln_lenL = param.min_aln_lenL;

    int NAA = options.NAA;
    int len_eff = len;
    param.ControlShortCoverage(len_eff, options);
    param.ComputeRequiredBases(options.NAA, 2, options);

    buf.EncodeWords(seq, options.NAA, false);

    // if minimal alignment length > len, return
    // I can not return earlier, because I need to calc the word_encodes etc
    if (options.min_control > len) return 0; // return flag=0

    // lookup_aan
    int aan_no = len - options.NAA + 1;
    // int M = frag_size ? table.frag_count : S;
    // buf.CountWords(aan_no,  qid,word_table,false, required_aan);
    table.CountWords(aan_no, qid, word_encodes, word_encodes_no, lookCounts, indexMapping, false, required_aan);
    // contained_in_old_lib()
    int len_upper_bound = param.len_upper_bound;
    int len_lower_bound = param.len_lower_bound;
    int band_left, band_right, best_score, band_width1, best_sum, len2, alnln, len_eff1;
    int tiden_no, band_center;
    float tiden_pc, distance = 0;
    int talign_info[5];
    int best1, sum;
    INTs *lookptr;
    char *seqj;
    int frg2 = frag_size ? (len - NAA + options.band_width) / frag_size + 1 + 1 : 0;
    int lens;
    int has_aa2 = 0;
    IndexCount *ic = lookCounts.items;
    ic = lookCounts.items;
    for (; ic->count; ic++) {
        if (!frag_size) {
            indexMapping[ic->index] = 0;
            if (ic->count < required_aan) continue;
        }
        // cerr<<"ic->index "<<ic->index<<endl;
        // Sequence *rep = table.sequences[ ic->index ];
        Sequence *rep = sequences[ic->index];
        if (rep->state & IS_REDUNDANT)continue;
        len2 = rep->size;
        if (len2 > len_upper_bound) continue;
        if (options.has2D && len2 < len_lower_bound) continue;

        param.ControlLongCoverage(len2, options);

        if (has_aa2 == 0) { // calculate AAP array
            buf.ComputeAAP(seqi, seq->size);
            has_aa2 = 1;
        }
        seqj = rep->data; // NR_seq[NR90_idx[j]];

        band_width1 = (options.band_width < len + len2 - 2) ? options.band_width : len + len2 - 2;
        diag_test_aapn(NAA1, seqj, len, len2, buf, best_sum, band_width1, band_left, band_center, band_right,
                       required_aa1);
        if (best_sum < required_aa2) continue;

        int rc = FAILED_FUNC;
#ifndef NO_AVX512
        if (options.print || aln_cover_flag) // return overlap region
            rc = rotation_band_align_AVX512(seqi, seqj, len, len2, mat, best_score, tiden_no, alnln, distance,
                                            talign_info, band_left, band_center, band_right, buf);
        else
            rc = rotation_band_align_AVX512(seqi, seqj, len, len2, mat, best_score, tiden_no, alnln, distance,
                                            talign_info, band_left, band_center, band_right, buf);
#else
        // auto t0 = std::chrono::high_resolution_clock::now();
        if (options.print || aln_cover_flag) // return overlap region
            rc = local_band_align(seqi, seqj, len, len2, mat, best_score, tiden_no, alnln, distance, talign_info,
                                  band_left, band_center, band_right, buf);
        else
            rc = local_band_align(seqi, seqj, len, len2, mat, best_score, tiden_no, alnln, distance, talign_info,
                                  band_left, band_center, band_right, buf);
#endif
        // auto t1 = std::chrono::high_resolution_clock::now();
        // param.local_align_total_time_ns += std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count();
        // param.local_align_call_count++;
        if (rc == FAILED_FUNC) continue;
        if (tiden_no < required_aa1) continue;
        lens = len;
        if (options.has2D && len > len2) lens = len2;
        len_eff1 = (options.global_identity == 0) ? alnln : (lens - talign_info[4]);
        tiden_pc = tiden_no / (float) len_eff1;
        if (options.useDistance) {
            if (distance > options.distance_thd) continue;
            if (distance >= seq->distance) continue; // existing distance
        } else {
            if (tiden_pc < options.cluster_thd) continue;
            if (tiden_pc <= seq->identity) continue; // existing iden_no
        }
        if (aln_cover_flag) {
            if (talign_info[3] - talign_info[2] + 1 < min_aln_lenL) continue;
            if (talign_info[1] - talign_info[0] + 1 < min_aln_lenS) continue;
        }
        flag = 1;
        if (rep->state & IS_REP) {
            seq->state |= IS_REDUNDANT;
            break;
        }
        buf.thread_edges.emplace_back(qid, ic->index);
    }
    return flag;
}
// int WorkingBuffer::CountWords(int aan_no, int qid,const std::vector<std::vector<std::pair<int,int>>>& word_table
// ,bool est, int min)
// {
// 	out_pairs.clear();
// 	visited.clear();
// 	auto add_count = [&](int tid, int add){
// 		if (counts[tid] == 0) visited.push_back(tid);
// 		counts[tid] += add;
// 	};
// 	for (int j0 = 0; j0 < aan_no; ++j0) {
// 		int bucket = word_encodes[j0];
// 		int qcnt   = word_encodes_no[j0];
// 		if (qcnt == 0) continue;

// 		const auto& hits = word_table[bucket];  // 已按 seq_id 降序
// 		int rest = aan_no - j0 + 1;

// 		// 只扫 id>qid 的前缀
// 		for (size_t k = 0; k < hits.size(); ++k) {
// 			int tid  = hits[k].first;
// 			if (tid >= qid) break;  // 及早终止
// 			int tcnt = hits[k].second;

// 			if (min > 0 && rest < min && counts[tid] == 0) continue;
// 			int add = (qcnt < tcnt) ? qcnt : tcnt;
// 			add_count(tid, add);
// 		}
// 	}
// 	out_pairs.reserve(visited.size());
// 	for (int tid : visited) {
// 		out_pairs.emplace_back(tid, counts[tid]);
// 		counts[tid] = 0;
// 	}
// 	return OK_FUNC;
// }

//--------------------------------------------------------
// Achieve efficient information output. -- By Guiliang Ma
// Used for fast conversion of uint32_t to char* / string
static inline char *fast_uint_to_str(char *buf, uint32_t val) {
    if (val == 0) {
        *buf++ = '0';
        return buf;
    }
    char tmp[12];
    int len = 0;
    while (val > 0) {
        tmp[len++] = val % 10 + '0';
        val /= 10;
    }
    for (int i = len - 1; i >= 0; i--) {
        *buf++ = tmp[i];
    }
    return buf;
}

// Used for fast conversion of float to char* / string with 2 decimal places
static inline char *fast_float_to_str_2dec(char *buf, float val)
{
    // 1. printf("%.2f") 会把 float 提升为 double
    double x = (double)val;

    // 2. 放大 100 倍，使用 nearbyint 做银行家舍入（ties-to-even）
    double scaled  = x * 100.0;
    long long iv   = (long long) nearbyint(scaled);  

    // 3. 拆分整数和两位小数
    unsigned int integer = (unsigned int)(iv / 100);
    unsigned int decimal = (unsigned int)(iv % 100);

    // 4. 输出整数部分
    buf = fast_uint_to_str(buf, integer);

    // 5. 输出 ".xx"
    *buf++ = '.';
    *buf++ = '0' + (decimal / 10);
    *buf++ = '0' + (decimal % 10);

    // 如你不需要百分号，就不要加
    *buf++ = '%';

    return buf;
}


// Used for fast counting the number of digits of a uint32_t
static inline int count_digits(uint32_t val) {
    if (val == 0) return 1;
    int count = 0;
    while (val > 0) {
        val /= 10;
        count++;
    }
    return count;
}

// Maintain partition information describing the data range each thread needs to process.
struct PartitionInfo {
    size_t start_idx;
    size_t end_idx;
    size_t file_offset;
    size_t estimated_size;

    // 构造函数
    PartitionInfo() : start_idx(0), end_idx(0), file_offset(0), estimated_size(0) {}
};

// Used for finding the boundaries of clusters in the sorted result.
std::vector<size_t> FindClusterBoundaries(const std::vector<MasterSeqInfo> &all_infos, int num_partitions) {
    std::vector<size_t> boundaries;
    boundaries.reserve(num_partitions + 1);
    boundaries.push_back(0);
    if (all_infos.empty() || num_partitions <= 0) {
        boundaries.push_back(all_infos.size());
        return boundaries;
    }
    size_t total_size = all_infos.size();
    size_t target_chunk_size = total_size / num_partitions;
    for (int p = 1; p < num_partitions; p++) {
        size_t ideal_pos = p * target_chunk_size;

        // Search boundary of ideal pos + half of target chunk size
        size_t search_start = ideal_pos;
        size_t search_end = min(ideal_pos + target_chunk_size / 2, total_size);

        bool found = false;
        for (size_t i = search_start; i < search_end; ++i) {
            if (all_infos[i].cluster_id != all_infos[i - 1].cluster_id) {
                boundaries.push_back(i);
                found = true;
                break;
            }
        }
        if (!found) {
            boundaries.push_back(search_end);
        }
    }
    boundaries.push_back(total_size);
    return boundaries;
}

// Accurately calculate the actual size of each block
std::vector<PartitionInfo> CalculatePartitionInfo(const std::vector<MasterSeqInfo> &all_infos,
                                                  const std::vector<size_t> &boundaries, const int num_threads) {
    omp_set_num_threads(num_threads);
    size_t num_partitions = boundaries.size() - 1;
    // cout << ">>>> num_partitions: " << num_partitions << endl;
    std::vector<PartitionInfo> partitions(num_partitions);

#pragma omp parallel for schedule(static)
    for (size_t p = 0; p < num_partitions; ++p) {
        PartitionInfo &info = partitions[p];
        info.start_idx = boundaries[p];
        info.end_idx = boundaries[p + 1];

        bool first_in_partition = true;
        size_t partition_size = 0;
        uint32_t current_cluster = all_infos[info.start_idx].cluster_id;
        int seq_idx = 0;

        for (size_t i = info.start_idx; i < info.end_idx; ++i) {
            const MasterSeqInfo &m = all_infos[i];

            if (m.cluster_id != current_cluster || first_in_partition) {
                // if (!first_in_partition) {
                //     seq_index = 0;
                // }
                // >Cluster [cluster_id]\n
                // 9:">Cluster "
                // count_digits(m.cluster_id): cluster_id
                // 1: "\n"
                first_in_partition = false;
                partition_size += 9 + count_digits(m.cluster_id) + 1;
                current_cluster = m.cluster_id;
                seq_idx = 0;
            }

            // [seq_idx]\t[size]aa, >[name]... at [identity]%\n
            // [seq_idx]\t[size]aa, >[name]... *\n

            // count_digits(seq_idx): seq_idx
            // 1: "\t"
            // count_digits(m.size): size
            // 5: "aa, >"
            // name.length(): name
            partition_size += count_digits(seq_idx) + 1;
            partition_size += count_digits(m.size) + 5;
            partition_size += m.name.length();


            if (m.identity == 0) {
                // 6:"... *\n"
                partition_size += 6;
            } else {
                // "... at [identity]%\n"
                float percent = m.identity * 100.0f;
                if(options.print){
                    // cmp_fout << m.coverage[0] << ":" << ccoverage[1] << ":" << ccoverage[2] << ":" << ccoverage[3] << "/";
                    partition_size += count_digits(m.coverage[0])+count_digits(m.coverage[1])+count_digits(m.coverage[2])+count_digits(m.coverage[3])+4;
                }
                // 7:"... at "
                partition_size += 7;
                // 7 or 6:
                //   7:"100.00%"
                //   6:"xx.xx%"
                partition_size += (percent >= 99.995f) ? 7 : 6;
                // 1:"\n"
                partition_size += 1;
            }
            seq_idx++;
        }

        info.estimated_size = partition_size;
    }

    // char *buff;
    // fast_float_to_str_2dec(buff, 100.00f);
    // exit(0);

    std::vector<size_t> sizes(num_partitions);
    for (size_t p = 0; p < num_partitions; ++p) {
        sizes[p] = partitions[p].estimated_size;
    }

    std::vector<size_t> offsets(num_partitions);
    std::partial_sum(sizes.begin(), sizes.end() - 1, offsets.begin() + 1);
    offsets[0] = 0;

    for (size_t p = 0; p < num_partitions; ++p) {
        partitions[p].file_offset = offsets[p];
    }

    return partitions;
}

// Write a partition to a memory-mapped file
void WritePartitionToMmap(const std::vector<MasterSeqInfo> &all_infos, const PartitionInfo &partition,
                          char *mapped_base) {
    char *ptr = mapped_base + partition.file_offset;
    uint32_t current_cluster = all_infos[partition.start_idx].cluster_id;
    int seq_index = 0;
    bool first_in_partition = true;
    for (size_t i = partition.start_idx; i < partition.end_idx; i++) {
        const MasterSeqInfo &m = all_infos[i];
        if (m.cluster_id != current_cluster || first_in_partition) {
            if (!first_in_partition) {
                seq_index = 0;
            }
            current_cluster = m.cluster_id;
            first_in_partition = false;
            memcpy(ptr, ">Cluster ", 9);
            ptr += 9;
            ptr = fast_uint_to_str(ptr, m.cluster_id);
            *ptr++ = '\n';
        }
        ptr = fast_uint_to_str(ptr, seq_index);
        *ptr++ = '\t';
        ptr = fast_uint_to_str(ptr, m.size);
        memcpy(ptr, "aa, >", 5);
        ptr += 5;
        size_t name_len = m.name.length();
        // if (i == 1670149) cout << ">>>> name: " << m.name << endl;
        memcpy(ptr, m.name.c_str(), name_len);
        ptr += name_len;
        if (m.identity == 0) {
            memcpy(ptr, "... *\n", 6);
            ptr += 6;
        } else {
            memcpy(ptr, "... at ", 7);
            ptr += 7;
            if(options.print){
                ptr = fast_uint_to_str(ptr,m.coverage[0]);
                *ptr++=':';
                ptr = fast_uint_to_str(ptr,m.coverage[1]);
                *ptr++=':';
                ptr = fast_uint_to_str(ptr,m.coverage[2]);
                *ptr++=':';
                ptr = fast_uint_to_str(ptr,m.coverage[3]);
                *ptr++='/';
            }
            ptr = fast_float_to_str_2dec(ptr, m.identity * 100.0f);
            *ptr++ = '\n';
        }
        seq_index++;
    }
}

void SequenceDB::SortAndWriteResult(vector<MasterSeqInfo> &all_infos, const Options &options) {
    tbb::global_control(tbb::global_control::max_allowed_parallelism, options.threads);
    int num_threads = options.threads;
    omp_set_num_threads(num_threads);
    double written_mb = 0.0;
    double written_time_naive = 0.0;
    double written_time_mmap = 0.0;

    // cout << "--------------------------------" << endl;
    // cout << "[1/6] Sorting and writing result..." << endl;
    // ips2ra::sort(all_infos.begin(), all_infos.end(), [](const MasterSeqInfo &info) { return info.cluster_id; });

    ips2ra::parallel::sort(all_infos.begin(), all_infos.end(),
                           [](const MasterSeqInfo &info) { return info.cluster_id; });
    // cout << "--------------------------------" << endl;

    // Naive methods
    auto start_time = std::chrono::high_resolution_clock::now();

    // 设置浮点数输出精度为2位小数
#ifdef OUTPUT_CMP
    string cmp_output = options.output + ".clstr.cmp";
    ofstream cmp_fout(cmp_output);

    cmp_fout << std::fixed << std::setprecision(2);

    if (all_infos.empty()) {
        cmp_fout.close();
        return;
    }

    uint32_t current_cluster = all_infos[0].cluster_id;
    int seq_index = 0;
    bool first_entry = true;

    for (size_t i = 0; i < all_infos.size(); ++i) {
        const MasterSeqInfo &m = all_infos[i];

        // 检测cluster变化，输出cluster标题
        if (m.cluster_id != current_cluster || first_entry) {
            if (!first_entry) {
                seq_index = 0; // 新cluster，重置索引
            }
            current_cluster = m.cluster_id;
            first_entry = false;

            cmp_fout << ">Cluster " << m.cluster_id << "\n";
        }

        // 输出序列信息行
        cmp_fout << seq_index << "\t" << m.size << "aa, >" << m.name << "...";

        if (m.identity == 0) {
            // 代表序列
            cmp_fout << " *\n";
        } else {
            // 其他序列，显示相似度
            if (options.print)
 			cmp_fout << m.coverage[0] << ":" << ccoverage[1] << ":" << ccoverage[2] << ":" << ccoverage[3] << "/";
            cmp_fout << " at " << (m.identity * 100.0f) << "%\n";
        }
        
        seq_index++;
    }

    cmp_fout.close();
#endif
    // auto end_time = std::chrono::high_resolution_clock::now();
    // auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    // cout << "Naive methods time: " << duration.count() << " ms" << endl;
    // written_time_naive = (double) duration.count() / 1000.0;
    // naive methods end

    // start_time = std::chrono::high_resolution_clock::now();

    // cout << "[2/6] Finding cluster boundaries..." << endl;
    if (all_infos.empty()) {
        return;
    }

    // cout << "--------------------------------" << endl;
    std::vector<size_t> boundaries = FindClusterBoundaries(all_infos, num_threads);
    // for (size_t i = 0; i < boundaries.size(); i++) {
    //     if (i != boundaries.size() - 1) {
    //         cout << boundaries[i] << " (";
    //         cout << all_infos[boundaries[i]].cluster_id << ", " << all_infos[boundaries[i + 1] - 1].cluster_id << ")
    //         "; cout << boundaries[i + 1] << "\n";
    //     }
    // }
    // cout << endl;
    // cout << "--------------------------------" << endl;

    // cout << "\tUsing [" << num_threads << "] threads, [" << boundaries.size() - 1 << "] partitions" << endl;

    // cout << "[3/6] Calculating partition sizes..." << endl;

    std::vector<PartitionInfo> partitions = CalculatePartitionInfo(all_infos, boundaries, num_threads);
    // for (const auto &p : partitions) {
    //     cout << ">>>> partition: " << p.start_idx << " " << p.end_idx << " " << all_infos[p.start_idx].cluster_id <<
    //     " "
    //          << all_infos[p.end_idx - 1].cluster_id << " " << p.estimated_size << endl;
    // }
    size_t total_size = 0;
    for (const auto &p : partitions) {
        total_size += p.estimated_size;
    }
    // double file_size = (double) total_size / (1024.0 * 1024.0);
    // cout << ">>>> file_size: " << file_size << " MB" << endl;
    size_t safe_size = total_size + (total_size * 10 / 100);

    // cout << "[4/6] Creating memory-map..." << endl;
    string clstr_output = options.output + ".clstr";
    // string clstr_output = "/home/bigssd/mgl_data/huge.clstr";
    int fd = open(clstr_output.c_str(), O_RDWR | O_CREAT | O_TRUNC, 0644);
    if (fd == -1) {
        cerr << "ERROR: Failed to create mmap file: " << strerror(errno) << endl;
        return;
    }
    if (ftruncate(fd, safe_size) == -1) {
        cerr << "ERROR: Failed to allocate memory space: " << strerror(errno) << endl;
        close(fd);
        return;
    }

    char *mapped = (char *) mmap(nullptr, safe_size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
    if (mapped == MAP_FAILED) {
        cerr << "ERROR: mmap failed: " << strerror(errno) << endl;
        close(fd);
        return;
    }

    char *tmp = mapped;

    // cout << "[5/6] Parallel Writing ..." << endl;

    // auto write_start = std::chrono::high_resolution_clock::now();
#pragma omp parallel for schedule(dynamic)
    for (size_t p = 0; p < partitions.size(); p++) {
        WritePartitionToMmap(all_infos, partitions[p], mapped);
    }

    // tbb::parallel_for(tbb::blocked_range<size_t>(0, partitions.size()), [&](const tbb::blocked_range<size_t> &range)
    // {
    //     for (size_t p = range.begin(); p != range.end(); ++p) {
    //         WritePartitionToMmap(all_infos, partitions[p], mapped);
    //     }
    // });

    // auto write_end = std::chrono::high_resolution_clock::now();
    // auto write_duration = std::chrono::duration_cast<std::chrono::milliseconds>(write_end - write_start);
    // cout << "      Completed in " << write_duration.count() << " ms" << endl;

    // start_time = std::chrono::high_resolution_clock::now();

    // cout << "[6/6] Finalizing..." << endl;
    msync(mapped, total_size + total_size * 3 / 100, MS_SYNC);

    // end_time = std::chrono::high_resolution_clock::now();
    // duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

    // cout << "msync time: " << duration.count() << " ms" << endl;
    // written_mb = (double) total_size / (1024.0 * 1024.0);
    // written_time_mmap = (double) duration.count() / 1000.0;
    // cout << "Naive Written Speed: " << written_mb / written_time_naive << " MB/s" << endl;
    // cout << "Mmap Written Speed: " << written_mb / written_time_mmap << " MB/s" << endl;

    munmap(mapped, safe_size);
    ftruncate(fd, total_size);
    close(fd);

    cout << "[DONE] Cluster Information written to: " << clstr_output << endl;
    // clstr_fout.close();
    return;
}

void SequenceDB::DoClustering_MPI(const Options &options, int my_rank, bool master, bool worker, int worker_rank,
                                  const char *output, MPI_Comm worker_comm) {
    int rank_size;
    MPI_Comm_size(MPI_COMM_WORLD, &rank_size);
    int source = 0;
    int i, j, jj, k;
    int T = options.threads;
    int NAA = options.NAA;
    double aa1_cutoff = options.cluster_thd;
    double aas_cutoff = 1 - (1 - options.cluster_thd) * 4;
    double aan_cutoff = 1 - (1 - options.cluster_thd) * options.NAA;
    int seq_no = total_num;
    int frag_no = seq_no;
    int frag_size = options.frag_size;
    int len, len_bound;
    int flag;
    string preprocess_output_dir = options.preprocess_dir;
    if (!preprocess_output_dir.empty() && preprocess_output_dir.back() != '/' && preprocess_output_dir.back() != '\\') {
        preprocess_output_dir += '/';
    }
    if (not options.isEST)
        cal_aax_cutoff(aa1_cutoff, aas_cutoff, aan_cutoff, options.cluster_thd, options.tolerance,
                       naa_stat_start_percent, naa_stat, NAA);
    Vector<WorkingParam> params(T);
    Vector<WorkingBuffer> buffers(T);

    for (i = 0; i < T; i++) {
        params[i].Set(aa1_cutoff, aas_cutoff, aan_cutoff);
        buffers[i].Set(frag_no, max_len, NAAN, options);
    }

    WordTable word_table(options.NAA, NAAN);

    // size_t mem_need = MinimalMemory(frag_no, buffers[0].total_bytes, T, options);

    // size_t mem_limit = MemoryLimit(mem_need, options);
    // size_t mem, mega = 50000;

    // int remaining = 0;
    // int local_remaining = 0;
    // Options opts(options);
    // size_t max_items = opts.max_entries;
    // size_t max_seqs = opts.max_sequences;

    long *info_buf;
    long *cluster_id_buf;
    long *seqs_suffix_buf;
    long *indexCount_buf;
    long long *prefix_buf;
    long long indexCount_buf_size;
    long prefix_size;
    int record_last = 0;
    int record = 0;
    const int NUM_LOCKS = 262144;
    const int LOCK_MASK = 0x3FFFF;
    PaddedLock *locks = new PaddedLock[NUM_LOCKS];

    omp_set_num_threads(T);
#pragma omp parallel for schedule(static)
    for (int i = 0; i < NUM_LOCKS; ++i) {
        omp_init_lock(&locks[i].lock);
    }

    if (master) {
        cerr << "chunks_num  " << chunks_num << endl;
        if (rank_size <= 1) {
            cerr << "no workers found" << endl;
            exit(0);
        }

        std::vector<std::vector<int>> neigh;
        // vector<vector<pair<int, int>>> all_wordtable(NAAN);

        // string clstr_output = options.output + '/' + std::to_string(my_rank) + ".clstr";
        // ofstream clstr_fout(clstr_output);

        // Modify by MGL: remove the rank directory
        // string rep_output = options.output + '/' + std::to_string(my_rank);
        string rep_output = options.output+ ".txt";
        ofstream fout(rep_output);

        std::vector<std::vector<string>> id_tables(rank_size - 1);
        vector<int> read_flag(chunks_num, 0);
        int output_index = 0;
        int file_index = 0;
        vector<gzFile> chunk_fp(rank_size - 1, nullptr);
        vector<kseq_t *> chunk_kseq(rank_size - 1, nullptr);
        int last_rep_index = 0;
        int ibcast_flag = 1;
        int C = 0;
        int last_seq_num = 0;
        int chunk_id = 0;
        int start_id = -1;
        int end_id = -1;
        int len = 0;
        long long now_bytes = 0;

        int target_worker = 1;
        int start_global_id = sequences.size();
        bool read_stop = 0;
        int send_file_index = 0;
        kseq_t *seq = nullptr;
        Clear();
        Sequence one;
        std::string file = preprocess_output_dir + "_proc" + std::to_string(file_index) + ".fa";

        send_file_index = file_index;
        if (chunk_kseq[file_index] == nullptr) {
            gzFile fp = gzopen(file.c_str(), "r");
            if (!fp) {
                fprintf(stderr, "Cannot open file: %s\n", file.c_str());
                exit(1);
            }
            chunk_fp[file_index] = fp;
            seq = kseq_init(fp);

        } else
            seq = chunk_kseq[file_index];
        char *id_ptr = new char[max_idf + 1];
        char *data_ptr = new char[max_len + 1];
        while ((len = kseq_read(seq)) >= 0) {
            memcpy(id_ptr, seq->name.s, seq->name.l);
            id_ptr[seq->name.l] = 0;
            one.identifier = id_ptr;

            memcpy(data_ptr, seq->seq.s, seq->seq.l);
            data_ptr[seq->seq.l] = 0;
            one.data = data_ptr;

            one.size = len;
            one.index = sequences.size() + start_global_id;
            one.master_flag = 1;
            sequences.Append(new Sequence(one));
            id_tables[file_index].emplace_back(seq->name.s, seq->name.l);
            now_bytes += len;
            sequences[sequences.size() - 1]->ConvertBases();

            if (sequences.size() >= first_chunk_size) {
                chunks_id.push_back(chunk_id);
                chunk_id++;
                chunk_kseq[file_index] = seq;
                now_bytes = 0;
                break;
            }
        }
        one.identifier = nullptr;
        one.data = nullptr;
        delete[] id_ptr;
        delete[] data_ptr;
        file_index = (file_index + 1) % (rank_size - 1);
        // vector<vector<string>> clusters_identifier;
        // vector<vector<float>> clusters_identity;
        // vector<vector<int>> clusters_size;
        // vector<vector<int>> clusters_coverage;
        // vector<string> rep_identifier_cur;
        // vector<int> rep_size_cur;
        // vector<string> rep_identifier_next;
        // vector<int> rep_size_next;
        omp_set_num_threads(T);

        for (i = 0; i < chunks_num; i++) {
            auto start = std::chrono::high_resolution_clock::now();
            int start_rep_suffix = rep_seqs.size();
            int N = sequences.size();
            neigh.clear();
            neigh.assign(N, {});
            int centers = 0;
            double tA = get_time();
            if (T == 1) {
                for (j = 0; j < N; j++) {
                    Sequence *seq = sequences[j];
                    if (seq->state & IS_REDUNDANT) continue;

                    ClusterOne_single(seq, j, word_table, params[0], buffers[0], options, centers);
                    if (seq->state & IS_REP) {
                        fout << ">" << seq->identifier << "\n";
                        fout << seq->true_data << "\n";

                        // clstr_fout << seq->cluster_id << "\t";
                        // clstr_fout << seq->size << "\t0\t>" << seq->identifier << "\n";
                    }
                }
            } else {
#pragma omp parallel for schedule(dynamic, 1)
                for (j = 0; j < N; j++) {
                    Sequence *seq = sequences[j];
                    if (seq->state & IS_REDUNDANT) continue;
                    int tid = omp_get_thread_num();
                    ClusterOne_worker(seq, j, word_table, params[tid], buffers[tid], options, locks, NUM_LOCKS,
                                      LOCK_MASK);
                }
#pragma omp parallel for schedule(dynamic)
                for (long long j = 0; j < (long long) word_table.NAAN; ++j) {
                    NVector<IndexCount> &row = word_table.indexCounts[j];

                    if (row.size < 2 || row.items == NULL) continue;
                   ips4o::sort(row.items, row.items + row.size,
                              [](const IndexCount &a, const IndexCount &b) { return a.index < b.index; });
                }

                double tA1 = get_time();
                std::vector<size_t> cnt(N, 0);
                int num_batches = 20;
                int batch_size = (N + num_batches - 1) / num_batches;
                for (int batch_idx = 0; batch_idx < num_batches; ++batch_idx)
                {
                    int start = batch_idx * batch_size;
                    int end = (start + batch_size > N) ? N : (start + batch_size);

                    if (start >= end)
                        break; 
#pragma omp parallel for schedule(dynamic, 1)
                    for (int j = start; j < end; j++)
                    {
                        int tid = omp_get_thread_num();
                        Sequence *seq = sequences[j];
                        if (seq->state & IS_REDUNDANT)
                            continue;
                        int flag = CheckOne_master(seq, j, word_table, params[tid], buffers[tid], options);
                        if (flag == 0)
                            seq->state |= IS_REP;
                    }

                    for (int t = 0; t < T; t++)
                    {
                        auto &vec = buffers[t].thread_edges;
                        for (auto &e : vec)
                            ++cnt[e.first];
                    }

                    #pragma omp parallel for schedule(static)
                    for (int u = start; u < end; ++u)
                        if (cnt[u])
                            neigh[u].reserve(neigh[u].size() + cnt[u]);

                    for (int t = 0; t < T; t++)
                    {
                        auto &vec = buffers[t].thread_edges;
                        for (auto &e : vec)
                            neigh[e.first].push_back(e.second);
                        vec.clear(); 
                    }

                    for (int j = start; j < end; ++j)
                    {
                        Sequence *seq = sequences[j];
                        if (seq->state & IS_REDUNDANT)
                            continue;

                        if (!neigh[j].empty())
                        {
                            for (int v : neigh[j])
                            {
                                if (sequences[v]->state & IS_REP)
                                {
                                    seq->state |= IS_REDUNDANT;
                                    break;
                                }
                            }
                        }
                        if (seq->state & IS_REDUNDANT)
                            continue;

                        int size = rep_seqs.size();
                        rep_seqs.Append(seq->index);
                        seq->cluster_id = size;
                        seq->identity = 0;
                        seq->state |= IS_REP;
                        seq->table_idx = centers;
                        ++centers;
                    }
                } // End of batch loop
                for (int j = 0; j < N; ++j)
                {
                    Sequence *seq = sequences[j];
                    if (seq->state & IS_REDUNDANT)
                        continue;
                    fout << ">" << seq->identifier << "\n";
                    fout << seq->true_data << "\n";
                    // clstr_fout << ...
                }

                double tA2 = get_time();
                std::cerr << "checkone time   : " << (tA2 - tA1) << " s\n";
            }

            printf("\n%9i  finished  %9i  clusters\n", start_global_id + sequences.size(), rep_seqs.size());
            int end_rep_suffix = rep_seqs.size();
            encode_WordTable(info_buf, i, start_rep_suffix, end_rep_suffix, cluster_id_buf, seqs_suffix_buf,
                             indexCount_buf, prefix_buf, indexCount_buf_size, prefix_size, send_file_index,
                             start_global_id);
            word_table.Clear();
            if (i > 0) {
                MPI_Request request[4] = {MPI_REQUEST_NULL};
                MPI_Ibcast(&ibcast_flag, 1, MPI_INT, source, MPI_COMM_WORLD, &request[0]);
                MPI_Ibcast((void *) info_buf, 7, MPI_LONG, source, MPI_COMM_WORLD, &request[1]);
                MPI_Ibcast((void *) cluster_id_buf, (int) info_buf[1], MPI_LONG, source, MPI_COMM_WORLD, &request[2]);
                MPI_Ibcast((void *) seqs_suffix_buf, (int) info_buf[1], MPI_LONG, source, MPI_COMM_WORLD, &request[3]);
                MPI_Waitall(4, request, MPI_STATUSES_IGNORE);
                cerr << "send over " << info_buf[0] << endl;
                double tB = get_time();
                std::cerr << "master time   : " << (tB - tA) << " s\n";
                //--------------------------------------------------
                // 			clusters_identifier.resize(C);
                // 			clusters_size.resize(C);
                // 			clusters_identity.resize(C);
                // 			clusters_coverage.resize(C);
                // 			std::vector<omp_lock_t> locks(C);
                // 			for (int c = 0; c < C; ++c)
                // 				omp_init_lock(&locks[c]);

                // #pragma omp parallel
                // #pragma omp single
                // 			for (int src = 1; src < rank_size; ++src)
                // 			{
                // 				int N;
                // 				MPI_Recv(&N, 1, MPI_INT, src, 101, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                // 				int *prefix_seq = (int *)malloc((C + 1) * sizeof(int));
                // 				int *flat_size = (int *)malloc(N * sizeof(int));
                // 				float *flat_identity = (float *)malloc(N * sizeof(float));
                // 				int *flat_coverage = (int *)malloc(N * 4 * sizeof(int));
                // 				char *flat_identifier = (char *)malloc((size_t)N * (max_idf + 1));

                // 				MPI_Recv(prefix_seq, C + 1, MPI_INT, src, 110, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                // 				MPI_Recv(flat_size, N, MPI_INT, src, 111, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                // 				MPI_Recv(flat_identity, N, MPI_FLOAT, src, 112, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                // 				MPI_Recv(flat_coverage, N * 4, MPI_INT, src, 113, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                // 				MPI_Recv(flat_identifier, N * (max_idf + 1), MPI_CHAR, src, 114, MPI_COMM_WORLD,
                // MPI_STATUS_IGNORE);

                // #pragma omp task firstprivate(N, prefix_seq, flat_size, flat_identity, flat_coverage,
                // flat_identifier)
                // 					{
                // 						for (int c = 0; c < C; ++c)
                // 						{
                // 							int begin = prefix_seq[c], end = prefix_seq[c + 1];
                // 							std::vector<std::string> ids;
                // 							std::vector<int> sizes, covs;
                // 							std::vector<float> idts;
                // 							ids.reserve(end - begin);
                // 							sizes.reserve(end - begin);
                // 							idts.reserve(end - begin);
                // 							covs.reserve((end - begin) * 4);

                // 							for (int i = begin; i < end; ++i)
                // 							{
                // 								const char *idptr = flat_identifier + (size_t)i * (max_idf + 1);
                // 								ids.emplace_back(idptr);
                // 								sizes.emplace_back(flat_size[i]);
                // 								idts.emplace_back(flat_identity[i]);
                // 								for (int k = 0; k < 4; ++k)
                // 									covs.emplace_back(flat_coverage[i * 4 + k]);
                // 							}

                // 							omp_set_lock(&locks[c]);
                // 							clusters_identifier[c].insert(clusters_identifier[c].end(), ids.begin(),
                // ids.end()); 							clusters_size[c].insert(clusters_size[c].end(), sizes.begin(),
                // sizes.end()); 							clusters_identity[c].insert(clusters_identity[c].end(),
                // idts.begin(), idts.end()); clusters_coverage[c].insert(clusters_coverage[c].end(), covs.begin(),
                // covs.end()); 							omp_unset_lock(&locks[c]);
                // 						}
                // 						free(prefix_seq);
                // 						free(flat_size);
                // 						free(flat_identity);
                // 						free(flat_coverage);
                // 						free(flat_identifier);
                // 					}
                // 				}
                // #pragma omp taskwait

                // 				for (int c = 0; c < C; ++c)
                // 					omp_destroy_lock(&locks[c]);

                // for (int kk = 0; kk < C; kk++)
                // {
                // 	clstr_fout << ">Cluster " << (output_index + kk) << "\n";
                // 	clstr_fout << 0 << '\t' << rep_size_cur[kk] << "aa, >" << rep_identifier_cur[kk] << "..." << " *"
                // <<"\n"; for (int kkk = 0; kkk < clusters_identifier[kk].size(); kkk++)
                // {
                // 	clstr_fout << kkk + 1 << '\t' << clusters_size[kk][kkk] << "aa, >" << clusters_identifier[kk][kkk]
                // << "..."; 	int *c = &clusters_coverage[kk][kkk * 4]; 	clstr_fout << " at "; 	if (options.print)
                // 		clstr_fout << c[0] << ":" << c[1] << ":" << c[2] << ":" << c[3] << "/";
                // 	clstr_fout << std::fixed << std::setprecision(2) << (clusters_identity[kk][kkk] * 100) << "%";

                // 	clstr_fout << "\n";
                // }
                // }

                // 				clusters_identifier.clear();

                // 				clusters_size.clear();
                // 				clusters_identity.clear();
                // 				clusters_coverage.clear();
                // 				rep_size_cur.clear();
                // 				rep_identifier_cur.clear();
                // 				read_flag[i] = 1;
            } else {
                MPI_Bcast((void *) info_buf, 7, MPI_LONG, source, MPI_COMM_WORLD);
                MPI_Bcast((void *) cluster_id_buf, (int) info_buf[1], MPI_LONG, source, MPI_COMM_WORLD);
                MPI_Bcast((void *) seqs_suffix_buf, (int) info_buf[1], MPI_LONG, source, MPI_COMM_WORLD);
            }
            output_index = last_rep_index;
            C = rep_seqs.size() - last_rep_index;
            // rep_identifier_cur.swap(rep_identifier_next);
            // rep_size_cur.swap(rep_size_next);
            // rep_identifier_next.clear();
            // rep_size_next.clear();
            // if (i == chunks_num - 1)
            // {
            // 				clusters_identifier.resize(C);
            // 				clusters_size.resize(C);
            // 				clusters_identity.resize(C);
            // 				clusters_coverage.resize(C);
            // 				std::vector<omp_lock_t> locks(C);
            // 				for (int c = 0; c < C; ++c)
            // 					omp_init_lock(&locks[c]);

            // #pragma omp parallel
            // #pragma omp single
            // 			for (int src = 1; src < rank_size; ++src)
            // 			{
            // 				int N;
            // 				MPI_Recv(&N, 1, MPI_INT, src, 101, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            // 				int *prefix_seq = (int *)malloc((C + 1) * sizeof(int));
            // 				int *flat_size = (int *)malloc(N * sizeof(int));
            // 				float *flat_identity = (float *)malloc(N * sizeof(float));
            // 				int *flat_coverage = (int *)malloc(N * 4 * sizeof(int));
            // 				char *flat_identifier = (char *)malloc((size_t)N * (max_idf + 1));

            // 				MPI_Recv(prefix_seq, C + 1, MPI_INT, src, 110, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            // 				MPI_Recv(flat_size, N, MPI_INT, src, 111, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            // 				MPI_Recv(flat_identity, N, MPI_FLOAT, src, 112, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            // 				MPI_Recv(flat_coverage, N * 4, MPI_INT, src, 113, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            // 				MPI_Recv(flat_identifier, N * (max_idf + 1), MPI_CHAR, src, 114, MPI_COMM_WORLD,
            // MPI_STATUS_IGNORE); #pragma omp task firstprivate(N, prefix_seq, flat_size, flat_identity, flat_coverage,
            // flat_identifier)
            // 					{
            // 						for (int c = 0; c < C; ++c)
            // 						{
            // 							int begin = prefix_seq[c], end = prefix_seq[c + 1];
            // 							std::vector<std::string> ids;
            // 							std::vector<int> sizes, covs;
            // 							std::vector<float> idts;
            // 							ids.reserve(end - begin);
            // 							sizes.reserve(end - begin);
            // 							idts.reserve(end - begin);
            // 							covs.reserve((end - begin) * 4);

            // 							for (int i = begin; i < end; ++i)
            // 							{
            // 								const char *idptr = flat_identifier + (size_t)i * (max_idf + 1);
            // 								ids.emplace_back(idptr);
            // 								sizes.emplace_back(flat_size[i]);
            // 								idts.emplace_back(flat_identity[i]);
            // 								for (int k = 0; k < 4; ++k)
            // 									covs.emplace_back(flat_coverage[i * 4 + k]);
            // 							}

            // 							omp_set_lock(&locks[c]);
            // 							clusters_identifier[c].insert(clusters_identifier[c].end(), ids.begin(),
            // ids.end()); 							clusters_size[c].insert(clusters_size[c].end(), sizes.begin(),
            // sizes.end()); 							clusters_identity[c].insert(clusters_identity[c].end(),
            // idts.begin(), idts.end()); clusters_coverage[c].insert(clusters_coverage[c].end(), covs.begin(),
            // covs.end()); 							omp_unset_lock(&locks[c]);
            // 						}
            // 						free(prefix_seq);
            // 						free(flat_size);
            // 						free(flat_identity);
            // 						free(flat_coverage);
            // 						free(flat_identifier);
            // 					}
            // 				}
            // #pragma omp taskwait

            // 				for (int c = 0; c < C; ++c)
            // 					omp_destroy_lock(&locks[c]);
            // for (int kk = 0; kk < C; kk++)
            // {
            // 	clstr_fout << ">Cluster " << (last_rep_index + kk) << "\n";
            // 	clstr_fout << 0 << '\t' << rep_size_cur[kk] << "aa, >" << rep_identifier_cur[kk] << "..." << " *"
            // <<"\n";
            // 	// for (int kkk = 0; kkk < clusters_identifier[kk].size(); kkk++)
            // 	// {
            // 	// 	clstr_fout << kkk + 1 << '\t' << clusters_size[kk][kkk] << "aa, >" << clusters_identifier[kk][kkk]
            // << "...";
            // 	// 	int *c = &clusters_coverage[kk][kkk * 4];
            // 	// 	clstr_fout << " at ";
            // 	// 	if (options.print)
            // 	// 		clstr_fout << c[0] << ":" << c[1] << ":" << c[2] << ":" << c[3] << "/";
            // 	// 	clstr_fout << std::fixed << std::setprecision(2) << (clusters_identity[kk][kkk] * 100) << "%";

            // 	// 	clstr_fout << "\n";
            // 	// }
            // }
            // 				clusters_identifier.clear();

            // 				clusters_size.clear();
            // 				clusters_identity.clear();
            // 				clusters_coverage.clear();
            // 				rep_size_cur.clear();
            // 				rep_identifier_cur.clear();
            // }
            if (i == chunks_num - 1) {
                double tb1 = get_time();
                vector<MasterSeqInfo> all_infos(total_num);
                size_t base = 0;
                for (int src = 1; src < rank_size; ++src) {
                    int C_remote = id_tables[src - 1].size();
                    std::vector<RedundantSeqInfoHeader> buf(C_remote);

                    int max_seqs_hdr = INT_MAX / (int) sizeof(RedundantSeqInfoHeader);
                    int CHUNK_SEQS = min(max_seqs_hdr, C_remote);

                    int received = 0;
                    while (received < C_remote) {
                        int n = min(CHUNK_SEQS, C_remote - received);

                        MPI_Recv(buf.data() + received, n * (int) sizeof(RedundantSeqInfoHeader), MPI_BYTE, src, 500,
                                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                        received += n;
                    }
#pragma omp parallel for
                    for (int j = 0; j < C_remote; ++j) {
                        const RedundantSeqInfoHeader &h = buf[j];
                        MasterSeqInfo &m = all_infos[base + j];
                        m.cluster_id = h.cluster_id;
                        m.size = h.size;
                        m.identity = h.identity;
                        m.coverage[0] = h.coverage[0];
                        m.coverage[1] = h.coverage[1];
                        m.coverage[2] = h.coverage[2];
                        m.coverage[3] = h.coverage[3];
                        m.name = id_tables[src - 1][j];
                    }
                    base += C_remote;
                }

                double tb2 = get_time();
                cerr << "total recv time " << tb2 - tb1 << endl;

                cout << "[DONE] Representative Information written to: " << rep_output << endl;
                SortAndWriteResult(all_infos, options);
                double tb3 = get_time();
                cerr << "total write time " << tb3 - tb2 << endl;
                break;
            }
            last_rep_index = rep_seqs.size();
            start_global_id += sequences.size();
            for (int i = 0; i < sequences.size(); i++) delete sequences[i];
            sequences.clear();
            Sequence one;
            std::string file = preprocess_output_dir + "_proc" + std::to_string(file_index) + ".fa";
            cerr << "master  file_index      " << file_index << endl;
            send_file_index = file_index;
            if (chunk_kseq[file_index] == nullptr) {
                gzFile fp = gzopen(file.c_str(), "r");
                if (!fp) {
                    fprintf(stderr, "Cannot open file: %s\n", file.c_str());
                    exit(1);
                }
                chunk_fp[file_index] = fp;
                seq = kseq_init(fp);
            } else
                seq = chunk_kseq[file_index];
            char *id_ptr = new char[max_idf + 1];
            char *data_ptr = new char[max_len + 1];
            while ((len = kseq_read(seq)) >= 0) {
                memcpy(id_ptr, seq->name.s, seq->name.l);
                id_ptr[seq->name.l] = 0;
                one.identifier = id_ptr;

                memcpy(data_ptr, seq->seq.s, seq->seq.l);
                data_ptr[seq->seq.l] = 0;
                one.data = data_ptr;
                one.size = len;
                one.index = sequences.size() + start_global_id;
                one.master_flag = 1;
                sequences.Append(new Sequence(one));
                id_tables[file_index].emplace_back(seq->name.s, seq->name.l);
                now_bytes += len;
                sequences[sequences.size() - 1]->ConvertBases();
                if (now_bytes > chunk_bytes || sequences.size() >= chunk_size) {
                    chunks_id.push_back(chunk_id);
                    chunk_id++;
                    chunk_kseq[file_index] = seq;
                    now_bytes = 0;
                    break;
                }
            }
            if (now_bytes > 0) {
                cerr << "first seq  " << sequences[0]->identifier << endl;
                cerr << "last seq    " << sequences[sequences.size() - 1]->identifier << endl;
                chunks_id.push_back(chunk_id);
            }
            file_index = (file_index + 1) % (rank_size - 1);
            one.identifier = nullptr;
            one.data = nullptr;
            delete[] id_ptr;
            delete[] data_ptr;
            int continue_size = 0;
            int size = sequences.size();
            int *rep_chunk = (int *) malloc(size * 2 * sizeof(int));
            if (rank_size == 2) target_worker = 0;
            cerr << "receive size " << size << endl;
            MPI_Recv(rep_chunk, size * 2, MPI_INT, target_worker + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            cerr << "Recevie chunk " << i + 1 << " by worker " << target_worker + 1 << endl;
            target_worker = (target_worker + 1) % (rank_size - 1);
#pragma omp parallel for 

            for (int j = 0; j < sequences.size(); j++) {
                int index = j * 2;
                Sequence *seq = sequences[j];
                if (rep_chunk[index] == 0) {
                    continue;
                }
                seq->state = (short) rep_chunk[index];
                seq->cluster_id = rep_chunk[index + 1];
            }

            free(rep_chunk);
            rep_chunk = NULL;
            word_table.Clear();

            // free momery
            free(info_buf);
            free(cluster_id_buf);
            free(seqs_suffix_buf);
            free(prefix_buf);
            free(indexCount_buf);
            info_buf = NULL;
            cluster_id_buf = NULL;
            seqs_suffix_buf = NULL;
            prefix_buf = NULL;
            indexCount_buf = NULL;
            auto end1 = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed1 = end1 - start;
            std::cerr << "-----chunk " << i << "  total time   " << elapsed1.count() << " 秒\n";
        }
        for (int i = 0; i < rank_size - 1; i++) {
            if (chunk_kseq[i]) kseq_destroy(chunk_kseq[i]);
            if (chunk_fp[i]) gzclose(chunk_fp[i]);
        }

        fout.close();
        // clstr_fout.close();
        free(info_buf);
        free(cluster_id_buf);
        free(seqs_suffix_buf);
        free(prefix_buf);
        free(indexCount_buf);
        info_buf = NULL;
        cluster_id_buf = NULL;
        seqs_suffix_buf = NULL;
        prefix_buf = NULL;
        indexCount_buf = NULL;
    }

    if (worker) {
        vector<int> read_chunk(chunks_num, 0);
        kseq_t *seq = nullptr;
        vector<gzFile> chunk_fp(rank_size - 1, nullptr);
        vector<kseq_t *> chunk_kseq(rank_size - 1, nullptr);
        std::vector<int> rep_chunk_buf;
        rep_chunk_buf.reserve(1024);
        int top = 0;
        int bottom = sub_chunks.size();
        int start = 0;
        int now_byte = 0;
        int *rep_chunk;
        int chunksNum = all_chunks.size() * (rank_size - 1) + 1;
        int remaining = 0;
        long start_id = 0;
        int read_stop = 0;
        Slot slots[2];
        int cur = 0, next = 1;
        int first_flag = 1;
        double total_time = 0;
        double total_wait_time = 0;
        post_ibcasts_for_this_block(slots[cur], source, MPI_COMM_WORLD);
        while (1) {
            int tb[2];
            int remote_idxs[3];
            int ibcast_flag = 0;
            int done_flag = 0;
            int over_flag = 0;
            int now_rank = 0;
            MPI_Request request = MPI_REQUEST_NULL;
            int send_flag = 0;
            int *prefix_seq, *flat_size, *flat_coverage;
            float *flat_identity;
            char *flat_identifier;
            int soure_chunk = slots[cur].info[0];
            int file_index = slots[cur].info[5];
            start_id = slots[cur].info[6];

            if (soure_chunk < chunks_num - 1) {
                MPI_Ibcast(&ibcast_flag, 1, MPI_INT, source, MPI_COMM_WORLD, &request);
            }
            if (read_chunk[soure_chunk] == 0) {
                int id;
                for (int i = 0; i < rep_sequences.size(); i++) delete rep_sequences[i];
                rep_sequences.clear();
                Sequence one;
                std::string file = preprocess_output_dir + "_proc" + std::to_string(file_index) + ".fa";
                cerr << "file_index      " << file_index << endl;
                if (chunk_kseq[file_index] == nullptr) {
                    gzFile fp = gzopen(file.c_str(), "r");
                    if (!fp) {
                        fprintf(stderr, "Cannot open file: %s\n", file.c_str());
                        exit(1);
                    }
                    chunk_fp[file_index] = fp;
                    seq = kseq_init(fp);
                } else
                    seq = chunk_kseq[file_index];
                char *id_ptr = new char[max_idf + 1];
                char *data_ptr = new char[max_len + 1];
                while ((len = kseq_read(seq)) >= 0) {
                    memcpy(id_ptr, seq->name.s, seq->name.l);
                    id_ptr[seq->name.l] = 0;
                    one.identifier = id_ptr;

                    memcpy(data_ptr, seq->seq.s, seq->seq.l);
                    data_ptr[seq->seq.l] = 0;
                    one.data = data_ptr;

                    one.size = len;
                    one.index = rep_sequences.size();
                    one.master_flag = 0;
                    rep_sequences.Append(new Sequence(one));
                    now_byte += len;

                    rep_sequences[rep_sequences.size() - 1]->ConvertBases();
                    if ((soure_chunk != 0 && now_byte > chunk_bytes) ||
                        (soure_chunk == 0 && rep_sequences.size() >= first_chunk_size) ||
                        (rep_sequences.size() >= chunk_size)) {
                        chunk_kseq[file_index] = seq;
                        now_byte = 0;
                        break;
                    }
                }
                read_chunk[soure_chunk] = 1;
                one.identifier = nullptr;
                one.data = nullptr;
                delete[] id_ptr;
                delete[] data_ptr;
            }

            record += slots[cur].info[1];
            decode_WordTable(word_table, start, slots[cur],options);

            double t141 = get_time();
            omp_set_num_threads(T);
#pragma omp parallel for schedule(dynamic, 1)
            for (int j = 0; j < rep_sequences.size(); j++) {
                int tid = omp_get_thread_num();
                Sequence *seq = rep_sequences[j];
                if (seq->state & IS_REP)
                    ClusterOne_worker(seq, seq->table_idx, word_table, params[tid], buffers[tid], options, locks,
                                      NUM_LOCKS, LOCK_MASK);
            }
            double t142 = get_time();
            cerr << "build word table time   " << t142 - t141 << " chunk  " << soure_chunk << "  by rank  " << my_rank
                 << endl;
            total_time += t142 - t141;
            int remain_chunks = my_chunks.size() - start;

            double t14 = get_time();
            if (options.stealing) {
#pragma omp parallel
                {
                    int tid = omp_get_thread_num();

                    for (int i = 0; i < remain_chunks; ++i) {
                        int idx = i + start;

                        for (int s = 0; s < SUB; ++s) {
                            static int l_shared = 0, r_shared = -1;
                            static int have_task = 0;
#pragma omp master
                            {
                                MPI_Win_sync(win_ctrl_);
                                top = ptr_ctrl_[0];
                                bottom = ptr_ctrl_[1];
                                // MPI_Get(tb, 2, MPI_INT, worker_rank, 0, 2, MPI_INT, win_ctrl_);
                                // // MPI_Get(&top, 1, MPI_INT, worker_rank, 0, 1, MPI_INT, win_ctrl_);
                                // // MPI_Get(&bottom, 1, MPI_INT, worker_rank, 1, 1, MPI_INT, win_ctrl_);
                                // MPI_Win_flush_local(worker_rank, win_ctrl_);
                                // top = tb[0];
                                // bottom = tb[1];
                                // cerr<<"bottom "<<bottom<<endl;
                                if (top <= bottom) {
                                    int dec = 1;
                                    // Task t;
                                    

                                    MPI_Fetch_and_op(&dec, &top, MPI_INT, worker_rank, 0, MPI_SUM, win_ctrl_);
                                    MPI_Win_flush(worker_rank, win_ctrl_);
                                    // MPI_Get(&t, sizeof(Task), MPI_BYTE, worker_rank, top, sizeof(Task), MPI_BYTE,
                                    //         win_tasks_);
                                    // MPI_Win_flush_local(worker_rank, win_tasks_);
                                    // l_shared = t.l;
                                    // r_shared = t.r;
                                    l_shared = sub_chunks[top].first;
                                    r_shared = sub_chunks[top].second;
                                    have_task = 1;
                                } else {
                                    l_shared = 1;
                                    r_shared = 0;
                                    have_task = 0;
                                }
                            }
#pragma omp barrier

                            if (!have_task) {
                                break;
                            }

#pragma omp for schedule(dynamic,1)
                            for (int k = l_shared; k <= r_shared; ++k) {
                                Sequence *seq = sequences[k];
                                auto &m = meta_[k];
                                if ((m.state & IS_REDUNDANT) || (m.state & IS_REP)) continue;
                                if (tid == 0 && !over_flag)
                                {
                                    if (!done_flag)
                                    {
                                        int flag_local = 0;
                                        MPI_Test(&request, &flag_local, MPI_STATUS_IGNORE);
                                        if (ibcast_flag)
                                        {
                                            post_ibcasts_for_next_block(slots[next], source, MPI_COMM_WORLD);
                                            ibcast_flag = 0;
                                            done_flag = 1;
                                        }
                                    }
                                    if (done_flag)
                                    {
                                        int flag1 = 0;
                                        int flag2 = 0;
                                        MPI_Test(&slots[next].reqs[0], &flag1, MPI_STATUS_IGNORE);
                                        MPI_Test(&slots[next].reqs[1], &flag2, MPI_STATUS_IGNORE);
                                        if (flag1 && flag2)
                                            over_flag = 1;
                                    }
                                }

                                CheckOne_worker(seq, word_table, params[tid], buffers[tid], options, k);
                            }

// #pragma omp master
//                             {
//                                 const int cnt = r_shared - l_shared + 1;
//                                 MPI_Put(&meta_[l_shared], cnt * sizeof(SeqMeta), MPI_BYTE, worker_rank, l_shared,
//                                         cnt * sizeof(SeqMeta), MPI_BYTE, win_meta_);
//                                 MPI_Win_flush(worker_rank, win_meta_);
//                             }
                        }

#pragma omp master
                        {
                            if (chunks_id[idx] == soure_chunk + 1) {
                                const int size = my_chunks[idx].second - my_chunks[idx].first + 1;
                                rep_chunk_buf.resize((size_t) size * 2); // 只在容量不够时才会扩容

                                for (int j = my_chunks[idx].first; j <= my_chunks[idx].second; ++j) {
                                    int k = (j - my_chunks[idx].first) * 2;
                                    // Sequence *seq = sequences[j];
                                    auto &m = meta_[j];
                                    if (m.state & IS_REDUNDANT) {
                                        rep_chunk_buf[k] = (int)m.state;
                                        rep_chunk_buf[k + 1] = -1;
                                    } else {
                                        rep_chunk_buf[k] = 0;
                                        rep_chunk_buf[k + 1] = -1;
                                    }
                                }
                                MPI_Send(rep_chunk_buf.data(), size * 2, MPI_INT, source, 0, MPI_COMM_WORLD);
                            }
                        }
#pragma omp barrier
                    }
                }
                                double t155 = get_time();
                cerr << "-----checkone1 time  " << t155 - t14 << "  by rank  " << my_rank << endl;
                // int now_top;
                // double t1555 = get_time();
                // MPI_Win_sync(win_meta_);
                // double t1556 = get_time();
                //  cerr << "-----out time  " << t1556 - t1555 << "  by rank  " << my_rank << endl;
                // int origin_top = ctrl_[2];
                // MPI_Win_sync(win_ctrl_);
                // int now_top = ptr_ctrl_[0];
                // // MPI_Get(&now_top, 1, MPI_INT, worker_rank, 0, 1, MPI_INT, win_ctrl_);
                // // MPI_Win_flush_local(worker_rank, win_ctrl_);
                // Task t1,t2;
                // if (now_top > origin_top)
                // {
                //     // t2.r = sub_chunks[now_top-1].second;
                //     // t1.l = sub_chunks[origin_top].first;
                //     // // cerr << "origin_top " << origin_top << "now_top " << now_top << endl;
                //     // // cerr << "t2.r  " << t2.r << "t1.l " << t1.l << endl;

                //     // const int cnt = t2.r - t1.l + 1;
                //     // cerr << "cnt  " << cnt << endl;
                //     // MPI_Put(&meta_[t1.l], cnt * sizeof(SeqMeta), MPI_BYTE, worker_rank, t1.l,
                //     //         cnt * sizeof(SeqMeta), MPI_BYTE, win_meta_);
                //     // MPI_Win_flush(worker_rank, win_meta_);
                //     MPI_Win_sync(win_meta_);
                // }


                progress_running = true;
                std::thread progress(mpi_progress_thread);
                int worker_size = rank_size - 1;
                for (int offset = -3; offset <= 3; ++offset) {
                    int tt = (worker_rank + offset + worker_size) % worker_size;
                    if (tt == worker_rank) continue;
                    while (true) {
                        // int stealing_top, stealing_bottom, origin_top;
                        MPI_Get(remote_idxs, 3, MPI_INT, tt, 0, 3, MPI_INT, win_ctrl_);
                        MPI_Win_flush_local(tt, win_ctrl_);
                        int stealing_top = remote_idxs[0];
                        int stealing_bottom = remote_idxs[1];
                        int origin_top = remote_idxs[2];
                        // cerr<<"stealing rank "<<tt<<endl;
                        // cerr<<"stealing_top "<<stealing_top<<endl;
                        // cerr<<"stealing_bottom "<<stealing_bottom<<endl;
                        // cerr<<"origin_top "<<origin_top<<endl;
                        // MPI_Get(&stealing_top, 1, MPI_INT, tt, 0, 1, MPI_INT, win_ctrl_);
                        // MPI_Get(&stealing_bottom, 1, MPI_INT, tt, 1, 1, MPI_INT, win_ctrl_);
                        // MPI_Get(&origin_top, 1, MPI_INT, tt, 2, 1, MPI_INT, win_ctrl_);
                        // MPI_Win_flush_local(tt, win_ctrl_);
                        if (stealing_bottom - stealing_top > 1&&stealing_bottom-origin_top>SUB) {
                            // cerr<<"stealing rank "<<tt<<endl;
                            int dec = -1;
                            MPI_Fetch_and_op(&dec, &stealing_bottom, MPI_INT, tt, 1, MPI_SUM, win_ctrl_);
                            MPI_Win_flush(tt, win_ctrl_);
                            Task t;
                            MPI_Get(&t, sizeof(Task), MPI_BYTE, tt, stealing_bottom, sizeof(Task), MPI_BYTE,
                                    win_tasks_);
                            MPI_Win_flush_local(tt, win_tasks_);
                            
                            const int cnt = t.r - t.l + 1;
                            std::vector<SeqMeta> metas(cnt);

                            MPI_Get(metas.data(), (int) (cnt * sizeof(SeqMeta)), MPI_BYTE, tt, /*disp = 元素下标*/ t.l,
                                    (int) (cnt * sizeof(SeqMeta)), MPI_BYTE, win_meta_);

                            MPI_Win_flush_local(tt, win_meta_);

                            size_t total_bytes = 0;
                            for (const auto &m : metas) total_bytes += (size_t) m.data_len;
                            std::vector<uint8_t> slab(total_bytes);
                            size_t cursor = metas[0].data_off;
                            MPI_Get(slab.data(), total_bytes, MPI_BYTE, tt, (MPI_Aint) cursor, total_bytes, MPI_BYTE,
                                    win_pool_d_);
                            MPI_Win_flush_local(tt, win_pool_d_);
// #pragma omp parallel num_threads(T)
//                             {
// #pragma omp for schedule(dynamic,1)
//                                 for (int ttt = 0; ttt < cnt; ttt++) {
//                                     int tid = omp_get_thread_num();
//                                     auto &m = metas[ttt];
//                                     if ((m.state & IS_REDUNDANT) || (m.state & IS_REP)) continue;
//                                     Sequence *seq = new Sequence();

//                                     char *buf = (char *) malloc((size_t) m.data_len + 1);
//                                     memcpy(buf, slab.data() + m.data_off - cursor, (size_t) m.data_len);
//                                     buf[m.data_len] = '\0';
//                                     seq->data = buf;
//                                     seq->size = m.size;
//                                     seq->state = m.state;
//                                     seq->cluster_id = m.cluster_id;
//                                     seq->identity = m.identity;
//                                     seq->distance = m.distance;
//                                     for (int tttt = 0; tttt < 4; ++tttt) seq->coverage[tttt] = m.coverage[tttt];
//                                     CheckOne_stealing(seq, word_table, params[tid], buffers[tid], options);

//                                     if ((seq->state & IS_REDUNDANT)) {
//                                         auto &m = metas[ttt];
//                                         m.state = seq->state;
//                                         m.cluster_id = seq->cluster_id;
//                                         m.identity = seq->identity;
//                                         m.distance = seq->distance;
//                                         for (int tttt = 0; tttt < 4; ++tttt) m.coverage[tttt] = seq->coverage[tttt];
//                                     }

//                                     delete seq;
//                                 }
//                             }
#pragma omp parallel
                            {
                                std::vector<char> local_buf;

#pragma omp for schedule(dynamic, 1) 
                                for (int ttt = 0; ttt < cnt; ttt++)
                                {
                                    auto &m = metas[ttt];

                                    if ((m.state & IS_REDUNDANT) || (m.state & IS_REP))
                                        continue;
                                    Sequence seq ;
                                    if (local_buf.size() < m.data_len + 1)
                                    {
                                        local_buf.resize(m.data_len + 1);
                                    }
                                    size_t slab_offset = m.data_off - cursor;
                                    memcpy(local_buf.data(), slab.data() + slab_offset, m.data_len);
                                    local_buf[m.data_len] = '\0';
                                    seq.data = local_buf.data();
                                    seq.size = m.size;
                                    seq.state = m.state;
                                    seq.cluster_id = m.cluster_id;
                                    seq.identity = m.identity;
                                    seq.distance = m.distance;
                                    for (int tttt = 0; tttt < 4; ++tttt) seq.coverage[tttt] = m.coverage[tttt];
                                    int tid = omp_get_thread_num();
                                    CheckOne_stealing(&seq, word_table, params[tid], buffers[tid], options);
                                    seq.data = nullptr;
                                    if (seq.state & IS_REDUNDANT)
                                    {
                                        m.state = seq.state;
                                        m.cluster_id = seq.cluster_id;
                                        m.identity = seq.identity;
                                        m.distance = seq.distance;
                                        for (int tttt = 0; tttt < 4; ++tttt) seq.coverage[tttt] = m.coverage[tttt];
                                    }
                                }
                            } // end parallel

                            // int task_flag = 0;
                            // MPI_Get(&task_flag, 1, MPI_INT, tt, stealing_bottom, 1, MPI_INT, win_tasks_flag_);
                            // MPI_Win_flush_local(tt, win_tasks_flag_);
                            // if (task_flag == 0) {
                                MPI_Put(metas.data(), cnt * sizeof(SeqMeta), MPI_BYTE, tt, t.l, cnt * sizeof(SeqMeta),
                                        MPI_BYTE, win_meta_);
                                MPI_Win_flush(tt, win_meta_);
                                // int finished_val = 1;

                                // MPI_Put(&finished_val, 1, MPI_INT, tt, stealing_bottom, 1, MPI_INT, win_tasks_flag_);
                                // MPI_Win_flush(tt, win_tasks_flag_);
                                // }

                        } else
                            break;
                    }
                }
                progress_running = false;
                progress.join();
            } else {
#pragma omp parallel
                {
                    int tid = omp_get_thread_num();

                    for (int i = 0; i < remain_chunks; ++i) {
                        int idx = i + start;

#pragma omp for schedule(dynamic, 1)
                        for (int j = my_chunks[idx].first; j <= my_chunks[idx].second; ++j) {
                            Sequence *seq = sequences[j];
                            if ((seq->state & IS_REDUNDANT) || (seq->state & IS_REP)) continue;

                            if (tid == 0 && !over_flag) {
                                if (!done_flag) {
                                    int flag_local = 0;
                                    MPI_Test(&request, &flag_local, MPI_STATUS_IGNORE);
                                    if (ibcast_flag) {
                                        post_ibcasts_for_next_block(slots[next], source, MPI_COMM_WORLD);
                                        ibcast_flag = 0;
                                        done_flag = 1;
                                    }
                                }
                                if (done_flag) {
                                    int flag1 = 0;
                                    int flag2 = 0;
                                    MPI_Test(&slots[next].reqs[0], &flag1, MPI_STATUS_IGNORE);
                                    MPI_Test(&slots[next].reqs[1], &flag2, MPI_STATUS_IGNORE);
                                    if (flag1 && flag2) over_flag = 1;
                                }
                            }
                            CheckOne(seq, word_table, params[tid], buffers[tid], options);
                        }

#pragma omp master
                        {
                            if (chunks_id[idx] == soure_chunk + 1) {
                                int size = my_chunks[idx].second - my_chunks[idx].first + 1;
                                rep_chunk_buf.resize((size_t) size * 2); // 只在容量不够时才会扩容

                                for (int j = my_chunks[idx].first; j <= my_chunks[idx].second; ++j) {
                                    int k = (j - my_chunks[idx].first) * 2;
                                    Sequence *seq = sequences[j];
                                    if (seq->state & IS_REDUNDANT) {
                                        rep_chunk_buf[k] = (int) seq->state;
                                        rep_chunk_buf[k + 1] = -1;
                                    } else {
                                        rep_chunk_buf[k] = 0;
                                        rep_chunk_buf[k + 1] = -1;
                                    }
                                }

                                MPI_Send(rep_chunk_buf.data(), size * 2, MPI_INT, source, 0, MPI_COMM_WORLD);
                            }
                        }
                    }
                }
            }
            double t15 = get_time();
            cerr << "-----checkone time  " << t15 - t14 << " chunk  " << soure_chunk << "  by rank  " << my_rank
                 << endl;

            word_table.Clear();
            slots[cur].release();
            double t16 = get_time();
            if (!done_flag && soure_chunk < chunks_num - 1) {
                MPI_Wait(&request, MPI_STATUS_IGNORE);

                post_ibcasts_for_next_block(slots[next], source, MPI_COMM_WORLD);
            }

            wait_all(slots[next]);
            double t17 = get_time();
            cerr << "-----wait time  " << t17 - t16 << "  by rank  " << my_rank << endl;
            total_wait_time += t17 - t16;
            std::swap(cur, next);
            if (options.stealing) {

                MPI_Barrier(worker_comm);
                // MPI_Win_sync(win_ctrl_);
        //         bottom = ptr_ctrl_[1];
        //         // MPI_Get(&bottom, 1, MPI_INT, worker_rank, 1, 1, MPI_INT, win_ctrl_);
        //         // MPI_Win_flush_local(worker_rank, win_ctrl_);
        //         if (bottom < sub_chunks.size() - 1)
        //         {
        //             int start_chunk_idx = bottom + 1;
        //             int end_chunk_idx = sub_chunks.size() - 1;

        //             int global_start = sub_chunks[start_chunk_idx].first;
        //             int global_end = sub_chunks[end_chunk_idx].second;

        //             size_t total_count = global_end - global_start + 1;

        //             std::vector<SeqMeta> all_metas(total_count);
        //             MPI_Win_sync(win_meta_);
        
        // #pragma omp parallel for schedule(static)
        // for (int i = 0; i < (int)total_count; ++i) {
        //     int seq_idx = global_start + i;
        //     const auto &m = meta_[seq_idx];
            
        //     Sequence *seq = sequences[seq_idx];

        //     if ((seq->state & IS_REDUNDANT) || (seq->state & IS_REP)) continue;

        //     if (m.state & IS_REDUNDANT) {
        //         seq->state = m.state;
        //         seq->cluster_id = m.cluster_id;
        //         seq->identity = m.identity;
        //         seq->distance = m.distance;
        //         for (int t = 0; t < 4; ++t) seq->coverage[t] = m.coverage[t];
        //     }
        // }
        //         }
            }
            // double tb1 = get_time();
            // int C = record - record_last;
            // 	vector<vector<string>> clusters_identifier(C);
            // 	vector<vector<float>> clusters_identity(C);
            // 	vector<vector<int>> clusters_size(C);
            // 	vector<vector<int>> clusters_coverage(C);
            // 	for (i = 0; i < remain_chunks; i++)

            // 	{
            // 		int idx ;
            // 		 idx = i + start;
            // 		for (j = my_chunks[idx].first; j <= my_chunks[idx].second; j++)
            // 		{
            // 			Sequence *seq = sequences[j];
            // 			if (seq->state & IS_REDUNDANT && seq->cluster_id != -1)
            // 			{

            // 				clusters_identifier[seq->cluster_id -
            // record_last].emplace_back(std::string(seq->identifier)); 				clusters_size[seq->cluster_id -
            // record_last].emplace_back(seq->size); 				clusters_identity[seq->cluster_id -
            // record_last].emplace_back(seq->identity); 				clusters_coverage[seq->cluster_id -
            // record_last].emplace_back(seq->coverage[0]); 				clusters_coverage[seq->cluster_id -
            // record_last].emplace_back(seq->coverage[1]); 				clusters_coverage[seq->cluster_id -
            // record_last].emplace_back(seq->coverage[2]); 				clusters_coverage[seq->cluster_id -
            // record_last].emplace_back(seq->coverage[3]); 				seq->cluster_id = -1;
            // 			}
            // 		}
            // 	}
            // 	for (int kk = 0; kk < C; ++kk)
            // 	{
            // 		for (size_t m = 0; m < clusters_identifier[kk].size(); ++m)
            // 		{
            // 			clstr_fout  << (record_last + kk) << "\n";
            // 			clstr_fout << (m + 1) << '\t'
            // 					   << clusters_size[kk][m] << "aa, >"
            // 					   << clusters_identifier[kk][m] << "...";

            // 			int *cov = &clusters_coverage[kk][m * 4];

            // 			clstr_fout << " at ";
            // 			if (options.print)
            // 				clstr_fout << cov[0] << ":" << cov[1] << ":" << cov[2] << ":" << cov[3] << "/";

            // 			clstr_fout << std::fixed << std::setprecision(2)
            // 					   << (clusters_identity[kk][m] * 100) << "%\n";
            // 		}
            // 	}

            // 	int  N;
            // 	int IDLEN=max_idf+1;
            // 	send_cluster(clusters_identifier,clusters_size,clusters_identity,clusters_coverage,prefix_seq,flat_size,flat_identity,flat_coverage,flat_identifier,C,N,IDLEN);
            // 	int CHAR_TOTAL = N * (max_idf+1);
            // 	MPI_Send(&N, 1, MPI_INT, 0, 101, MPI_COMM_WORLD);
            // 	MPI_Send(prefix_seq, C + 1, MPI_INT, 0, 110, MPI_COMM_WORLD);
            // 	MPI_Send(flat_size, N, MPI_INT, 0, 111, MPI_COMM_WORLD);
            // 	MPI_Send(flat_identity, N, MPI_FLOAT, 0, 112, MPI_COMM_WORLD);
            // 	MPI_Send(flat_coverage, N * 4, MPI_INT, 0, 113, MPI_COMM_WORLD);
            // 	MPI_Send(flat_identifier, CHAR_TOTAL, MPI_CHAR, 0, 114, MPI_COMM_WORLD);
            // 	double tb2 = get_time();
            // 	cerr << "-----send time  " << tb2 - tb1 << "  by rank  " << my_rank << endl;
            // 	cerr<<"send cluster   "<<soure_chunk<<endl;
            // 	free(prefix_seq);
            // 	free(flat_size);
            // 	free(flat_identity);
            // 	free(flat_coverage);
            // 	free(flat_identifier);
            // 	prefix_seq = nullptr;
            // 	flat_size = nullptr;
            // 	flat_identity = nullptr;
            // 	flat_coverage = nullptr;
            // 	flat_identifier = nullptr;
            record_last = record;
            if (soure_chunk == chunks_num - 1) {
                cerr << "total build time " << total_time << endl;
                cerr << "total wait time " << total_wait_time << endl;
                double tb1 = get_time();
                int C = (int) sequences.size();
                vector<RedundantSeqInfoHeader> info_headers(C);
                if (options.stealing)
                {
#pragma omp parallel for
                    for (int ii = 0; ii < (int)my_chunks.size(); ++ii)
                    {
                        for (int j = my_chunks[ii].first; j <= my_chunks[ii].second; ++j)
                        {
                            auto &m = meta_[j];
                            RedundantSeqInfoHeader &h = info_headers[j];
                            h.cluster_id = m.cluster_id;
                            h.size = m.size;
                            h.identity = m.identity;
                            h.coverage[0] = m.coverage[0];
                            h.coverage[1] = m.coverage[1];
                            h.coverage[2] = m.coverage[2];
                            h.coverage[3] = m.coverage[3];
                        }
                    }
                }
                else
                {
#pragma omp parallel for
                    for (int ii = 0; ii < (int)my_chunks.size(); ++ii)
                    {
                        for (int j = my_chunks[ii].first; j <= my_chunks[ii].second; ++j)
                        {
                            Sequence *seq = sequences[j];
                            RedundantSeqInfoHeader &h = info_headers[j];
                            h.cluster_id = seq->cluster_id;
                            h.size = seq->size;
                            h.identity = seq->identity;
                            h.coverage[0] = seq->coverage[0];
                            h.coverage[1] = seq->coverage[1];
                            h.coverage[2] = seq->coverage[2];
                            h.coverage[3] = seq->coverage[3];
                        }
                    }
                }

                std::vector<MPI_Request> reqs;
                int max_seqs_hdr = INT_MAX / (int) sizeof(RedundantSeqInfoHeader);
                int CHUNK_SEQS = max_seqs_hdr;
                for (int offset = 0; offset < C; offset += CHUNK_SEQS) {
                    int n = min(CHUNK_SEQS, C - offset);
                    MPI_Request req;
                    MPI_Isend(info_headers.data() + offset, n * (int) sizeof(RedundantSeqInfoHeader), MPI_BYTE, 0, 500,
                              MPI_COMM_WORLD, &req);
                    reqs.push_back(req);
                }

                sequences.clear();
                MPI_Waitall(reqs.size(), reqs.data(), MPI_STATUSES_IGNORE);
                double tb2 = get_time();
                cerr << "total send time " << tb2 - tb1 << endl;
                //             #pragma omp parallel for num_threads(T)
                // for (int c = 0; c < C; ++c) {
                //     int base = prefix_seq[c];
                //     for (int j = 0; j < seq_cnt[c]; ++j) {
                //         int g = base + j;

                //         // identifier 拷贝 + 补零
                //         const std::string &id = clusters_identifier[c][j];
                //         char *dst = flat_identifier + g * IDLEN;
                //         memcpy(dst, id.data(), id.size());
                //         memset(dst + id.size(), 0, IDLEN - id.size());

                //         flat_size[g] = clusters_size[c][j];
                //         flat_identity[g] = clusters_identity[c][j];

                //         for (int k = 0; k < 4; ++k) flat_coverage[g * 4 + k] = clusters_coverage[c][j * 4 + k];
                //     }
                // }
                // vector<vector<string>> clusters_identifier(C);
                // vector<vector<float>> clusters_identity(C);
                // vector<vector<int>> clusters_size(C);
                // vector<vector<int>> clusters_coverage(C);
                // for (int ii = 0; ii < my_chunks.size(); ii++)

                // {
                //     for (j = my_chunks[ii].first; j <= my_chunks[ii].second; j++) {
                //         Sequence *seq = sequences[j];
                //         if (seq->state & IS_REDUNDANT) {
                //             clusters_identifier[seq->cluster_id].emplace_back(std::string(seq->identifier));
                //             clusters_size[seq->cluster_id].emplace_back(seq->size);
                //             clusters_identity[seq->cluster_id].emplace_back(seq->identity);
                //             clusters_coverage[seq->cluster_id].emplace_back(seq->coverage[0]);
                //             clusters_coverage[seq->cluster_id].emplace_back(seq->coverage[1]);
                //             clusters_coverage[seq->cluster_id].emplace_back(seq->coverage[2]);
                //             clusters_coverage[seq->cluster_id].emplace_back(seq->coverage[3]);
                //         }
                //     }
                // }
                // for (int kk = 0; kk < C; ++kk) {
                //     for (size_t m = 0; m < clusters_identifier[kk].size(); ++m) {
                //         clstr_fout << kk << "\t";
                //         clstr_fout << clusters_size[kk][m] << "\t";
                //         clstr_fout << clusters_identity[kk][m] << "\t";
                //         clstr_fout << ">" << clusters_identifier[kk][m] << "\t";

                //         int *cov = &clusters_coverage[kk][m * 4];
                //         if (options.print) clstr_fout << cov[0] << "\t" << cov[1] << "\t" << cov[2] << "\t" <<
                //         cov[3]; clstr_fout << "\n";
                //     }
                // }
                // double tb2 = get_time();
                // cerr << "write time " << tb2 - tb1 << endl;
                // clstr_fout.close();

                break;
            }
            if (chunks_id[start] == soure_chunk) {
                start++;
            }
            if (options.stealing) {
                // cerr<<"ptr_ctrl_[1]; "<<ptr_ctrl_[1]<<endl;
                top = start * SUB;
                bottom = sub_chunks.size() - 1;
                ctrl_[0] = top;
                ctrl_[1] = bottom;
                ctrl_[2] = top;
                
                MPI_Put(ctrl_, 3, MPI_INT, worker_rank, 0, 3, MPI_INT, win_ctrl_);
                // // tasks_flag.assign(sub_chunks.size(), 0);
                // // MPI_Put(tasks_flag.data(), sub_chunks.size() * sizeof(int), MPI_INT, worker_rank, 0,
                // //         sub_chunks.size() * sizeof(int), MPI_INT, win_tasks_flag_);

                MPI_Win_flush(worker_rank, win_ctrl_);
                // MPI_Win_flush(worker_rank, win_tasks_flag_);

                // MPI_Barrier(worker_comm);
            }
        }
        for (int i = 0; i < rank_size - 1; i++) {
            if (chunk_kseq[i]) kseq_destroy(chunk_kseq[i]);
            if (chunk_fp[i]) gzclose(chunk_fp[i]);
        }
    }
#pragma omp parallel for schedule(static)
    for (int i = 0; i < NUM_LOCKS; ++i) {
        omp_destroy_lock(&locks[i].lock);
    }
    delete[] locks;
    if (worker && options.stealing)
    {
        // Ensure all workers have finished RMA operations before tearing down windows.
        MPI_Barrier(worker_comm);

        auto cleanup_win = [](MPI_Win &w)
        {
            if (w != MPI_WIN_NULL)
            {
                MPI_Win_flush_all(w);
                MPI_Win_unlock_all(w);
                MPI_Win_free(&w);
            }
        };

        cleanup_win(win_pool_d_);
        cleanup_win(win_meta_);
        cleanup_win(win_ctrl_);
        cleanup_win(win_tasks_);
        cleanup_win(win_tasks_flag_);
        ptr_ctrl_ = nullptr;
    }
    MPI_Barrier(MPI_COMM_WORLD);

    if (master) {
        printf("\n%9lli  finished  %9i  clusters\n", total_num, rep_seqs.size());
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

int SequenceDB::CheckOne_worker(Sequence *seq, WordTable &table, WorkingParam &param, WorkingBuffer &buf,
                                const Options &options, int id) {
    int len = seq->size;
    // cerr<<seq->data<<endl;
    param.len_upper_bound = upper_bound_length_rep(len, options);
    if (options.isEST) return CheckOneEST(seq, table, param, buf, options);
    return CheckOneAA_worker(seq, table, param, buf, options, id);
}

int SequenceDB::CheckOne(Sequence *seq, WordTable &table, WorkingParam &param, WorkingBuffer &buf,
                         const Options &options) {
    int len = seq->size;
    // cerr<<seq->data<<endl;
    param.len_upper_bound = upper_bound_length_rep(len, options);
    if (options.isEST) return CheckOneEST(seq, table, param, buf, options);
    return CheckOneAA(seq, table, param, buf, options);
}
int SequenceDB::CheckOne_stealing(Sequence *seq, WordTable &table, WorkingParam &param, WorkingBuffer &buf,
                         const Options &options) {
    int len = seq->size;
    // cerr<<seq->data<<endl;
    param.len_upper_bound = upper_bound_length_rep(len, options);
    if (options.isEST) return CheckOneEST(seq, table, param, buf, options);
    return CheckOneAA_stealing(seq, table, param, buf, options);
}

int SequenceDB::CheckOneAA_worker(Sequence *seq, WordTable &table, WorkingParam &param, WorkingBuffer &buf,
                                  const Options &options, int id) {
    NVector<IndexCount> &lookCounts = buf.lookCounts;
    NVector<uint32_t> &indexMapping = buf.indexMapping;
    Vector<INTs> &word_encodes_no = buf.word_encodes_no;
    Vector<INTs> &aap_list = buf.aap_list;
    Vector<INTs> &aap_begin = buf.aap_begin;
    Vector<int> &word_encodes = buf.word_encodes;
    Vector<int> &taap = buf.taap;
    double aa1_cutoff = param.aa1_cutoff;
    double aa2_cutoff = param.aas_cutoff;
    double aan_cutoff = param.aan_cutoff;
    // cerr<<seq->data<<endl;
    char *seqi = seq->data;
    int j, k, j1, len = seq->size;
    int flag = 0;
    int frag_size = options.frag_size;
    int &aln_cover_flag = param.aln_cover_flag;
    int &required_aa1 = param.required_aa1;
    int &required_aa2 = param.required_aas;
    int &required_aan = param.required_aan;
    int &min_aln_lenS = param.min_aln_lenS;
    int &min_aln_lenL = param.min_aln_lenL;

    int NAA = options.NAA;
    int S = table.sequences.size();
    int len_eff = len;

    if (S) {
        int min = table.sequences[S - 1]->size;
        if (min < len) {
            if (len * options.diff_cutoff2 > min) min = (int) (len * options.diff_cutoff2);
            if ((len - options.diff_cutoff_aa2) > min) min = len - options.diff_cutoff_aa2;
            len_eff = min;
        }
    }

    // liwz 2016 01, seq is too short for the shortest (longer) seq in word_table to satisfy -aL option
    // longer seqeunce * -aL -band_width
    if (S) {
        int min = table.sequences[S - 1]->size;
        int min_red = min * options.long_coverage - options.band_width;
        if (len < min_red) return 0; // return flag=0
    }

    param.ControlShortCoverage(len_eff, options);
    param.ComputeRequiredBases(options.NAA, 2, options);

    buf.EncodeWords(seq, options.NAA, false);

    // if minimal alignment length > len, return
    // I can not return earlier, because I need to calc the word_encodes etc
    if (options.min_control > len) return 0; // return flag=0

    // lookup_aan
    int aan_no = len - options.NAA + 1;
    int M = frag_size ? table.frag_count : S;
    table.CountWords(aan_no, word_encodes, word_encodes_no, lookCounts, indexMapping, false, required_aan);

    // contained_in_old_lib()
    int len_upper_bound = param.len_upper_bound;
    int len_lower_bound = param.len_lower_bound;
    int band_left, band_right, best_score, band_width1, best_sum, len2, alnln, len_eff1;
    int tiden_no, band_center;
    float tiden_pc, distance = 0;
    int talign_info[5];
    int best1, sum;
    INTs *lookptr;
    char *seqj;
    int frg2 = frag_size ? (len - NAA + options.band_width) / frag_size + 1 + 1 : 0;
    int lens;
    int has_aa2 = 0;

    IndexCount *ic = lookCounts.items;
    ic = lookCounts.items;
    for (; ic->count; ic++) {
        if (!frag_size) {
            indexMapping[ic->index] = 0;
            if (ic->count < required_aan) continue;
        }

        Sequence *rep = table.sequences[ic->index];
        // if(my_rank==3)
        // cerr<<"error   "<<rep->data<<endl;
        len2 = rep->size;
        if (len2 > len_upper_bound || len2 < len) continue;
        if (options.has2D && len2 < len_lower_bound) continue;
        if (frag_size) {
            uint32_t *ims = &indexMapping[ic->index];
            int count = ic->count;
            k = (len2 - NAA) / frag_size + 1;
            sum = 0;
            for (j1 = 0; j1 < frg2; j1++) {
                uint32_t im = ims[j1];
                if (im) sum += lookCounts[im - 1].count;
            }
            count = sum;
            for (j1 = frg2; j1 < k; j1++) {
                uint32_t im1 = ims[j1];
                uint32_t im2 = ims[j1 - frg2];
                if (im1) sum += lookCounts[im1 - 1].count;
                if (im2) sum -= lookCounts[im2 - 1].count;
                if (sum > count) count = sum;
            }
            if (count < required_aan) continue;
        }

        param.ControlLongCoverage(len2, options);

        if (has_aa2 == 0) { // calculate AAP array
            buf.ComputeAAP(seqi, seq->size);
            has_aa2 = 1;
        }
        seqj = rep->data; // NR_seq[NR90_idx[j]];

        band_width1 = (options.band_width < len + len2 - 2) ? options.band_width : len + len2 - 2;

        diag_test_aapn(NAA1, seqj, len, len2, buf, best_sum, band_width1, band_left, band_center, band_right,
                       required_aa1);
        if (best_sum < required_aa2) continue;

        int rc = FAILED_FUNC;
#ifndef NO_AVX512
        if (options.print || aln_cover_flag) // return overlap region
            rc = rotation_band_align_AVX512(seqi, seqj, len, len2, mat, best_score, tiden_no, alnln, distance,
                                            talign_info, band_left, band_center, band_right, buf);
        else
            rc = rotation_band_align_AVX512(seqi, seqj, len, len2, mat, best_score, tiden_no, alnln, distance,
                                            talign_info, band_left, band_center, band_right, buf);
#else
        // auto t0 = std::chrono::high_resolution_clock::now();
        if (options.print || aln_cover_flag) // return overlap region
            rc = local_band_align(seqi, seqj, len, len2, mat, best_score, tiden_no, alnln, distance, talign_info,
                                  band_left, band_center, band_right, buf);
        else
            rc = local_band_align(seqi, seqj, len, len2, mat, best_score, tiden_no, alnln, distance, talign_info,
                                  band_left, band_center, band_right, buf);
#endif

        if (rc == FAILED_FUNC) continue;
        if (tiden_no < required_aa1) continue;
        lens = len;
        if (options.has2D && len > len2) lens = len2;
        len_eff1 = (options.global_identity == 0) ? alnln : (lens - talign_info[4]);
        tiden_pc = tiden_no / (float) len_eff1;
        if (options.useDistance) {
            if (distance > options.distance_thd) continue;
            // if (distance >= seq->distance) continue; // existing distance
            if (distance >= meta_[id].distance) continue; // existing distance
        } else {
            if (tiden_pc < options.cluster_thd) continue;
            // if (tiden_pc <= seq->identity) continue; // existing iden_no
            if (tiden_pc <= meta_[id].identity) continue; // existing iden_no
        }
        if (aln_cover_flag) {
            if (talign_info[3] - talign_info[2] + 1 < min_aln_lenL) continue;
            if (talign_info[1] - talign_info[0] + 1 < min_aln_lenS) continue;
        }
        // if (options.has2D) seq->state |= IS_REDUNDANT;
        if (options.has2D) meta_[id].state = IS_REDUNDANT;
        flag = 1;
        // seq->identity = tiden_pc;
        // seq->cluster_id = rep->cluster_id;
        // seq->distance = distance;
        // seq->coverage[0] = talign_info[0] + 1;
        // seq->coverage[1] = talign_info[1] + 1;
        // seq->coverage[2] = talign_info[2] + 1;
        // seq->coverage[3] = talign_info[3] + 1;
        meta_[id].cluster_id = rep->cluster_id;
        meta_[id].identity = tiden_pc;
        meta_[id].distance = distance;
        for (int t = 0; t < 4; ++t) meta_[id].coverage[t] = talign_info[t] + 1;

        if (not options.cluster_best) break;
        update_aax_cutoff(aa1_cutoff, aa2_cutoff, aan_cutoff, options.tolerance, naa_stat_start_percent, naa_stat, NAA,
                          tiden_pc);
        param.ComputeRequiredBases(options.NAA, 2, options);
    }
    if (frag_size) ic = lookCounts.items;
    while (ic->count) {
        indexMapping[ic->index] = 0;
        ic += 1;
    }
    lookCounts.size = 0;
    if (flag == 1) { // if similar to old one delete it

        if (!options.cluster_best) {
            seq->Clear();
            meta_[id].state = IS_REDUNDANT;
            // seq->state |= IS_REDUNDANT;
        }
    }
    return flag;
}
int SequenceDB::CheckOneAA(Sequence *seq, WordTable &table, WorkingParam &param, WorkingBuffer &buf,
                           const Options &options) {
    NVector<IndexCount> &lookCounts = buf.lookCounts;
    NVector<uint32_t> &indexMapping = buf.indexMapping;
    Vector<INTs> &word_encodes_no = buf.word_encodes_no;
    Vector<INTs> &aap_list = buf.aap_list;
    Vector<INTs> &aap_begin = buf.aap_begin;
    Vector<int> &word_encodes = buf.word_encodes;
    Vector<int> &taap = buf.taap;
    double aa1_cutoff = param.aa1_cutoff;
    double aa2_cutoff = param.aas_cutoff;
    double aan_cutoff = param.aan_cutoff;
    // cerr<<seq->data<<endl;
    char *seqi = seq->data;
    int j, k, j1, len = seq->size;
    int flag = 0;
    int frag_size = options.frag_size;
    int &aln_cover_flag = param.aln_cover_flag;
    int &required_aa1 = param.required_aa1;
    int &required_aa2 = param.required_aas;
    int &required_aan = param.required_aan;
    int &min_aln_lenS = param.min_aln_lenS;
    int &min_aln_lenL = param.min_aln_lenL;

    int NAA = options.NAA;
    int S = table.sequences.size();
    int len_eff = len;

    if (S) {
        int min = table.sequences[S - 1]->size;
        if (min < len) {
            if (len * options.diff_cutoff2 > min) min = (int) (len * options.diff_cutoff2);
            if ((len - options.diff_cutoff_aa2) > min) min = len - options.diff_cutoff_aa2;
            len_eff = min;
        }
    }

    // liwz 2016 01, seq is too short for the shortest (longer) seq in word_table to satisfy -aL option
    // longer seqeunce * -aL -band_width
    if (S) {
        int min = table.sequences[S - 1]->size;
        int min_red = min * options.long_coverage - options.band_width;
        if (len < min_red) return 0; // return flag=0
    }

    param.ControlShortCoverage(len_eff, options);
    param.ComputeRequiredBases(options.NAA, 2, options);

    buf.EncodeWords(seq, options.NAA, false);

    // if minimal alignment length > len, return
    // I can not return earlier, because I need to calc the word_encodes etc
    if (options.min_control > len) return 0; // return flag=0

    // lookup_aan
    int aan_no = len - options.NAA + 1;
    int M = frag_size ? table.frag_count : S;
    table.CountWords(aan_no, word_encodes, word_encodes_no, lookCounts, indexMapping, false, required_aan);

    // contained_in_old_lib()
    int len_upper_bound = param.len_upper_bound;
    int len_lower_bound = param.len_lower_bound;
    int band_left, band_right, best_score, band_width1, best_sum, len2, alnln, len_eff1;
    int tiden_no, band_center;
    float tiden_pc, distance = 0;
    int talign_info[5];
    int best1, sum;
    INTs *lookptr;
    char *seqj;
    int frg2 = frag_size ? (len - NAA + options.band_width) / frag_size + 1 + 1 : 0;
    int lens;
    int has_aa2 = 0;

    IndexCount *ic = lookCounts.items;
    ic = lookCounts.items;
    for (; ic->count; ic++) {
        if (!frag_size) {
            indexMapping[ic->index] = 0;
            if (ic->count < required_aan) continue;
        }

        Sequence *rep = table.sequences[ic->index];
        // if(my_rank==3)
        // cerr<<"error   "<<rep->data<<endl;
        len2 = rep->size;
        if (len2 > len_upper_bound || len2 < len) continue;
        if (options.has2D && len2 < len_lower_bound) continue;
        if (frag_size) {
            uint32_t *ims = &indexMapping[ic->index];
            int count = ic->count;
            k = (len2 - NAA) / frag_size + 1;
            sum = 0;
            for (j1 = 0; j1 < frg2; j1++) {
                uint32_t im = ims[j1];
                if (im) sum += lookCounts[im - 1].count;
            }
            count = sum;
            for (j1 = frg2; j1 < k; j1++) {
                uint32_t im1 = ims[j1];
                uint32_t im2 = ims[j1 - frg2];
                if (im1) sum += lookCounts[im1 - 1].count;
                if (im2) sum -= lookCounts[im2 - 1].count;
                if (sum > count) count = sum;
            }
            if (count < required_aan) continue;
        }

        param.ControlLongCoverage(len2, options);

        if (has_aa2 == 0) { // calculate AAP array
            buf.ComputeAAP(seqi, seq->size);
            has_aa2 = 1;
        }
        seqj = rep->data; // NR_seq[NR90_idx[j]];

        band_width1 = (options.band_width < len + len2 - 2) ? options.band_width : len + len2 - 2;

        diag_test_aapn(NAA1, seqj, len, len2, buf, best_sum, band_width1, band_left, band_center, band_right,
                       required_aa1);
        if (best_sum < required_aa2) continue;

        int rc = FAILED_FUNC;
#ifndef NO_AVX512
        if (options.print || aln_cover_flag) // return overlap region
            rc = rotation_band_align_AVX512(seqi, seqj, len, len2, mat, best_score, tiden_no, alnln, distance,
                                            talign_info, band_left, band_center, band_right, buf);
        else
            rc = rotation_band_align_AVX512(seqi, seqj, len, len2, mat, best_score, tiden_no, alnln, distance,
                                            talign_info, band_left, band_center, band_right, buf);
#else
        // auto t0 = std::chrono::high_resolution_clock::now();
        if (options.print || aln_cover_flag) // return overlap region
            rc = local_band_align(seqi, seqj, len, len2, mat, best_score, tiden_no, alnln, distance, talign_info,
                                  band_left, band_center, band_right, buf);
        else
            rc = local_band_align(seqi, seqj, len, len2, mat, best_score, tiden_no, alnln, distance, talign_info,
                                  band_left, band_center, band_right, buf);
#endif

        if (rc == FAILED_FUNC) continue;
        if (tiden_no < required_aa1) continue;
        lens = len;
        if (options.has2D && len > len2) lens = len2;
        len_eff1 = (options.global_identity == 0) ? alnln : (lens - talign_info[4]);
        tiden_pc = tiden_no / (float) len_eff1;
        if (options.useDistance) {
            if (distance > options.distance_thd) continue;
            if (distance >= seq->distance) continue; // existing distance
        } else {
            if (tiden_pc < options.cluster_thd) continue;
            if (tiden_pc <= seq->identity) continue; // existing iden_no
        }
        if (aln_cover_flag) {
            if (talign_info[3] - talign_info[2] + 1 < min_aln_lenL) continue;
            if (talign_info[1] - talign_info[0] + 1 < min_aln_lenS) continue;
        }
        if (options.has2D) seq->state |= IS_REDUNDANT;
        flag = 1;
        seq->identity = tiden_pc;
        seq->cluster_id = rep->cluster_id;
        seq->distance = distance;
        seq->coverage[0] = talign_info[0] + 1;
        seq->coverage[1] = talign_info[1] + 1;
        seq->coverage[2] = talign_info[2] + 1;
        seq->coverage[3] = talign_info[3] + 1;
        if (not options.cluster_best) break;
        update_aax_cutoff(aa1_cutoff, aa2_cutoff, aan_cutoff, options.tolerance, naa_stat_start_percent, naa_stat, NAA,
                          tiden_pc);
        param.ComputeRequiredBases(options.NAA, 2, options);
    }
    if (frag_size) ic = lookCounts.items;
    while (ic->count) {
        indexMapping[ic->index] = 0;
        ic += 1;
    }
    lookCounts.size = 0;
    if (flag == 1) { // if similar to old one delete it

        if (!options.cluster_best) {
            seq->Clear();
            seq->state |= IS_REDUNDANT;
        }
    }
    return flag;
}
int SequenceDB::CheckOneAA_stealing(Sequence *seq, WordTable &table, WorkingParam &param, WorkingBuffer &buf,
                           const Options &options) {
    NVector<IndexCount> &lookCounts = buf.lookCounts;
    NVector<uint32_t> &indexMapping = buf.indexMapping;
    Vector<INTs> &word_encodes_no = buf.word_encodes_no;
    Vector<INTs> &aap_list = buf.aap_list;
    Vector<INTs> &aap_begin = buf.aap_begin;
    Vector<int> &word_encodes = buf.word_encodes;
    Vector<int> &taap = buf.taap;
    double aa1_cutoff = param.aa1_cutoff;
    double aa2_cutoff = param.aas_cutoff;
    double aan_cutoff = param.aan_cutoff;
    // cerr<<seq->data<<endl;
    char *seqi = seq->data;
    int j, k, j1, len = seq->size;
    int flag = 0;
    int frag_size = options.frag_size;
    int &aln_cover_flag = param.aln_cover_flag;
    int &required_aa1 = param.required_aa1;
    int &required_aa2 = param.required_aas;
    int &required_aan = param.required_aan;
    int &min_aln_lenS = param.min_aln_lenS;
    int &min_aln_lenL = param.min_aln_lenL;

    int NAA = options.NAA;
    int S = table.sequences.size();
    int len_eff = len;

    if (S) {
        int min = table.sequences[S - 1]->size;
        if (min < len) {
            if (len * options.diff_cutoff2 > min) min = (int) (len * options.diff_cutoff2);
            if ((len - options.diff_cutoff_aa2) > min) min = len - options.diff_cutoff_aa2;
            len_eff = min;
        }
    }

    // liwz 2016 01, seq is too short for the shortest (longer) seq in word_table to satisfy -aL option
    // longer seqeunce * -aL -band_width
    if (S) {
        int min = table.sequences[S - 1]->size;
        int min_red = min * options.long_coverage - options.band_width;
        if (len < min_red) return 0; // return flag=0
    }

    param.ControlShortCoverage(len_eff, options);
    param.ComputeRequiredBases(options.NAA, 2, options);

    buf.EncodeWords(seq, options.NAA, false);

    // if minimal alignment length > len, return
    // I can not return earlier, because I need to calc the word_encodes etc
    if (options.min_control > len) return 0; // return flag=0

    // lookup_aan
    int aan_no = len - options.NAA + 1;
    int M = frag_size ? table.frag_count : S;
    table.CountWords(aan_no, word_encodes, word_encodes_no, lookCounts, indexMapping, false, required_aan);

    // contained_in_old_lib()
    int len_upper_bound = param.len_upper_bound;
    int len_lower_bound = param.len_lower_bound;
    int band_left, band_right, best_score, band_width1, best_sum, len2, alnln, len_eff1;
    int tiden_no, band_center;
    float tiden_pc, distance = 0;
    int talign_info[5];
    int best1, sum;
    INTs *lookptr;
    char *seqj;
    int frg2 = frag_size ? (len - NAA + options.band_width) / frag_size + 1 + 1 : 0;
    int lens;
    int has_aa2 = 0;

    IndexCount *ic = lookCounts.items;
    ic = lookCounts.items;
    for (; ic->count; ic++) {
        if (!frag_size) {
            indexMapping[ic->index] = 0;
            if (ic->count < required_aan) continue;
        }

        Sequence *rep = table.sequences[ic->index];
        // if(my_rank==3)
        // cerr<<"error   "<<rep->data<<endl;
        len2 = rep->size;
        if (len2 > len_upper_bound || len2 < len) continue;
        if (options.has2D && len2 < len_lower_bound) continue;
        if (frag_size) {
            uint32_t *ims = &indexMapping[ic->index];
            int count = ic->count;
            k = (len2 - NAA) / frag_size + 1;
            sum = 0;
            for (j1 = 0; j1 < frg2; j1++) {
                uint32_t im = ims[j1];
                if (im) sum += lookCounts[im - 1].count;
            }
            count = sum;
            for (j1 = frg2; j1 < k; j1++) {
                uint32_t im1 = ims[j1];
                uint32_t im2 = ims[j1 - frg2];
                if (im1) sum += lookCounts[im1 - 1].count;
                if (im2) sum -= lookCounts[im2 - 1].count;
                if (sum > count) count = sum;
            }
            if (count < required_aan) continue;
        }

        param.ControlLongCoverage(len2, options);

        if (has_aa2 == 0) { // calculate AAP array
            buf.ComputeAAP(seqi, seq->size);
            has_aa2 = 1;
        }
        seqj = rep->data; // NR_seq[NR90_idx[j]];

        band_width1 = (options.band_width < len + len2 - 2) ? options.band_width : len + len2 - 2;

        diag_test_aapn(NAA1, seqj, len, len2, buf, best_sum, band_width1, band_left, band_center, band_right,
                       required_aa1);
        if (best_sum < required_aa2) continue;

        int rc = FAILED_FUNC;
#ifndef NO_AVX512
        if (options.print || aln_cover_flag) // return overlap region
            rc = rotation_band_align_AVX512(seqi, seqj, len, len2, mat, best_score, tiden_no, alnln, distance,
                                            talign_info, band_left, band_center, band_right, buf);
        else
            rc = rotation_band_align_AVX512(seqi, seqj, len, len2, mat, best_score, tiden_no, alnln, distance,
                                            talign_info, band_left, band_center, band_right, buf);
#else
        // auto t0 = std::chrono::high_resolution_clock::now();
        if (options.print || aln_cover_flag) // return overlap region
            rc = local_band_align(seqi, seqj, len, len2, mat, best_score, tiden_no, alnln, distance, talign_info,
                                  band_left, band_center, band_right, buf);
        else
            rc = local_band_align(seqi, seqj, len, len2, mat, best_score, tiden_no, alnln, distance, talign_info,
                                  band_left, band_center, band_right, buf);
#endif

        if (rc == FAILED_FUNC) continue;
        if (tiden_no < required_aa1) continue;
        lens = len;
        if (options.has2D && len > len2) lens = len2;
        len_eff1 = (options.global_identity == 0) ? alnln : (lens - talign_info[4]);
        tiden_pc = tiden_no / (float) len_eff1;
        if (options.useDistance) {
            if (distance > options.distance_thd) continue;
            if (distance >= seq->distance) continue; // existing distance
        } else {
            if (tiden_pc < options.cluster_thd) continue;
            if (tiden_pc <= seq->identity) continue; // existing iden_no
        }
        if (aln_cover_flag) {
            if (talign_info[3] - talign_info[2] + 1 < min_aln_lenL) continue;
            if (talign_info[1] - talign_info[0] + 1 < min_aln_lenS) continue;
        }
        if (options.has2D) seq->state |= IS_REDUNDANT;
        flag = 1;
        seq->identity = tiden_pc;
        seq->cluster_id = rep->cluster_id;
        seq->distance = distance;
        seq->coverage[0] = talign_info[0] + 1;
        seq->coverage[1] = talign_info[1] + 1;
        seq->coverage[2] = talign_info[2] + 1;
        seq->coverage[3] = talign_info[3] + 1;
        if (not options.cluster_best) break;
        update_aax_cutoff(aa1_cutoff, aa2_cutoff, aan_cutoff, options.tolerance, naa_stat_start_percent, naa_stat, NAA,
                          tiden_pc);
        param.ComputeRequiredBases(options.NAA, 2, options);
    }
    if (frag_size) ic = lookCounts.items;
    while (ic->count) {
        indexMapping[ic->index] = 0;
        ic += 1;
    }
    lookCounts.size = 0;
    if (flag == 1) { // if similar to old one delete it

        if (!options.cluster_best) {
            // seq->Clear();
            seq->state |= IS_REDUNDANT;
        }
    }
    return flag;
}
int SequenceDB::CheckOneEST(Sequence *seq, WordTable &table, WorkingParam &param, WorkingBuffer &buf,
                            const Options &options) {
    NVector<IndexCount> &lookCounts = buf.lookCounts;
    NVector<uint32_t> &indexMapping = buf.indexMapping;
    Vector<INTs> &word_encodes_no = buf.word_encodes_no;
    Vector<INTs> &aap_list = buf.aap_list;
    Vector<INTs> &aap_begin = buf.aap_begin;
    Vector<int> &word_encodes = buf.word_encodes;
    Vector<int> &taap = buf.taap;
    Vector<int> &aan_list_comp = buf.aan_list_comp;
    char *seqi_comp = &buf.seqi_comp[0];

    int &aln_cover_flag = param.aln_cover_flag;
    int &required_aa1 = param.required_aa1;
    int &required_aas = param.required_aas;
    int &required_aan = param.required_aan;
    int &min_aln_lenS = param.min_aln_lenS;
    int &min_aln_lenL = param.min_aln_lenL;

    char *seqi = seq->data;
    int j, len = seq->size;
    int flag = 0;
    int S = table.sequences.size();
    int len_eff = len;
    if (S) {
        int min = table.sequences[S - 1]->size;
        if (min < len) {
            if (len * options.diff_cutoff2 > min) min = (int) (len * options.diff_cutoff2);
            if ((len - options.diff_cutoff_aa2) > min) min = len - options.diff_cutoff_aa2;
            len_eff = min;
        }
    }

    // liwz 2016 01, seq is too short for the shortest (longer) seq in word_table to satisfy -aL option
    // longer seqeunce * -aL -band_width
    if (S) {
        int min = table.sequences[S - 1]->size;
        int min_red = min * options.long_coverage - options.band_width;
        if (len < min_red) return 0; // return flag=0
    }

    param.ControlShortCoverage(len_eff, options);
    param.ComputeRequiredBases(options.NAA, 4, options);
    int skip = buf.EncodeWords(seq, options.NAA, true);
    required_aan -= skip;
    required_aas -= skip;
    required_aa1 -= skip;
    if (required_aan <= 0) required_aan = 1;
    if (required_aas <= 0) required_aas = 1;
    if (required_aa1 <= 0) required_aa1 = 1;

    // if minimal alignment length > len, return
    // I can not return earlier, because I need to calc the word_encodes etc
    if (options.min_control > len) return 0; // return flag=0

    int aan_no = len - options.NAA + 1;

    // contained_in_old_lib()
    int len_upper_bound = param.len_upper_bound;
    int len_lower_bound = param.len_lower_bound;
    int band_left, band_right, best_score, band_width1, best_sum, len2, alnln, len_eff1;
    int tiden_no, band_center;
    float tiden_pc, distance = 0;
    int talign_info[5];
    int j0, comp, lens;
    char *seqj;

    for (comp = 0; comp < 2; comp++) {
        if (comp) {
            for (j0 = 0; j0 < aan_no; j0++) {
                j = word_encodes[j0];
                if (j < 0)
                    aan_list_comp[j0] = j;
                else
                    aan_list_comp[j0] = Comp_AAN_idx[j];
            }
            make_comp_iseq(len, seqi_comp, seqi);
            seqi = seqi_comp;
        }
        int has_aas = 0;

        if (comp) {
            table.CountWords(aan_no, aan_list_comp, word_encodes_no, lookCounts, indexMapping, true, required_aan);
        } else {
            table.CountWords(aan_no, word_encodes, word_encodes_no, lookCounts, indexMapping, true, required_aan);
        }

        IndexCount *ic = lookCounts.items;
        ic = lookCounts.items;
        for (; ic->count; ic++) {
            indexMapping[ic->index] = 0;
            if (ic->count < required_aan) continue;
            Sequence *rep = table.sequences[ic->index];

            len2 = rep->size;
            if (len2 > len_upper_bound) continue;
            if (options.has2D && len2 < len_lower_bound) continue;

            seqj = rep->data;

            param.ControlLongCoverage(len2, options);

            if (has_aas == 0) { // calculate AAP array
                buf.ComputeAAP2(seqi, seq->size);
                has_aas = 1;
            }

            band_width1 = (options.band_width < len + len2 - 2) ? options.band_width : len + len2 - 2;
            diag_test_aapn_est(NAA1, seqj, len, len2, buf, best_sum, band_width1, band_left, band_center, band_right,
                               required_aa1);
            if (best_sum < required_aas) continue;
            // if( comp and flag and (not options.cluster_best) and j > rep->cluster_id ) goto Break;

            int rc = FAILED_FUNC;
            if (options.print || aln_cover_flag) { // return overlap region
                rc = local_band_align(seqi, seqj, len, len2, mat, best_score, tiden_no, alnln, distance, talign_info,
                                      band_left, band_center, band_right, buf);
                if (comp) {
                    talign_info[0] = len - talign_info[0] - 1;
                    talign_info[1] = len - talign_info[1] - 1;
                }
            } else {
                // printf( "%5i %5i %5i %5i\n", band_width1, band_right-band_left, band_left, band_right );
                rc = local_band_align(seqi, seqj, len, len2, mat, best_score, tiden_no, alnln, distance, talign_info,
                                      band_left, band_center, band_right, buf);
            }
            if (rc == FAILED_FUNC) continue;
            // printf( "%i  %i  %i\n", best_score, tiden_no, required_aa1 );
            if (tiden_no < required_aa1) continue;
            if (options.is454) {
                if (talign_info[2] != talign_info[0]) continue; // same start
                if (talign_info[0] > 1) continue;               // one mismatch allowed at beginning
                if ((len - talign_info[1]) > 2) continue;       // one mismatch allowed at end
            }

            lens = len;
            if (options.has2D && len > len2) lens = len2;
            len_eff1 = (options.global_identity == 0) ? alnln : (lens - talign_info[4]);
            tiden_pc = tiden_no / (float) len_eff1;
            // printf( "%i %f\n", tiden_no, tiden_pc );
            if (options.useDistance) {
                if (distance > options.distance_thd) continue;
                if (options.cluster_best and distance >= seq->distance) continue; // existing distance
            } else {
                if (tiden_pc < options.cluster_thd) continue;
                if (options.cluster_best and tiden_pc < seq->identity) continue; // existing iden_no
            }
            if (aln_cover_flag) {
                if (talign_info[3] - talign_info[2] + 1 < min_aln_lenL) continue;
                if (comp) {
                    if (talign_info[0] - talign_info[1] + 1 < min_aln_lenS) continue;
                } else {
                    if (talign_info[1] - talign_info[0] + 1 < min_aln_lenS) continue;
                }
            }
            if (options.cluster_best and fabs(tiden_pc - seq->identity) < 1E-9 and rep->cluster_id >= seq->cluster_id)
                continue;
            if ((not options.cluster_best) and flag != 0 and rep->cluster_id >= seq->cluster_id) continue;
            flag = comp ? -1 : 1;
            seq->identity = tiden_pc;
            seq->distance = distance;
            seq->cluster_id = rep->cluster_id;
            seq->coverage[0] = talign_info[0] + 1;
            seq->coverage[1] = talign_info[1] + 1;
            seq->coverage[2] = talign_info[2] + 1;
            seq->coverage[3] = talign_info[3] + 1;
            if (not options.cluster_best) break;
        }
        while (ic->count) {
            indexMapping[ic->index] = 0;
            ic += 1;
        }
        lookCounts.size = 0;
        if (not options.option_r) break;
    }
    if ((flag == 1) || (flag == -1)) { // if similar to old one delete it
        if (!options.cluster_best) {
            seq->Clear();
            seq->state |= IS_REDUNDANT;
        }
        if (flag == -1)
            seq->state |= IS_MINUS_STRAND;
        else
            seq->state &= ~IS_MINUS_STRAND;
    }
    return flag;
}
void SequenceDB::ComputeDistance(const Options &options) {
    int i, j, N = sequences.size();
    int best_score, best_sum;
    int band_width1, band_left, band_center, band_right, required_aa1;
    int tiden_no, alnln;
    int talign_info[5];
    float distance;
    WorkingBuffer buf(N, max_len, options);

    Vector<NVector<float>> dists(N, NVector<float>(N));

    Sequence comseq(*sequences[0]);

    for (i = 0; i < N; i++) {
        Sequence *seq = sequences[i];
        char *seqi = seq->data;
        int len = seq->size;
        buf.EncodeWords(seq, options.NAA, false);
        buf.ComputeAAP2(seqi, seq->size);
        dists[i][i] = 0.0;
        if ((i + 1) % 1000 == 0) printf("%9i\n", (i + 1));
        for (j = 0; j < i; j++) {
            Sequence *rep = sequences[j];
            char *seqj = rep->data;
            int len2 = rep->size;
            band_width1 = (options.band_width < len + len2 - 2) ? options.band_width : len + len2 - 2;
            diag_test_aapn_est(NAA1, seqj, len, len2, buf, best_sum, band_width1, band_left, band_center, band_right,
                               0);
            local_band_align(seqi, seqj, len, len2, mat, best_score, tiden_no, alnln, distance, talign_info, band_left,
                             band_center, band_right, buf);
            dists[seq->index][rep->index] = dists[rep->index][seq->index] = distance;
        }
        if (not options.option_r) break;
        comseq.index = seq->index;
        comseq.size = len;
        for (j = 0; j < len; j++) comseq.data[i] = seq->data[len - i - 1];
        seqi = comseq.data;
        buf.EncodeWords(&comseq, options.NAA, false);
        buf.ComputeAAP2(seqi, seq->size);
        for (j = 0; j < i; j++) {
            Sequence *rep = sequences[j];
            char *seqj = rep->data;
            int len2 = rep->size;
            band_width1 = (options.band_width < len + len2 - 2) ? options.band_width : len + len2 - 2;
            diag_test_aapn_est(NAA1, seqj, len, len2, buf, best_sum, band_width1, band_left, band_center, band_right,
                               0);
            local_band_align(seqi, seqj, len, len2, mat, best_score, tiden_no, alnln, distance, talign_info, band_left,
                             band_center, band_right, buf);
            if (distance < dists[seq->index][rep->index])
                dists[seq->index][rep->index] = dists[rep->index][seq->index] = distance;
        }
    }
    std::string output = options.output + ".dist";
    FILE *fout = fopen(output.c_str(), "w+");
    fprintf(fout, "1");
    for (i = 1; i < N; i++) fprintf(fout, "\t%i", i + 1);
    fprintf(fout, "\n");
    for (i = 0; i < N; i++) {
        fprintf(fout, "%g", dists[i][0]);
        for (j = 1; j < N; j++) fprintf(fout, "\t%g", dists[i][j]);
        fprintf(fout, "\n");
    }
    fclose(fout);
}

int calc_ann_list(int len, char *seqi, int NAA, int &aan_no, Vector<int> &aan_list, Vector<INTs> &aan_list_no,
                  bool est) {
    int i, j, k, i0, i1, k1;

    // check_aan_list
    aan_no = len - NAA + 1;
    for (j = 0; j < aan_no; j++) {
        aan_list[j] = 0;
        for (k = 0, k1 = NAA - 1; k < NAA; k++, k1--) aan_list[j] += seqi[j + k] * NAAN_array[k1];
    }
    if (est) {
        // for the short word containing 'N', mask it to '-1'
        for (j = 0; j < len; j++) {
            if (seqi[j] >= 4) { // here N is 4
                i0 = (j - NAA + 1 > 0) ? j - NAA + 1 : 0;
                i1 = j < aan_no ? j : aan_no - 1;
                for (i = i0; i <= i1; i++) aan_list[i] = -1;
            }
        }
    }

    std::sort(aan_list.begin(), aan_list.begin() + aan_no);
    for (j = 0; j < aan_no; j++) aan_list_no[j] = 1;
    for (j = aan_no - 1; j; j--) {
        if (aan_list[j] == aan_list[j - 1]) {
            aan_list_no[j - 1] += aan_list_no[j];
            aan_list_no[j] = 0;
        }
    }
    return OK_FUNC;
} // END calc_ann_list

void make_comp_short_word_index(int NAA, int *NAAN_array, Vector<int> &Comp_AAN_idx) {
    int i, j, k, icomp, k1;
    int c[4] = {3, 2, 1, 0};
    unsigned char short_word[32]; // short_word[12] is enough

    int NAA1 = NAAN_array[1];
    int NAAN = NAAN_array[NAA];

    for (i = 0; i < NAAN; i++) {
        // decompose i back to short_word
        for (k = i, j = 0; j < NAA; j++) {
            short_word[j] = (unsigned char) (k % NAA1);
            k = k / NAA1;
        }

        // calc_comp_aan_list
        icomp = 0;
        for (k = 0, k1 = NAA - 1; k < NAA; k++, k1--) icomp += c[short_word[k1]] * NAAN_array[k];

        Comp_AAN_idx[i] = icomp;
    }
} // make_comp_short_word_index

// quick_sort_idx calling (a, idx, 0, no-1)
// sort a with another array idx
// so that idx rearranged
int quick_sort_idx(int *a, int *idx, int lo0, int hi0) {
    int lo = lo0;
    int hi = hi0;
    int mid;
    int tmp;

    if (hi0 > lo0) {
        mid = a[(lo0 + hi0) / 2];

        while (lo <= hi) {
            while ((lo < hi0) && (a[lo] < mid)) lo++;
            while ((hi > lo0) && (a[hi] > mid)) hi--;
            if (lo <= hi) {
                tmp = a[lo];
                a[lo] = a[hi];
                a[hi] = tmp;
                tmp = idx[lo];
                idx[lo] = idx[hi];
                idx[hi] = tmp;
                lo++;
                hi--;
            }
        } // while

        if (lo0 < hi) quick_sort_idx(a, idx, lo0, hi);
        if (lo < hi0) quick_sort_idx(a, idx, lo, hi0);
    } // if ( hi0 > lo0)
    return 0;
} // quick_sort_idx

// decreasing can not use reverse of quick_sort_idx due to tie
// quick_sort_idxr calling (a, idx, 0, no-1)
// sort a with another array idx
// so that idx rearranged
int quick_sort_idxr(int *a, int *idx, int lo0, int hi0) {
    int lo = lo0;
    int hi = hi0;
    int mid;
    int tmp;

    if (hi0 > lo0) {
        mid = a[(lo0 + hi0) / 2];

        while (lo <= hi) {
            while ((lo < hi0) && (a[lo] > mid)) lo++;
            while ((hi > lo0) && (a[hi] < mid)) hi--;
            if (lo <= hi) {
                tmp = a[lo];
                a[lo] = a[hi];
                a[hi] = tmp;
                tmp = idx[lo];
                idx[lo] = idx[hi];
                idx[hi] = tmp;
                lo++;
                hi--;
            }
        } // while

        if (lo0 < hi) quick_sort_idxr(a, idx, lo0, hi);
        if (lo < hi0) quick_sort_idxr(a, idx, lo, hi0);
    } // if ( hi0 > lo0)
    return 0;
} // quick_sort_idxr

/////////////////////////// END ALL ////////////////////////

int naa_stat_start_percent = 40;
int naa_stat[5][61][4] = {

    // cover 0.99
    {
        // N=5   N=4   N=3   N=2
        {
            0,
            0,
            0,
            7,
        }, // 40%
        {
            0,
            0,
            0,
            8,
        }, // 41%
        {
            0,
            0,
            0,
            9,
        }, // 42%
        {
            0,
            0,
            0,
            9,
        }, // 43%
        {
            0,
            0,
            1,
            10,
        }, // 44%
        {
            0,
            0,
            1,
            11,
        }, // 45%
        {
            0,
            0,
            1,
            12,
        }, // 46%
        {
            0,
            0,
            2,
            13,
        }, // 47%
        {
            0,
            0,
            2,
            14,
        }, // 48%
        {
            0,
            0,
            4,
            16,
        }, // 49%
        {
            0,
            0,
            4,
            16,
        }, // 50%
        {
            0,
            0,
            5,
            17,
        }, // 51%
        {
            0,
            0,
            5,
            18,
        }, // 52%
        {
            0,
            0,
            7,
            20,
        }, // 53%
        {
            0,
            1,
            7,
            21,
        }, // 54%
        {
            0,
            1,
            7,
            21,
        }, // 55%
        {
            0,
            2,
            8,
            23,
        }, // 56%
        {
            0,
            2,
            8,
            25,
        }, // 57%
        {
            0,
            2,
            10,
            25,
        }, // 58%
        {
            0,
            3,
            10,
            26,
        }, // 59%
        {
            0,
            4,
            13,
            28,
        }, // 60%
        {
            0,
            5,
            13,
            30,
        }, // 61%
        {
            0,
            5,
            14,
            30,
        }, // 62%
        {
            1,
            6,
            15,
            33,
        }, // 63%
        {
            2,
            7,
            17,
            34,
        }, // 64%
        {
            2,
            7,
            17,
            35,
        }, // 65%
        {
            2,
            9,
            20,
            37,
        }, // 66%
        {
            4,
            10,
            20,
            37,
        }, // 67%
        {
            4,
            11,
            22,
            40,
        }, // 68%
        {
            5,
            12,
            24,
            41,
        }, // 69%
        {
            5,
            12,
            25,
            42,
        }, // 70%
        {
            6,
            16,
            27,
            43,
        }, // 71%
        {
            8,
            16,
            27,
            45,
        }, // 72%
        {
            9,
            17,
            29,
            47,
        }, // 73%
        {
            10,
            18,
            31,
            47,
        }, // 74%
        {
            10,
            20,
            32,
            50,
        }, // 75%
        {
            12,
            20,
            32,
            51,
        }, // 76%
        {
            14,
            22,
            36,
            54,
        }, // 77%
        {
            15,
            24,
            37,
            55,
        }, // 78%
        {
            17,
            26,
            41,
            58,
        }, // 79%
        {
            18,
            29,
            41,
            59,
        }, // 80%
        {
            20,
            30,
            45,
            60,
        }, // 81%
        {
            24,
            35,
            48,
            62,
        }, // 82%
        {
            26,
            36,
            48,
            64,
        }, // 83%
        {
            27,
            38,
            51,
            65,
        }, // 84%
        {
            31,
            43,
            54,
            68,
        }, // 85%
        {
            35,
            43,
            55,
            70,
        }, // 86%
        {
            36,
            48,
            60,
            71,
        }, // 87%
        {
            36,
            50,
            61,
            73,
        }, // 88%
        {
            40,
            50,
            61,
            75,
        }, // 89%
        {
            45,
            54,
            65,
            75,
        }, // 90%
        {
            52,
            60,
            70,
            79,
        }, // 91%
        {
            53,
            62,
            71,
            81,
        }, // 92%
        {
            57,
            66,
            75,
            84,
        }, // 93%
        {
            57,
            66,
            76,
            85,
        }, // 94%
        {
            64,
            71,
            78,
            85,
        }, // 95%
        {
            70,
            75,
            82,
            89,
        }, // 96%
        {
            77,
            81,
            86,
            92,
        }, // 97%
        {
            82,
            86,
            90,
            94,
        }, // 98%
        {
            83,
            87,
            91,
            95,
        }, // 99%
        {
            91,
            93,
            95,
            97,
        }, // 100%
    },
    // cover 0.95
    {
        // N=5   N=4   N=3   N=2
        {
            0,
            0,
            1,
            9,
        }, // 40%
        {
            0,
            0,
            2,
            10,
        }, // 41%
        {
            0,
            0,
            2,
            11,
        }, // 42%
        {
            0,
            0,
            3,
            12,
        }, // 43%
        {
            0,
            0,
            3,
            12,
        }, // 44%
        {
            0,
            0,
            4,
            14,
        }, // 45%
        {
            0,
            0,
            4,
            14,
        }, // 46%
        {
            0,
            1,
            5,
            16,
        }, // 47%
        {
            0,
            1,
            6,
            17,
        }, // 48%
        {
            0,
            2,
            7,
            19,
        }, // 49%
        {
            0,
            2,
            8,
            19,
        }, // 50%
        {
            0,
            2,
            8,
            20,
        }, // 51%
        {
            0,
            2,
            9,
            21,
        }, // 52%
        {
            0,
            4,
            10,
            23,
        }, // 53%
        {
            1,
            4,
            11,
            24,
        }, // 54%
        {
            1,
            4,
            11,
            24,
        }, // 55%
        {
            1,
            5,
            13,
            26,
        }, // 56%
        {
            2,
            5,
            13,
            27,
        }, // 57%
        {
            2,
            6,
            15,
            29,
        }, // 58%
        {
            2,
            7,
            15,
            30,
        }, // 59%
        {
            3,
            8,
            16,
            31,
        }, // 60%
        {
            4,
            8,
            18,
            32,
        }, // 61%
        {
            4,
            9,
            18,
            33,
        }, // 62%
        {
            5,
            11,
            20,
            36,
        }, // 63%
        {
            6,
            12,
            22,
            37,
        }, // 64%
        {
            6,
            12,
            22,
            38,
        }, // 65%
        {
            8,
            14,
            24,
            40,
        }, // 66%
        {
            8,
            15,
            25,
            41,
        }, // 67%
        {
            10,
            16,
            27,
            42,
        }, // 68%
        {
            10,
            18,
            28,
            45,
        }, // 69%
        {
            11,
            18,
            29,
            45,
        }, // 70%
        {
            14,
            21,
            31,
            47,
        }, // 71%
        {
            14,
            22,
            32,
            48,
        }, // 72%
        {
            14,
            22,
            33,
            50,
        }, // 73%
        {
            17,
            24,
            36,
            52,
        }, // 74%
        {
            17,
            25,
            36,
            52,
        }, // 75%
        {
            18,
            27,
            39,
            54,
        }, // 76%
        {
            20,
            29,
            41,
            56,
        }, // 77%
        {
            21,
            31,
            42,
            58,
        }, // 78%
        {
            21,
            31,
            46,
            60,
        }, // 79%
        {
            27,
            35,
            46,
            60,
        }, // 80%
        {
            28,
            37,
            50,
            63,
        }, // 81%
        {
            31,
            38,
            50,
            64,
        }, // 82%
        {
            34,
            43,
            53,
            66,
        }, // 83%
        {
            36,
            45,
            54,
            67,
        }, // 84%
        {
            41,
            50,
            60,
            70,
        }, // 85%
        {
            43,
            51,
            60,
            71,
        }, // 86%
        {
            45,
            54,
            63,
            74,
        }, // 87%
        {
            48,
            55,
            64,
            75,
        }, // 88%
        {
            54,
            60,
            68,
            78,
        }, // 89%
        {
            55,
            62,
            71,
            80,
        }, // 90%
        {
            56,
            63,
            71,
            80,
        }, // 91%
        {
            64,
            70,
            76,
            84,
        }, // 92%
        {
            69,
            74,
            80,
            86,
        }, // 93%
        {
            73,
            78,
            83,
            88,
        }, // 94%
        {
            74,
            78,
            84,
            89,
        }, // 95%
        {
            80,
            84,
            87,
            91,
        }, // 96%
        {
            83,
            86,
            90,
            93,
        }, // 97%
        {
            86,
            89,
            92,
            95,
        }, // 98%
        {
            91,
            93,
            95,
            97,
        }, // 99%
        {
            92,
            93,
            95,
            97,
        }, // 100%
    },
    // cover 0.9
    {
        // N=5   N=4   N=3   N=2
        {
            0,
            0,
            2,
            11,
        }, // 40%
        {
            0,
            0,
            3,
            12,
        }, // 41%
        {
            0,
            0,
            3,
            12,
        }, // 42%
        {
            0,
            1,
            4,
            13,
        }, // 43%
        {
            0,
            1,
            5,
            14,
        }, // 44%
        {
            0,
            1,
            5,
            15,
        }, // 45%
        {
            0,
            1,
            6,
            16,
        }, // 46%
        {
            0,
            2,
            7,
            18,
        }, // 47%
        {
            0,
            2,
            7,
            18,
        }, // 48%
        {
            0,
            3,
            9,
            20,
        }, // 49%
        {
            1,
            4,
            9,
            20,
        }, // 50%
        {
            1,
            4,
            10,
            21,
        }, // 51%
        {
            1,
            4,
            11,
            23,
        }, // 52%
        {
            2,
            5,
            12,
            24,
        }, // 53%
        {
            2,
            5,
            12,
            25,
        }, // 54%
        {
            2,
            6,
            13,
            26,
        }, // 55%
        {
            3,
            7,
            14,
            28,
        }, // 56%
        {
            3,
            7,
            15,
            28,
        }, // 57%
        {
            4,
            8,
            16,
            30,
        }, // 58%
        {
            5,
            9,
            17,
            31,
        }, // 59%
        {
            5,
            10,
            18,
            32,
        }, // 60%
        {
            6,
            11,
            20,
            35,
        }, // 61%
        {
            6,
            11,
            20,
            35,
        }, // 62%
        {
            7,
            13,
            22,
            38,
        }, // 63%
        {
            8,
            14,
            23,
            39,
        }, // 64%
        {
            8,
            15,
            24,
            39,
        }, // 65%
        {
            10,
            16,
            26,
            42,
        }, // 66%
        {
            10,
            17,
            27,
            42,
        }, // 67%
        {
            12,
            19,
            29,
            44,
        }, // 68%
        {
            13,
            20,
            30,
            46,
        }, // 69%
        {
            13,
            21,
            31,
            47,
        }, // 70%
        {
            16,
            23,
            33,
            48,
        }, // 71%
        {
            18,
            25,
            34,
            50,
        }, // 72%
        {
            18,
            26,
            36,
            51,
        }, // 73%
        {
            19,
            28,
            38,
            53,
        }, // 74%
        {
            20,
            29,
            38,
            53,
        }, // 75%
        {
            23,
            30,
            41,
            56,
        }, // 76%
        {
            24,
            33,
            43,
            57,
        }, // 77%
        {
            26,
            34,
            45,
            59,
        }, // 78%
        {
            28,
            37,
            48,
            61,
        }, // 79%
        {
            30,
            37,
            48,
            62,
        }, // 80%
        {
            33,
            42,
            52,
            64,
        }, // 81%
        {
            35,
            43,
            53,
            65,
        }, // 82%
        {
            38,
            47,
            56,
            68,
        }, // 83%
        {
            40,
            47,
            56,
            68,
        }, // 84%
        {
            44,
            53,
            61,
            71,
        }, // 85%
        {
            45,
            53,
            62,
            73,
        }, // 86%
        {
            50,
            58,
            66,
            75,
        }, // 87%
        {
            51,
            58,
            66,
            76,
        }, // 88%
        {
            57,
            63,
            71,
            79,
        }, // 89%
        {
            60,
            66,
            72,
            81,
        }, // 90%
        {
            62,
            68,
            75,
            83,
        }, // 91%
        {
            70,
            74,
            80,
            85,
        }, // 92%
        {
            74,
            78,
            82,
            88,
        }, // 93%
        {
            85,
            87,
            90,
            92,
        }, // 94%
        {
            86,
            88,
            90,
            92,
        }, // 95%
        {
            87,
            89,
            91,
            93,
        }, // 96%
        {
            87,
            89,
            92,
            94,
        }, // 97%
        {
            89,
            91,
            93,
            96,
        }, // 98%
        {
            93,
            94,
            96,
            97,
        }, // 99%
        {
            94,
            95,
            97,
            98,
        }, // 100%
    },
    // cover 0.8
    {
        // N=5   N=4   N=3   N=2
        {
            0,
            1,
            4,
            13,
        }, // 40%
        {
            0,
            1,
            5,
            13,
        }, // 41%
        {
            0,
            1,
            5,
            14,
        }, // 42%
        {
            0,
            2,
            6,
            15,
        }, // 43%
        {
            0,
            2,
            6,
            16,
        }, // 44%
        {
            0,
            2,
            7,
            17,
        }, // 45%
        {
            1,
            3,
            8,
            18,
        }, // 46%
        {
            1,
            4,
            9,
            20,
        }, // 47%
        {
            1,
            4,
            9,
            20,
        }, // 48%
        {
            2,
            5,
            11,
            22,
        }, // 49%
        {
            2,
            5,
            11,
            22,
        }, // 50%
        {
            2,
            6,
            12,
            24,
        }, // 51%
        {
            3,
            6,
            13,
            25,
        }, // 52%
        {
            3,
            7,
            14,
            26,
        }, // 53%
        {
            4,
            8,
            14,
            27,
        }, // 54%
        {
            4,
            8,
            15,
            28,
        }, // 55%
        {
            5,
            9,
            17,
            30,
        }, // 56%
        {
            5,
            9,
            17,
            30,
        }, // 57%
        {
            6,
            11,
            19,
            32,
        }, // 58%
        {
            7,
            12,
            20,
            34,
        }, // 59%
        {
            8,
            12,
            20,
            34,
        }, // 60%
        {
            9,
            14,
            22,
            37,
        }, // 61%
        {
            9,
            14,
            23,
            37,
        }, // 62%
        {
            10,
            16,
            25,
            39,
        }, // 63%
        {
            11,
            17,
            26,
            41,
        }, // 64%
        {
            12,
            18,
            27,
            41,
        }, // 65%
        {
            13,
            20,
            28,
            43,
        }, // 66%
        {
            14,
            21,
            30,
            45,
        }, // 67%
        {
            15,
            22,
            31,
            46,
        }, // 68%
        {
            17,
            24,
            33,
            48,
        }, // 69%
        {
            17,
            24,
            34,
            48,
        }, // 70%
        {
            19,
            26,
            36,
            50,
        }, // 71%
        {
            20,
            27,
            37,
            51,
        }, // 72%
        {
            21,
            29,
            39,
            53,
        }, // 73%
        {
            23,
            31,
            41,
            55,
        }, // 74%
        {
            23,
            31,
            41,
            55,
        }, // 75%
        {
            26,
            34,
            44,
            58,
        }, // 76%
        {
            28,
            36,
            46,
            59,
        }, // 77%
        {
            29,
            37,
            47,
            60,
        }, // 78%
        {
            34,
            41,
            50,
            62,
        }, // 79%
        {
            34,
            42,
            51,
            63,
        }, // 80%
        {
            38,
            45,
            55,
            66,
        }, // 81%
        {
            39,
            46,
            55,
            67,
        }, // 82%
        {
            44,
            51,
            60,
            70,
        }, // 83%
        {
            44,
            51,
            60,
            70,
        }, // 84%
        {
            49,
            56,
            64,
            73,
        }, // 85%
        {
            50,
            57,
            64,
            74,
        }, // 86%
        {
            57,
            63,
            69,
            77,
        }, // 87%
        {
            58,
            64,
            70,
            78,
        }, // 88%
        {
            68,
            71,
            76,
            82,
        }, // 89%
        {
            68,
            72,
            77,
            83,
        }, // 90%
        {
            75,
            79,
            81,
            85,
        }, // 91%
        {
            86,
            87,
            89,
            90,
        }, // 92%
        {
            88,
            89,
            90,
            92,
        }, // 93%
        {
            90,
            91,
            92,
            93,
        }, // 94%
        {
            91,
            92,
            93,
            94,
        }, // 95%
        {
            92,
            94,
            94,
            95,
        }, // 96%
        {
            93,
            94,
            95,
            96,
        }, // 97%
        {
            94,
            95,
            95,
            96,
        }, // 98%
        {
            94,
            95,
            96,
            98,
        }, // 99%
        {
            95,
            96,
            97,
            98,
        }, // 100%
    },
    // cover 0.6
    {
        // N=5   N=4   N=3   N=2
        {
            1,
            2,
            6,
            15,
        }, // 40%
        {
            1,
            3,
            7,
            16,
        }, // 41%
        {
            1,
            3,
            8,
            17,
        }, // 42%
        {
            2,
            4,
            9,
            18,
        }, // 43%
        {
            2,
            4,
            9,
            19,
        }, // 44%
        {
            2,
            5,
            10,
            20,
        }, // 45%
        {
            3,
            5,
            10,
            21,
        }, // 46%
        {
            3,
            6,
            12,
            22,
        }, // 47%
        {
            3,
            6,
            12,
            23,
        }, // 48%
        {
            4,
            8,
            14,
            25,
        }, // 49%
        {
            4,
            8,
            14,
            25,
        }, // 50%
        {
            5,
            8,
            15,
            26,
        }, // 51%
        {
            5,
            9,
            16,
            27,
        }, // 52%
        {
            6,
            10,
            17,
            29,
        }, // 53%
        {
            6,
            11,
            18,
            30,
        }, // 54%
        {
            7,
            11,
            18,
            31,
        }, // 55%
        {
            8,
            12,
            20,
            32,
        }, // 56%
        {
            8,
            13,
            20,
            33,
        }, // 57%
        {
            10,
            14,
            22,
            35,
        }, // 58%
        {
            10,
            15,
            23,
            37,
        }, // 59%
        {
            11,
            16,
            24,
            37,
        }, // 60%
        {
            12,
            18,
            26,
            39,
        }, // 61%
        {
            13,
            18,
            26,
            40,
        }, // 62%
        {
            14,
            20,
            28,
            42,
        }, // 63%
        {
            16,
            22,
            30,
            43,
        }, // 64%
        {
            16,
            22,
            31,
            44,
        }, // 65%
        {
            17,
            23,
            32,
            45,
        }, // 66%
        {
            18,
            25,
            33,
            47,
        }, // 67%
        {
            19,
            26,
            35,
            48,
        }, // 68%
        {
            21,
            27,
            36,
            50,
        }, // 69%
        {
            22,
            29,
            37,
            51,
        }, // 70%
        {
            24,
            30,
            39,
            52,
        }, // 71%
        {
            25,
            32,
            41,
            53,
        }, // 72%
        {
            26,
            33,
            42,
            55,
        }, // 73%
        {
            29,
            35,
            44,
            57,
        }, // 74%
        {
            29,
            36,
            45,
            57,
        }, // 75%
        {
            32,
            39,
            48,
            60,
        }, // 76%
        {
            34,
            41,
            50,
            61,
        }, // 77%
        {
            36,
            43,
            51,
            62,
        }, // 78%
        {
            40,
            46,
            54,
            65,
        }, // 79%
        {
            40,
            46,
            54,
            65,
        }, // 80%
        {
            46,
            52,
            59,
            68,
        }, // 81%
        {
            46,
            52,
            60,
            69,
        }, // 82%
        {
            53,
            59,
            65,
            73,
        }, // 83%
        {
            54,
            60,
            66,
            73,
        }, // 84%
        {
            63,
            67,
            73,
            78,
        }, // 85%
        {
            68,
            71,
            75,
            79,
        }, // 86%
        {
            78,
            80,
            82,
            85,
        }, // 87%
        {
            79,
            81,
            83,
            85,
        }, // 88%
        {
            83,
            85,
            86,
            87,
        }, // 89%
        {
            85,
            86,
            87,
            89,
        }, // 90%
        {
            86,
            88,
            89,
            90,
        }, // 91%
        {
            88,
            89,
            90,
            91,
        }, // 92%
        {
            90,
            90,
            91,
            92,
        }, // 93%
        {
            91,
            92,
            92,
            93,
        }, // 94%
        {
            92,
            93,
            94,
            94,
        }, // 95%
        {
            94,
            94,
            95,
            95,
        }, // 96%
        {
            95,
            95,
            96,
            96,
        }, // 97%
        {
            95,
            96,
            97,
            97,
        }, // 98%
        {
            96,
            96,
            97,
            98,
        }, // 99%
        {
            97,
            98,
            98,
            99,
        }, // 100%
    },
};
