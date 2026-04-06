#include "cdhit-common.h"

#include <algorithm>
#include <chrono>
#include <iostream>

#ifndef NO_OPENMP
#include <omp.h>
#else
#define omp_set_num_threads(T) (T = T)
#define omp_get_thread_num() 0
#endif

void SequenceDB::DoClustering_MasterKernel_Test(const Options &options, const char *output) {
    const int T = options.threads;
    const int NAA = options.NAA;
    double aa1_cutoff = options.cluster_thd;
    double aas_cutoff = 1 - (1 - options.cluster_thd) * 4;
    double aan_cutoff = 1 - (1 - options.cluster_thd) * options.NAA;
    const int frag_no = total_num;
    string preprocess_output_dir = options.preprocess_dir;
    if (!preprocess_output_dir.empty() && preprocess_output_dir.back() != '/' && preprocess_output_dir.back() != '\\') {
        preprocess_output_dir += '/';
    }

    if (!options.isEST) {
        cal_aax_cutoff(aa1_cutoff, aas_cutoff, aan_cutoff, options.cluster_thd, options.tolerance,
                       naa_stat_start_percent, naa_stat, NAA);
    }

    Vector<WorkingParam> params(T);
    Vector<WorkingBuffer> buffers(T);
    for (int i = 0; i < T; ++i) {
        params[i].Set(aa1_cutoff, aas_cutoff, aan_cutoff);
        buffers[i].Set(frag_no, max_len, NAAN, options);
    }
    WordTable word_table(options.NAA, NAAN);

    const int NUM_LOCKS = 262144;
    const int LOCK_MASK = 0x3FFFF;
    PaddedLock *locks = new PaddedLock[NUM_LOCKS];
    omp_set_num_threads(T);
#pragma omp parallel for schedule(static)
    for (int i = 0; i < NUM_LOCKS; ++i) {
        omp_init_lock(&locks[i].lock);
    }

    string rep_output = string(output) + ".txt";
    ofstream fout(rep_output);
    if (!fout) {
        cerr << "Cannot open output file: " << rep_output << endl;
        for (int i = 0; i < NUM_LOCKS; ++i) omp_destroy_lock(&locks[i].lock);
        delete[] locks;
        return;
    }

    const int file_slots = total_mpi_num - 1;
    if (file_slots <= 0) {
        cerr << "no workers found" << endl;
        exit(0);
    }
    vector<gzFile> chunk_fp(file_slots, nullptr);
    vector<kseq_t *> chunk_kseq(file_slots, nullptr);
    int file_index = 0;
    long long start_global_id = sequences.size();
    Clear();
    int chunk_id = 0;
    long long now_bytes = 0;
    int round_id = 0;
    kseq_t *seq = nullptr;
    Sequence one;
    std::string file = preprocess_output_dir + "_proc" + std::to_string(file_index) + ".fa";
    if (chunk_kseq[file_index] == nullptr) {
        gzFile fp = gzopen(file.c_str(), "r");
        if (!fp) {
            fprintf(stderr, "Cannot open file: %s\n", file.c_str());
            exit(1);
        }
        chunk_fp[file_index] = fp;
        seq = kseq_init(fp);
    } else {
        seq = chunk_kseq[file_index];
    }
    int len = -1;
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
        now_bytes += len;
        sequences[sequences.size() - 1]->ConvertBases();

        if (sequences.size() >= first_chunk_size) {
            chunks_id.push_back(chunk_id);
            ++chunk_id;
            chunk_kseq[file_index] = seq;
            now_bytes = 0;
            break;
        }
    }
    one.identifier = nullptr;
    one.data = nullptr;
    delete[] id_ptr;
    delete[] data_ptr;
    file_index = (file_index + 1) % file_slots;
    for (int i = 0; i < chunks_num; ++i) {

        auto start = std::chrono::high_resolution_clock::now();
        int rep_before = rep_seqs.size();
        const int N = sequences.size();
        std::vector<std::vector<int>> neigh;
        neigh.assign(N, {});
        int centers = 0;

#pragma omp parallel for schedule(dynamic, 1)
        for (int j = 0; j < N; ++j) {
            Sequence *seq = sequences[j];
            if (seq->state & IS_REDUNDANT) continue;
            int tid = omp_get_thread_num();
            ClusterOne_worker(seq, j, word_table, params[tid], buffers[tid], options, locks, NUM_LOCKS, LOCK_MASK);
        }

#pragma omp parallel for schedule(dynamic)
        for (long long j = 0; j < (long long) word_table.NAAN; ++j) {
            NVector<IndexCount> &row = word_table.indexCounts[j];
            if (row.size < 2 || row.items == NULL) continue;
            std::sort(row.items, row.items + row.size,
                      [](const IndexCount &a, const IndexCount &b) { return a.index < b.index; });
        }

        std::vector<size_t> cnt(N, 0);
        int num_batches = 20;
        int batch_size = (N + num_batches - 1) / num_batches;
        for (int batch_idx = 0; batch_idx < num_batches; ++batch_idx) {
            int start_j = batch_idx * batch_size;
            int end_j = (start_j + batch_size > N) ? N : (start_j + batch_size);
            if (start_j >= end_j) break;

#pragma omp parallel for schedule(dynamic, 1)
            for (int j = start_j; j < end_j; ++j) {
                int tid = omp_get_thread_num();
                Sequence *seq = sequences[j];
                if (seq->state & IS_REDUNDANT) continue;
                int flag = CheckOne_master(seq, j, word_table, params[tid], buffers[tid], options);
                if (flag == 0) seq->state |= IS_REP;
            }

            std::fill(cnt.begin() + start_j, cnt.begin() + end_j, 0);
            for (int t = 0; t < T; ++t) {
                auto &vec = buffers[t].thread_edges;
                for (auto &e : vec) {
                    if (e.first >= start_j && e.first < end_j) ++cnt[e.first];
                }
            }

#pragma omp parallel for schedule(static)
            for (int u = start_j; u < end_j; ++u) {
                if (cnt[u]) neigh[u].reserve(neigh[u].size() + cnt[u]);
            }

            for (int t = 0; t < T; ++t) {
                auto &vec = buffers[t].thread_edges;
                for (auto &e : vec) {
                    if (e.first >= start_j && e.first < end_j) neigh[e.first].push_back(e.second);
                }
                vec.clear();
            }

            for (int j = start_j; j < end_j; ++j) {
                Sequence *seq = sequences[j];
                if (seq->state & IS_REDUNDANT) continue;

                if (!neigh[j].empty()) {
                    for (int v : neigh[j]) {
                        if (sequences[v]->state & IS_REP) {
                            seq->state |= IS_REDUNDANT;
                            break;
                        }
                    }
                }
                if (seq->state & IS_REDUNDANT) continue;

                int size = rep_seqs.size();
                rep_seqs.Append(seq->index);
                seq->cluster_id = size;
                seq->identity = 0;
                seq->state |= IS_REP;
                seq->table_idx = centers;
                ++centers;
            }
        }

        for (int j = 0; j < N; ++j) {
            Sequence *seq = sequences[j];
            if (seq->state & IS_REDUNDANT) continue;
            fout << ">" << seq->identifier << "\n";
            fout << seq->true_data << "\n";
        }

        word_table.Clear();
        start_global_id += N;
        int rep_after = rep_seqs.size();
        int rep_added = rep_after - rep_before;

        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        cerr << "[MasterKernelTest] chunk " << round_id << " seqs_read=" << N << " clusters_added=" << rep_added
             << " total_clusters=" << rep_after << " done in " << elapsed.count() << " s" << endl;

        ++round_id;
        if (i == chunks_num - 1) break;
        for (int i = 0; i < sequences.size(); ++i) delete sequences[i];
        sequences.clear();
        Sequence one;
        std::string file = preprocess_output_dir + "_proc" + std::to_string(file_index) + ".fa";
        if (chunk_kseq[file_index] == nullptr) {
            gzFile fp = gzopen(file.c_str(), "r");
            if (!fp) {
                fprintf(stderr, "Cannot open file: %s\n", file.c_str());
                exit(1);
            }
            chunk_fp[file_index] = fp;
            seq = kseq_init(fp);
        } else {
            seq = chunk_kseq[file_index];
        }
        int len = -1;
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
            now_bytes += len;
            sequences[sequences.size() - 1]->ConvertBases();
            if (now_bytes > chunk_bytes || sequences.size() >= chunk_size) {
                chunks_id.push_back(chunk_id);
                ++chunk_id;
                chunk_kseq[file_index] = seq;
                now_bytes = 0;
                break;
            }
        }
        if (now_bytes > 0) {
            chunks_id.push_back(chunk_id);
        }
        one.identifier = nullptr;
        one.data = nullptr;
        delete[] id_ptr;
        delete[] data_ptr;
        file_index = (file_index + 1) % file_slots;
    }

    for (int i = 0; i < file_slots; ++i) {
        if (chunk_kseq[i]) kseq_destroy(chunk_kseq[i]);
        if (chunk_fp[i]) gzclose(chunk_fp[i]);
    }

    fout.close();

    for (int i = 0; i < NUM_LOCKS; ++i) omp_destroy_lock(&locks[i].lock);
    delete[] locks;

    cout << "[DONE] Master kernel test output: " << rep_output << endl;
}
