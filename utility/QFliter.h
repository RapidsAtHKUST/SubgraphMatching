//
// Created by ssunah on 11/1/19.
//

#ifndef SUBGRAPHMATCHING_QFLITER_H
#define SUBGRAPHMATCHING_QFLITER_H

#pragma once

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>

#include <omp.h> //OpenMP

#include <iostream>
#include <vector>
#include <fstream>
#include <set>
#include <atomic> //CAS
#include <chrono>

#include "han/intersection_algos.hpp"


using namespace std;
using namespace std::chrono;

struct BSRSet {
    int *base_ = nullptr;
    int *states_ = nullptr;
    int size_ = 0;

    BSRSet() = default;
};

struct BSRGraph {
    vector<BSRSet> bsrs;
    int max_d_ = 0;

    BSRGraph() = default;

    template<typename OFF, typename T>
    void load(size_t num_vertices, OFF &off_beg, OFF &node_off_end, T &adj) {
        bsrs.resize(num_vertices);

        int max_d = 0;
        {
            if (node_off_end != off_beg) {
                cout << "err" << endl;
                for (auto u = 0u; u < num_vertices; u++) {
                    node_off_end[u + 1] = static_cast<uint32_t>(
                            lower_bound(adj + off_beg[u], adj + off_beg[u + 1], u) - adj);
                }
            }
            for (auto u = 0u; u < num_vertices; u++) {
                max_d = max<int>(max_d, node_off_end[u + 1] - off_beg[u]);
            }
            int *tmp_base = new int[max_d];
            int *tmp_state = new int[max_d];

            for (int i = 0; i < num_vertices; i++) {

                auto degree = node_off_end[i + 1] - off_beg[i];
                auto tmp_size = offline_uint_trans_bsr(reinterpret_cast<int *>(adj) + off_beg[i], degree, tmp_base,
                                                       tmp_state);
                assert(tmp_size<= degree);
                bsrs[i].base_ = new int[degree + 16];
                bsrs[i].states_ = new int[degree + 16];
                bsrs[i].size_ = tmp_size;
                memcpy(bsrs[i].base_, tmp_base, static_cast<size_t>(tmp_size) * sizeof(int));
                memcpy(bsrs[i].states_, tmp_state, static_cast<size_t>(tmp_size) * sizeof(int));
            }

        }
        max_d_ = max_d;
    }
};

#endif //SUBGRAPHMATCHING_QFLITER_H
