//
// Created by ssunah on 11/20/18.
//

#include "FilterVertices.h"
#include "GenerateFilteringPlan.h"
#include <memory.h>
#include <utility/graphoperations.h>
#include <vector>
#include <algorithm>
#define INVALID_VERTEX_ID 100000000

bool
FilterVertices::LDFFilter(const Graph *data_graph, const Graph *query_graph, ui **&candidates, ui *&candidates_count) {
    allocateBuffer(data_graph, query_graph, candidates, candidates_count);

    for (ui i = 0; i < query_graph->getVerticesCount(); ++i) {
        LabelID label = query_graph->getVertexLabel(i);
        ui degree = query_graph->getVertexDegree(i);

        ui data_vertex_num;
        const ui* data_vertices = data_graph->getVerticesByLabel(label, data_vertex_num);

        for (ui j = 0; j < data_vertex_num; ++j) {
            ui data_vertex = data_vertices[j];
            if (data_graph->getVertexDegree(data_vertex) >= degree) {
                candidates[i][candidates_count[i]++] = data_vertex;
            }
        }

        if (candidates_count[i] == 0) {
            return false;
        }
    }

    return true;
}

bool
FilterVertices::NLFFilter(const Graph *data_graph, const Graph *query_graph, ui **&candidates, ui *&candidates_count) {
    allocateBuffer(data_graph, query_graph, candidates, candidates_count);

    for (ui i = 0; i < query_graph->getVerticesCount(); ++i) {
        VertexID query_vertex = i;
        computeCandidateWithNLF(data_graph, query_graph, query_vertex, candidates_count[query_vertex], candidates[query_vertex]);

        if (candidates_count[query_vertex] == 0) {
            return false;
        }
    }

    return true;
}

bool
FilterVertices::GQLFilter(const Graph *data_graph, const Graph *query_graph, ui **&candidates, ui *&candidates_count) {
    // Local refinement.
    if (!NLFFilter(data_graph, query_graph, candidates, candidates_count))
        return false;

    // Allocate buffer.
    ui query_vertex_num = query_graph->getVerticesCount();
    ui data_vertex_num = data_graph->getVerticesCount();

    bool** valid_candidates = new bool*[query_vertex_num];
    for (ui i = 0; i < query_vertex_num; ++i) {
        valid_candidates[i] = new bool[data_vertex_num];
        memset(valid_candidates[i], 0, sizeof(bool) * data_vertex_num);
    }

    ui query_graph_max_degree = query_graph->getGraphMaxDegree();
    ui data_graph_max_degree = data_graph->getGraphMaxDegree();

    int* left_to_right_offset = new int[query_graph_max_degree + 1];
    int* left_to_right_edges = new int[query_graph_max_degree * data_graph_max_degree];
    int* left_to_right_match = new int[query_graph_max_degree];
    int* right_to_left_match = new int[data_graph_max_degree];
    int* match_visited = new int[data_graph_max_degree + 1];
    int* match_queue = new int[query_vertex_num];
    int* match_previous = new int[data_graph_max_degree + 1];

    // Record valid candidate vertices for each query vertex.
    for (ui i = 0; i < query_vertex_num; ++i) {
        VertexID query_vertex = i;
        for (ui j = 0; j < candidates_count[query_vertex]; ++j) {
            VertexID data_vertex = candidates[query_vertex][j];
            valid_candidates[query_vertex][data_vertex] = true;
        }
    }

    // Global refinement.
    for (ui l = 0; l < 2; ++l) {
        for (ui i = 0; i < query_vertex_num; ++i) {
            VertexID query_vertex = i;
            for (ui j = 0; j < candidates_count[query_vertex]; ++j) {
                VertexID data_vertex = candidates[query_vertex][j];

                if (data_vertex == INVALID_VERTEX_ID)
                    continue;

                if (!verifyExactTwigIso(data_graph, query_graph, data_vertex, query_vertex, valid_candidates,
                                        left_to_right_offset, left_to_right_edges, left_to_right_match,
                                        right_to_left_match, match_visited, match_queue, match_previous)) {
                    candidates[query_vertex][j] = INVALID_VERTEX_ID;
                    valid_candidates[query_vertex][data_vertex] = false;
                }
            }
        }
    }

    // Compact candidates.
    compactCandidates(candidates, candidates_count, query_vertex_num);

    // Release memory.
    for (ui i = 0; i < query_vertex_num; ++i) {
        delete[] valid_candidates[i];
    }
    delete[] valid_candidates;
    delete[] left_to_right_offset;
    delete[] left_to_right_edges;
    delete[] left_to_right_match;
    delete[] right_to_left_match;
    delete[] match_visited;
    delete[] match_queue;
    delete[] match_previous;

    return isCandidateSetValid(candidates, candidates_count, query_vertex_num);
}

bool
FilterVertices::TSOFilter(const Graph *data_graph, const Graph *query_graph, ui **&candidates, ui *&candidates_count,
                          ui *&order, TreeNode *&tree) {
    allocateBuffer(data_graph, query_graph, candidates, candidates_count);
    GenerateFilteringPlan::generateTSOFilterPlan(data_graph, query_graph, tree, order);

    ui query_vertex_num = query_graph->getVerticesCount();

    // Get the candidates of the start vertex.
    VertexID start_vertex = order[0];
    computeCandidateWithNLF(data_graph, query_graph, start_vertex, candidates_count[start_vertex], candidates[start_vertex]);

    ui* updated_flag = new ui[data_graph->getVerticesCount()];
    ui* flag = new ui[data_graph->getVerticesCount()];
    std::fill(flag, flag + data_graph->getVerticesCount(), 0);

    for (ui i = 1; i < query_vertex_num; ++i) {
        VertexID query_vertex = order[i];
        TreeNode& node = tree[query_vertex];
        generateCandidates(data_graph, query_graph, query_vertex, &node.parent_, 1, candidates, candidates_count, flag, updated_flag);
    }

    for (int i = query_vertex_num - 1; i >= 0; --i) {
        VertexID query_vertex = order[i];
        TreeNode& node = tree[query_vertex];
        if (node.children_count_ > 0) {
            pruneCandidates(data_graph, query_graph, query_vertex, node.children_, node.children_count_, candidates, candidates_count, flag, updated_flag);
        }
    }

    compactCandidates(candidates, candidates_count, query_vertex_num);

    delete[] updated_flag;
    delete[] flag;
    return isCandidateSetValid(candidates, candidates_count, query_vertex_num);
}

bool
FilterVertices::CFLFilter(const Graph *data_graph, const Graph *query_graph, ui **&candidates, ui *&candidates_count,
                          ui *&order, TreeNode *&tree) {
    allocateBuffer(data_graph, query_graph, candidates, candidates_count);
    int level_count;
    ui* level_offset;
    GenerateFilteringPlan::generateCFLFilterPlan(data_graph, query_graph, tree, order, level_count, level_offset);

    VertexID start_vertex = order[0];
    computeCandidateWithNLF(data_graph, query_graph, start_vertex, candidates_count[start_vertex], candidates[start_vertex]);

    ui* updated_flag = new ui[data_graph->getVerticesCount()];
    ui* flag = new ui[data_graph->getVerticesCount()];
    std::fill(flag, flag + data_graph->getVerticesCount(), 0);

    // Top-down generation.
    for (int i = 1; i < level_count; ++i) {
        // Forward generation.
        for (int j = level_offset[i]; j < level_offset[i + 1]; ++j) {
            VertexID query_vertex = order[j];
            TreeNode& node = tree[query_vertex];
            generateCandidates(data_graph, query_graph, query_vertex, node.bn_, node.bn_count_, candidates, candidates_count, flag, updated_flag);
        }

        // Backward prune.
        for (int j = level_offset[i + 1] - 1; j >= level_offset[i]; --j) {
            VertexID query_vertex = order[j];
            TreeNode& node = tree[query_vertex];

            if (node.fn_count_ > 0) {
                pruneCandidates(data_graph, query_graph, query_vertex, node.fn_, node.fn_count_, candidates, candidates_count, flag, updated_flag);
            }
        }
    }

    // Bottom-up refinement.
    for (int i = level_count - 2; i >= 0; --i) {
        for (int j = level_offset[i]; j < level_offset[i + 1]; ++j) {
            VertexID query_vertex = order[j];
            TreeNode& node = tree[query_vertex];

            if (node.under_level_count_ > 0) {
                pruneCandidates(data_graph, query_graph, query_vertex, node.under_level_, node.under_level_count_, candidates, candidates_count, flag, updated_flag);
            }
        }
    }


    compactCandidates(candidates, candidates_count, query_graph->getVerticesCount());

    delete[] updated_flag;
    delete[] flag;
    return isCandidateSetValid(candidates, candidates_count, query_graph->getVerticesCount());
}

bool
FilterVertices::DPisoFilter(const Graph *data_graph, const Graph *query_graph, ui **&candidates, ui *&candidates_count,
                            ui *&order, TreeNode *&tree) {
    if (!LDFFilter(data_graph, query_graph, candidates, candidates_count))
        return false;

    GenerateFilteringPlan::generateDPisoFilterPlan(data_graph, query_graph, tree, order);

    ui query_vertices_num = query_graph->getVerticesCount();
    ui* updated_flag = new ui[data_graph->getVerticesCount()];
    ui* flag = new ui[data_graph->getVerticesCount()];
    std::fill(flag, flag + data_graph->getVerticesCount(), 0);

    // The number of refinement is k. According to the original paper, we set k as 3.
    for (ui k = 0; k < 3; ++k) {
        if (k % 2 == 0) {
            for (int i = 1; i < query_vertices_num; ++i) {
                VertexID query_vertex = order[i];
                TreeNode& node = tree[query_vertex];
                pruneCandidates(data_graph, query_graph, query_vertex, node.bn_, node.bn_count_, candidates, candidates_count, flag, updated_flag);
            }
        }
        else {
            for (int i = query_vertices_num - 2; i >= 0; --i) {
                VertexID query_vertex = order[i];
                TreeNode& node = tree[query_vertex];
                pruneCandidates(data_graph, query_graph, query_vertex, node.fn_, node.fn_count_, candidates, candidates_count, flag, updated_flag);
            }
        }
    }

    compactCandidates(candidates, candidates_count, query_graph->getVerticesCount());

    delete[] updated_flag;
    delete[] flag;
    return isCandidateSetValid(candidates, candidates_count, query_graph->getVerticesCount());
}

bool
FilterVertices::CECIFilter(const Graph *data_graph, const Graph *query_graph, ui **&candidates, ui *&candidates_count,
                           ui *&order, TreeNode *&tree,  std::vector<std::unordered_map<VertexID, std::vector<VertexID >>> &TE_Candidates,
                           std::vector<std::vector<std::unordered_map<VertexID, std::vector<VertexID>>>> &NTE_Candidates) {
    GenerateFilteringPlan::generateCECIFilterPlan(data_graph, query_graph, tree, order);

    allocateBuffer(data_graph, query_graph, candidates, candidates_count);

    ui query_vertices_count = query_graph->getVerticesCount();
    ui data_vertices_count = data_graph->getVerticesCount();
    // Find the pivots.
    VertexID root = order[0];
    computeCandidateWithNLF(data_graph, query_graph, root, candidates_count[root], candidates[root]);

    if (candidates_count[root] == 0)
        return false;

    // TE_Candidates construction and filtering.
    std::vector<ui> updated_flag(data_vertices_count);
    std::vector<ui> flag(data_vertices_count);
    std::fill(flag.begin(), flag.end(), 0);
    std::vector<bool> visited_query_vertex(query_vertices_count);
    std::fill(visited_query_vertex.begin(), visited_query_vertex.end(), false);

    visited_query_vertex[root] = true;

    TE_Candidates.resize(query_vertices_count);

    for (ui i = 1; i < query_vertices_count; ++i) {
        VertexID u = order[i];
        TreeNode& u_node = tree[u];
        VertexID u_p = tree[u].parent_;

        ui u_label = query_graph->getVertexLabel(u);
        ui u_degree = query_graph->getVertexDegree(u);
#if OPTIMIZED_LABELED_GRAPH == 1
        const std::unordered_map<LabelID, ui>* u_nlf = query_graph->getVertexNLF(u);
#endif
        candidates_count[u] = 0;

        visited_query_vertex[u] = true;
        VertexID* frontiers = candidates[u_p];
        ui frontiers_count = candidates_count[u_p];

        for (ui j = 0; j < frontiers_count; ++j) {
            VertexID v_f = frontiers[j];

            if (v_f == INVALID_VERTEX_ID)
                continue;

            ui nbrs_cnt;
            const VertexID* nbrs = data_graph->getVertexNeighbors(v_f, nbrs_cnt);

            auto iter_pair = TE_Candidates[u].emplace(v_f, std::vector<VertexID>());
            for (ui k = 0; k < nbrs_cnt; ++k) {
                VertexID v = nbrs[k];

                if (data_graph->getVertexLabel(v) == u_label && data_graph->getVertexDegree(v) >= u_degree) {

                    // NLF check

#if OPTIMIZED_LABELED_GRAPH == 1
                    const std::unordered_map<LabelID, ui>* v_nlf = data_graph->getVertexNLF(v);

                    if (v_nlf->size() >= u_nlf->size()) {
                        bool is_valid = true;

                        for (auto element : *u_nlf) {
                            auto iter = v_nlf->find(element.first);
                            if (iter == v_nlf->end() || iter->second < element.second) {
                                is_valid = false;
                                break;
                            }
                        }

                        if (is_valid) {
                            iter_pair.first->second.push_back(v);
                            if (flag[v] == 0) {
                                flag[v] = 1;
                                candidates[u][candidates_count[u]++] = v;
                            }
                        }
                    }
#else
                    iter_pair.first->second.push_back(v);
                    if (flag[v] == 0) {
                        flag[v] = 1;
                        candidates[u][candidates_count[u]++] = v;
                    }
#endif
                }
            }

            if (iter_pair.first->second.empty()) {
                frontiers[j] = INVALID_VERTEX_ID;
                for (ui k = 0; k < tree[u_p].children_count_; ++k) {
                    VertexID u_c = tree[u_p].children_[k];
                    if (visited_query_vertex[u_c]) {
                        TE_Candidates[u_c].erase(v_f);
                    }
                }
            }
        }

        if (candidates_count[u] == 0)
            return false;

        for (ui j = 0; j < candidates_count[u]; ++j) {
            VertexID v = candidates[u][j];
            flag[v] = 0;
        }
    }

    // NTE_Candidates construction and filtering.
    NTE_Candidates.resize(query_vertices_count);
    for (auto& value : NTE_Candidates) {
        value.resize(query_vertices_count);
    }

    for (ui i = 1; i < query_vertices_count; ++i) {
        VertexID u = order[i];
        TreeNode &u_node = tree[u];

        ui u_label = query_graph->getVertexLabel(u);
        ui u_degree = query_graph->getVertexDegree(u);
#if OPTIMIZED_LABELED_GRAPH == 1
        const std::unordered_map<LabelID, ui> *u_nlf = query_graph->getVertexNLF(u);
#endif
        for (ui l = 0; l < u_node.bn_count_; ++l) {
            VertexID u_p = u_node.bn_[l];
            VertexID *frontiers = candidates[u_p];
            ui frontiers_count = candidates_count[u_p];

            for (ui j = 0; j < frontiers_count; ++j) {
                VertexID v_f = frontiers[j];

                if (v_f == INVALID_VERTEX_ID)
                    continue;

                ui nbrs_cnt;
                const VertexID *nbrs = data_graph->getVertexNeighbors(v_f, nbrs_cnt);

                auto iter_pair = NTE_Candidates[u][u_p].emplace(v_f, std::vector<VertexID>());
                for (ui k = 0; k < nbrs_cnt; ++k) {
                    VertexID v = nbrs[k];

                    if (data_graph->getVertexLabel(v) == u_label && data_graph->getVertexDegree(v) >= u_degree) {

                        // NLF check
#if OPTIMIZED_LABELED_GRAPH == 1
                        const std::unordered_map<LabelID, ui> *v_nlf = data_graph->getVertexNLF(v);

                        if (v_nlf->size() >= u_nlf->size()) {
                            bool is_valid = true;

                            for (auto element : *u_nlf) {
                                auto iter = v_nlf->find(element.first);
                                if (iter == v_nlf->end() || iter->second < element.second) {
                                    is_valid = false;
                                    break;
                                }
                            }

                            if (is_valid) {
                                iter_pair.first->second.push_back(v);
                            }
                        }
#else
                        iter_pair.first->second.push_back(v);
#endif
                    }
                }

                if (iter_pair.first->second.empty()) {
                    frontiers[j] = INVALID_VERTEX_ID;
                    for (ui k = 0; k < tree[u_p].children_count_; ++k) {
                        VertexID u_c = tree[u_p].children_[k];
                        TE_Candidates[u_c].erase(v_f);
                    }
                }
            }
        }
    }

    // Reverse BFS refine.
    std::vector<std::vector<ui>> cardinality(query_vertices_count);
    for (ui i = 0; i < query_vertices_count; ++i) {
        cardinality[i].resize(candidates_count[i], 1);
    }

    std::vector<ui> local_cardinality(data_vertices_count);
    std::fill(local_cardinality.begin(), local_cardinality.end(), 0);

    for (int i = query_vertices_count - 1; i >= 0; --i) {
        VertexID u = order[i];
        TreeNode& u_node = tree[u];

        ui flag_num = 0;
        ui updated_flag_count = 0;

        // Compute the intersection of TE_Candidates and NTE_Candidates.
        for (ui j = 0; j < candidates_count[u]; ++j) {
            VertexID v = candidates[u][j];

            if (v == INVALID_VERTEX_ID)
                continue;

            if (flag[v] == flag_num) {
                flag[v] += 1;
                updated_flag[updated_flag_count++] = v;
            }
        }

        for (ui j = 0; j < u_node.bn_count_; ++j) {
            VertexID u_bn = u_node.bn_[j];
            flag_num += 1;
            for (auto iter = NTE_Candidates[u][u_bn].begin(); iter != NTE_Candidates[u][u_bn].end(); ++iter) {
                for (auto v : iter->second) {
                    if (flag[v] == flag_num) {
                        flag[v] += 1;
                    }
                }
            }
        }

        flag_num += 1;

        // Get the cardinality of the candidates of u.
        for (ui j = 0; j < candidates_count[u]; ++j) {
            VertexID v = candidates[u][j];
            if (v != INVALID_VERTEX_ID && flag[v] == flag_num) {
                local_cardinality[v] = cardinality[u][j];
            }
            else {
                cardinality[u][j] = 0;
            }
        }

        VertexID u_p = u_node.parent_;
        VertexID* frontiers = candidates[u_p];
        ui frontiers_count = candidates_count[u_p];

        // Loop over TE_Candidates.
        for (ui j = 0; j < frontiers_count; ++j) {
            VertexID v_f = frontiers[j];

            if (v_f == INVALID_VERTEX_ID) {
                cardinality[u_p][j] = 0;
                continue;
            }

            ui temp_score = 0;
            for (auto iter = TE_Candidates[u][v_f].begin(); iter != TE_Candidates[u][v_f].end();) {
                VertexID v = *iter;
                temp_score += local_cardinality[v];
                if (local_cardinality[v] == 0) {
                    iter = TE_Candidates[u][v_f].erase(iter);
                    for (ui k = 0; k < u_node.children_count_; ++k) {
                        VertexID u_c = u_node.children_[k];
                        TE_Candidates[u_c].erase(v);
                    }

                    for (ui k = 0; k < u_node.fn_count_; ++k) {
                        VertexID u_c = u_node.fn_[k];
                        NTE_Candidates[u_c][u].erase(v);
                    }
                }
                else {
                    ++iter;
                }
            }

            cardinality[u_p][j] *= temp_score;
        }

        // Clear updated flag.
        for (ui j = 0; j < updated_flag_count; ++j) {
            flag[updated_flag[j]] = 0;
            local_cardinality[updated_flag[j]] = 0;
        }
    }

    compactCandidates(candidates, candidates_count, query_vertices_count);
    sortCandidates(candidates, candidates_count, query_vertices_count);


    for (ui i = 0; i < query_vertices_count; ++i) {
        if (candidates_count[i] == 0) {
            return false;
        }
    }

    for (ui i = 1; i < query_vertices_count; ++i) {
        VertexID u = order[i];
        TreeNode& u_node = tree[u];

        // Clear TE_Candidates.
        {
            VertexID u_p = u_node.parent_;
            auto iter = TE_Candidates[u].begin();
            while (iter != TE_Candidates[u].end()) {
                VertexID v_f = iter->first;
                if (!std::binary_search(candidates[u_p], candidates[u_p] + candidates_count[u_p], v_f)) {
                    iter = TE_Candidates[u].erase(iter);
                }
                else {
                    std::sort(iter->second.begin(), iter->second.end());
                    iter++;
                }
            }
        }

        // Clear NTE_Candidates.
        {
            for (ui j = 0; j < u_node.bn_count_; ++j) {
                VertexID u_p = u_node.bn_[j];
                auto iter = NTE_Candidates[u][u_p].end();
                while (iter != NTE_Candidates[u][u_p].end()) {
                    VertexID v_f = iter->first;
                    if (!std::binary_search(candidates[u_p], candidates[u_p] + candidates_count[u_p], v_f)) {
                        iter = NTE_Candidates[u][u_p].erase(iter);
                    }
                    else {
                        std::sort(iter->second.begin(), iter->second.end());
                        iter++;
                    }
                }
            }
        }
    }

//    for (ui i = 0; i < query_vertices_count; ++i) {
//        VertexID u = i;
//        std::cout << u << ':';
//        for (ui j = 0; j < candidates_count[u]; ++j) {
//            std::cout << candidates[u][j] << ' ';
//        }
//        std::cout << std::endl;
//    }
//
//    for (ui i = 1; i < query_vertices_count; ++i) {
//        VertexID u = order[i];
//        // TE_Candidates
//        std::cout << "TE_Candidates: " << u << ',' << tree[u].parent_ << std::endl;
//        for (auto iter = TE_Candidates[u].begin(); iter != TE_Candidates[u].end(); ++iter) {
//            std::cout << iter->first << ": ";
//            for (auto v : iter->second) {
//                std::cout << v << ' ';
//                if (!data_graph->checkEdgeExistence(iter->first, v)) {
//                    std::cout << "Edge does not exist" << std::endl;
//                }
//            }
//            std::cout << std::endl;
//        }
//        std::cout << "-----" << std::endl;
//        for (ui j = 0; j < tree[u].bn_count_; ++j) {
//            VertexID u_bn = tree[u].bn_[j];
//            std::cout << "NTE_Candidates: " << u << ',' << u_bn << std::endl;
//            for (auto iter = NTE_Candidates[u][u_bn].begin(); iter != NTE_Candidates[u][u_bn].end(); ++iter) {
//                std::cout << iter->first << ": ";
//                for (auto v : iter->second) {
//                    std::cout << v << ' ';
//                    if (!data_graph->checkEdgeExistence(iter->first, v)) {
//                        std::cout << "Edge does not exist" << std::endl;
//                    }
//                }
//                std::cout << std::endl;
//            }
//            std::cout << "-----" << std::endl;
//        }
//    }

    return true;
}

void FilterVertices::allocateBuffer(const Graph *data_graph, const Graph *query_graph, ui **&candidates,
                                    ui *&candidates_count) {
    ui query_vertex_num = query_graph->getVerticesCount();
    ui candidates_max_num = data_graph->getGraphMaxLabelFrequency();

    candidates_count = new ui[query_vertex_num];
    memset(candidates_count, 0, sizeof(ui) * query_vertex_num);

    candidates = new ui*[query_vertex_num];

    for (ui i = 0; i < query_vertex_num; ++i) {
        candidates[i] = new ui[candidates_max_num];
    }
}

bool
FilterVertices::verifyExactTwigIso(const Graph *data_graph, const Graph *query_graph, ui data_vertex, ui query_vertex,
                                   bool **valid_candidates, int *left_to_right_offset, int *left_to_right_edges,
                                   int *left_to_right_match, int *right_to_left_match, int* match_visited,
                                   int* match_queue, int* match_previous) {
    // Construct the bipartite graph between N(query_vertex) and N(data_vertex)
    ui left_partition_size;
    ui right_partition_size;
    const VertexID* query_vertex_neighbors = query_graph->getVertexNeighbors(query_vertex, left_partition_size);
    const VertexID* data_vertex_neighbors = data_graph->getVertexNeighbors(data_vertex, right_partition_size);

    ui edge_count = 0;
    for (int i = 0; i < left_partition_size; ++i) {
        VertexID query_vertex_neighbor = query_vertex_neighbors[i];
        left_to_right_offset[i] = edge_count;

        for (int j = 0; j < right_partition_size; ++j) {
            VertexID data_vertex_neighbor = data_vertex_neighbors[j];

            if (valid_candidates[query_vertex_neighbor][data_vertex_neighbor]) {
                left_to_right_edges[edge_count++] = j;
            }
        }
    }
    left_to_right_offset[left_partition_size] = edge_count;

    memset(left_to_right_match, -1, left_partition_size * sizeof(int));
    memset(right_to_left_match, -1, right_partition_size * sizeof(int));

    GraphOperations::match_bfs(left_to_right_offset, left_to_right_edges, left_to_right_match, right_to_left_match,
                               match_visited, match_queue, match_previous, left_partition_size, right_partition_size);
    for (int i = 0; i < left_partition_size; ++i) {
        if (left_to_right_match[i] == -1)
            return false;
    }

    return true;
}

void FilterVertices::compactCandidates(ui **&candidates, ui *&candidates_count, ui query_vertex_num) {
    for (ui i = 0; i < query_vertex_num; ++i) {
        VertexID query_vertex = i;
        ui next_position = 0;
        for (ui j = 0; j < candidates_count[query_vertex]; ++j) {
            VertexID data_vertex = candidates[query_vertex][j];

            if (data_vertex != INVALID_VERTEX_ID) {
                candidates[query_vertex][next_position++] = data_vertex;
            }
        }

        candidates_count[query_vertex] = next_position;
    }
}

bool FilterVertices::isCandidateSetValid(ui **&candidates, ui *&candidates_count, ui query_vertex_num) {
    for (ui i = 0; i < query_vertex_num; ++i) {
        if (candidates_count[i] == 0)
            return false;
    }
    return true;
}

void
FilterVertices::computeCandidateWithNLF(const Graph *data_graph, const Graph *query_graph, VertexID query_vertex,
                                               ui &count, ui *buffer) {
    LabelID label = query_graph->getVertexLabel(query_vertex);
    ui degree = query_graph->getVertexDegree(query_vertex);
#if OPTIMIZED_LABELED_GRAPH == 1
    const std::unordered_map<LabelID, ui>* query_vertex_nlf = query_graph->getVertexNLF(query_vertex);
#endif
    ui data_vertex_num;
    const ui* data_vertices = data_graph->getVerticesByLabel(label, data_vertex_num);
    count = 0;
    for (ui j = 0; j < data_vertex_num; ++j) {
        ui data_vertex = data_vertices[j];
        if (data_graph->getVertexDegree(data_vertex) >= degree) {

            // NFL check
#if OPTIMIZED_LABELED_GRAPH == 1
            const std::unordered_map<LabelID, ui>* data_vertex_nlf = data_graph->getVertexNLF(data_vertex);

            if (data_vertex_nlf->size() >= query_vertex_nlf->size()) {
                bool is_valid = true;

                for (auto element : *query_vertex_nlf) {
                    auto iter = data_vertex_nlf->find(element.first);
                    if (iter == data_vertex_nlf->end() || iter->second < element.second) {
                        is_valid = false;
                        break;
                    }
                }

                if (is_valid) {
                    if (buffer != NULL) {
                        buffer[count] = data_vertex;
                    }
                    count += 1;
                }
            }
#else
            if (buffer != NULL) {
                buffer[count] = data_vertex;
            }
            count += 1;
#endif
        }
    }

}

void FilterVertices::computeCandidateWithLDF(const Graph *data_graph, const Graph *query_graph, VertexID query_vertex,
                                             ui &count, ui *buffer) {
    LabelID label = query_graph->getVertexLabel(query_vertex);
    ui degree = query_graph->getVertexDegree(query_vertex);
    count = 0;
    ui data_vertex_num;
    const ui* data_vertices = data_graph->getVerticesByLabel(label, data_vertex_num);

    if (buffer == NULL) {
        for (ui i = 0; i < data_vertex_num; ++i) {
            VertexID v = data_vertices[i];
            if (data_graph->getVertexDegree(v) >= degree) {
                count += 1;
            }
        }
    }
    else {
        for (ui i = 0; i < data_vertex_num; ++i) {
            VertexID v = data_vertices[i];
            if (data_graph->getVertexDegree(v) >= degree) {
                buffer[count++] = v;
            }
        }
    }
}

void FilterVertices::generateCandidates(const Graph *data_graph, const Graph *query_graph, VertexID query_vertex,
                                       VertexID *pivot_vertices, ui pivot_vertices_count, VertexID **candidates,
                                       ui *candidates_count, ui *flag, ui *updated_flag) {
    LabelID query_vertex_label = query_graph->getVertexLabel(query_vertex);
    ui query_vertex_degree = query_graph->getVertexDegree(query_vertex);
#if OPTIMIZED_LABELED_GRAPH == 1
    const std::unordered_map<LabelID , ui>* query_vertex_nlf = query_graph->getVertexNLF(query_vertex);
#endif
    ui count = 0;
    ui updated_flag_count = 0;
    for (ui i = 0; i < pivot_vertices_count; ++i) {
        VertexID pivot_vertex = pivot_vertices[i];

        for (ui j = 0; j < candidates_count[pivot_vertex]; ++j) {
            VertexID v = candidates[pivot_vertex][j];

            if (v == INVALID_VERTEX_ID)
                continue;
            ui v_nbrs_count;
            const VertexID* v_nbrs = data_graph->getVertexNeighbors(v, v_nbrs_count);

            for (ui k = 0; k < v_nbrs_count; ++k) {
                VertexID v_nbr = v_nbrs[k];
                LabelID v_nbr_label = data_graph->getVertexLabel(v_nbr);
                ui v_nbr_degree = data_graph->getVertexDegree(v_nbr);

                if (flag[v_nbr] == count && v_nbr_label == query_vertex_label && v_nbr_degree >= query_vertex_degree) {
                    flag[v_nbr] += 1;

                    if (count == 0) {
                        updated_flag[updated_flag_count++] = v_nbr;
                    }
                }
            }
        }

        count += 1;
    }

    for (ui i = 0; i < updated_flag_count; ++i) {
        VertexID v = updated_flag[i];
        if (flag[v] == count) {
            // NLF filter.
#if OPTIMIZED_LABELED_GRAPH == 1
            const std::unordered_map<LabelID, ui>* data_vertex_nlf = data_graph->getVertexNLF(v);

            if (data_vertex_nlf->size() >= query_vertex_nlf->size()) {
                bool is_valid = true;

                for (auto element : *query_vertex_nlf) {
                    auto iter = data_vertex_nlf->find(element.first);
                    if (iter == data_vertex_nlf->end() || iter->second < element.second) {
                        is_valid = false;
                        break;
                    }
                }

                if (is_valid) {
                    candidates[query_vertex][candidates_count[query_vertex]++] = v;
                }
            }
#else
            candidates[query_vertex][candidates_count[query_vertex]++] = v;
#endif
        }
    }

    for (ui i = 0; i < updated_flag_count; ++i) {
        ui v = updated_flag[i];
        flag[v] = 0;
    }
}

void FilterVertices::pruneCandidates(const Graph *data_graph, const Graph *query_graph, VertexID query_vertex,
                                    VertexID *pivot_vertices, ui pivot_vertices_count, VertexID **candidates,
                                    ui *candidates_count, ui *flag, ui *updated_flag) {
    LabelID query_vertex_label = query_graph->getVertexLabel(query_vertex);
    ui query_vertex_degree = query_graph->getVertexDegree(query_vertex);

    ui count = 0;
    ui updated_flag_count = 0;
    for (ui i = 0; i < pivot_vertices_count; ++i) {
        VertexID pivot_vertex = pivot_vertices[i];

        for (ui j = 0; j < candidates_count[pivot_vertex]; ++j) {
            VertexID v = candidates[pivot_vertex][j];

            if (v == INVALID_VERTEX_ID)
                continue;
            ui v_nbrs_count;
            const VertexID* v_nbrs = data_graph->getVertexNeighbors(v, v_nbrs_count);

            for (ui k = 0; k < v_nbrs_count; ++k) {
                VertexID v_nbr = v_nbrs[k];
                LabelID v_nbr_label = data_graph->getVertexLabel(v_nbr);
                ui v_nbr_degree = data_graph->getVertexDegree(v_nbr);

                if (flag[v_nbr] == count && v_nbr_label == query_vertex_label && v_nbr_degree >= query_vertex_degree) {
                    flag[v_nbr] += 1;

                    if (count == 0) {
                        updated_flag[updated_flag_count++] = v_nbr;
                    }
                }
            }
        }

        count += 1;
    }

    for (ui i = 0; i < candidates_count[query_vertex]; ++i) {
        ui v = candidates[query_vertex][i];
        if (v == INVALID_VERTEX_ID)
            continue;

        if (flag[v] != count) {
            candidates[query_vertex][i] = INVALID_VERTEX_ID;
        }
    }

    for (ui i = 0; i < updated_flag_count; ++i) {
        ui v = updated_flag[i];
        flag[v] = 0;
    }
}

void FilterVertices::printCandidatesInfo(const Graph *query_graph, ui *candidates_count, std::vector<ui> &optimal_candidates_count) {
    std::vector<std::pair<VertexID, ui>> core_vertices;
    std::vector<std::pair<VertexID, ui>> tree_vertices;
    std::vector<std::pair<VertexID, ui>> leaf_vertices;

    ui query_vertices_num = query_graph->getVerticesCount();
    double sum = 0;
    double optimal_sum = 0;
    for (ui i = 0; i < query_vertices_num; ++i) {
        VertexID cur_vertex = i;
        ui count = candidates_count[cur_vertex];
        sum += count;
        optimal_sum += optimal_candidates_count[cur_vertex];

        if (query_graph->getCoreValue(cur_vertex) > 1) {
            core_vertices.emplace_back(std::make_pair(cur_vertex, count));
        }
        else {
            if (query_graph->getVertexDegree(cur_vertex) > 1) {
                tree_vertices.emplace_back(std::make_pair(cur_vertex, count));
            }
            else {
                leaf_vertices.emplace_back(std::make_pair(cur_vertex, count));
            }
        }
    }

    printf("#Candidate Information: CoreVertex(%zu), TreeVertex(%zu), LeafVertex(%zu)\n", core_vertices.size(), tree_vertices.size(), leaf_vertices.size());

    for (auto candidate_info : core_vertices) {
        printf("CoreVertex %u: %u, %u \n", candidate_info.first, candidate_info.second, optimal_candidates_count[candidate_info.first]);
    }

    for (auto candidate_info : tree_vertices) {
        printf("TreeVertex %u: %u, %u\n", candidate_info.first, candidate_info.second, optimal_candidates_count[candidate_info.first]);
    }

    for (auto candidate_info : leaf_vertices) {
        printf("LeafVertex %u: %u, %u\n", candidate_info.first, candidate_info.second, optimal_candidates_count[candidate_info.first]);
    }

    printf("Total #Candidates: %.1lf, %.1lf\n", sum, optimal_sum);
}

void FilterVertices::sortCandidates(ui **candidates, ui *candidates_count, ui num) {
    for (ui i = 0; i < num; ++i) {
        std::sort(candidates[i], candidates[i] + candidates_count[i]);
    }
}

double
FilterVertices::computeCandidatesFalsePositiveRatio(const Graph *data_graph, const Graph *query_graph, ui **candidates,
                                                    ui *candidates_count, std::vector<ui> &optimal_candidates_count) {
    ui query_vertices_count = query_graph->getVerticesCount();
    ui data_vertices_count = data_graph->getVerticesCount();

    std::vector<std::vector<ui>> candidates_copy(query_vertices_count);
    for (ui i = 0; i < query_vertices_count; ++i) {
        candidates_copy[i].resize(candidates_count[i]);
        std::copy(candidates[i], candidates[i] + candidates_count[i], candidates_copy[i].begin());
    }

    std::vector<int> flag(data_vertices_count, 0);
    std::vector<ui> updated_flag;
    std::vector<double> per_query_vertex_false_positive_ratio(query_vertices_count);
    optimal_candidates_count.resize(query_vertices_count);

    bool is_steady = false;
    while (!is_steady) {
        is_steady = true;
        for (ui i = 0; i < query_vertices_count; ++i) {
            ui u = i;

            ui u_nbr_cnt;
            const ui *u_nbrs = query_graph->getVertexNeighbors(u, u_nbr_cnt);

            ui valid_flag = 0;
            for (ui j = 0; j < u_nbr_cnt; ++j) {
                ui u_nbr = u_nbrs[j];

                for (ui k = 0; k < candidates_count[u_nbr]; ++k) {
                    ui v = candidates_copy[u_nbr][k];

                    if (v == INVALID_VERTEX_ID)
                        continue;

                    ui v_nbr_cnt;
                    const ui *v_nbrs = data_graph->getVertexNeighbors(v, v_nbr_cnt);

                    for (ui l = 0; l < v_nbr_cnt; ++l) {
                        ui v_nbr = v_nbrs[l];

                        if (flag[v_nbr] == valid_flag) {
                            flag[v_nbr] += 1;

                            if (valid_flag == 0) {
                                updated_flag.push_back(v_nbr);
                            }
                        }
                    }
                }
                valid_flag += 1;
            }

            for (ui j = 0; j < candidates_count[u]; ++j) {
                ui v = candidates_copy[u][j];

                if (v == INVALID_VERTEX_ID)
                    continue;

                if (flag[v] != valid_flag) {
                    candidates_copy[u][j] = INVALID_VERTEX_ID;
                    is_steady = false;
                }
            }

            for (auto v : updated_flag) {
                flag[v] = 0;
            }
            updated_flag.clear();
        }
    }

    double sum = 0;
    for (ui i = 0; i < query_vertices_count; ++i) {
        ui u = i;
        ui negative_count = 0;
        for (ui j = 0; j < candidates_count[u]; ++j) {
            ui v = candidates_copy[u][j];

            if (v == INVALID_VERTEX_ID)
                negative_count += 1;
        }

        per_query_vertex_false_positive_ratio[u] =
                (negative_count) / (double) candidates_count[u];
        sum += per_query_vertex_false_positive_ratio[u];
        optimal_candidates_count[u] = candidates_count[u] - negative_count;
    }

    return sum / query_vertices_count;
}
