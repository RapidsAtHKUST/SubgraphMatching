//
// Created by ssunah on 11/21/18.
//

#include "GenerateFilteringPlan.h"
#include "FilterVertices.h"
#include <queue>
#include <utility/graphoperations.h>

void GenerateFilteringPlan::generateTSOFilterPlan(const Graph *data_graph, const Graph *query_graph, TreeNode *&tree,
                                                   VertexID *&order) {
    VertexID start_vertex = selectTSOFilterStartVertex(data_graph, query_graph);
    VertexID* bfs_order;
    GraphOperations::bfsTraversal(query_graph, start_vertex, tree, bfs_order);
    GraphOperations::dfsTraversal(tree, start_vertex, query_graph->getVerticesCount(), order);
    delete[] bfs_order;
}

void GenerateFilteringPlan::generateCFLFilterPlan(const Graph *data_graph, const Graph *query_graph, TreeNode *&tree,
                                                  VertexID *&order, int &level_count, ui *&level_offset) {
    ui query_vertices_num = query_graph->getVerticesCount();
    VertexID start_vertex = selectCFLFilterStartVertex(data_graph, query_graph);
    GraphOperations::bfsTraversal(query_graph, start_vertex, tree, order);

    std::vector<ui> order_index(query_vertices_num);
    for (ui i = 0; i < query_vertices_num; ++i) {
        VertexID query_vertex = order[i];
        order_index[query_vertex] = i;
    }

    level_count = -1;
    level_offset = new ui[query_vertices_num + 1];

    for (ui i = 0; i < query_vertices_num; ++i) {
        VertexID u = order[i];
        tree[u].under_level_count_ = 0;
        tree[u].bn_count_ = 0;
        tree[u].fn_count_ = 0;

        if (tree[u].level_ != level_count) {
            level_count += 1;
            level_offset[level_count] = 0;
        }

        level_offset[level_count] += 1;

        ui u_nbrs_count;
        const VertexID* u_nbrs = query_graph->getVertexNeighbors(u, u_nbrs_count);
        for (ui j = 0; j < u_nbrs_count; ++j) {
            VertexID u_nbr = u_nbrs[j];

            if (tree[u].level_ == tree[u_nbr].level_) {
                if (order_index[u_nbr] < order_index[u]) {
                    tree[u].bn_[tree[u].bn_count_++] = u_nbr;
                }
                else {
                    tree[u].fn_[tree[u].fn_count_++] = u_nbr;
                }
            }
            else if (tree[u].level_ > tree[u_nbr].level_) {
                tree[u].bn_[tree[u].bn_count_++] = u_nbr;
            }
            else {
                tree[u].under_level_[tree[u].under_level_count_++] = u_nbr;
            }
        }
    }

    level_count += 1;

    ui prev_value = 0;
    for (ui i = 1; i <= level_count; ++i) {
        ui temp = level_offset[i];
        level_offset[i] = level_offset[i - 1] + prev_value;
        prev_value = temp;
    }
    level_offset[0] = 0;
}

void GenerateFilteringPlan::generateDPisoFilterPlan(const Graph *data_graph, const Graph *query_graph, TreeNode *&tree,
                                                    VertexID *&order) {
    VertexID start_vertex = selectDPisoStartVertex(data_graph, query_graph);
    GraphOperations::bfsTraversal(query_graph, start_vertex, tree, order);

    ui query_vertices_num = query_graph->getVerticesCount();
    std::vector<ui> order_index(query_vertices_num);
    for (ui i = 0; i < query_vertices_num; ++i) {
        VertexID query_vertex = order[i];
        order_index[query_vertex] = i;
    }

    for (ui i = 0; i < query_vertices_num; ++i) {
        VertexID u = order[i];
        tree[u].under_level_count_ = 0;
        tree[u].bn_count_ = 0;
        tree[u].fn_count_ = 0;

        ui u_nbrs_count;
        const VertexID* u_nbrs = query_graph->getVertexNeighbors(u, u_nbrs_count);
        for (ui j = 0; j < u_nbrs_count; ++j) {
            VertexID u_nbr = u_nbrs[j];
            if (order_index[u_nbr] < order_index[u]) {
                tree[u].bn_[tree[u].bn_count_++] = u_nbr;
            }
            else {
                tree[u].fn_[tree[u].fn_count_++] = u_nbr;
            }
        }
    }
}

void GenerateFilteringPlan::generateCECIFilterPlan(const Graph *data_graph, const Graph *query_graph, TreeNode *&tree,
                                                   VertexID *&order) {
    VertexID start_vertex = selectCECIStartVertex(data_graph, query_graph);
    GraphOperations::bfsTraversal(query_graph, start_vertex, tree, order);

    ui query_vertices_num = query_graph->getVerticesCount();
    std::vector<ui> order_index(query_vertices_num);
    for (ui i = 0; i < query_vertices_num; ++i) {
        VertexID query_vertex = order[i];
        order_index[query_vertex] = i;
    }

    for (ui i = 0; i < query_vertices_num; ++i) {
        VertexID u = order[i];
        tree[u].under_level_count_ = 0;
        tree[u].bn_count_ = 0;
        tree[u].fn_count_ = 0;

        ui u_nbrs_count;
        const VertexID* u_nbrs = query_graph->getVertexNeighbors(u, u_nbrs_count);
        for (ui j = 0; j < u_nbrs_count; ++j) {
            VertexID u_nbr = u_nbrs[j];
            if (u_nbr != tree[u].parent_ && order_index[u_nbr] < order_index[u]) {
                tree[u].bn_[tree[u].bn_count_++] = u_nbr;
                tree[u_nbr].fn_[tree[u_nbr].fn_count_++] = u;
            }
        }
    }
}


VertexID GenerateFilteringPlan::selectTSOFilterStartVertex(const Graph *data_graph, const Graph *query_graph) {
    auto rank_compare = [](std::pair<VertexID, double> l, std::pair<VertexID, double> r) {
        return l.second < r.second;
    };
    // Maximum priority queue.
    std::priority_queue<std::pair<VertexID, double>, std::vector<std::pair<VertexID, double>>, decltype(rank_compare)> rank_queue(rank_compare);

    // Compute the ranking.
    for (ui i = 0; i < query_graph->getVerticesCount(); ++i) {
        VertexID query_vertex = i;
        LabelID label = query_graph->getVertexLabel(query_vertex);
        ui degree = query_graph->getVertexDegree(query_vertex);
        ui frequency = data_graph->getLabelsFrequency(label);
        double rank = frequency / (double)degree;
        rank_queue.push(std::make_pair(query_vertex, rank));
    }

    // Keep the top-3.
    while (rank_queue.size() > 3) {
        rank_queue.pop();
    }

    // Pick the one with the smallest number of candidates.
    VertexID start_vertex = 0;
    ui min_candidates_num = data_graph->getGraphMaxLabelFrequency() + 1;
    while (!rank_queue.empty()) {
        VertexID query_vertex = rank_queue.top().first;

        if (rank_queue.size() == 1) {
            ui count;
            FilterVertices::computeCandidateWithNLF(data_graph, query_graph, query_vertex, count);
            if (count < min_candidates_num) {
                start_vertex = query_vertex;
            }
        }
        else {
            LabelID label = query_graph->getVertexLabel(query_vertex);
            ui frequency = data_graph->getLabelsFrequency(label);
            if (frequency / (double)data_graph->getVerticesCount() <= 0.05) {
                ui count;
                FilterVertices::computeCandidateWithNLF(data_graph, query_graph, query_vertex, count);
                if (count < min_candidates_num) {
                    start_vertex = query_vertex;
                    min_candidates_num = count;
                }
            }
        }
        rank_queue.pop();
    }

    return start_vertex;
}

VertexID GenerateFilteringPlan::selectCFLFilterStartVertex(const Graph *data_graph, const Graph *query_graph) {
    auto rank_compare = [](std::pair<VertexID, double> l, std::pair<VertexID, double> r) {
        return l.second < r.second;
    };

    std::priority_queue<std::pair<VertexID, double>, std::vector<std::pair<VertexID, double>>, decltype(rank_compare)> rank_queue(rank_compare);

    // Compute the ranking.
    for (ui i = 0; i < query_graph->getVerticesCount(); ++i) {
        VertexID query_vertex = i;

        if (query_graph->get2CoreSize() == 0 || query_graph->getCoreValue(query_vertex) > 1) {
            LabelID label = query_graph->getVertexLabel(query_vertex);
            ui degree = query_graph->getVertexDegree(query_vertex);
            ui frequency = data_graph->getLabelsFrequency(label);
            double rank = frequency / (double) degree;
            rank_queue.push(std::make_pair(query_vertex, rank));
        }
    }

    // Keep the top-3.
    while (rank_queue.size() > 3) {
        rank_queue.pop();
    }

    VertexID start_vertex = 0;
    double min_score = data_graph->getGraphMaxLabelFrequency() + 1;

    while (!rank_queue.empty()) {
        VertexID query_vertex = rank_queue.top().first;
        ui count;
        FilterVertices::computeCandidateWithNLF(data_graph, query_graph, query_vertex, count);
        double cur_score = count / (double) query_graph->getVertexDegree(query_vertex);

        if (cur_score < min_score) {
            start_vertex = query_vertex;
            min_score = cur_score;
        }
        rank_queue.pop();
    }

    return start_vertex;
}

VertexID GenerateFilteringPlan::selectDPisoStartVertex(const Graph *data_graph, const Graph *query_graph) {
    double min_score = data_graph->getVerticesCount();
    VertexID start_vertex = 0;

    for (ui i = 0; i < query_graph->getVerticesCount(); ++i) {
        ui degree = query_graph->getVertexDegree(i);
        if (degree <= 1)
            continue;

        ui count = 0;
        FilterVertices::computeCandidateWithLDF(data_graph, query_graph, i, count);
        double cur_score = count / (double)degree;
        if (cur_score < min_score) {
            min_score = cur_score;
            start_vertex = i;
        }
    }

    return start_vertex;
}

VertexID GenerateFilteringPlan::selectCECIStartVertex(const Graph *data_graph, const Graph *query_graph) {
    double min_score = data_graph->getVerticesCount();
    VertexID start_vertex = 0;

    for (ui i = 0; i < query_graph->getVerticesCount(); ++i) {
        ui degree = query_graph->getVertexDegree(i);
        ui count = 0;
        FilterVertices::computeCandidateWithNLF(data_graph, query_graph, i, count);
        double cur_score = count / (double)degree;
        if (cur_score < min_score) {
            min_score = cur_score;
            start_vertex = i;
        }
    }

    return start_vertex;
}