//
// Created by ssunah on 11/20/18.
//

#include "BuildTable.h"
#include <vector>
#include <algorithm>

#if ENABLE_QFLITER == 1
BSRGraph*** BuildTable::qfliter_bsr_graph_;
#endif

void BuildTable::buildTables(const Graph *data_graph, const Graph *query_graph, ui **candidates, ui *candidates_count,
                             Edges ***edge_matrix) {
    ui query_vertices_num = query_graph->getVerticesCount();
    ui* flag = new ui[data_graph->getVerticesCount()];
    ui* updated_flag = new ui[data_graph->getVerticesCount()];
    std::fill(flag, flag + data_graph->getVerticesCount(), 0);

    for (ui i = 0; i < query_vertices_num; ++i) {
        for (ui j = 0; j < query_vertices_num; ++j) {
            edge_matrix[i][j] = NULL;
        }
    }

    std::vector<VertexID> build_table_order(query_vertices_num);
    for (ui i = 0; i < query_vertices_num; ++i) {
        build_table_order[i] = i;
    }

    std::sort(build_table_order.begin(), build_table_order.end(), [query_graph](VertexID l, VertexID r) {
        if (query_graph->getVertexDegree(l) == query_graph->getVertexDegree(r)) {
            return l < r;
        }
        return query_graph->getVertexDegree(l) > query_graph->getVertexDegree(r);
    });

    std::vector<ui> temp_edges(data_graph->getEdgesCount() * 2);

    for (auto u : build_table_order) {
        ui u_nbrs_count;
        const VertexID* u_nbrs = query_graph->getVertexNeighbors(u, u_nbrs_count);

        ui updated_flag_count = 0;

        for (ui i = 0; i < u_nbrs_count; ++i) {
            VertexID u_nbr = u_nbrs[i];

            if (edge_matrix[u][u_nbr] != NULL)
                continue;

            if (updated_flag_count == 0) {
                for (ui j = 0; j < candidates_count[u]; ++j) {
                    VertexID v = candidates[u][j];
                    flag[v] = j + 1;
                    updated_flag[updated_flag_count++] = v;
                }
            }

            edge_matrix[u_nbr][u] = new Edges;
            edge_matrix[u_nbr][u]->vertex_count_ = candidates_count[u_nbr];
            edge_matrix[u_nbr][u]->offset_ = new ui[candidates_count[u_nbr] + 1];

            edge_matrix[u][u_nbr] = new Edges;
            edge_matrix[u][u_nbr]->vertex_count_ = candidates_count[u];
            edge_matrix[u][u_nbr]->offset_ = new ui[candidates_count[u] + 1];
            std::fill(edge_matrix[u][u_nbr]->offset_, edge_matrix[u][u_nbr]->offset_ + candidates_count[u] + 1, 0);

            ui local_edge_count = 0;
            ui local_max_degree = 0;

            for (ui j = 0; j < candidates_count[u_nbr]; ++j) {
                VertexID v = candidates[u_nbr][j];
                edge_matrix[u_nbr][u]->offset_[j] = local_edge_count;

                ui v_nbrs_count;
                const VertexID* v_nbrs = data_graph->getVertexNeighbors(v, v_nbrs_count);

                ui local_degree = 0;

                for (ui k = 0; k < v_nbrs_count; ++k) {
                    VertexID v_nbr = v_nbrs[k];

                    if (flag[v_nbr] != 0) {
                        ui position = flag[v_nbr] - 1;
                        temp_edges[local_edge_count++] = position;
                        edge_matrix[u][u_nbr]->offset_[position + 1] += 1;
                        local_degree += 1;
                    }
                }

                if (local_degree > local_max_degree) {
                    local_max_degree = local_degree;
                }
            }

            edge_matrix[u_nbr][u]->offset_[candidates_count[u_nbr]] = local_edge_count;
            edge_matrix[u_nbr][u]->max_degree_ = local_max_degree;
            edge_matrix[u_nbr][u]->edge_count_ = local_edge_count;
            edge_matrix[u_nbr][u]->edge_ = new ui[local_edge_count];
            std::copy(temp_edges.begin(), temp_edges.begin() + local_edge_count, edge_matrix[u_nbr][u]->edge_);

            edge_matrix[u][u_nbr]->edge_count_ = local_edge_count;
            edge_matrix[u][u_nbr]->edge_ = new ui[local_edge_count];

            local_max_degree = 0;
            for (ui j = 1; j <= candidates_count[u]; ++j) {
                if (edge_matrix[u][u_nbr]->offset_[j] > local_max_degree) {
                    local_max_degree = edge_matrix[u][u_nbr]->offset_[j];
                }
                edge_matrix[u][u_nbr]->offset_[j] += edge_matrix[u][u_nbr]->offset_[j - 1];
            }

            edge_matrix[u][u_nbr]->max_degree_ = local_max_degree;

            for (ui j = 0; j < candidates_count[u_nbr]; ++j) {
                ui begin = j;
                for (ui k = edge_matrix[u_nbr][u]->offset_[begin]; k < edge_matrix[u_nbr][u]->offset_[begin + 1]; ++k) {
                    ui end = edge_matrix[u_nbr][u]->edge_[k];

                    edge_matrix[u][u_nbr]->edge_[edge_matrix[u][u_nbr]->offset_[end]++] = begin;
                }
            }

            for (ui j = candidates_count[u]; j >= 1; --j) {
                edge_matrix[u][u_nbr]->offset_[j] = edge_matrix[u][u_nbr]->offset_[j - 1];
            }
            edge_matrix[u][u_nbr]->offset_[0] = 0;
        }

        for (ui i = 0; i < updated_flag_count; ++i) {
            VertexID v = updated_flag[i];
            flag[v] = 0;
        }
    }

#if ENABLE_QFLITER == 1
    qfliter_bsr_graph_ = new BSRGraph**[query_vertices_num];
    for (ui i = 0; i < query_vertices_num; ++i) {
        qfliter_bsr_graph_[i] = new BSRGraph*[query_vertices_num];
        for (ui j = 0; j < query_vertices_num; ++j) {

            qfliter_bsr_graph_[i][j] = new BSRGraph[query_vertices_num];

            if (edge_matrix[i][j] != NULL) {
                qfliter_bsr_graph_[i][j]->load(edge_matrix[i][j]->vertex_count_,
                                               edge_matrix[i][j]->offset_, edge_matrix[i][j]->offset_,
                                               edge_matrix[i][j]->edge_);
            }
        }
    }
#endif
}

void BuildTable::printTableCardinality(const Graph *query_graph, Edges ***edge_matrix) {
    std::vector<std::pair<std::pair<VertexID, VertexID >, ui>> core_edges;
    std::vector<std::pair<std::pair<VertexID, VertexID >, ui>> tree_edges;
    std::vector<std::pair<std::pair<VertexID, VertexID >, ui>> leaf_edges;

    ui query_vertices_num = query_graph->getVerticesCount();

    double sum = 0;
    for (ui i = 0; i < query_vertices_num; ++i) {
        VertexID begin_vertex = i;

        for (ui j = i + 1; j < query_vertices_num; ++j) {
            VertexID end_vertex = j;

            if (query_graph->checkEdgeExistence(begin_vertex, end_vertex)) {
                ui cardinality = (*edge_matrix[begin_vertex][end_vertex]).edge_count_;
                sum += cardinality;
                if (query_graph->getCoreValue(begin_vertex) > 1 && query_graph->getCoreValue(end_vertex) > 1) {
                    core_edges.emplace_back(std::make_pair(std::make_pair(begin_vertex, end_vertex), cardinality));
                }
                else if (query_graph->getVertexDegree(begin_vertex) == 1 || query_graph->getVertexDegree(end_vertex) == 1) {
                    leaf_edges.emplace_back(std::make_pair(std::make_pair(begin_vertex, end_vertex), cardinality));
                }
                else {
                    tree_edges.emplace_back(std::make_pair(std::make_pair(begin_vertex, end_vertex), cardinality));
                }
            }
        }
    }

    printf("Index Info: CoreTable(%zu), TreeTable(%zu), LeafTable(%zu)\n", core_edges.size(), tree_edges.size(), leaf_edges.size());

    for (auto table_info : core_edges) {
        printf("CoreTable %u-%u: %u\n", table_info.first.first, table_info.first.second, table_info.second);
    }

    for (auto table_info : tree_edges) {
        printf("TreeTable %u-%u: %u\n", table_info.first.first, table_info.first.second, table_info.second);
    }

    for (auto table_info : leaf_edges) {
        printf("LeafTable %u-%u: %d\n", table_info.first.first, table_info.first.second, table_info.second);
    }

    printf("Total Cardinality: %.1lf\n", sum);
}

void BuildTable::printTableCardinality(const Graph *query_graph, TreeNode *tree, ui *order,
                                      std::vector<std::unordered_map<VertexID, std::vector<VertexID >>> &TE_Candidates,
                                      std::vector<std::vector<std::unordered_map<VertexID, std::vector<VertexID>>>> &NTE_Candidates) {
    std::vector<std::pair<std::pair<VertexID, VertexID >, ui>> core_edges;
    std::vector<std::pair<std::pair<VertexID, VertexID >, ui>> tree_edges;
    std::vector<std::pair<std::pair<VertexID, VertexID >, ui>> leaf_edges;

    ui query_vertices_num = query_graph->getVerticesCount();

    double sum = 0;
    for (ui i = 1; i < query_vertices_num; ++i) {
        VertexID begin_vertex = order[i];
        TreeNode& node = tree[begin_vertex];

        {
            VertexID end_vertex = node.parent_;
            ui cardinality = 0;
            for (auto iter = TE_Candidates[begin_vertex].begin(); iter != TE_Candidates[begin_vertex].end(); ++iter) {
                cardinality += iter->second.size();
            }

            sum += cardinality;
            if (query_graph->getCoreValue(begin_vertex) > 1 && query_graph->getCoreValue(end_vertex) > 1) {
                core_edges.emplace_back(std::make_pair(std::make_pair(begin_vertex, end_vertex), cardinality));
            }
            else if (query_graph->getVertexDegree(begin_vertex) == 1 || query_graph->getVertexDegree(end_vertex) == 1) {
                leaf_edges.emplace_back(std::make_pair(std::make_pair(begin_vertex, end_vertex), cardinality));
            }
            else {
                tree_edges.emplace_back(std::make_pair(std::make_pair(begin_vertex, end_vertex), cardinality));
            }
        }
        {
            for (ui j = 0; j < node.bn_count_; ++j) {
                VertexID end_vertex = node.bn_[j];

                ui cardinality = 0;
                for (auto iter = NTE_Candidates[begin_vertex][end_vertex].begin(); iter != NTE_Candidates[begin_vertex][end_vertex].end(); ++iter) {
                    cardinality += iter->second.size();
                }
                sum += cardinality;
                if (query_graph->getCoreValue(begin_vertex) > 1 && query_graph->getCoreValue(end_vertex) > 1) {
                    core_edges.emplace_back(std::make_pair(std::make_pair(begin_vertex, end_vertex), cardinality));
                }
                else if (query_graph->getVertexDegree(begin_vertex) == 1 || query_graph->getVertexDegree(end_vertex) == 1) {
                    leaf_edges.emplace_back(std::make_pair(std::make_pair(begin_vertex, end_vertex), cardinality));
                }
                else {
                    tree_edges.emplace_back(std::make_pair(std::make_pair(begin_vertex, end_vertex), cardinality));
                }
            }
        }
    }

    printf("Index Info: CoreTable(%zu), TreeTable(%zu), LeafTable(%zu)\n", core_edges.size(), tree_edges.size(), leaf_edges.size());

    for (auto table_info : core_edges) {
        printf("CoreTable %u-%u: %u\n", table_info.first.first, table_info.first.second, table_info.second);
    }

    for (auto table_info : tree_edges) {
        printf("TreeTable %u-%u: %u\n", table_info.first.first, table_info.first.second, table_info.second);
    }

    for (auto table_info : leaf_edges) {
        printf("LeafTable %u-%u: %d\n", table_info.first.first, table_info.first.second, table_info.second);
    }

    printf("Total Cardinality: %.1lf\n", sum);
}

void BuildTable::printTableInfo(const Graph *query_graph, Edges ***edge_matrix) {
    ui query_vertices_num = query_graph->getVerticesCount();
    printf("Index Info:\n");

    for (ui i = 0; i < query_vertices_num; ++i) {
        VertexID begin_vertex = i;
        for (ui j = 0; j < query_vertices_num; ++j) {
            VertexID end_vertex = j;

            if (begin_vertex < end_vertex && query_graph->checkEdgeExistence(begin_vertex, end_vertex)) {
                printTableInfo(begin_vertex, end_vertex, edge_matrix);
                printTableInfo(end_vertex, begin_vertex, edge_matrix);
            }

        }
    }
}

void BuildTable::printTableInfo(VertexID begin_vertex, VertexID end_vertex, Edges ***edge_matrix) {
    Edges& edge = *edge_matrix[begin_vertex][end_vertex];
    ui vertex_count = edge.vertex_count_;
    ui edge_count = edge.edge_count_;
    ui max_degree = edge.max_degree_;
    double average = (double) edge_count /vertex_count;

    printf("R(%u, %u): Edge Count: %u, Vertex Count: %u, Max Degree: %u, Avg Degree: %.4lf\n",
           begin_vertex, end_vertex, edge_count, vertex_count, max_degree, average);
}

size_t BuildTable::computeMemoryCostInBytes(const Graph *query_graph, ui *candidates_count, Edges ***edge_matrix) {
    size_t memory_cost_in_bytes = 0;
    size_t per_element_size = sizeof(ui);

    ui query_vertices_num = query_graph->getVerticesCount();
    for (ui i = 0; i < query_vertices_num; ++i) {
        memory_cost_in_bytes += candidates_count[i] * per_element_size;
    }

    for (ui i = 0; i < query_vertices_num; ++i) {
        VertexID begin_vertex = i;
        for (ui j = 0; j < query_vertices_num; ++j) {
            VertexID end_vertex = j;

            if (begin_vertex < end_vertex && query_graph->checkEdgeExistence(begin_vertex, end_vertex)) {
                Edges& edge = *edge_matrix[begin_vertex][end_vertex];
                memory_cost_in_bytes += edge.edge_count_ * per_element_size + edge.vertex_count_ * per_element_size;

                Edges& reverse_edge = *edge_matrix[end_vertex][begin_vertex];
                memory_cost_in_bytes += reverse_edge.edge_count_ * per_element_size + reverse_edge.vertex_count_ * per_element_size;
            }

        }
    }

    return memory_cost_in_bytes;
}

size_t BuildTable::computeMemoryCostInBytes(const Graph *query_graph, ui *candidates_count, ui *order, TreeNode *tree,
                                            std::vector<std::unordered_map<VertexID, std::vector<VertexID >>> &TE_Candidates,
                                            std::vector<std::vector<std::unordered_map<VertexID, std::vector<VertexID>>>> &NTE_Candidates) {
    size_t memory_cost_in_bytes = 0;
    size_t per_element_size = sizeof(ui);

    ui query_vertices_num = query_graph->getVerticesCount();
    for (ui i = 0; i < query_vertices_num; ++i) {
        memory_cost_in_bytes += candidates_count[i] * per_element_size;
    }

    for (ui i = 1; i < query_vertices_num; ++i) {
        VertexID u = order[i];
        TreeNode& u_node = tree[u];

        // NTE_Candidates
        for (ui j = 0; j < u_node.bn_count_; ++j) {
            VertexID u_bn = u_node.bn_[j];
            memory_cost_in_bytes += NTE_Candidates[u][u_bn].size() * per_element_size;
            for (auto iter = NTE_Candidates[u][u_bn].begin(); iter != NTE_Candidates[u][u_bn].end(); ++iter) {
                memory_cost_in_bytes += iter->second.size() * per_element_size;
            }
        }

        // TE_Candidates
        memory_cost_in_bytes += TE_Candidates[u].size() * per_element_size;
        for (auto iter = TE_Candidates[u].begin(); iter != TE_Candidates[u].end(); ++iter) {
            memory_cost_in_bytes += iter->second.size() * per_element_size;
        }
    }

    return memory_cost_in_bytes;
}


