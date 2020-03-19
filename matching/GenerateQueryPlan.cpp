//
// Created by ssunah on 11/20/18.
//

#include "GenerateQueryPlan.h"
#include "utility/computesetintersection.h"
#include <vector>
#include <limits>
#include <cassert>
#include <algorithm>
#include <utility/graphoperations.h>

void GenerateQueryPlan::generateGQLQueryPlan(const Graph *data_graph, const Graph *query_graph, ui *candidates_count,
                                             ui *&order, ui *&pivot) {
     /**
      * Select the vertex v such that (1) v is adjacent to the selected vertices; and (2) v has the minimum number of candidates.
      */
     std::vector<bool> visited_vertices(query_graph->getVerticesCount(), false);
     std::vector<bool> adjacent_vertices(query_graph->getVerticesCount(), false);
     order = new ui[query_graph->getVerticesCount()];
     pivot = new ui[query_graph->getVerticesCount()];

     VertexID start_vertex = selectGQLStartVertex(query_graph, candidates_count);
     order[0] = start_vertex;
     updateValidVertices(query_graph, start_vertex, visited_vertices, adjacent_vertices);

     for (ui i = 1; i < query_graph->getVerticesCount(); ++i) {
          VertexID next_vertex;
          ui min_value = data_graph->getVerticesCount() + 1;
          for (ui j = 0; j < query_graph->getVerticesCount(); ++j) {
               VertexID cur_vertex = j;

               if (!visited_vertices[cur_vertex] && adjacent_vertices[cur_vertex]) {
                    if (candidates_count[cur_vertex] < min_value) {
                         min_value = candidates_count[cur_vertex];
                         next_vertex = cur_vertex;
                    }
                    else if (candidates_count[cur_vertex] == min_value && query_graph->getVertexDegree(cur_vertex) > query_graph->getVertexDegree(next_vertex)) {
                         next_vertex = cur_vertex;
                    }
               }
          }
          updateValidVertices(query_graph, next_vertex, visited_vertices, adjacent_vertices);
          order[i] = next_vertex;
     }

     // Pick a pivot randomly.
     for (ui i = 1; i < query_graph->getVerticesCount(); ++i) {
         VertexID u = order[i];
         for (ui j = 0; j < i; ++j) {
             VertexID cur_vertex = order[j];
             if (query_graph->checkEdgeExistence(u, cur_vertex)) {
                 pivot[i] = cur_vertex;
                 break;
             }
         }
     }
}

void GenerateQueryPlan::generateQSIQueryPlan(const Graph *data_graph, const Graph *query_graph, Edges ***edge_matrix,
                                             ui *&order, ui *&pivot) {
     /**
      * Generate a minimum spanning tree.
      */
     std::vector<bool> visited_vertices(query_graph->getVerticesCount(), false);
     std::vector<bool> adjacent_vertices(query_graph->getVerticesCount(), false);
     order = new ui[query_graph->getVerticesCount()];
     pivot = new ui[query_graph->getVerticesCount()];

     std::pair<VertexID, VertexID> start_edge = selectQSIStartEdge(query_graph, edge_matrix);
     order[0] = start_edge.first;
     order[1] = start_edge.second;
     pivot[1] = start_edge.first;
     updateValidVertices(query_graph, start_edge.first, visited_vertices, adjacent_vertices);
     updateValidVertices(query_graph, start_edge.second, visited_vertices, adjacent_vertices);

     for (ui l = 2; l < query_graph->getVerticesCount(); ++l) {
          ui min_value = std::numeric_limits<ui>::max();
          ui max_degree = 0;
          ui max_adjacent_selected_vertices = 0;
          std::pair<VertexID, VertexID> selected_edge;
          for (ui i = 0; i < query_graph->getVerticesCount(); ++i) {
               VertexID begin_vertex = i;
               if (visited_vertices[begin_vertex]) {
                    ui nbrs_cnt;
                    const VertexID *nbrs = query_graph->getVertexNeighbors(begin_vertex, nbrs_cnt);
                    for (ui j = 0; j < nbrs_cnt; ++j) {
                         VertexID end_vertex = nbrs[j];

                         if (!visited_vertices[end_vertex]) {
                              ui cur_value = (*edge_matrix[begin_vertex][end_vertex]).edge_count_;
                              ui cur_degree = query_graph->getVertexDegree(end_vertex);
                              ui adjacent_selected_vertices = 0;
                              ui end_vertex_nbrs_count;
                              const VertexID *end_vertex_nbrs = query_graph->getVertexNeighbors(end_vertex,
                                                                                                end_vertex_nbrs_count);

                              for (ui k = 0; k < end_vertex_nbrs_count; ++k) {
                                   VertexID end_vertex_nbr = end_vertex_nbrs[k];

                                   if (visited_vertices[end_vertex_nbr]) {
                                        adjacent_selected_vertices += 1;
                                   }
                              }

                              if (cur_value < min_value || (cur_value == min_value && adjacent_selected_vertices < max_adjacent_selected_vertices)
                                      || (cur_value == min_value && adjacent_selected_vertices == max_adjacent_selected_vertices && cur_degree > max_degree)) {
                                   selected_edge = std::make_pair(begin_vertex, end_vertex);
                                   min_value = cur_value;
                                   max_degree = cur_degree;
                                   max_adjacent_selected_vertices = adjacent_selected_vertices;
                              }
                         }
                    }
               }
          }

          order[l] = selected_edge.second;
          pivot[l] = selected_edge.first;
          updateValidVertices(query_graph, selected_edge.second, visited_vertices, adjacent_vertices);
     }
}

VertexID GenerateQueryPlan::selectGQLStartVertex(const Graph *query_graph, ui *candidates_count) {
    /**
     * Select the vertex with the minimum number of candidates as the start vertex.
     * Tie Handling:
     *  1. degree
     *  2. label id
     */

     ui start_vertex = 0;

     for (ui i = 1; i < query_graph->getVerticesCount(); ++i) {
          VertexID cur_vertex = i;

          if (candidates_count[cur_vertex] < candidates_count[start_vertex]) {
               start_vertex = cur_vertex;
          }
          else if (candidates_count[cur_vertex] == candidates_count[start_vertex]
                   && query_graph->getVertexDegree(cur_vertex) > query_graph->getVertexDegree(start_vertex)) {
               start_vertex = cur_vertex;
          }
     }

     return start_vertex;
}

void GenerateQueryPlan::updateValidVertices(const Graph *query_graph, VertexID query_vertex, std::vector<bool> &visited,
                                            std::vector<bool> &adjacent) {
     visited[query_vertex] = true;
     ui nbr_cnt;
     const ui* nbrs = query_graph->getVertexNeighbors(query_vertex, nbr_cnt);

     for (ui i = 0; i < nbr_cnt; ++i) {
          ui nbr = nbrs[i];
          adjacent[nbr] = true;
     }
}

std::pair<VertexID, VertexID> GenerateQueryPlan::selectQSIStartEdge(const Graph *query_graph, Edges ***edge_matrix) {
     /**
      * Select the edge with the minimum number of candidate edges.
      * Tie Handling:
      *  1. the sum of the degree values of the end points of an edge (prioritize the edge with smaller degree).
      *  2. label id
      */
     ui min_value = std::numeric_limits<ui>::max();
     ui min_degree_sum = std::numeric_limits<ui>::max();

     std::pair<VertexID, VertexID> start_edge = std::make_pair(0, 1);
     for (ui i = 0; i < query_graph->getVerticesCount(); ++i) {
          VertexID begin_vertex = i;
          ui nbrs_cnt;
          const ui* nbrs = query_graph->getVertexNeighbors(begin_vertex, nbrs_cnt);

          for (ui j = 0; j < nbrs_cnt; ++j) {
               VertexID end_vertex = nbrs[j];
               ui cur_value = (*edge_matrix[begin_vertex][end_vertex]).edge_count_;
               ui cur_degree_sum = query_graph->getVertexDegree(begin_vertex) + query_graph->getVertexDegree(end_vertex);

               if (cur_value < min_value || (cur_value == min_value
                   && (cur_degree_sum < min_degree_sum))) {
                    min_value = cur_value;
                    min_degree_sum = cur_degree_sum;

                    start_edge = query_graph->getVertexDegree(begin_vertex) < query_graph->getVertexDegree(end_vertex) ?
                            std::make_pair(end_vertex, begin_vertex) : std::make_pair(begin_vertex, end_vertex);
               }
          }
     }

     return start_edge;
}

void GenerateQueryPlan::generateTSOQueryPlan(const Graph *query_graph, Edges ***edge_matrix, ui *&order, ui *&pivot,
                                             TreeNode *tree, ui *dfs_order) {
    /**
     * Order the root to leaf paths according to the estimated number of embeddings and generate the matching order.
     */
    ui query_vertices_num = query_graph->getVerticesCount();
    std::vector<std::vector<ui>> paths;
    paths.reserve(query_vertices_num);

    std::vector<ui> single_path;
    single_path.reserve(query_vertices_num);

    generateRootToLeafPaths(tree, dfs_order[0], single_path, paths);

    std::vector<std::pair<double, std::vector<ui>*>> path_orders;

    for (std::vector<ui>& path : paths) {

        std::vector<size_t> estimated_embeddings_num;

        ui non_tree_edges_count = generateNoneTreeEdgesCount(query_graph, tree, path);
        estimatePathEmbeddsingsNum(path, edge_matrix, estimated_embeddings_num);
        double score = estimated_embeddings_num[0] / (double) (non_tree_edges_count + 1);
        path_orders.emplace_back(std::make_pair(score, &path));
    }

    std::sort(path_orders.begin(), path_orders.end(), [](std::pair<double, std::vector<ui>*> l, std::pair<double, std::vector<ui>*> r)
    { return l.first < r.first; });

    std::vector<bool> visited_vertices(query_vertices_num, false);
    order = new ui[query_vertices_num];
    pivot = new ui[query_vertices_num];

    ui count = 0;
    for (auto& path : path_orders) {
        for (ui i = 0; i < path.second->size(); ++i) {
            VertexID vertex = path.second->at(i);
            if (!visited_vertices[vertex]) {
                order[count] = vertex;
                if (i != 0) {
                    pivot[count] = path.second->at(i - 1);
                }
                count += 1;
                visited_vertices[vertex] = true;
            }
        }
    }
}

void GenerateQueryPlan::generateCFLQueryPlan(const Graph *data_graph, const Graph *query_graph, Edges ***edge_matrix,
                                             ui *&order, ui *&pivot, TreeNode *tree, ui *bfs_order, ui *candidates_count) {
    ui query_vertices_num = query_graph->getVerticesCount();
    VertexID root_vertex = bfs_order[0];
    order = new ui[query_vertices_num];
    pivot = new ui[query_vertices_num];
    std::vector<bool> visited_vertices(query_vertices_num, false);

    std::vector<std::vector<ui>> core_paths;
    std::vector<std::vector<std::vector<ui>>> forests;
    std::vector<ui> leaves;

    generateLeaves(query_graph, leaves);
    if (query_graph->getCoreValue(root_vertex) > 1) {
        std::vector<ui> temp_core_path;
        generateCorePaths(query_graph, tree, root_vertex, temp_core_path, core_paths);
        for (ui i = 0; i < query_vertices_num; ++i) {
            VertexID cur_vertex = i;
            if (query_graph->getCoreValue(cur_vertex) > 1) {
                std::vector<std::vector<ui>> temp_tree_paths;
                std::vector<ui> temp_tree_path;
                generateTreePaths(query_graph, tree, cur_vertex, temp_tree_path, temp_tree_paths);
                if (!temp_tree_paths.empty()) {
                    forests.emplace_back(temp_tree_paths);
                }
            }
        }
    }
    else {
        std::vector<std::vector<ui>> temp_tree_paths;
        std::vector<ui> temp_tree_path;
        generateTreePaths(query_graph, tree, root_vertex, temp_tree_path, temp_tree_paths);
        if (!temp_tree_paths.empty()) {
            forests.emplace_back(temp_tree_paths);
        }
    }

    // Order core paths.
    ui selected_vertices_count = 0;
    order[selected_vertices_count++] = root_vertex;
    visited_vertices[root_vertex] = true;

    if (!core_paths.empty()) {
        std::vector<std::vector<size_t>> paths_embededdings_num;
        std::vector<ui> paths_non_tree_edge_num;
        for (auto& path : core_paths) {
            ui non_tree_edge_num = generateNoneTreeEdgesCount(query_graph, tree, path);
            paths_non_tree_edge_num.push_back(non_tree_edge_num + 1);

            std::vector<size_t> path_embeddings_num;
            estimatePathEmbeddsingsNum(path, edge_matrix, path_embeddings_num);
            paths_embededdings_num.emplace_back(path_embeddings_num);
        }

        // Select the start path.
        double min_value = std::numeric_limits<double>::max();
        ui selected_path_index = 0;

        for (ui i = 0; i < core_paths.size(); ++i) {
            double cur_value = paths_embededdings_num[i][0] / (double) paths_non_tree_edge_num[i];

            if (cur_value < min_value) {
                min_value = cur_value;
                selected_path_index = i;
            }
        }


        for (ui i = 1; i < core_paths[selected_path_index].size(); ++i) {
            order[selected_vertices_count] = core_paths[selected_path_index][i];
            pivot[selected_vertices_count] = core_paths[selected_path_index][i - 1];
            selected_vertices_count += 1;
            visited_vertices[core_paths[selected_path_index][i]] = true;
        }

        core_paths.erase(core_paths.begin() + selected_path_index);
        paths_embededdings_num.erase(paths_embededdings_num.begin() + selected_path_index);
        paths_non_tree_edge_num.erase(paths_non_tree_edge_num.begin() + selected_path_index);

        while (!core_paths.empty()) {
            min_value = std::numeric_limits<double>::max();
            selected_path_index = 0;

            for (ui i = 0; i < core_paths.size(); ++i) {
                ui path_root_vertex_idx = 0;
                for (ui j = 0; j < core_paths[i].size(); ++j) {
                    VertexID cur_vertex = core_paths[i][j];

                    if (visited_vertices[cur_vertex])
                        continue;

                    path_root_vertex_idx = j - 1;
                    break;
                }

                double cur_value = paths_embededdings_num[i][path_root_vertex_idx] / (double)candidates_count[core_paths[i][path_root_vertex_idx]];
                if (cur_value < min_value) {
                    min_value = cur_value;
                    selected_path_index = i;
                }
            }

            for (ui i = 1; i < core_paths[selected_path_index].size(); ++i) {
                if (visited_vertices[core_paths[selected_path_index][i]])
                    continue;

                order[selected_vertices_count] = core_paths[selected_path_index][i];
                pivot[selected_vertices_count] = core_paths[selected_path_index][i - 1];
                selected_vertices_count += 1;
                visited_vertices[core_paths[selected_path_index][i]] = true;
            }

            core_paths.erase(core_paths.begin() + selected_path_index);
            paths_embededdings_num.erase(paths_embededdings_num.begin() + selected_path_index);
        }
    }

    // Order tree paths.
    for (auto& tree_paths : forests) {
        std::vector<std::vector<size_t>> paths_embededdings_num;
        for (auto& path : tree_paths) {
            std::vector<size_t> path_embeddings_num;
            estimatePathEmbeddsingsNum(path, edge_matrix, path_embeddings_num);
            paths_embededdings_num.emplace_back(path_embeddings_num);
        }

        while (!tree_paths.empty()) {
            double min_value = std::numeric_limits<double>::max();
            ui selected_path_index = 0;

            for (ui i = 0; i < tree_paths.size(); ++i) {
                ui path_root_vertex_idx = 0;
                for (ui j = 0; j < tree_paths[i].size(); ++j) {
                    VertexID cur_vertex = tree_paths[i][j];

                    if (visited_vertices[cur_vertex])
                        continue;

                    path_root_vertex_idx = j == 0 ? j : j - 1;
                    break;
                }

                double cur_value = paths_embededdings_num[i][path_root_vertex_idx] / (double)candidates_count[tree_paths[i][path_root_vertex_idx]];
                if (cur_value < min_value) {
                    min_value = cur_value;
                    selected_path_index = i;
                }
            }

            for (ui i = 0; i < tree_paths[selected_path_index].size(); ++i) {
                if (visited_vertices[tree_paths[selected_path_index][i]])
                    continue;

                order[selected_vertices_count] = tree_paths[selected_path_index][i];
                pivot[selected_vertices_count] = tree_paths[selected_path_index][i - 1];
                selected_vertices_count += 1;
                visited_vertices[tree_paths[selected_path_index][i]] = true;
            }

            tree_paths.erase(tree_paths.begin() + selected_path_index);
            paths_embededdings_num.erase(paths_embededdings_num.begin() + selected_path_index);
        }
    }

    // Order the leaves.
    while (!leaves.empty()) {
        double min_value = std::numeric_limits<double>::max();
        ui selected_leaf_index = 0;

        for (ui i = 0; i < leaves.size(); ++i) {
            VertexID vertex = leaves[i];
            double cur_value = candidates_count[vertex];

            if (cur_value < min_value) {
                min_value = cur_value;
                selected_leaf_index = i;
            }
        }

        if (!visited_vertices[leaves[selected_leaf_index]]) {
            order[selected_vertices_count] = leaves[selected_leaf_index];
            pivot[selected_vertices_count] = tree[leaves[selected_leaf_index]].parent_;
            selected_vertices_count += 1;
            visited_vertices[leaves[selected_leaf_index]] = true;
        }
        leaves.erase(leaves.begin() + selected_leaf_index);
    }
}

void GenerateQueryPlan::generateRootToLeafPaths(TreeNode *tree_node, VertexID cur_vertex, std::vector<ui> &cur_path,
                                                std::vector<std::vector<ui>> &paths) {
    TreeNode& cur_node = tree_node[cur_vertex];
    cur_path.push_back(cur_vertex);

    if (cur_node.children_count_ == 0) {
        paths.emplace_back(cur_path);
    }
    else {
        for (ui i = 0; i < cur_node.children_count_; ++i) {
            VertexID next_vertex = cur_node.children_[i];
            generateRootToLeafPaths(tree_node, next_vertex, cur_path, paths);
        }
    }

    cur_path.pop_back();
}

void GenerateQueryPlan::estimatePathEmbeddsingsNum(std::vector<ui> &path, Edges ***edge_matrix,
                                                   std::vector<size_t> &estimated_embeddings_num) {
    assert(path.size() > 1);
    std::vector<size_t> parent;
    std::vector<size_t> children;

    estimated_embeddings_num.resize(path.size() - 1);
    Edges& last_edge = *edge_matrix[path[path.size() - 2]][path[path.size() - 1]];
    children.resize(last_edge.vertex_count_);

    size_t sum = 0;
    for (ui i = 0; i < last_edge.vertex_count_; ++i) {
        children[i] = last_edge.offset_[i + 1] - last_edge.offset_[i];
        sum += children[i];
    }

    estimated_embeddings_num[path.size() - 2] = sum;

    for (int i = path.size() - 2; i >= 1; --i) {
        ui begin = path[i - 1];
        ui end = path[i];

        Edges& edge = *edge_matrix[begin][end];
        parent.resize(edge.vertex_count_);

        sum = 0;
        for (ui j = 0; j < edge.vertex_count_; ++j) {

            size_t local_sum = 0;
            for (ui k = edge.offset_[j]; k < edge.offset_[j + 1]; ++k) {
                ui nbr = edge.edge_[k];
                local_sum += children[nbr];
            }

            parent[j] = local_sum;
            sum += local_sum;
        }

        estimated_embeddings_num[i - 1] = sum;
        parent.swap(children);
    }
}

ui GenerateQueryPlan::generateNoneTreeEdgesCount(const Graph *query_graph, TreeNode *tree_node, std::vector<ui> &path) {
    ui non_tree_edge_count = query_graph->getVertexDegree(path[0]) - tree_node[path[0]].children_count_;

    for (ui i = 1; i < path.size(); ++i) {
        VertexID vertex = path[i];
        non_tree_edge_count += query_graph->getVertexDegree(vertex) - tree_node[vertex].children_count_ - 1;
    }

    return non_tree_edge_count;
}

void GenerateQueryPlan::generateCorePaths(const Graph *query_graph, TreeNode *tree_node, VertexID cur_vertex,
                                          std::vector<ui> &cur_core_path, std::vector<std::vector<ui>> &core_paths) {
    TreeNode& node = tree_node[cur_vertex];
    cur_core_path.push_back(cur_vertex);

    bool is_core_leaf = true;
    for (ui i = 0; i < node.children_count_; ++i) {
        VertexID child = node.children_[i];
        if (query_graph->getCoreValue(child) > 1) {
            generateCorePaths(query_graph, tree_node, child, cur_core_path, core_paths);
            is_core_leaf = false;
        }
    }

    if (is_core_leaf) {
        core_paths.emplace_back(cur_core_path);
    }
    cur_core_path.pop_back();
}

void GenerateQueryPlan::generateTreePaths(const Graph *query_graph, TreeNode *tree_node, VertexID cur_vertex,
                                          std::vector<ui> &cur_tree_path, std::vector<std::vector<ui>> &tree_paths) {
    TreeNode& node = tree_node[cur_vertex];
    cur_tree_path.push_back(cur_vertex);

    bool is_tree_leaf = true;
    for (ui i = 0; i < node.children_count_; ++i) {
        VertexID child = node.children_[i];
        if (query_graph->getVertexDegree(child) > 1) {
            generateTreePaths(query_graph, tree_node, child, cur_tree_path, tree_paths);
            is_tree_leaf = false;
        }
    }

    if (is_tree_leaf && cur_tree_path.size() > 1) {
        tree_paths.emplace_back(cur_tree_path);
    }
    cur_tree_path.pop_back();
}

void GenerateQueryPlan::generateLeaves(const Graph *query_graph, std::vector<ui> &leaves) {
    for (ui i = 0; i < query_graph->getVerticesCount(); ++i) {
        VertexID cur_vertex = i;
        if (query_graph->getVertexDegree(cur_vertex) == 1) {
            leaves.push_back(cur_vertex);
        }
    }
}

void GenerateQueryPlan::checkQueryPlanCorrectness(const Graph *query_graph, ui *order, ui *pivot) {
    ui query_vertices_num = query_graph->getVerticesCount();
    std::vector<bool> visited_vertices(query_vertices_num, false);
    // Check whether each query vertex is in the order.
    for (ui i = 0; i < query_vertices_num; ++i) {
        VertexID vertex = order[i];
        assert(vertex < query_vertices_num && vertex >= 0);

        visited_vertices[vertex] = true;
    }

    for (ui i = 0; i < query_vertices_num; ++i) {
        VertexID vertex = i;
        assert(visited_vertices[vertex]);
    }

    // Check whether the order is connected.

    std::fill(visited_vertices.begin(), visited_vertices.end(), false);
    visited_vertices[order[0]] = true;
    for (ui i = 1; i < query_vertices_num; ++i) {
        VertexID vertex = order[i];
        VertexID pivot_vertex = pivot[i];
        assert(query_graph->checkEdgeExistence(vertex, pivot_vertex));
        assert(visited_vertices[pivot_vertex]);
        visited_vertices[vertex] = true;
    }
}

void GenerateQueryPlan::printQueryPlan(const Graph *query_graph, ui *order) {
    ui query_vertices_num = query_graph->getVerticesCount();
    printf("Query Plan: ");
    for (ui i = 0; i < query_vertices_num; ++i) {
        printf("%u, ", order[i]);
    }
    printf("\n");

    printf("%u: N/A\n", order[0]);
    for (ui i = 1; i < query_vertices_num; ++i) {
        VertexID end_vertex = order[i];
        printf("%u: ", end_vertex);
        for (ui j = 0; j < i; ++j) {
            VertexID begin_vertex = order[j];
            if (query_graph->checkEdgeExistence(begin_vertex, end_vertex)) {
                printf("R(%u, %u), ", begin_vertex, end_vertex);
            }
        }
        printf("\n");
    }
}

void GenerateQueryPlan::printSimplifiedQueryPlan(const Graph *query_graph, ui *order) {
    ui query_vertices_num = query_graph->getVerticesCount();
    printf("Query Plan: ");
    for (ui i = 0; i < query_vertices_num; ++i) {
        printf("%u ", order[i]);
    }
    printf("\n");
}

void GenerateQueryPlan::generateDSPisoQueryPlan(const Graph *query_graph, Edges ***edge_matrix, ui *&order, ui *&pivot,
                                                TreeNode *tree, ui *bfs_order, ui *candidates_count, ui **&weight_array) {
    ui query_vertices_num = query_graph->getVerticesCount();
    order = new ui[query_vertices_num];
    pivot = new ui[query_vertices_num];

    for (ui i = 0; i < query_vertices_num; ++i) {
        order[i] = bfs_order[i];
    }

    for (ui i = 1; i < query_vertices_num; ++i) {
        pivot[i] = tree[order[i]].parent_;
    }

    // Compute weight array.
    weight_array = new ui*[query_vertices_num];
    for (ui i = 0; i < query_vertices_num; ++i) {
        weight_array[i] = new ui[candidates_count[i]];
        std::fill(weight_array[i], weight_array[i] + candidates_count[i], std::numeric_limits<ui>::max());
    }

    for (int i = query_vertices_num - 1; i >= 0; --i) {
        VertexID vertex = order[i];
        TreeNode& node = tree[vertex];
        bool set_to_one = true;

        for (ui j = 0; j < node.fn_count_; ++j) {
            VertexID child = node.fn_[j];
            TreeNode& child_node = tree[child];

            if (child_node.bn_count_ == 1) {
                set_to_one = false;
                Edges& cur_edge = *edge_matrix[vertex][child];
                for (ui k = 0; k < candidates_count[vertex]; ++k) {
                    ui cur_candidates_count = cur_edge.offset_[k + 1] - cur_edge.offset_[k];
                    ui* cur_candidates = cur_edge.edge_ + cur_edge.offset_[k];

                    ui weight = 0;

                    for (ui l = 0; l < cur_candidates_count; ++l) {
                        ui candidates = cur_candidates[l];
                        weight += weight_array[child][candidates];
                    }

                    if (weight < weight_array[vertex][k])
                        weight_array[vertex][k] = weight;
                }
            }
        }

        if (set_to_one) {
            std::fill(weight_array[vertex], weight_array[vertex] + candidates_count[vertex], 1);
        }
    }
}

void GenerateQueryPlan::generateCECIQueryPlan(const Graph *query_graph, TreeNode *tree, ui *bfs_order, ui *&order,
                                              ui *&pivot) {
    ui query_vertices_num = query_graph->getVerticesCount();
    order = new ui[query_vertices_num];
    pivot = new ui[query_vertices_num];

    for (ui i = 0; i < query_vertices_num; ++i) {
        order[i] = bfs_order[i];
    }

    for (ui i = 1; i < query_vertices_num; ++i) {
        pivot[i] = tree[order[i]].parent_;
    }
}

void GenerateQueryPlan::generateRIQueryPlan(const Graph *data_graph, const Graph *query_graph, ui *&order, ui *&pivot) {
    ui query_vertices_num = query_graph->getVerticesCount();
    order = new ui[query_vertices_num];
    pivot = new ui[query_vertices_num];

    std::vector<bool> visited(query_vertices_num, false);
    // Select the vertex with the maximum degree as the start vertex.
    order[0] = 0;
    for (ui i = 1; i < query_vertices_num; ++i) {
        if (query_graph->getVertexDegree(i) > query_graph->getVertexDegree(order[0])) {
            order[0] = i;
        }
    }
    visited[order[0]] = true;
    // Order vertices.
    std::vector<ui> tie_vertices;
    std::vector<ui> temp;

    for (ui i = 1; i < query_vertices_num; ++i) {
        // Select the vertices with the maximum number of backward neighbors.
        ui max_bn = 0;
        for (ui u = 0; u < query_vertices_num; ++u) {
            if (!visited[u]) {
                // Compute the number of backward neighbors of u.
                ui cur_bn = 0;
                for (ui j = 0; j < i; ++j) {
                    ui uu = order[j];
                    if (query_graph->checkEdgeExistence(u, uu)) {
                        cur_bn += 1;
                    }
                }

                // Update the vertices under consideration.
                if (cur_bn > max_bn) {
                    max_bn = cur_bn;
                    tie_vertices.clear();
                    tie_vertices.push_back(u);
                } else if (cur_bn == max_bn) {
                    tie_vertices.push_back(u);
                }
            }
        }

        if (tie_vertices.size() != 1) {
            temp.swap(tie_vertices);
            tie_vertices.clear();

            ui count = 0;
            std::vector<ui> u_fn;
            for (auto u : temp) {
                // Compute the number of vertices in the matching order that has at least one vertex not in the matching order && connected with u.

                // Get the neighbors of u that are not in the matching order.
                ui un_count;
                const ui* un = query_graph->getVertexNeighbors(u, un_count);
                for (ui j = 0; j < un_count; ++j) {
                    if (!visited[un[j]]) {
                        u_fn.push_back(un[j]);
                    }
                }

                // Compute the valid number of vertices.
                ui cur_count = 0;
                for (ui j = 0; j < i; ++j) {
                    ui uu = order[j];
                    ui uun_count;
                    const ui* uun = query_graph->getVertexNeighbors(uu, uun_count);
                    ui common_neighbor_count = 0;
                    ComputeSetIntersection::ComputeCandidates(uun, uun_count, u_fn.data(), (ui)u_fn.size(), common_neighbor_count);
                    if (common_neighbor_count > 0) {
                        cur_count += 1;
                    }
                }

                u_fn.clear();

                // Update the vertices under consideration.
                if (cur_count > count) {
                    count = cur_count;
                    tie_vertices.clear();
                    tie_vertices.push_back(u);
                }
                else if (cur_count == count){
                    tie_vertices.push_back(u);
                }
            }
        }

        if (tie_vertices.size() != 1) {
            temp.swap(tie_vertices);
            tie_vertices.clear();

            ui count = 0;
            std::vector<ui> u_fn;
            for (auto u : temp) {
                // Compute the number of vertices not in the matching order && not the neighbor of vertices in the matching order, but is connected with u.

                // Get the neighbors of u that are not in the matching order.
                ui un_count;
                const ui* un = query_graph->getVertexNeighbors(u, un_count);
                for (ui j = 0; j < un_count; ++j) {
                    if (!visited[un[j]]) {
                        u_fn.push_back(un[j]);
                    }
                }

                // Compute the valid number of vertices.
                ui cur_count = 0;
                for (auto uu : u_fn) {
                    bool valid = true;

                    for (ui j = 0; j < i; ++j) {
                        if (query_graph->checkEdgeExistence(uu, order[j])) {
                            valid = false;
                            break;
                        }
                    }

                    if (valid) {
                        cur_count += 1;
                    }
                }

                u_fn.clear();

                // Update the vertices under consideration.
                if (cur_count > count) {
                    count = cur_count;
                    tie_vertices.clear();
                    tie_vertices.push_back(u);
                }
                else if (cur_count == count){
                    tie_vertices.push_back(u);
                }
            }
        }

        order[i] = tie_vertices[0];

        visited[order[i]] = true;
        for (ui j = 0; j < i; ++j) {
            if (query_graph->checkEdgeExistence(order[i], order[j])) {
                pivot[i] = order[j];
                break;
            }
        }

        tie_vertices.clear();
        temp.clear();
    }
}

void
GenerateQueryPlan::generateVF2PPQueryPlan(const Graph *data_graph, const Graph *query_graph, ui *&order, ui *&pivot) {
    ui query_vertices_num = query_graph->getVerticesCount();
    order = new ui[query_vertices_num];
    pivot = new ui[query_vertices_num];

    ui property_count = 0;
    std::vector<std::vector<ui>> properties(query_vertices_num);
    std::vector<bool> order_type(query_vertices_num, true);     // True: Ascending, False: Descending.
    std::vector<ui> vertices;
    std::vector<bool> in_matching_order(query_vertices_num, false);

    for (ui i = 0; i < query_vertices_num; ++i) {
        properties[i].resize(3);
    }

    // Select the root vertex with the rarest node labels and the largest degree.
    property_count = 2;
    order_type[0] = true;
    order_type[1] = false;

    for (ui u = 0; u < query_vertices_num; ++u) {
        vertices.push_back(u);
        properties[u][0] = data_graph->getLabelsFrequency(query_graph->getVertexLabel(u));
        properties[u][1] = query_graph->getVertexDegree(u);
    }

    auto order_lambda = [&properties, &order_type, property_count](ui l, ui r) -> bool {
        for (ui x = 0; x < property_count; ++x) {
            if (properties[l][x] == properties[r][x])
                continue;

            if (order_type[0]) {
                return properties[l][x] < properties[r][x];
            }
            else {
                return properties[l][x] > properties[r][x];
            }
        }

        return l < r;
    };

    std::stable_sort(vertices.begin(), vertices.end(), order_lambda);
    order[0] = vertices[0];
    in_matching_order[order[0]] = true;

    vertices.clear();
    TreeNode* tree;
    ui* bfs_order;
    GraphOperations::bfsTraversal(query_graph, order[0], tree, bfs_order);


    property_count = 3;
    order_type[0] = false;
    order_type[1] = false;
    order_type[2] = true;

    ui level = 1;
    ui count = 1;
    ui num_vertices_in_matching_order = 1;
    while (num_vertices_in_matching_order < query_vertices_num) {
        // Get the vertices in current level.
        while (count < query_vertices_num) {
            ui u = bfs_order[count];
            if (tree[u].level_ == level) {
                vertices.push_back(u);
                count += 1;
            }
            else {
                level += 1;
                break;
            }
        }

        // Process a level in the BFS tree.
        while(!vertices.empty()) {
            // Set property.
            for (auto u : vertices) {
                ui un_count;
                const ui* un = query_graph->getVertexNeighbors(u, un_count);

                ui bn_count = 0;
                for (ui i = 0; i < un_count; ++i) {
                    ui uu = un[i];
                    if (in_matching_order[uu]) {
                        bn_count += 1;
                    }
                }

                properties[u][0] = bn_count;
                properties[u][1] = query_graph->getVertexDegree(u);
                properties[u][2] = data_graph->getLabelsFrequency(query_graph->getVertexLabel(u));
            }


            std::sort(vertices.begin(), vertices.end(), order_lambda);
            pivot[num_vertices_in_matching_order] = tree[vertices[0]].parent_;
            order[num_vertices_in_matching_order++] = vertices[0];

            in_matching_order[vertices[0]] = true;

            vertices.erase(vertices.begin());
        }
    }


    delete[] tree;
    delete[] bfs_order;
}

void GenerateQueryPlan::generateOrderSpectrum(const Graph *query_graph, std::vector<std::vector<ui>> &spectrum,
                                              ui num_spectrum_limit) {
    ui query_vertices_num = query_graph->getVerticesCount();

    std::vector<bool> visited(query_vertices_num, false);
    std::vector<ui> matching_order;

    // Get the core vertices and the non-core vertices.
    std::vector<ui> core_vertices;
    std::vector<ui> noncore_vertices;
    for (ui u = 0; u < query_vertices_num; ++u) {
        if (query_graph->getCoreValue(u) >= 2) {
            core_vertices.push_back(u);
        }
        else {
            noncore_vertices.push_back(u);
        }
    }

    // Sort the vertices by the core value and the degree value.
    auto compare = [query_graph](ui l, ui r) -> bool {
        if (query_graph->getCoreValue(l) != query_graph->getCoreValue(r)) {
            return query_graph->getCoreValue(l) > query_graph->getCoreValue(r);
        }

        if (query_graph->getVertexDegree(l) != query_graph->getVertexDegree(r)) {
            return query_graph->getVertexDegree(l) > query_graph->getVertexDegree(r);
        }

        return l < r;
    };

    std::sort(core_vertices.begin(), core_vertices.end(), compare);
    std::sort(noncore_vertices.begin(), noncore_vertices.end(), compare);
    // Permutate the vertices. Keep the matching order is connected.

    std::vector<ui> cur_index(query_vertices_num, 0);

    std::vector<ui>& candidates = core_vertices.empty() ? noncore_vertices : core_vertices;

    for (ui i = 0; i < candidates.size(); ++i) {
        ui u = candidates[i];

        matching_order.push_back(u);
        visited[u] = true;

        int cur_depth = 1;
        cur_index[cur_depth] = 0;

        while (true) {
            while (cur_index[cur_depth] < candidates.size()) {
                u = candidates[cur_index[cur_depth]++];

                // If u has been in the matching order, then skip it.
                if (visited[u])
                    continue;

                // Determine whether u is connected to the vertices in the matching order.
                bool valid = false;

                for (auto v : matching_order) {
                    if (query_graph->checkEdgeExistence(u, v)) {
                        valid = true;
                        break;
                    }
                }

                if (!valid)
                    continue;

                // Update.
                matching_order.push_back(u);
                visited[u] = true;

                // All the core vertices are in the matching order.
                if (matching_order.size() == candidates.size()) {
                    spectrum.emplace_back(matching_order);

                    visited[matching_order.back()] = false;
                    matching_order.pop_back();

                    if (spectrum.size() >= num_spectrum_limit) {
                        goto EXIT;
                    }
                    else {
                        break;
                    }
                }
                else {
                    cur_depth += 1;
                    cur_index[cur_depth] = 0;
                }
            }

            visited[matching_order.back()] = false;
            matching_order.pop_back();
            cur_depth -= 1;

            if (cur_depth == 0)
                break;
        }
    }
    EXIT:
    // Order the non-core vertices. All matching orders in the spectrum have the same order of non-core vertices.
    if (!core_vertices.empty()) {
        matching_order.clear();

        for (auto u : core_vertices) {
            visited[u] = true;
        }
        while (matching_order.size() < noncore_vertices.size()) {
            for (auto u : noncore_vertices) {
                if (visited[u])
                    continue;

                bool valid = false;
                for (ui v = 0; v < query_vertices_num; ++v) {
                    if (visited[v] && query_graph->checkEdgeExistence(u, v)) {
                        valid = true;
                        break;
                    }
                }

                if (valid) {
                    matching_order.push_back(u);
                    visited[u] = true;
                    break;
                }
            }
        }

        for (auto &order : spectrum) {
            order.insert(order.end(), matching_order.begin(), matching_order.end());
        }
    }

    // Check the correctness of each matching order.
    for (auto& order : spectrum) {
        checkQueryPlanCorrectness(query_graph, order.data());
    }
}

void GenerateQueryPlan::checkQueryPlanCorrectness(const Graph *query_graph, ui *order) {
    ui query_vertices_num = query_graph->getVerticesCount();
    std::vector<bool> visited_vertices(query_vertices_num, false);
    // Check whether each query vertex is in the order.
    for (ui i = 0; i < query_vertices_num; ++i) {
        VertexID vertex = order[i];
        assert(vertex < query_vertices_num && vertex >= 0);

        visited_vertices[vertex] = true;
    }

    for (ui i = 0; i < query_vertices_num; ++i) {
        VertexID vertex = i;
        assert(visited_vertices[vertex]);
    }

    // Check whether the order is connected.

    std::fill(visited_vertices.begin(), visited_vertices.end(), false);
    visited_vertices[order[0]] = true;
    for (ui i = 1; i < query_vertices_num; ++i) {
        VertexID u = order[i];

        bool valid = false;
        for (ui j = 0; j < i; ++j) {
            VertexID v = order[j];
            if (query_graph->checkEdgeExistence(u, v)) {
                valid = true;
                break;
            }
        }

        assert(valid);
        visited_vertices[u] = true;
    }
}
