//
// Created by ssunah on 11/20/18.
//

#ifndef SUBGRAPHMATCHING_BUILDTABLE_H
#define SUBGRAPHMATCHING_BUILDTABLE_H

#include "graph/graph.h"
#include "utility/QFliter.h"
#include <vector>
class BuildTable {
public:
    static void buildTables(const Graph* data_graph, const Graph* query_graph, ui** candidates, ui* candidates_count,
                            Edges*** edge_matrix);

    static void printTableCardinality(const Graph* query_graph, Edges*** edge_matrix);
    static void printTableCardinality(const Graph *query_graph, TreeNode *tree, ui *order,
                                         std::vector<std::unordered_map<VertexID, std::vector<VertexID >>> &TE_Candidates,
                                         std::vector<std::vector<std::unordered_map<VertexID, std::vector<VertexID>>>> &NTE_Candidates);
    static void printTableInfo(const Graph* query_graph, Edges*** edge_matrix);
    static void printTableInfo(VertexID begin_vertex, VertexID end_vertex, Edges*** edge_matrix);
    static size_t computeMemoryCostInBytes(const Graph* query_graph, ui* candidates_count, Edges*** edge_matrix);
    static size_t computeMemoryCostInBytes(const Graph *query_graph, ui *candidates_count, ui *order, TreeNode *tree,
                                               std::vector<std::unordered_map<VertexID, std::vector<VertexID >>> &TE_Candidates,
                                               std::vector<std::vector<std::unordered_map<VertexID, std::vector<VertexID>>>> &NTE_Candidates);

#if ENABLE_QFLITER == 1
    static BSRGraph*** qfliter_bsr_graph_;
#endif
};


#endif //SUBGRAPHMATCHING_BUILDTABLE_H
