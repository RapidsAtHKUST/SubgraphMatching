//
// Created by ssunah on 11/21/18.
//

#ifndef SUBGRAPHMATCHING_GENERATEFILTERINGPLAN_H
#define SUBGRAPHMATCHING_GENERATEFILTERINGPLAN_H


#include "graph/graph.h"
#include "configuration/types.h"

class GenerateFilteringPlan {
public:
    static void generateTSOFilterPlan(const Graph *data_graph, const Graph *query_graph, TreeNode *&tree,
                                           VertexID *&order);
    static void generateCFLFilterPlan(const Graph *data_graph, const Graph *query_graph, TreeNode *&tree,
                                      VertexID *&order, int &level_count, ui *&level_offset);
    static void generateDPisoFilterPlan(const Graph *data_graph, const Graph *query_graph, TreeNode *&tree,
                                        VertexID *&order);
    static void generateCECIFilterPlan(const Graph *data_graph, const Graph *query_graph, TreeNode *&tree,
                                       VertexID *&order);
private:
    static VertexID selectTSOFilterStartVertex(const Graph *data_graph, const Graph *query_graph);
    static VertexID selectCFLFilterStartVertex(const Graph *data_graph, const Graph *query_graph);
    static VertexID selectDPisoStartVertex(const Graph *data_graph, const Graph *query_graph);
    static VertexID selectCECIStartVertex(const Graph *data_graph, const Graph *query_graph);
};


#endif //SUBGRAPHMATCHING_GENERATEFILTERINGPLAN_H
