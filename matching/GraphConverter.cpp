//
// Created by ssunah on 10/29/19.
//

#include "graph/graph.h"

int main(int argc, char** argv) {
    std::string input_src_file_path(argv[1]);
    std::string output_dst_file_path(argv[2]);

    Graph graph(false);
    graph.loadGraphFromFile(input_src_file_path);

    std::string output_dst_degree_file_path = output_dst_file_path + "_deg.bin";
    std::string output_dst_adj_file_path = output_dst_file_path + "_adj.bin";
    std::string output_dst_label_file_path = output_dst_file_path + "_label.bin";

    graph.storeComparessedGraph(output_dst_degree_file_path, output_dst_adj_file_path, output_dst_label_file_path);
    return 0;
}

