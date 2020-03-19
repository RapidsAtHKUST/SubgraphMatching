//
// Created by ssunah on 10/31/19.
//

//
// Created by ssunah on 10/31/19.
//

#include <cassert>
#include <algorithm>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <chrono>
#include <tuple>
#include <random>
using namespace std;

void write_csr(string& deg_output_file, string &adj_output_file, vector<pair<int, int>> &lines, unsigned int vertex_num) {
    auto edge_num = lines.size();
    vector<int> degree_arr(vertex_num, 0);
    vector<vector<int>> matrix(vertex_num);

    ofstream deg_ofs(deg_output_file, ios::binary);

    for (const auto &line : lines) {
        int src, dst;
        std::tie(src, dst) = line;
        degree_arr[src]++;
        degree_arr[dst]++;
        matrix[src].emplace_back(dst);
        matrix[dst].emplace_back(src);
    }

    cout << "begin write" << endl;
    cout << "Input |V| = " << vertex_num << ", |E| = " << edge_num << endl;
    int int_size = sizeof(int);
    deg_ofs.write(reinterpret_cast<const char *>(&int_size), 4);
    deg_ofs.write(reinterpret_cast<const char *>(&vertex_num), 4);
    deg_ofs.write(reinterpret_cast<const char *>(&edge_num), 4);
    deg_ofs.write(reinterpret_cast<const char *>(&degree_arr.front()), degree_arr.size() * 4);

    cout << "finish degree write..." << endl;
    ofstream adj_ofs(adj_output_file, ios::binary);
    for (auto &adj_arr: matrix) {
        adj_ofs.write(reinterpret_cast<const char *>(&adj_arr.front()), adj_arr.size() * 4);
    }
    cout << "finish edge write..." << endl;
}

void read_csr(string& degree_bin_path, string& adj_bin_path, vector<unsigned int>& offsets, vector<unsigned int>& neighbors) {
    std::ifstream deg_file(degree_bin_path, std::ios::binary);

    if (deg_file.is_open()) {
        std::cout << "Open degree file " << degree_bin_path << " successfully." << std::endl;
    }
    else {
        std::cerr << "Cannot open degree file " << degree_bin_path << " ." << std::endl;
        exit(-1);
    }

    int int_size;
    unsigned int vertex_count;
    unsigned int edge_count;
    deg_file.read(reinterpret_cast<char *>(&int_size), 4);
    deg_file.read(reinterpret_cast<char *>(&vertex_count), 4);
    deg_file.read(reinterpret_cast<char *>(&edge_count), 4);

    cout << "Input |V| = " << vertex_count << ", |E| = " << edge_count << endl;
    offsets.resize(vertex_count + 1);
    vector<unsigned int> degrees(vertex_count);

    deg_file.read(reinterpret_cast<char *>(&degrees.front()), sizeof(int) * vertex_count);

    deg_file.close();
    deg_file.clear();


    std::ifstream adj_file(adj_bin_path, std::ios::binary);

    if (adj_file.is_open()) {
        std::cout << "Open edge file " << adj_bin_path << " successfully." << std::endl;
    }
    else {
        std::cerr << "Cannot open edge file " << adj_bin_path << " ." << std::endl;
        exit(-1);
    }

    size_t neighbors_count = (size_t)edge_count * 2;
    neighbors.resize(neighbors_count);

    offsets[0] = 0;
    for (unsigned int i = 1; i <= vertex_count; ++i) {
        offsets[i] = offsets[i - 1] + degrees[i - 1];
    }

    for (unsigned int i = 0; i < vertex_count; ++i) {
        if (degrees[i] > 0) {
            adj_file.read(reinterpret_cast<char *>(&neighbors.data()[offsets[i]]), degrees[i] * sizeof(int));
        }
    }

    adj_file.close();
    adj_file.clear();
}

int main(int argc, char *argv[]) {
    string input_edge_prob(argv[1]);
    string input_degree_bin_path(argv[2]);
    string input_adj_bin_path(argv[3]);
    string output_degree_bin_path(argv[4]);
    string output_adj_bin_path(argv[5]);

    double prob_threshold = std::stod(input_edge_prob);

    cout << "start read..." << endl;

    vector<unsigned int> offset;
    vector<unsigned int> adj_bin;
    read_csr(input_degree_bin_path, input_adj_bin_path, offset, adj_bin);

    vector<pair<int, int>> edge_list;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0, 1);

    for (unsigned int i = 0; i < offset.size() - 1; ++i) {
        int src = i;
        for (unsigned int j = offset[i]; j < offset[i + 1]; ++j) {
            int dst = adj_bin[j];

            if (src < dst) {
                double prob = dis(gen);
                if (prob < prob_threshold)
                    edge_list.emplace_back(make_pair(src, dst));
            }
        }
    }


    cout << "start write..." << endl;
    write_csr(output_degree_bin_path, output_adj_bin_path, edge_list, offset.size() - 1);
    cout << "finish..." << endl;

    return 0;
}
