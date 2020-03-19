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

int main(int argc, char *argv[]) {
    string input_vertex_num(argv[1]);
    string input_label_num(argv[2]);
    string output_file_path(argv[3]);

    size_t vertex_num = std::stoul(input_vertex_num);
    int label_num = std::stoi(input_label_num);

    cout << "start assign..." << endl;
    std::vector<int> labels(vertex_num);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, label_num - 1);

    for (int i = 0; i < vertex_num; ++i) {
        labels[i] = dis(gen);
    }

    cout << "start write..." << endl;
    ofstream label_ofs(output_file_path, ios::binary);
    label_ofs.write(reinterpret_cast<const char *>(&labels.front()), labels.size() * 4);
    cout << "finish..." << endl;

    return 0;
}
