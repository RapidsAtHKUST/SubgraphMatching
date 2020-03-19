#include "util.hpp"

void quit()
{
    system("pause");
    exit(0);
}

void align_malloc(void **memptr, size_t alignment, size_t size)
{
    int malloc_flag = posix_memalign(memptr, alignment, size);
    if (malloc_flag) {
        std::cerr << "posix_memalign: " << strerror(malloc_flag) << std::endl;
        quit();
    }
}

std::string extract_filename(const std::string full_filename)
{
    int pos = full_filename.find_last_of('.');
    return full_filename.substr(0, pos);
}

int arg_pos(char *str, int argc, char **argv)
{
  int a;
  for (a = 1; a < argc; a++) if (!strcmp(str, argv[a])) {
    if (a == argc - 1) {
      printf("Argument missing for %s\n", str);
      quit();
    }
    return a;
  }
  return -1;
}

EdgeVector load_graph(const std::string path)
{
    EdgeVector edge_vec;
    FILE *fp = fopen(path.c_str(), "r");
    if (fp == NULL) {
        std::cout << "fail to open " << path << std::endl;
        quit();
    }

    char line[512];
    while (fgets(line, 512, fp) != NULL) {
        if (line[0] == '#') continue;
        int u = 0, v = 0;
        const char *c = line;
        while (isdigit(*c))
            u = (u << 1) + (u << 3) + (*c++ - 48);
        c++;
        while (isdigit(*c))
            v = (v << 1) + (v << 3) + (*c++ - 48);
        edge_vec.push_back(std::make_pair(u, v));
    }    
    fclose(fp);
    
    return edge_vec;
}

std::vector<int> load_vertex_order(const std::string path)
{    
    FILE *fp = fopen(path.c_str(), "r");
    if (fp == NULL) {
        std::cout << "fail to open " << path << std::endl;
        quit();
    }

    EdgeVector id_pair;
    char line[512];
    while (fgets(line, 512, fp) != NULL) {
        if (line[0] == '#') continue;
        int u = 0, v = 0;
        const char *c = line;
        while (isdigit(*c))
            u = (u << 1) + (u << 3) + (*c++ - 48);
        c++;
        while (isdigit(*c))
            v = (v << 1) + (v << 3) + (*c++ - 48);
        id_pair.push_back(std::make_pair(u, v));
    }  
    fclose(fp);

    std::vector<int> order(id_pair.size());
    for (auto& p : id_pair)
        order[p.first] = p.second;

    return order;
}

void save_graph(const std::string path, const EdgeVector& edge_vec)
{
    FILE *fp = fopen(path.c_str(), "w");
    if (fp == NULL) {
        std::cout << "fail to create " << path << std::endl;
        quit();
    }

    for (auto& e : edge_vec) {
        fprintf(fp, "%d %d\n", e.first, e.second);
    }
    fclose(fp);
}

void save_newid(const std::string path, std::vector<int> org2newid)
{
    FILE *fp = fopen(path.c_str(), "w");
    if (fp == NULL) {
        std::cout << "fail to create " << path << std::endl;
        quit();
    }

    for (int i = 0; i < (int)org2newid.size(); ++i)
        fprintf(fp, "%d %d\n", i, org2newid[i]);
    fclose(fp);
}

bool edge_idpair_cmp(const Edge& a, const Edge& b)
{
    if(a.first == b.first) return a.second < b.second;
    else return a.first < b.first;
}