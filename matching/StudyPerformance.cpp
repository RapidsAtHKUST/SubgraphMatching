//
// Created by ssunah on 12/3/18.
//

#include <chrono>
#include <future>
#include <thread>
#include <fstream>

#include "matchingcommand.h"
#include "graph/graph.h"
#include "GenerateFilteringPlan.h"
#include "FilterVertices.h"
#include "BuildTable.h"
#include "GenerateQueryPlan.h"
#include "EvaluateQuery.h"

#define NANOSECTOSEC(elapsed_time) ((elapsed_time)/(double)1000000000)
#define BYTESTOMB(memory_cost) ((memory_cost)/(double)(1024 * 1024))

size_t enumerate(Graph* data_graph, Graph* query_graph, Edges*** edge_matrix, ui** candidates, ui* candidates_count,
                ui* matching_order, size_t output_limit) {
    static ui order_id = 0;

    order_id += 1;

    auto start = std::chrono::high_resolution_clock::now();
    size_t call_count = 0;
    size_t embedding_count = EvaluateQuery::LFTJ(data_graph, query_graph, edge_matrix, candidates, candidates_count,
                               matching_order, output_limit, call_count);

    auto end = std::chrono::high_resolution_clock::now();
    double enumeration_time_in_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
#ifdef SPECTRUM
    if (EvaluateQuery::exit_) {
        printf("Spectrum Order %u status: Timeout\n", order_id);
    }
    else {
        printf("Spectrum Order %u status: Complete\n", order_id);
    }
#endif
    printf("Spectrum Order %u Enumerate time (seconds): %.4lf\n", order_id, NANOSECTOSEC(enumeration_time_in_ns));
    printf("Spectrum Order %u #Embeddings: %zu\n", order_id, embedding_count);
    printf("Spectrum Order %u Call Count: %zu\n", order_id, call_count);
    printf("Spectrum Order %u Per Call Count Time (nanoseconds): %.4lf\n", order_id, enumeration_time_in_ns / (call_count == 0 ? 1 : call_count));

    return embedding_count;
}

void spectrum_analysis(Graph* data_graph, Graph* query_graph, Edges*** edge_matrix, ui** candidates, ui* candidates_count,
                       size_t output_limit, std::vector<std::vector<ui>>& spectrum, size_t time_limit_in_sec) {

    for (auto& order : spectrum) {
        std::cout << "----------------------------" << std::endl;
        ui* matching_order = order.data();
        GenerateQueryPlan::printSimplifiedQueryPlan(query_graph, matching_order);

        std::future<size_t> future = std::async(std::launch::async, [data_graph, query_graph, edge_matrix, candidates, candidates_count,
                                                                     matching_order, output_limit](){
            return enumerate(data_graph, query_graph, edge_matrix, candidates, candidates_count, matching_order, output_limit);
        });

        std::cout << "execute...\n";
        std::future_status status;
        do {
            status = future.wait_for(std::chrono::seconds(time_limit_in_sec));
            if (status == std::future_status::deferred) {
                std::cout << "Deferred\n";
                exit(-1);
            } else if (status == std::future_status::timeout) {
#ifdef SPECTRUM
                EvaluateQuery::exit_ = true;
#endif
            }
        } while (status != std::future_status::ready);
    }
}

int main(int argc, char** argv) {
    MatchingCommand command(argc, argv);
    std::string input_query_graph_file = command.getQueryGraphFilePath();
    std::string input_data_graph_file = command.getDataGraphFilePath();
    std::string input_filter_type = command.getFilterType();
    std::string input_order_type = command.getOrderType();
    std::string input_engine_type = command.getEngineType();
    std::string input_max_embedding_num = command.getMaximumEmbeddingNum();
    std::string input_time_limit = command.getTimeLimit();
    std::string input_order_num = command.getOrderNum();
    std::string input_distribution_file_path = command.getDistributionFilePath();
    std::string input_csr_file_path = command.getCSRFilePath();
    /**
     * Output the command line information.
     */
    std::cout << "Command Line:" << std::endl;
    std::cout << "\tData Graph CSR: " << input_csr_file_path << std::endl;
    std::cout << "\tData Graph: " << input_data_graph_file << std::endl;
    std::cout << "\tQuery Graph: " << input_query_graph_file << std::endl;
    std::cout << "\tFilter Type: " << input_filter_type << std::endl;
    std::cout << "\tOrder Type: " << input_order_type << std::endl;
    std::cout << "\tEngine Type: " << input_engine_type << std::endl;
    std::cout << "\tOutput Limit: " << input_max_embedding_num << std::endl;
    std::cout << "\tTime Limit (seconds): " << input_time_limit << std::endl;
    std::cout << "\tOrder Num: " << input_order_num << std::endl;
    std::cout << "\tDistribution File Path: " << input_distribution_file_path << std::endl;
    std::cout << "--------------------------------------------------------------------" << std::endl;

    /**
     * Load input graphs.
     */
    std::cout << "Load graphs..." << std::endl;

    auto start = std::chrono::high_resolution_clock::now();

    Graph* query_graph = new Graph(true);
    query_graph->loadGraphFromFile(input_query_graph_file);
    query_graph->buildCoreTable();

    Graph* data_graph = new Graph(true);

    if (input_csr_file_path.empty()) {
        data_graph->loadGraphFromFile(input_data_graph_file);
    }
    else {
        std::string degree_file_path = input_csr_file_path + "_deg.bin";
        std::string edge_file_path = input_csr_file_path + "_adj.bin";
        std::string label_file_path = input_csr_file_path + "_label.bin";
        data_graph->loadGraphFromFileCompressed(degree_file_path, edge_file_path, label_file_path);
    }

    auto end = std::chrono::high_resolution_clock::now();

    double load_graphs_time_in_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

    std::cout << "-----" << std::endl;
    std::cout << "Query Graph Meta Information" << std::endl;
    query_graph->printGraphMetaData();
    std::cout << "-----" << std::endl;
    data_graph->printGraphMetaData();

    std::cout << "--------------------------------------------------------------------" << std::endl;

    /**
     * Start queries.
     */

    std::cout << "Start queries..." << std::endl;
    std::cout << "-----" << std::endl;
    std::cout << "Filter candidates..." << std::endl;

    start = std::chrono::high_resolution_clock::now();

    ui** candidates = NULL;
    ui* candidates_count = NULL;
    ui* tso_order = NULL;
    TreeNode* tso_tree = NULL;
    ui* cfl_order = NULL;
    TreeNode* cfl_tree = NULL;
    ui* dpiso_order = NULL;
    TreeNode* dpiso_tree = NULL;
    TreeNode* ceci_tree = NULL;
    ui* ceci_order = NULL;
    std::vector<std::unordered_map<VertexID, std::vector<VertexID >>> TE_Candidates;
    std::vector<std::vector<std::unordered_map<VertexID, std::vector<VertexID>>>> NTE_Candidates;
    if (input_filter_type == "LDF") {
        FilterVertices::LDFFilter(data_graph, query_graph, candidates, candidates_count);
    } else if (input_filter_type == "NLF") {
        FilterVertices::NLFFilter(data_graph, query_graph, candidates, candidates_count);
    } else if (input_filter_type == "GQL") {
        FilterVertices::GQLFilter(data_graph, query_graph, candidates, candidates_count);
    } else if (input_filter_type == "TSO") {
        FilterVertices::TSOFilter(data_graph, query_graph, candidates, candidates_count, tso_order, tso_tree);
    } else if (input_filter_type == "CFL") {
        FilterVertices::CFLFilter(data_graph, query_graph, candidates, candidates_count, cfl_order, cfl_tree);
    } else if (input_filter_type == "DPiso") {
        FilterVertices::DPisoFilter(data_graph, query_graph, candidates, candidates_count, dpiso_order, dpiso_tree);
    } else if (input_filter_type == "CECI") {
        FilterVertices::CECIFilter(data_graph, query_graph, candidates, candidates_count, ceci_order, ceci_tree, TE_Candidates, NTE_Candidates);
    }  else {
        std::cout << "The specified filter type '" << input_filter_type << "' is not supported." << std::endl;
        exit(-1);
    }

    // Sort the candidates to support the set intersections
    if (input_filter_type != "CECI")
        FilterVertices::sortCandidates(candidates, candidates_count, query_graph->getVerticesCount());

    end = std::chrono::high_resolution_clock::now();
    double filter_vertices_time_in_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

    // Compute the candidates false positive ratio.
#ifdef OPTIMAL_CANDIDATES
    std::vector<ui> optimal_candidates_count;
    double avg_false_positive_ratio = FilterVertices::computeCandidatesFalsePositiveRatio(data_graph, query_graph, candidates,
                                                                                          candidates_count, optimal_candidates_count);
    FilterVertices::printCandidatesInfo(query_graph, candidates_count, optimal_candidates_count);
#endif
    std::cout << "-----" << std::endl;
    std::cout << "Build indices..." << std::endl;

    start = std::chrono::high_resolution_clock::now();

    Edges ***edge_matrix = NULL;
    if (input_filter_type != "CECI") {
        edge_matrix = new Edges **[query_graph->getVerticesCount()];
        for (ui i = 0; i < query_graph->getVerticesCount(); ++i) {
            edge_matrix[i] = new Edges *[query_graph->getVerticesCount()];
        }

        BuildTable::buildTables(data_graph, query_graph, candidates, candidates_count, edge_matrix);
    }

    end = std::chrono::high_resolution_clock::now();
    double build_table_time_in_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

    size_t memory_cost_in_bytes = 0;
    if (input_filter_type != "CECI") {
        memory_cost_in_bytes = BuildTable::computeMemoryCostInBytes(query_graph, candidates_count, edge_matrix);
        BuildTable::printTableCardinality(query_graph, edge_matrix);
    }
    else {
        memory_cost_in_bytes = BuildTable::computeMemoryCostInBytes(query_graph, candidates_count, ceci_order, ceci_tree,
                TE_Candidates, NTE_Candidates);
        BuildTable::printTableCardinality(query_graph, ceci_tree, ceci_order, TE_Candidates, NTE_Candidates);
    }
    std::cout << "-----" << std::endl;
    std::cout << "Generate a matching order..." << std::endl;

    start = std::chrono::high_resolution_clock::now();

    ui* matching_order = NULL;
    ui* pivots = NULL;
    ui** weight_array = NULL;

    size_t order_num = 0;
    sscanf(input_order_num.c_str(), "%zu", &order_num);

    std::vector<std::vector<ui>> spectrum;
    if (input_order_type == "QSI") {
        GenerateQueryPlan::generateQSIQueryPlan(data_graph, query_graph, edge_matrix, matching_order, pivots);
    } else if (input_order_type == "GQL") {
        GenerateQueryPlan::generateGQLQueryPlan(data_graph, query_graph, candidates_count, matching_order, pivots);
    } else if (input_order_type == "TSO") {
        if (tso_tree == NULL) {
            GenerateFilteringPlan::generateTSOFilterPlan(data_graph, query_graph, tso_tree, tso_order);
        }
        GenerateQueryPlan::generateTSOQueryPlan(query_graph, edge_matrix, matching_order, pivots, tso_tree, tso_order);
    } else if (input_order_type == "CFL") {
        if (cfl_tree == NULL) {
            int level_count;
            ui* level_offset;
            GenerateFilteringPlan::generateCFLFilterPlan(data_graph, query_graph, cfl_tree, cfl_order, level_count, level_offset);
            delete[] level_offset;
        }
        GenerateQueryPlan::generateCFLQueryPlan(data_graph, query_graph, edge_matrix, matching_order, pivots, cfl_tree, cfl_order, candidates_count);
    } else if (input_order_type == "DPiso") {
        if (dpiso_tree == NULL) {
            GenerateFilteringPlan::generateDPisoFilterPlan(data_graph, query_graph, dpiso_tree, dpiso_order);
        }

        GenerateQueryPlan::generateDSPisoQueryPlan(query_graph, edge_matrix, matching_order, pivots, dpiso_tree, dpiso_order,
                                                    candidates_count, weight_array);
    }
    else if (input_order_type == "CECI") {
        GenerateQueryPlan::generateCECIQueryPlan(query_graph, ceci_tree, ceci_order, matching_order, pivots);
    }
    else if (input_order_type == "RI") {
        GenerateQueryPlan::generateRIQueryPlan(data_graph, query_graph, matching_order, pivots);
    }
    else if (input_order_type == "VF2PP") {
        GenerateQueryPlan::generateVF2PPQueryPlan(data_graph, query_graph, matching_order, pivots);
    }
    else if (input_order_type == "Spectrum") {
        GenerateQueryPlan::generateOrderSpectrum(query_graph, spectrum, order_num);
    }
    else {
        std::cout << "The specified order type '" << input_order_type << "' is not supported." << std::endl;
    }

    end = std::chrono::high_resolution_clock::now();
    double generate_query_plan_time_in_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

    if (input_order_type != "Spectrum") {
        GenerateQueryPlan::checkQueryPlanCorrectness(query_graph, matching_order, pivots);
        GenerateQueryPlan::printSimplifiedQueryPlan(query_graph, matching_order);
    }
    else {
        std::cout << "Generate " << spectrum.size() << " matching orders." << std::endl;
    }

    std::cout << "-----" << std::endl;
    std::cout << "Enumerate..." << std::endl;
    size_t output_limit = 0;
    size_t embedding_count = 0;
    if (input_max_embedding_num == "MAX") {
        output_limit = std::numeric_limits<size_t>::max();
    }
    else {
        sscanf(input_max_embedding_num.c_str(), "%zu", &output_limit);
    }


#if ENABLE_QFLITER == 1
    EvaluateQuery::qfliter_bsr_graph_ = BuildTable::qfliter_bsr_graph_;
#endif

    size_t call_count = 0;
    size_t time_limit = 0;
    sscanf(input_time_limit.c_str(), "%zu", &time_limit);

    start = std::chrono::high_resolution_clock::now();

    if (input_engine_type == "EXPLORE") {
        embedding_count = EvaluateQuery::exploreGraph(data_graph, query_graph, edge_matrix, candidates,
                                                      candidates_count, matching_order, pivots, output_limit, call_count);
    } else if (input_engine_type == "LFTJ") {
        embedding_count = EvaluateQuery::LFTJ(data_graph, query_graph, edge_matrix, candidates, candidates_count,
                                              matching_order, output_limit, call_count);
    } else if (input_engine_type == "GQL") {
        embedding_count = EvaluateQuery::exploreGraphQLStyle(data_graph, query_graph, candidates, candidates_count,
                                                             matching_order, output_limit, call_count);
    } else if (input_engine_type == "QSI") {
        embedding_count = EvaluateQuery::exploreQuickSIStyle(data_graph, query_graph, candidates, candidates_count,
                                                             matching_order, pivots, output_limit, call_count);
    }
    else if (input_engine_type == "DPiso") {
        embedding_count = EvaluateQuery::exploreDPisoStyle(data_graph, query_graph, dpiso_tree,
                                                           edge_matrix, candidates, candidates_count,
                                                           weight_array, dpiso_order, output_limit,
                                                           call_count);
//        embedding_count = EvaluateQuery::exploreDPisoRecursiveStyle(data_graph, query_graph, dpiso_tree,
//                                                           edge_matrix, candidates, candidates_count,
//                                                           weight_array, dpiso_order, output_limit,
//                                                           call_count);
    }
    else if (input_engine_type == "Spectrum") {
        spectrum_analysis(data_graph, query_graph, edge_matrix, candidates, candidates_count, output_limit, spectrum, time_limit);
    }
    else if (input_engine_type == "CECI") {
        embedding_count = EvaluateQuery::exploreCECIStyle(data_graph, query_graph, ceci_tree, candidates, candidates_count, TE_Candidates,
                NTE_Candidates, ceci_order, output_limit, call_count);
    }
    else {
        std::cout << "The specified engine type '" << input_engine_type << "' is not supported." << std::endl;
        exit(-1);
    }

    end = std::chrono::high_resolution_clock::now();
    double enumeration_time_in_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

#ifdef DISTRIBUTION
    std::ofstream outfile (input_distribution_file_path , std::ofstream::binary);
    outfile.write((char*)EvaluateQuery::distribution_count_, sizeof(size_t) * data_graph->getVerticesCount());
    delete[] EvaluateQuery::distribution_count_;
#endif

    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cout << "Release memories..." << std::endl;
    /**
     * Release the allocated memories.
     */
    delete[] candidates_count;
    delete[] tso_order;
    delete[] tso_tree;
    delete[] cfl_order;
    delete[] cfl_tree;
    delete[] dpiso_order;
    delete[] dpiso_tree;
    delete[] ceci_order;
    delete[] ceci_tree;
    delete[] matching_order;
    delete[] pivots;
    for (ui i = 0; i < query_graph->getVerticesCount(); ++i) {
        delete[] candidates[i];
    }
    delete[] candidates;

    if (edge_matrix != NULL) {
        for (ui i = 0; i < query_graph->getVerticesCount(); ++i) {
            for (ui j = 0; j < query_graph->getVerticesCount(); ++j) {
                delete edge_matrix[i][j];
            }
            delete[] edge_matrix[i];
        }
        delete[] edge_matrix;
    }
    if (weight_array != NULL) {
        for (ui i = 0; i < query_graph->getVerticesCount(); ++i) {
            delete[] weight_array[i];
        }
        delete[] weight_array;
    }

    delete query_graph;
    delete data_graph;

    /**
     * End.
     */
    std::cout << "--------------------------------------------------------------------" << std::endl;
    double preprocessing_time_in_ns = filter_vertices_time_in_ns + build_table_time_in_ns + generate_query_plan_time_in_ns;
    double total_time_in_ns = preprocessing_time_in_ns + enumeration_time_in_ns;

    printf("Load graphs time (seconds): %.4lf\n", NANOSECTOSEC(load_graphs_time_in_ns));
    printf("Filter vertices time (seconds): %.4lf\n", NANOSECTOSEC(filter_vertices_time_in_ns));
    printf("Build table time (seconds): %.4lf\n", NANOSECTOSEC(build_table_time_in_ns));
    printf("Generate query plan time (seconds): %.4lf\n", NANOSECTOSEC(generate_query_plan_time_in_ns));
    printf("Enumerate time (seconds): %.4lf\n", NANOSECTOSEC(enumeration_time_in_ns));
    printf("Preprocessing time (seconds): %.4lf\n", NANOSECTOSEC(preprocessing_time_in_ns));
    printf("Total time (seconds): %.4lf\n", NANOSECTOSEC(total_time_in_ns));
    printf("Memory cost (MB): %.4lf\n", BYTESTOMB(memory_cost_in_bytes));
    printf("#Embeddings: %zu\n", embedding_count);
    printf("Call Count: %zu\n", call_count);
    printf("Per Call Count Time (nanoseconds): %.4lf\n", enumeration_time_in_ns / (call_count == 0 ? 1 : call_count));
    std::cout << "End." << std::endl;
    return 0;
}