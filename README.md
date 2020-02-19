# SubgraphMatching
## Introduction
We study the performance of eight representative in-memory
subgraph matching algorithms. Specifically, we put QuickSI,
GraphQL, CFL, CECI, DP-iso, RI and VF2++ in a common
framework to compare them on the following four aspects:
(1) method of filtering candidate vertices in the data graph;
(2) method of ordering query vertices; (3) method of enumer-
ating partial results; and (4) other optimization techniques.
Then, we compare the overall performance of these algo-
rithms with Glasgow, an algorithm based on the constraint
programming. Through experiments, we find that (1) the fil-
tering method of GraphQL is competitive to that of the latest
algorithms CFL, CECI and DP-iso in terms of pruning power;
(2) the ordering methods in GraphQL and RI are usually the
most effective; (3) the set intersection based local candidate
computation in CECI and DP-iso performs the best in the
enumeration; and (4) the failing sets pruning in DP-iso can
significantly improve the performance when queries become
large. Based on these new results, we recommend users to
adopt specific techniques depending on the data graph den-
sity and query size.

For the details, please refer to our SIGMOD'2020 paper
"In-Memory Subgraph Matching: an In-depth Study"
by [Dr. Shixuan Sun](https://github.com/shixuansun) and [Prof. Qiong Luo](http://www.cse.ust.hk/~luo/).
If you have any further questions, please feel free to contact us.

Please cite our paper, if you use our source code.

* "Shixuan Sun and Qiong Luo. In-Memory Subgraph Matching: an In-depth Study. SIGMOD 2020."



## Compile
Under the root directory of the project, execute the following commands to compile the source code.

```zsh
mkdir build
cd build
cmake ..
make
```

## Test
Execute the following commands to test the correctness of the binary file.

```zsh
cd test
python test.py ../build/matching/SubgraphMatching.out
```

## Execute
After compiling the source code, you can find the binary file 'SubgraphMatching.out' under the 'build/matching' directory. 
Execute the binary with the following command ./SubgraphMatching.out -d data_graphs -q query_graphs
-filter method_of_filtering_candidate_vertices -order method_of_ordering_query_vertices -engine method_of_enumerating_partial_results -num number_of_embeddings,
in which -d specifies the input of the data graphs and -q specifies the input of the query graphs.
The -filter parameter gives the filtering method, the -order specifies the ordering method, and the -engine
sets the enumeration method. The -num parameter sets the maximum number of embeddings that you would like to find.
If the number of embeddings enumerated reaches the limit or all the results have been found, then the program will terminate.
Set -num as 'MAX' to find all results.

Example (Use the filtering and ordering methods of GraphQL to generate the candidate vertex sets and the matching order respectively.
Enumerate results with the set-intersection based local candidate computation method):


```zsh
./SubgraphMatching.out -d ../../test/sample_dataset/test_case_1.graph -q ../../test/sample_dataset/query1_positive.graph -filter GQL -order GQL -engine LFTJ -num MAX
```

## Input
Both the input query graph and data graph are vertex-labeled.
Each graph starts with 't N M' where N is the number of vertices and M is the number of edges. A vertex and an edge are formatted
as 'v VertexID LabelId Degree' and 'e VertexId VertexId' respectively. Note that we require that the vertex
id is started from 0 and the range is [0,N - 1] where V is the vertex set. The following
is an input sample. You can also find sample data sets and query sets under the test folder.

Example:

```zsh
t 5 6
v 0 0 2
v 1 1 3
v 2 2 3
v 3 1 2
v 4 2 2
e 0 1
e 0 2
e 1 2
e 1 3
e 2 4
e 3 4
```

## Techniques Supported

The filtering methods that generate candidate vertex sets. Note that (1) every filtering method enables LDF by default; and (2)
We optimize the filtering method of TurboIso with the technique of CFL.

|Parameter of Command Line (-filter) | Description |
| :-----------------------------------: | :-------------: |
|LDF| the label and degree filter |
|NLF| the neighborhood label frequency filter |
|GQL| the filtering method of GraphQL |
|TSO| the filtering method of TurboIso |
|CFL| the filtering method of CFL|
|DPiso| the filtering method of DP-iso |
|CECI| the filtering method of CECI |

The ordering methods that generate matching order.

|Parameter of Command Line (-order) | Description |
| :-----------------------------------: | :-------------: |
|QSI| the ordering method of QuickSI |
|GQL| the ordering method of GraphQL |
|TSO| the ordering method of TurboIso |
|CFL| the ordering method of CFL |
|DPiso| the ordering method of DP-iso |
|CECI| the ordering method of CECI|
|RI| the ordering method of RI |
|VF2++| the ordering method of VF2++ |

The enumeration methods that find all results.

|Parameter of Command Line (-engine) | Description |
| :-----------------------------------: | :-------------: |
|QSI| the enumeration method of QuickSI and RI (Algorithm 2) |
|GQL| the enumeration method of GraphQL (Algorithm 3) |
|EXPLORE| the enumeration method of CFL (Algorithm 4) |
|DPiso| the enumeration method of DP-iso (Algorithm 5) |
|CECI| the enumeration method of CECI (Algorithm 5) |
|LFTJ| the set-intersection based local candidates computation |

We can execute a query by specifying the filtering method, the ordering method, and the enumeration method.
For example, we evaluate the query with the GraphQL algorithm in the following command.

```zsh
./SubgraphMatching.out -d ../../test/sample_dataset/test_case_1.graph -q ../../test/sample_dataset/query1_positive.graph -filter GQL -order GQL -engine GQL -num MAX
```

Apart from evaluating queries with a specific algorithm, we can also execute queries by integrating techniques from
different algorithms. For example, we generate the candidates with the label and degree filter, obtain the matching order
with the ordering method of CFL, and finally enumerate the results with the set-intersection based local candidates
computation method. Note that the source code cannot support the filtering/ordering/enumeration orders of CECI to integrate with
other algorithms.

```zsh
./SubgraphMatching.out -d ../../test/sample_dataset/test_case_1.graph -q ../../test/sample_dataset/query1_positive.graph -filter LDF -order CFL -engine LFTJ -num MAX
```

In addition to setting the filtering, ordering and enumeration methods, you can also configure the set intersection
algorithms and optimization techniques by defining macros in 'configuration/config.h'. Please see the comments
in 'configuration/config.h' for the detail. We briefly introduce these methods in the following.


| Macro | Description |
| :-----------------------------------: | :-------------: |
|HYBRID 0| a hybrid method handling the cardinality skew by integrating the merge-based method with the galloping-based method  |
|HYBRID 1| the merge-based set intersection |
|SI 0 | Accelerate the set intersection with AVX2 |
|SI 1 | Accelerate the set intersection with AVX512 |
|SI 2 | Scalar set intersection |
|ENABLE_QFLITER | the [QFilter](https://dl.acm.org/doi/10.1145/3183713.3196924) set intersection algorithm |
|ENABLE_FAILING_SET | the failing set pruning technique in DP-iso |

Furthermore, our source code can support the spectrum analysis on a query. Specifically, we evaluate a query by
generate a number of matching orders by performing the permutation of query vertices. For each order, we set a time limit for evaluating the query.
In the following command line, we generate 1000 matching orders and set the time limit for each order as 60 seconds.
Note that you need to first define SPECTRUM in 'configuration/config.h' file before executing the command.

```zsh
./SubgraphMatching.out -d ../../test/sample_dataset/test_case_1.graph -q ../../test/sample_dataset/query1_positive.graph -filter LDF -order Spectrum -engine Spectrum -order_num 100 -time_limit 60 -num 100000
```

## Recommendation
We recommend you to use GQL, CFL or DPiso as the filtering method and use
GQL or RI as the ordering method. We always recommend you to use LFTJ as the enumeration method.
For the large queries, we recommend you to enable the failing set pruning.
For the dense data graphs, we recommend you to enable QFilter. Otherwise, use the hybrid set intersection method with AVX2.




## Experiment Datasets
The real world datasets and the corresponding query sets used in our paper can be downloaded here.
As the synthetic datasets are large, we do not upload them. You can easily generate them by following the instruction in our paper.

