## IncHL+
This is an implementation of the paper<br/>

Muhammad Farhan and Qing Wang, **[Efficient Maintenance of Distance Labelling for Incremental Updates in Large Dynamic Graphs.](https://arxiv.org/abs/2102.08529)**

## Sample data format
The format of the dataset text file is as follows: <br/>
Line 1 : |V| |E| <br/>
Line 2 : vertex_u deg_u v1 ... vn, where v1 to vn are neighbors of u. Note that the vertex id are from 0 to (V-1), where V is the number of vertices. There are no self loops in the graph, i.e., no edge from any vertex to itself. 

Please use the following commands to test IncHL_Plus.

## 1 - Compile source files

$ g++ -O3 -std=c++11 main.cpp -o run

## 2 - Construct Labelling
./run construct_labelling @1 @2 @3<br/>
@1: name of the dataset<br/>
@2: number of landmarks<br/>
@3: file to store labelling

Example:<br/>
./run construct_labelling graph.txt 20 graph_labelling

## 3 - Update Labelling
./run update_labelling @1 @2 @3 @4<br/>
@1: name of the dataset<br/>
@2: number of landmarks<br/>
@3: file to load the labelling from<br/>
@4: file containing updates<br/>

Example:<br/>
./run update_labelling graph.txt 20 graph_labelling updates.txt 0

## 4 - Perform distance queries
./run query-dis @1 @2 @3 @4 @5<br/>
@1: name of the dataset<br/>
@2: number of landmarks<br/>
@3: file to load the labelling from<br/>
@4: file containing query pairs<br/>
@5: file to write query results<br/>

Example:<br/>
./run query_labelling graph.txt 20 graph_labelling query_pairs.txt query_results.txt
