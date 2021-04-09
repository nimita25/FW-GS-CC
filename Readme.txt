This code generates clusters for arbitrarily weighted incomplete graph with
edges that have two labels. If you want to generate clusters for a random
graph, skip to Step 2. Otherwise, if you have a graph, in the first
step, any input graph is translated to the appropriate graph using Jaccard
coefficient for each pair of nodes i and j. Please refer to the paper for
more .
Step 1: If you have the '.mat' file that has
(a) a sparse matrix Problem.W1 contains edges whose weights indicate
dissimilarity between nodes, and
(b) a sparse matrix Problem.W2 contains edges whose weights indicate
similarity between nodes, then store in the directory: 'CCGraphs' and go to
Step 2.
Else Store the input graph ('.mat' file) in the directory 'Datasets' and run
the following in the command prompt:
CreateCCGraph(filename)
For example, filename = 'email-Eu-core.mat'
This will create a new dataset in the directory 'CCGraphs' with the name
['CC-',filename]. Depending on the size of the graph, this step might take
a while.

Step 2: Execute the following:
(a) If you want to solve CC on a random graph
FWGS('R', [V,degree], e1, max_time, MemLog)
where V number of nodes and average degree of each node equal to 'degree'

(b) If you want to solve CC on the input graph
FWGS('S',filename, e1, max_time, MemLog)
where filename is the '.mat' that consists of graph information in the format
given in Step 1.

Set e1 = epsilon, the relative error to generate the solution to SDP.
Optionally, set
(i) MemLog = 1 if you want to track memory usage
(ii) max_time = max time to run the algorithm (in seconds)
