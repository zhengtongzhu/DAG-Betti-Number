# DAG Betti Number Computation
Maximal betti number computation for DAG (directed acyclic graph). The main function and the comparisons with the general algorithm are provided.<br />
For the comparison, I generate random graphs with a specified structure, the edges are randomly sampled from a base graph (fully connected) without replacement. The density represents the proportion of edges selected from the base graph.<br />
For each density value ranging from 0.1 to 0.9 (increased by 0.1 each), 200 realizations are done and the total time spent is recorded.<br />
Graph 1: 2 layers with 10 nodes in each layer.<br />
Graph 2: 3 layers with 10 nodes in each layer.<br />
Graph 3: 4 layers with 4 nodes in the first layer and 10 in other layers.
