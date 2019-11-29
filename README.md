# Graph similarity and matching

## Motivation
Biological networks are one of most complex systems. They exhibit relationships between between molecules, entities and represent structural properties. Studying these networks is crucial to understand life systems and processes. Graph based data structures can efficiently capture the dynamics of such networks. A variety of static and dynamic graph based methods have been applied for such studies. Our project is inspired by the paper "Top-k Similar Graph Matching Using TraM in Biological Networks" https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=6226354

## Problem Statement
Given a graph database and a query graph, perform approximate graph matching to find similarity scores between the two graphs and substructures from the data graph that match the query graph. Approximate graph matching algorithms are an efficient alternative to exact matching. Such a matching algorithm can serve applications in querying in PPI networks, genome data graphs and for searching, comparing, matching and retrieving relevant information embedded in graphs

## Complexity
Since the algorithmic complexity for beta signature is O(beta*n*n) due to matrix-vector multiplication, a large biological graph containing above 10000 nodes will take a long time to be processed and evaluated. Therefore biological graphs within this size limit were used. Graphs for gene expression and protein-protein interaction in the form of egde list were obtained for worm (C.elegans), fly (D.melanogaster) and humans (H.sapiens) from WormNet (https://www.inetbio.org/wormnet/downloadnetwork.php)

![Data Source](/datasrc.png)

## Algorithm Design
The algorithm mentioned in the paper is followed. A distance value is obtained based on the conceptual likeness and similarity wrt topology. A random walk score is calculated for every vertex which represents global structure. An associated reset probability is defined using beta. pt+1 = (1 - β) * pt + β * p0, where pt represents a vector whose ith element is the probability of being at node vi at time step t. Next, a beta signature is evaluated at multiple betas. From this, a beta similarity is obtained between 2 graphs. This score is used to perform graph matching. Candidate subgraphs are generated and pruned. These are matched with the query graph and k related graphs are obtained.
 

## Findings

## Conclusions
These experiments indicate that the similarity score clearly explain the similarity in the species. 

## Code running instructions



