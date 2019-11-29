# Graph similarity and matching

## Motivation
Biological networks are one of most complex systems. They exhibit relationships between between molecules, entities and represent structural properties. Studying these networks is crucial to understand life. Graph based data structures can efficiently capture the dynamics of such networks. A variety of static and dynamic graph based methods have been applied for such studies

## Problem Statement
Given a graph database and a query graph, perform approximate graph matching to find similarity scores between the two graphs and substructures from the data graph that match the query graph. Approximate graph matching algorithms are an efficient alternative to exact matching. A random walk score is calculated for every vertex which represents global structure. An associated reset probability is defined using beta Next, a beta signature is evaluated at multiple betas. From this, a beta similarity is obtained between 2 graphs. Using this score, a matching is performed. Such a matching algorithm can serve applications such as querying in PPI networks, Genome data graphs and for searching, comparing, matching and retrieving relevant information embedded in graphs

## Complexity
Since the algorithmic complexity for beta signature is O(beta*n*n) due to matrix-vector multiplication, a large biological graph containing above 10000 nodes will take a long time to be processed and evaluated. Therefore biological graphs within this size limit were used. Graphs for gene expression and protein-protein interaction in the form of egde list were obtained for worm (C.elegans), fly (D.melanogaster) and humans (H.sapiens) from WormNet (https://www.inetbio.org/wormnet/downloadnetwork.php)

/home/iiitd/Desktop/Screenshot from 2019-11-30 02-08-46.png
