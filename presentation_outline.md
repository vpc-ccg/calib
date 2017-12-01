# Barcoded Read Error Correction Using Localing Sensitive Hashing and MinHashing
## Problem
- Illumina sequencing errors can be resolved using barcode tagging
- Error correction requires correct clustering of tagged barcodes
- But barcode tags are subject to sequencing errors themselves
### Formulation
***Cluster barcoded reads using their barcode and sequence similarity, and then correct each cluster's reads using the consensus of each cluster***

## Approach
### Preliminaries
#### Barcode similarity
Two barcodes are considered similar if they have a ***Hamming distance*** of less than a user defined threshold, *error_tolerance*. 
#### Minimizers
A minimizer of a given sequence is lexigraphicallly the smallest k-mer that exists in a the sequence, for a given k value. In our pipeline we pick k to be 8. We also split each mate's sequence (after trimming the barcode) of a read into 3 equal parts, and extract a minimizer from each part. Thus, each read will have a total of 6 minimizers, 3 from the left mate and 3 from the right mate.
### Graph formulation
We define a graph, $G$, that is composed of nodes, $V$. Each node is a set of reads that have identical barcodes and each of their 6 minimizers identical. An edge exists between any pair of nodes iff:

- They have simlar barcodes, and
- 1/3 of their left minimizers match, and 1/3 of their right minimizers match.

#### Read clusters
The final clusters of reads are defined to be the connected components of the graph.

### Error correction 
Error correction on the bases is done column by column by weighted vote. Weight are determined by quality score of a base. The voting winner quality score is calculated using likelihood of the error given repeated observations and their own quailty scores.

## Results
### Efficient graph building
To effiently find neighbors in the graph, we use locality sensitive hashing, allowing us to find all nodes neighbors in linear time (given that the barcodes are small). After that, edges that do not satisfy the minimizers condition are removed.
### Simulated dataset results
#### Simulated dataset
- 10K barcodes
- 1M molecules
- 60 PCR duplicated products
- Random subsampling from 60\*1M = 60M possible reads down to 6M reads
#### Accuracy
- 1 error
- 2 errors
#### Time and space
- 1 error
- 2 errors
## Things to do/consider
- Varying k-mer size
- Varying number of k-mers
- Varying the minimum number of matching k-mers for the condition
