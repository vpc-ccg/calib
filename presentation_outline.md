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
A minimizer of a given sequence is lexigraphicallly the smallest k-mer that exists in a the sequence, for a given k value. In our pipeline we pick k to be 8. We also split a read's sequence (the part of the read after trimming the barcode) into 3 equal parts, and extract a minimizer from each part. 

We define a graph, $G$, that is composed of nodes, $V$. Each connected componenet 
