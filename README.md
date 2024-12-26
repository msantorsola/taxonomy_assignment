# clustering_bioimaging

The project concerned the quality controls of radioactive substances for PET scan with a fast decay and produced in situ at regional hospitals.

A bioinformatics pipeline was designed to perform the full-length 16S ribosomal RNA (or 16S rRNA) amplicon analysis from High-Throughput Sequencing (HTS).

Firstly, ssu-align software package was used for identifying, aligning and masking bacterial 16S rRNA sequences.
Next, the multi-threaded version of RAxML under the GTRCAT approximation model was used for phylogenetic analysis. A phylogenetic tree including 12,474 bacteria taxa from SILVA 128 ribosomal RNA (rRNA) database was used as reference.
Best phylogeny and the taxonomic assignment of query sequences were released in a tabular file format by using an in-house developed script, including the probability value of belonging to a certain taxon.


