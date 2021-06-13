# BI290 Course Project Proposal

[toc]

## Requirements

- Your team should choose one project from the following proposals, of course you can also device one.
- The submitted materials should include:
  - a **written report** comprised of
    -  the **motivation** (describing the problem), **background**
    - the **solution** (describing the algorithmic details with pseudocode),
    - the **experiment** (conducting the algorithm on the simulated or real-world data),
    - the **conclusion** and
    - your **perspective** over the further improvement
  - a **package** which contains
    - the **source** codes written in C
    - the **Makefile** specifying the compilation procedure,
    - the **man**ual,
    - the **test** data, and 
    - the **demo**（运行过程）
  - a **presentation (e.g., ppt)** to present the project
- All the materials should be uploaded to [Github](https://github.com) or [Gitee](https://gitee.com/education) through Git through the cooperation of all the members in your team.



## Proposals



### 1. Biological sequence compression algorithms

The storage of the sequencing data has become an emergent issue in the field of high-throughput biological science. Most of the generated data are in the format of FASTQ, which is a raw-text format. Here we require you to design a set of lossless algorithms from following choices, and compare their compression rates for a set of data. 

#### Reference

- https://doi.org/10.1016/j.ygeno.2014.08.007
- https://doi.org/10.1155/2016/3528406

### 2. All-pairs suffix-prefix matching algorithms

Determining the suffix-prefix overlap between reads is a key step in assembling reads into contigs, which is the core algorithm in the field of high-throughput sequencing biotechnology.

#### Reference

- https://doi.org/10.5339/qfarc.2016.HBPP2066
- https://doi.org/10.1016/j.tcs.2017.07.013
- https://doi.org/10.1109/ICOSEC49089.2020.9215350
- https://link.springer.com/chapter/10.1007/978-3-642-13509-5_8
- https://drops.dagstuhl.de/opus/volltexte/2017/7617/pdf/LIPIcs-SEA-2017-14.pdf



### 3. All-pairs shortest path in a protein-protein interaction network

The shortest path between two proteins in a protein-protein interaction network may work as the indication of their biological relatedness. And a protein-protein interaction dataset can be modeled as a undirected network, either weighted or unweighted.  Various of algorithms are available for computing the shortest paths in a graph:

- **Floyd-Warshall Algorithm**
- **Dijkstra** algorithm
- **Bellman-Ford** algorithm

#### Reference

- Stark C, Breitkreutz B-J, Reguly T, Boucher L, Breitkreutz A, Tyers M. BioGRID: a general repository for interaction datasets. *Nucleic acids research.* 2006;34(suppl 1):D535–D539.
- Wilson N. Human Protein Reference Database. *Nature Reviews Molecular Cell Biology.* 2004;5(1):4–4.



### 4. Finding lowest common ancestor for given go terms in the GO graph

Finding the lowest common ancestor (LCA) is an important step in calculating the semantic similarity between two given genes, based on their annotations as well as the gene ontology structure, which is a complex graph to describe the relationships of the go terms.

#### Reference

- http://geneontology.org/
- https://www.nature.com/articles/s41598-018-33219-y
- https://doi.org/10.1155/2013/793091



### 5. Finding the exact matches (e.g., transposons) in a set of sequences

Finding the exact matches of short elements is an essential task in bioinformatics. One of the solution is to build the index for the set of sequences using either hashing, suffix tree, suffix array, or through Burrows-wheeler transform.

#### Reference

- https://link.springer.com/article/10.1007%2Fs007780200064

- https://dx.doi.org/10.1186%2Fs12859-018-2094-5
- https://doi.org/10.1371/journal.pone.0090581
- http://dx.doi.org/10.1093/bioinformatics/btaa911
- https://doi.org/10.1587/transinf.E94.D.2557



### 6. TrimReads

The high-throughput sequencing has brought us a vast amount of genetic sequences. However, some of the sequences are of low quality. Therefore, we need to develop a tool to trim the resulting sequences by removing the low-quality segments which should be located at either end of the reads.

- You can trim the sequence base-by-base according to the quality threshold.
- You can also trim the sequence window-by-window based on the average quality score of a given window size.
- These settings should be given by a set of options like: `--base-threshold 25 --window-threshold 20`



---





