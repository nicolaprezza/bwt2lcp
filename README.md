#BWT2LCP

### Overview

We offer two fast and memory-efficient internal-memory algorithms to induce the LCP array from the BWT of a DNA text o text collection. The algorithms simulate a navigation of the nodes and leaves of the generalized suffix tree of the collection using negligible space (O(log n) words) on top of the inpt (BWT) and output (LCP/DA). 

**bwt2lcp**: Induce the LCP array from the BWT of a text or a text collection over DNA alphabet {A,C,G,N,T,#}, where # is the string terminator (can be changed with option -t). If no 'N's are present, the algorithms switch automatically to more efficient data structures on alphabet {A,C,G,T,#}.

**merge_bwt**: Merges the (e)BWTs of two DNA string collections. Optionally, computes also the Document Array (DA) and LCP array of the merged collection. 

 Based on an extension (to BWTs of collections) of the suffix-tree navigation algorithm described in  the paper "*Linear time construction of compressed text indices in compact space*" by Djamal Belazzougui. The extension includes navigation of the suffix tree leaves and several optimizations to reduce the number of visited leaves.

### RAM usage

Let B be the number of Bytes needed to represent an LCP value. B can be set to 1, 2, 4, or 8 using option -l. When option 'N's are not present (alphabet = A,C,G,T,#) these are the RAM requirements:

- To merge two eBWTs, of size m and n: **(m+n)*0.625 Bytes**.
- To merge two eBWTs, of size m and n, and induce the LCP of the merged collection: **(m+n)*(B+0.625) Bytes**.
- To induce the LCP from the BWT, of size n, of a string collection: **n*(B+0.5) Bytes**.

If 'N's are present (i.e. alphabet is A,C,G,N,T,#), then the algorithms use **0.05*N** additional Bytes, where N is the total number of processed bases. 

The string terminator can be changed with option -t. In any case, the terminator is considered to be smaller than A,C,G,N,T (i.e. no matter the chosen character for the terminator).

### Running time

**Short answer**: O(n), where n is the total collection(s) size. 

**Long answer**: the library is based on a new packed DNA string (4 bits/base) with very efficient rank/access operations: a parallel rank for all letters A,C,G,T,# takes O(1) time and causes just one cache miss. Also a packed DNA string on alphabet A,C,G,T,N,# is provided. This structure uses 4.38 bits/base and still offers 1 cache miss per parallel rank on all alphabet's letters. 

There are at most n leaves and n internal suffix tree nodes. For each internal node, we perform (n_child + 1) ranks (for all letters), which overall causes at most 3n cache misses. Leaves navigation, on the other hand, generate at most 2n cache misses. To this, we need to add at most n cache misses to write LCP values. 

Overall, **bwt2lcp** causes therefore at most 6n cache misses. **merge_bwt**, on the other hand, generates at most 12n cache misses since the input BWTs are separate and we also need to compute the DA array. Also in practice, it can be seen that merge_bwt is twice as slow as induce_lcp.

**In practice**: Several optimizations have been introduced to reduce the number of visited leaves. These optimizations, together with the fact that we often write several adjacent LCP/DA entries and that several nodes fit within a cache line (because they have short intervals), bring down the above numbers to around **1.5n  cache-misses for bwt2lcp** and **3n cache misses for merge_bwt** (these numbers have been computed experimentally with a memory profiler).

### Publications

*Nicola Prezza and Giovanna Rosone, 2019. Space-Efficient Computation of the LCP Array from the Burrows-Wheeler Transform. Proceedings of the 30th Annual Symposium on Combinatorial Pattern Matching (CPM).*

### Funding

Supported by the project Italian MIUR-SIR CMACBioSeq ("Combinatorial methods for analysis and compression of biological sequences") grant n.~RBSI146R5L, PI: Giovanna Rosone

### Install

~~~~
git clone https://github.com/nicolaprezza/bwt2lcp
cd bwt2lcp
mkdir build
cd build
cmake ..
make
~~~~

### Run

Input BWTs must be ASCII files containing only A,C,G,N,T, and #. (terminator # can be changed with option -t)

- To induce the LCP from the eBWT of a string collection using B = 1,2,4,8 bytes for each LCP entry. Output is stored in out.lcp.
~~~~
bwt2lcp -i bwt -o out.lcp -l B
~~~~ 
- To merge two eBWTs of DNA collections:
~~~~
merge_bwt -1 bwt1 -2 bwt2 -o out
~~~~ 
- To merge two eBWTs and induce the LCP of the merged collection using B = 1,2,4,8 bytes for each LCP entry. Output is stored in out.lcp.
~~~~
merge_bwt -1 bwt1 -2 bwt2 -o out -l B
~~~~ 
- To merge two eBWTs and output the Document Array (ASCII file containing '0' and '1' characters) of the merged collections. Output is stored in out.lcp and out.da.
~~~~
merge_bwt -1 bwt1 -2 bwt2 -o out -d
~~~~ 






