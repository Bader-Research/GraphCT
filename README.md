# GraphCT

David A. Bader

GraphCT is a parallel toolkit for analyzing massive graphs on the massively multithreaded Cray XMT. Designed from the ground up to take advantage of the large shared memory and fine-grained synchronization of the Cray XMT, GraphCT provides an optimized library of functions including clustering coefficients, connected components, betweenness centrality, graph diameter, and distribution statistics. Functions are provided to read and write large data files in parallel using the shared memory. A proof-of-concept scripting and command-line interface is available for non-programmers. GraphCT uses a single common graph representation for all analysis kernels, and a straightforward API makes it easy to extend the library by implementing your own custom routines.


