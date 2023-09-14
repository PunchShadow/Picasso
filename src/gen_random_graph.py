import networkx as nx
import scipy as sp
import scipy.io

n = 2 ** 10 
m = 2 ** 18
G = nx.dense_gnm_random_graph(n,m)
A = nx.to_scipy_sparse_array(G)

fileName = "/scratch/gilbreth/sferdou/dense-random/dense_gnm1.mtx"
sp.io.mmwrite(fileName,A)
