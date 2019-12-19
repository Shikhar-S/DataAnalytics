import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

#read graph
G=nx.read_gml('./../data/dolphins/dolphins.gml')
#get noralized laplacian
A = nx.adjacency_matrix(G).todense()
deg = A.sum(axis=1)
deg =np.squeeze(np.array(deg))
D_hf_inv=np.diag(deg**-0.5)
D=np.diag(deg)
L=D-A
L_norm=np.matmul(np.matmul(D_hf_inv,L),D_hf_inv)

#compute eigenvecs,eigenvalues
vals,vec=np.linalg.eig(L_norm)

# sort
vec = vec[:,np.argsort(vals)]
vals = vals[np.argsort(vals)]

# use Fiedler value to find best cut to separate data based on sign
community_a = vec[:,1] > 0

#obtain list of node names
V=list(G.nodes)
cluster_a=[]
cluster_b=[]
for i,x in enumerate(community_a):
    if x:
        cluster_a.append(V[i])
    else:
        cluster_b.append(V[i])


pos=nx.spring_layout(G) # positions for all nodes
# draw nodes
nx.draw_networkx_nodes(G,pos,
                       nodelist=cluster_a,
                       node_color='r',
                       node_size=100)
nx.draw_networkx_nodes(G,pos,
                       nodelist=cluster_b,
                       node_color='b',
                       node_size=100)
#draw edges
nx.draw_networkx_edges(G,pos,width=1.0,alpha=0.5)
#mark labels
labels={v:v for k,v in enumerate(V)}
nx.draw_networkx_labels(G,pos,labels,font_size=8)

plt.axis('off')
plt.show() 