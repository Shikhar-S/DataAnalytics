import networkx as nx
import community as comm
import matplotlib.pyplot as plt

G=nx.read_gml('../data/dolphins/dolphins.gml')
partition=comm.best_partition(G)
communities=list(set(partition.values()))

pos=nx.spring_layout(G)
color=['b','g','r','c','m','y','k']


for community in communities:
    list_nodes=[nodes for nodes in partition.keys() if partition[nodes]==community]
    nx.draw_networkx_nodes(G,pos,list_nodes,node_size=100,node_color=color[community])
    
nx.draw_networkx_edges(G,pos,alpha=0.5)

V=list(G.nodes)
labels={v:v for k,v in enumerate(V)}

nx.draw_networkx_labels(G,pos,labels,font_size=8)

plt.show()