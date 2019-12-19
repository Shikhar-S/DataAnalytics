import numpy as np
import pandas as pd
#read data
df=pd.read_csv('./../data/Circle.csv')
#get np vecs
X=df['x'].to_numpy()
Y=df['y'].to_numpy()
points=np.vstack((X, Y)).T
#returns edge weight between point i and j
def get_edge_weight(i,j):
    dis=np.sum((points[i]-points[j])**2)
    sigma_sq=0.005 #based on experimentation
    return np.exp(-(dis/sigma_sq))

A=np.zeros(shape=(points.shape[0],points.shape[0]))
#initialise adjacencey matrix
for i in range(points.shape[0]):
    for j in range(points.shape[0]):
        A[i][j]=get_edge_weight(i,j)

#compute laplacian
D = np.diag(np.squeeze(np.array(A.sum(axis=1))))
L = D-A
#compute eigenvecs and eigenvalues
vals,vec=np.linalg.eig(L)
# sort
vec = vec[:,np.argsort(vals)]
vals = vals[np.argsort(vals)]
#get first two eigenvecs for kmeans clustering
data=vec[:,:2]
#perform k means
from sklearn.cluster import KMeans
kmeans = KMeans(n_clusters=2).fit(data)
circle_labels=kmeans.labels_
#plot with color
import matplotlib.pyplot as plt
X=points[circle_labels==1]
Y=points[circle_labels==0]
plt.scatter(X[:,0],X[:,1],color='r')
plt.scatter(Y[:,0],Y[:,1],color='b')
plt.xlabel('x')
plt.ylabel('y')
plt.show()