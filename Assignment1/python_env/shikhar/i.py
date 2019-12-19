import pandas as pd
data=pd.read_csv('../data/datafile.csv')
#Remove wrong data
data=data[data['Error.In.Data']==0]
print("Original data shape",data.shape)

#getting the number of matches in remaining data
num_matches=len(set(data['Match']))
print("Number of matches:",num_matches)

#filter data by first innings and project only relevant columns
filtered_data=data[data['Innings']==1]
filtered_data=filtered_data[['Match','Innings','Over','Total.Runs','Total.Out','Innings.Total.Runs']]
print("Filtered Data shape",filtered_data.shape)

'''
u: overs to go
w: wickets in hand
R: runs scored from that point onwards
z0_L: parameters of model
'''
import numpy as np
#z0_L denotes 11 parameters : Z(1),...,Z(10) and L in this order
z0_L=np.random.rand(11)

#returns run scoring potential given u overs to go,w wickets in hand and 11 parameters
def Z(u,w,z0_L):
    z_indx=w-1#to correctly index onto the right parameter acc to wicket in hand
    return z0_L[z_indx]*(1-np.exp(-(z0_L[10]*u)/z0_L[z_indx])) # Implementation of model foumula Z(u,w;theta)=Z(W)(1-exp(-Lu/Z(w)))

#returns y-Z(u,w,theta) or the residual. Numpy automatically squares this to fit least square sum.
def residual(z0_L,u,w,R):
    return Z(u,w,z0_L)-R

u_train=50-filtered_data['Over'].values #overs to go
w_train=10-filtered_data['Total.Out'].values #wickets in hand
R_train=filtered_data['Innings.Total.Runs'].values-filtered_data['Total.Runs'].values #Runs scored from data point onwards

from scipy.optimize import least_squares
res_lsq = least_squares(residual, z0_L, args=(u_train,w_train,R_train)) #library fn to minimize least sq error sum

print("Initial params (Z1 to Z10 and L):",z0_L)
print("Trained params (Z1 to Z10 and L):",res_lsq.x) 
print("Cost at solution:",res_lsq.cost)

# %matplotlib inline
import matplotlib.pyplot as plt
fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 12
fig_size[1] = 9
plt.rcParams["figure.figsize"] = fig_size
for wickets_in_hand in range(1,11):
    x=np.arange(50)
    y=(Z(np.arange(50),wickets_in_hand,res_lsq.x)/Z(50,10,res_lsq.x))*100
    plt.plot(x,y,label=str(wickets_in_hand))
plt.xlabel('Overs In Hand')
plt.ylabel('Resource Percentage Remaining')
plt.legend()
plt.yticks(np.arange(0, 101, 10)) 
plt.grid()
plt.show()