import pandas as pd
df=pd.read_csv('./../data/01_data_mars_opposition.csv')
from math import sin, cos
import numpy as np
#this method takes in np arrays of degree,minute,seconds and zodiac index and returns np array of radians
def get_as_radian(degree,minute,second,zodiac_index):
    base_angle=zodiac_index*30
    deg_angle=base_angle+degree+(minute/60.0)+(second/3600.0)
    return np.radians(deg_angle)

alpha=get_as_radian(df['Degree'],df['Minute'],df['Second'],df['ZodiacIndex'])
beta=get_as_radian(df['DegreeMean'],df['MinuteMean'],df['SecondMean'],df['ZodiacIndexAverageSun'])


#radius is the distance from center, alpha= mars_sun angle, beta= mars_avsun angle, x, y are params as in slides
def radius(alpha,beta,x,y):
    #A,B,C are angles of triangle formed by joining vertices sun, av sun and mars
    A=beta-y
    B=np.pi-(alpha-y)
    C=np.pi-A-B
    a=(np.sin(A)/np.sin(C))*(1+x) #used sine rule to get the side opposite to average sun vertex of triangle
    return np.sqrt(1+a**2-2*a*np.cos(B)) #used cosine rule to get dis of mars from center.

def theta(alpha,beta,x,y):
    A=beta-y
    B=np.pi-(alpha-y)
    C=np.pi-A-B
    a=(np.sin(A)/np.sin(C))*(1+x) #used sine rule to get the side opposite to average sun vertex of triangle
    return y+np.arcsin((np.sin(B)*a)/radius(alpha,beta,x,y))

def geometric_mean(values):
    a = np.log(values)
    return np.exp(a.sum()/len(a))

def objective(x_y,alpha,beta):
    r=radius(alpha,beta,x_y[0],x_y[1])
    return np.log(np.mean(r))-np.log(geometric_mean(r))

from scipy.optimize import minimize,Bounds,shgo
# x_y=np.array([0.11442556, 0.86202794])
# [0.56917875 0.82095843]
bounds = Bounds([0, 0], [np.inf, 2*np.pi])

#Grid search for finding best x and y guesses
divisions=5 #controls granularity of grid
best_result=None
print('Starting grid search on ',divisions*divisions,' pairs')
for x in np.arange(0,10,10.0/divisions):
    for y in np.arange(0,2*np.pi,2*np.pi/divisions):
        x_y=(x,y)
        print('Checking optima for x and y guess:=',x_y)
        res=minimize(objective,x_y,args=(alpha,beta),bounds=bounds) #perform regression with this guess of x and y
        if best_result is None or (res.success and res.fun < best_result.fun):
            best_result=res
print('-----------------------------------')
print('-----------BEST RESULT-------------')
print (best_result)
print('-----------------------------------')

print ('Final x and y values',best_result.x)
print ('Objective loss',best_result.fun)
#computing radius and theta
r=radius(alpha,beta,best_result.x[0],best_result.x[1])
th=theta(alpha,beta,best_result.x[0],best_result.x[1])

print ('AM radius',np.mean(r))
print ('GM radius',geometric_mean(r))
print ('Objective',best_result.fun)
print('-----------------------------------')

for R,T in zip(r,th):
    print ('(r,theta)=',(R, T))
