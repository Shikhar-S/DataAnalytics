def get_as_radian(degree,minute,second=0,zodiac_index=0):
    base_angle=zodiac_index*30
    deg_angle=base_angle+degree+(minute/60.0)+(second/3600.0)
    return np.radians(deg_angle)
#returns x y coordinates of mars from triangulation
def get_coordinates(A,B,C,D):
    x=(np.tan(D)*np.cos(B)-np.tan(C)*np.cos(A)-np.sin(B)+np.sin(A))/(np.tan(D)-np.tan(C))
    y=np.sin(A)+np.tan(C)*(x-np.cos(A))
    return (x,y)



def circle_fit_objective(r,x,y):
    return np.sum((np.sqrt(x**2+y**2)-r)**2)

def ellipse_fit_objective(params,u,v):
    a=params[0]
    eccentricity=params[1]
    b=np.sqrt(a*a*(1-eccentricity**2))
    delta=params[2]
    X=u*np.cos(delta)+v*np.sin(delta)
    Y=v*np.cos(delta)-u*np.sin(delta)+np.sqrt(a**2-b**2)
    f=(X/a)**2+(Y/b)**2-1
    return np.sum(f**2)

#returns mars heliocentric lats from geocentric lats and dis of mars from sun
def get_mars_helio_lat(mars_lat_geo,radius):
    y= (radius-1)*np.tan(mars_lat_geo)
    return np.arctan(y/radius)

#returns mars' 3d coordinates from latitude(alpha) and longitude(beta)
def get_mars_3d(alpha,beta):
    X=-np.cos(alpha)*np.sin(beta)
    Y=np.cos(alpha)*np.cos(beta)
    Z=np.sin(alpha)
    return (X,Y,Z)

def cost(x,y,z,normal_dirn): #squared sum of distances of (x,y,z) from plane with normal_dirn as normal
    dis=(x*normal_dirn[0]+y*normal_dirn[1]+z*normal_dirn[2])/np.sqrt(np.sum(normal_dirn**2))
    return np.sum(dis**2)
#check pdf for computation
def grad_sig_f_sq(x,y,z,normal_dirn):
    S=np.sum(normal_dirn**2)
    A,B,C=normal_dirn
    V=(A*x+B*y+C*z)
    df_da=(1/S)*(2*np.sum(V*x))-(2*A/(S*S))*(np.sum(V**2))
    df_db=(1/S)*(2*np.sum(V*y))-(2*B/(S*S))*(np.sum(V**2))
    df_dc=(1/S)*(2*np.sum(V*z))-(2*C/(S*S))*(np.sum(V**2))
    return np.array([df_da,df_db,df_dc])

#exact line search to get best alpha using dichotomous search
def dichotomous(x,y,z,theta,grad,left,right,tolerance=1e-5):
    while(right-left>tolerance):
        mid=(left+right)/2
        mid_l=mid-tolerance/4
        mid_r=mid+tolerance/4
        cost_l=cost(x,y,z,theta-grad*mid_l)
        cost_r=cost(x,y,z,theta-grad*mid_r)
        if cost_l>cost_r:
            left=mid_l
        else:
            right=mid_r
    return (left+right)/2

#steepest descent using exact line search.
def grad_desc(x,y,z):
    theta=np.random.uniform(low=-1, high=1, size=(3,))
    eps=1e-5
    gf=grad_sig_f_sq(x,y,z,theta)
    mod_gf=np.sum(gf**2)
    i=0
    while(mod_gf>eps):
        i+=1
        mod_gf=np.sum(gf**2)
        gf=grad_sig_f_sq(x,y,z,theta)
        alpha_min=dichotomous(x,y,z,theta,gf,0,1)
        theta=theta-alpha_min*gf
    return theta,cost(x,y,z,theta)

#returns inclination of a vec from xy plane (ecliptic plane)
def get_inclination(plane):
    sin_inclination=np.sqrt(np.sum(np.cross(plane,np.array([0,0,1]))**2))/(np.sqrt(np.sum(plane**2))) #001 is normal to xy plane.(horizontal plane)
    inclination=np.arcsin(sin_inclination)
    return np.degrees(inclination)

#returns mars' 3d coordinates from equation of plane and x,y coordinates (also the projection on ecliptic)
def get_mars_plane_projection(normal_dirn,x,y):
    z=(normal_dirn[0]*x+normal_dirn[1]*y)/normal_dirn[2]
    return (x,y,z)

#given n, n>1 3d points lying on a plane and normal to the plane this method returns the points in 2d coordinates u-v on plane.
def coordinate_transform(X,Y,Z,normal_dirn):
    u=np.array([X[0]-X[1],Y[0]-Y[1],Z[0]-Z[1]]) #vector on plane
    v=np.cross(normal_dirn,u)
    u=u/np.sqrt(np.sum(u**2)) #unit vector
    v=v/np.sqrt(np.sum(v**2)) 
    #let U,V be the respective coordinates in this new system then,
    U=np.dot(u,(X,Y,Z))
    V=np.dot(v,(X,Y,Z))
    return (U,V)

#returns uv coordinates for given point with parameter t on ellipse with eccentricity,semi major axis, rotation and focus on origin
def get_uv(t,delta,A,ecc):
    B=A*np.sqrt(1-ecc**2)
    u_=A*np.cos(t)
    v_=B*np.sin(t)-A*ecc
    x=u_*np.cos(delta)-v_*np.sin(delta)
    y=u_*np.sin(delta)+v_*np.cos(delta)
    return (x,y)

import pandas as pd
import numpy as np
df=pd.read_csv('../data/01_data_mars_triangulation.csv')

#ques 2
A=[] #heliocentric angle of earth
B=[] #paired heliocentric angle of earth
C=[] #geocentric angle of mars
D=[] #paired geocentric angle of mars

for i in range(1,6):
    rows=df[df['PairIndex']==i]
    A.append(get_as_radian(rows.iloc[0]['DegreeEarthLocationHelioCentric'],rows.iloc[0]['MinuteEarthLocationHelioCentric']))
    B.append(get_as_radian(rows.iloc[1]['DegreeEarthLocationHelioCentric'],rows.iloc[1]['MinuteEarthLocationHelioCentric']))
    C.append(get_as_radian(rows.iloc[0]['DegreeMarsLocationGeoCentric'],rows.iloc[0]['MinuteMarsLocationGeoCentric']))
    D.append(get_as_radian(rows.iloc[1]['DegreeMarsLocationGeoCentric'],rows.iloc[1]['MinuteMarsLocationGeoCentric']))
A=np.array(A)
B=np.array(B)
C=np.array(C)
D=np.array(D)

x,y=get_coordinates(A,B,C,D)

from scipy.optimize import minimize,Bounds
r=np.random.rand(1)
res=minimize(circle_fit_objective,r,args=(x,y))
Radius=res.x
print("Radius of Circle on Ecleptic plae",Radius)


#ques 3
df_opp=pd.read_csv('../data/01_data_mars_opposition.csv')
mars_lat_geo=get_as_radian(df_opp['LatDegree'],df_opp['LatMinute'])
mars_lat_helio=get_mars_helio_lat(mars_lat_geo,Radius)
print("Mars' heliocentric latitudes(radians)",mars_lat_helio)
mars_long_helio=get_as_radian(df_opp['Degree'],df_opp['Minute'],df_opp['Second'],df_opp['ZodiacIndex'])

m_x,m_y,m_z=get_mars_3d(mars_lat_helio,mars_long_helio)
res=grad_desc(m_x,m_y,m_z)
inclination=get_inclination(res[0])
print("Inclination with ecliptic",int(inclination),"degrees", (inclination-int(inclination))*60,"minutes")
print("Cost of regression",res[1])
mars_plane=(res[0]/np.sqrt(np.sum(res[0]**2)))


#ques 4
X,Y,Z=get_mars_plane_projection(mars_plane,x,y)
print("Mars' five locations on a plane (X,Y,Z) arrays",(X,Y,Z))
U,V=coordinate_transform(X,Y,Z,mars_plane)
#fitting circle
r_mc=np.random.rand(1)
res_mc=minimize(circle_fit_objective,r_mc,args=(U,V))
print("Sum of losses for circle on mars plane",res_mc.fun)
#fitting ellipse
r_me=np.random.rand(3)
bounds=Bounds([0,0,0],[np.inf,1,2*np.pi])
res_me=minimize(ellipse_fit_objective,r_me,args=(U,V),bounds=bounds)
print("Sum of losses for ellipse on mars' plane",res_me.fun)

#plots

import matplotlib.pyplot as plt

plt.plot(0,0,'ro')

t=np.linspace(0,2*np.pi,100)

circ_x=Radius*np.cos(t)
circ_y=Radius*np.sin(t)
plt.plot(circ_x,circ_y,'y')

A,ecc,delta=res_me.x
ellip_x,ellip_y=get_uv(t,delta,A,ecc)
plt.plot(ellip_x,ellip_y)
M_x,M_y,M_z=get_mars_plane_projection(mars_plane,m_x,m_y)
M_u,M_v=coordinate_transform(M_x,M_y,M_z,mars_plane)

T=np.arctan(M_v/M_u)
m12_x,m12_y=get_uv(T,delta,A,ecc)
plt.plot(m12_x,m12_y,'o')
plt.show()