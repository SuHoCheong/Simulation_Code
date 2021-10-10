import numpy as np 
from matplotlib import pyplot as plt 
from scipy.integrate import odeint
import random

def inside_boundary(x,y,r):
    if (x**2 + y**2 < r**2):
        return True
    else:
        return False

def newlocation(x0,y0,x1, y1, a, b,r):
    m = (y1 - y0)/(x1 - x0)
    #print("m",m)
    d = (y1) - (m*x1)
    #print("d",d)
    interx0 = (a + b*m - d*m + np.sqrt(r**2*(1 + m**2) - (b-m*a-d)**2))/(1 + m**2)
    interx1 = (a + b*m - d*m - np.sqrt(r**2*(1 + m**2) - (b-m*a-d)**2))/(1 + m**2)

    intery0 = (d + a*m + b*m**2 + m*np.sqrt(r**2*(1 + m**2) - (b-m*a-d)**2))/(1 + m**2)
    intery1 = (d + a*m + b*m**2 - m*np.sqrt(r**2*(1 + m**2) - (b-m*a-d)**2))/(1 + m**2)

    if np.sqrt((interx0 - x1)**2+(intery0 - y1)**2) > np.sqrt((interx1 - x1)**2+(intery1 - y1)**2):
        interx = interx1
        intery = intery1
    else:
        interx = interx0
        intery = intery0

    tangent_m = (interx - a)/(b- intery)
    #print("tangent_m:",tangent_m)
    tangent_d = (intery) - (tangent_m*interx)
    #print("tangent_d:",tangent_d)

    normal_m = -(1)/(tangent_m)
    #print("normal_m:",normal_m)
    normal_d = (intery) - (normal_m*interx)
    #print("normal_d:",normal_d)

    LENGTH = np.sqrt(1**2 + normal_m**2)

    unit_normal_vector = [1/LENGTH,normal_m/LENGTH]
    #print("normal:",unit_normal_vector)

    newx = x1 - 2*np.dot(unit_normal_vector,(np.dot(unit_normal_vector,(x1-interx))))
    newy = y1 - 2*np.dot(unit_normal_vector,(np.dot(unit_normal_vector,(y1-intery))))

    return[newx, newy]

def plot(x,y,a,b,r):

    theta = np.linspace(0, 2*np.pi, 100)

    r = np.sqrt(9)

    x1 = r*np.cos(theta)
    x2 = r*np.sin(theta)

    fig, ax = plt.subplots(1)

    ax.plot(x1, x2)
    ax.set_aspect(1)

    plt.xlim(-3.5,3.5)
    plt.ylim(-3.5,3.5)
    plt.plot(x,y)
    plt.show()

    plt.grid(linestyle='--')

def Dr(k_b,T,mu,R):
    return ((k_b*T)/(8*np.pi*mu*R**3))

def Dt(k_b,T,mu,R):
    return ((k_b*T)/(6*np.pi*mu*R))

def Particle(t, tau):
    v = 20 #0, #2 ,6, 20
    T = 300
    mu = 0.001
    R = 1
    k_b = 1.38064852e-5
    Omega = 0
    x0 = 0
    y0 = 0
    r = 3
    theta0 = 0
    D_r = Dr(k_b,T,mu,R)
    D_t = Dt(k_b,T,mu,R)

    x = [x0]
    y = [y0]
    theta = [theta0]

    for i in range(len(t)):
        theta.append(theta[i] + Omega*tau + np.sqrt(2*D_r*tau)*np.random.normal(0,1))
        x.append(x[i] + v*np.cos(theta[i])*tau + np.sqrt(2*D_t*tau)*np.random.normal(0,1))
        y.append(y[i] + v*np.sin(theta[i])*tau + np.sqrt(2*D_t*tau)*np.random.normal(0,1))

        if (inside_boundary(x[i+1],y[i+1],r)):
            continue
        else:
            new_location = newlocation(x[i],y[i],x[i+1],y[i+1],x0,y0,r)
            x[i+1] = new_location[0]
            y[i+1] = new_location[1]

    plot(x,y,0,0,r)

def MSD(T, TAU):
    v = 5
    Temp = 300
    mu = 0.001
    R = 1
    k_b = 1.38064852e-5
    Omega = 0
    x0 = 0
    y0 = 0
    r = 3
    theta0 = 0
    D_r = Dr(k_b,Temp,mu,R)
    print(1/D_r)
    D_t = Dt(k_b,Temp,mu,R)
    MSD = (4*D_t + 2*v**2*(1/D_r))*T + (2*(v**2*(1/D_r)**2))*(np.exp(-(2*T)/(1/D_r)) - 1)

    plt.plot(T,MSD)
    plt.show()

t_start = 0
t_end = 20
t_res = 10000
tau = (t_start + t_end)/(t_res)
t = np.linspace(t_start, t_end, t_res)
Particle(t, tau)

T_START = 10**(-2)
T_END = 10**(2)
T_RES = 100000
TAU = (T_START + T_END)/(T_RES)
T = np.linspace(T_START, T_END, T_RES)
#MSD(T, TAU)
