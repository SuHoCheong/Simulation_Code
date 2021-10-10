import numpy as np
import particle_class as particle
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from celluloid import Camera
import random

def Dr(k_b,T,mu,R):
    return ((k_b*T)/(8*np.pi*mu*R**3))

def Dt(k_b,T,mu,R):
    return ((k_b*T)/(6*np.pi*mu*R))

def inside_radius(xi,yi,xj,yj,r,mod):
    if np.sqrt((xi%mod - xj%mod)**2 + (yi%mod - yj%mod)**2) < 2*r:
        return True
    else:
        return False

def Particle(t, tau,camera):
    v = 5 #0, #2 ,6, 20
    T = 100
    mu = 0.001
    R = 1.5
    k_b = 1.38064852e-5
    Omega = 0
    x0 = 0
    y0 = 0
    r = 5
    length = 2
    theta0 = 0
    mod = 100
    D_r = Dr(k_b,T,mu,R)
    D_t = Dt(k_b,T,mu,R)

    x0array = np.linspace(2,98,25)
    #x0array = [0.4,0.6]
    y0array = np.linspace(2,98,25)
    #y0array = [0.4,0.6]

    FR = 0.01 #Flocking Radius
    N = 625
    particle_array = []
    xplot = []
    yplot = []
    v_array = []
    
    ax = plt.gca()
    ax.cla()

    for x in x0array:
        for y in y0array:
            particle_array.append(particle.particle(x,y,random.uniform(-np.pi,np.pi)))

    for i in range(N):
        xplot.append(particle_array[i].getxcoords())
        yplot.append(particle_array[i].getycoords())
        #print(particle_array[i].gettheta())

    for tstep in t:
        print("######################timestep#####################",tstep)
        #ax = plt.gca()
        #ax.cla()
        #plt.scatter(xplot,yplot,s=r,facecolors='none', edgecolors='r')
        for i in range(len(xplot)):
            circle1 = plt.Circle((xplot[i],yplot[i]),R,color=(0.1,0.2,0.3),fill=False)
            ax.add_patch(circle1)

            plt.arrow(xplot[i], yplot[i], 2*np.cos(particle_array[i].gettheta()), 2*np.sin(particle_array[i].gettheta()))


        #plt.axis('square')
        ax.set_aspect('equal', adjustable='box')
        ax.set_xlim((0, mod))
        ax.set_ylim((0, mod))
        camera.snap()

        xplot = []
        yplot = []

        for i in range(N):

            particle_array[i].setNEWtheta(particle_array[i].gettheta() + np.sqrt(2*D_r*tau)*np.random.normal(0,1))
            particle_array[i].setNEWx(particle_array[i].getxcoords() + v*np.cos(particle_array[i].getNEWtheta())*tau) 
            particle_array[i].setNEWy(particle_array[i].getycoords() + v*np.sin(particle_array[i].getNEWtheta())*tau) 


        for u in range(N):

            xu = particle_array[u].getxcoords()
            yu = particle_array[u].getycoords()
            #print("xu",xu,"yu",yu)

            for j in range(N):
                if u != j:
                    xj = particle_array[j].getxcoords()
                    yj = particle_array[j].getycoords()
                    #print("xj",xj,"yj",yj)

                    if inside_radius(xu,yu,xj,yj,R,mod):
                        lolr = np.sqrt((xu%mod - xj%mod)**2 + (yu%mod - yj%mod)**2)
                        overlap = 2*R - lolr

                        x_comp = (xu%mod - xj%mod)
                        y_comp = (yu%mod - yj%mod)

                        unit_dir = [x_comp, y_comp]/lolr
                        unit_dir = unit_dir*overlap

                        particle_array[u].setNEWx(particle_array[u].getxcoords() + (unit_dir[0])) 
                        particle_array[u].setNEWy(particle_array[u].getycoords() + (unit_dir[1])) 

                        particle_array[j].setNEWx(particle_array[j].getxcoords() - (unit_dir[0])) 
                        particle_array[j].setNEWy(particle_array[j].getycoords() - (unit_dir[1])) 

                        #print("yes")
                        #print("xu",xu,"yu",yu,"xj",xj,"yj",yj)
                        #print("newxu",particle_array[u].getNEWx(),"newyu",particle_array[u].getNEWy(),"newxj",particle_array[j].getNEWx(),"newyj",particle_array[j].getNEWy())
                        #print("distance",lolr)


                particle_array[u].setxcoords(particle_array[u].getNEWx())
                particle_array[u].setycoords(particle_array[u].getNEWy())
                particle_array[u].settheta(particle_array[u].getNEWtheta())

                particle_array[j].setxcoords(particle_array[j].getNEWx())
                particle_array[j].setycoords(particle_array[j].getNEWy())
                particle_array[j].settheta(particle_array[j].getNEWtheta())

            xplot.append(particle_array[u].getxcoords()%mod)
            yplot.append(particle_array[u].getycoords()%mod)

    print( np.sqrt(2*D_r*tau))
    return [t,v_array,xplot,yplot]

t_start = 0
t_end = 10
t_res = 200
tau = (t_start + t_end)/(t_res)
t = np.linspace(t_start, t_end, t_res)
camera = Camera(plt.figure())
solution = Particle(t, tau,camera)
anim = camera.animate(blit=True)
plt.show()

'''
plt.figure()
plt.title("Total Alignment of the system")
plt.ylabel("&\mu&")
plt.xlabel("x")
plt.plot(solution[0],solution[1])
plt.show()
'''