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

def inside_radius(xi,yi,xj,yj,FR):
    if abs(xi%1 - xj%1) < FR and abs(yi%1 - yj%1) < FR:
        return True
    else:
        return False

def Particle(t, tau,camera):
    v = 0.5 #0, #2 ,6, 20
    T = 100
    mu = 0.001
    R = 1
    k_b = 1.38064852e-5
    Omega = 0
    x0 = 0
    y0 = 0
    r = 3
    theta0 = 0
    D_r = Dr(k_b,T,mu,R)
    print(D_r)
    D_t = Dt(k_b,T,mu,R)

    x0array = np.linspace(0.1,0.9,15)
    #x0array = [0.1,0.4]
    y0array = np.linspace(0.1,0.9,15)
    #y0array = [0.1,0.4]

    FR = 0.02
    N = 225
    particle_array = []
    xplot = []
    yplot = []
    v_array = []

    for x in x0array:
        for y in y0array:
            particle_array.append(particle.particle(x,y,random.uniform(-np.pi,np.pi)))

    for i in range(N):
        xplot.append(particle_array[i].getxcoords())
        yplot.append(particle_array[i].getycoords())
        #print(particle_array[i].gettheta())

    for tstep in t:
        print("######################timestep#####################",tstep)

        plt.scatter(xplot,yplot,s=r,color='black')
        #if tstep == 1.2662662662662663:
        #plt.show()
        camera.snap()

        xplot = []
        yplot = []

        for i in range(N):
            thetai = particle_array[i].gettheta()
            #print("thetai",thetai)

            xi = particle_array[i].getxcoords()
            yi = particle_array[i].getycoords()

            thetam = [particle_array[i].gettheta()]
            #print(thetam)
            average = 0
            for j in range(N):
                if i != j:
                    xj = particle_array[j].getxcoords()
                    yj = particle_array[j].getycoords()

                    if inside_radius(xi,yi,xj,yj,FR):
                        thetam.append(particle_array[j].gettheta())
            
            #print(thetam)
            for p in range(len(thetam)):
                #print(thetam[p])
                average = average + thetam[p]
            #print("average",average)
            average = (average)/(len(thetam))
            #print("average after",average)

            particle_array[i].setNEWtheta(average + np.sqrt(2*D_r*tau)*np.random.normal(0,np.pi))
            particle_array[i].setNEWx(particle_array[i].getxcoords() + v*np.cos(particle_array[i].getNEWtheta())*tau) #+ np.sqrt(2*D_t*tau)*np.random.normal(0,1))
            particle_array[i].setNEWy(particle_array[i].getycoords() + v*np.sin(particle_array[i].getNEWtheta())*tau) #+ np.sqrt(2*D_t*tau)*np.random.normal(0,1))
            #print(particle_array[i].getNEWtheta())

        sintotal = 0
        costotal = 0
        
        for u in range(N):
            particle_array[u].setxcoords(particle_array[u].getNEWx())
            particle_array[u].setycoords(particle_array[u].getNEWy())
            particle_array[u].settheta(particle_array[u].getNEWtheta())
            xplot.append(particle_array[u].getxcoords()%1)
            yplot.append(particle_array[u].getycoords()%1)
            #print(particle_array[u].gettheta())
            sintotal = sintotal + np.sin(particle_array[u].gettheta())
            costotal = costotal + np.cos(particle_array[u].gettheta())

        v = (1/N)*np.sqrt(sintotal**2 + costotal**2)
        v_array.append(v)

    print(D_r)
    return [t,v_array,xplot,yplot]

t_start = 0
t_end = 5
t_res = 500
tau = (t_start + t_end)/(t_res)
t = np.linspace(t_start, t_end, t_res)
camera = Camera(plt.figure())
solution = Particle(t, tau,camera)
anim = camera.animate(blit=True)
plt.show()
plt.figure()
plt.scatter(solution[2],solution[3],s=3,color='black')
plt.show()


plt.figure()
plt.title("Total Alignment of the system")
plt.ylabel("alignment index")
plt.xlabel("t")
plt.plot(solution[0],solution[1])
plt.show()
