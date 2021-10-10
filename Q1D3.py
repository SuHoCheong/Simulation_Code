import numpy as np
import new_particle_class as new_particle
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation, PillowWriter 
from matplotlib.animation import FuncAnimation, writers
Writer = animation.writers['pillow']
writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800) 
from celluloid import Camera
import random
import math
from numba import njit

def inside_boundary(y0,L_y):
    if y0 < L_y/2 and y0 > -L_y/2:
        return False
    else:
        return True

def new_y_position(y1,L_y,i,particle_array):
    if y1 > L_y/2:
        dif = y1 - L_y/2 
        particle_array[i].setnewycoords(L_y/2 - dif)
    elif y1 < -L_y/2:
        dif = y1 + L_y/2
        particle_array[i].setnewycoords(-L_y/2 - dif)

def new_x_position(x1,L_x,i,particle_array):
    if x1 > L_x/2:
        dif = x1 - L_x/2
        particle_array[i].setnewxcoords(-L_x/2 + dif)
    elif x1 < -L_x/2:
        dif = x1 + L_x/2
        particle_array[i].setnewxcoords(L_x/2 + dif)

def effective_charge(M,magnetic_p):
    return M*np.sqrt(magnetic_p/4*np.pi)

def pClosest(points, K,R):
    points.sort(key = lambda K: np.sqrt((K[0]-R[0])**2 + (K[1]-R[1])**2))
    return points[:K]

def Average(lst):
    return sum(lst) / len(lst)

def Repulsive_Dipole_interaction(Qi,Qj,xi,yi,xj,yj,r,flag,L_x):
    if flag == 1:
        return ((Qi*Qj*3)/((np.sqrt((xi-xj)**2+(yi-yj)**2))**4)*((np.array([xi-xj,yi-yj]))/(np.sqrt((xi-xj)**2+(yi-yj)**2))))
    if flag == 2:
        return ((Qi*Qj*3)/((np.sqrt((xi-(xj-L_x))**2+(yi-yj)**2))**4)*((np.array([(xi)-(xj-L_x),yi-yj]))/(np.sqrt((xi-(xj-L_x))**2+(yi-yj)**2))))
    if flag == 3:
        return ((Qi*Qj*3)/((np.sqrt((xi-(xj+L_x))**2+(yi-yj)**2))**4)*((np.array([(xi)-(xj+L_x),yi-yj]))/(np.sqrt((xi-(xj+L_x))**2+(yi-yj)**2))))

def simulation(N_A,N_B,N,N_M,L_x,L_y,F_drx,r,t,tau,kb,T,Q_A,Q_B,Q_D,D_A,D_B,rho,mag_p):
    particle_array = []
    xplot = []
    yplot = []
    TotalAvA = []
    TotalAvB = []

    N_Bx = [-4.5,-4.1,-3.5,-2.9,-2.3,-1.8,-1.1,-0.5,0.1,0.7,1.3,1.8,2.4,3.0,3.6,4.2,4.7,-3.7,-2.5,-1.5,-0.3,1.3,2.5,3.9,-4.7,-4.3,-3.8,-2.8,-2.1,-1.5,-0.8,0.2,0.9,1.5,2.1,2.7,3.2,3.9,4.5,-0.21]
    N_By = [0.59,0.59,0.59,0.59,0.59,0.59,0.59,0.59,0.59,0.59,0.59,0.59,0.59,0.59,0.59,0.59,0.59,0.1,-0.1,0.1,-0.1,0.1,-0.1,0.1,-0.59,-0.59,-0.59,-0.59,-0.59,-0.59,-0.59,-0.59,-0.59,-0.59,-0.59,-0.59,-0.59,-0.59,-0.59,-0.56]

    #N_Ax = [-4.7,-4.4,-4.5,-4.3,-4.1,-3.3,-3.25,-3.15,3.07,-3.05,-3.0,-2.2,-2.1,-1.8,-1.7,-1.0,0.9,-0.7,-0.6,0.2,0.3,0.5,0.6,0.9,1.0,1.5,1.7,1.9,2.0,2.2,2.8,2.9,3.3,3.4,3.5,3.6,4.4,4.5,4.6,4.7]
    #N_Ay = [-0.3,0.2,-0.05,0.16,-0.2,-0.4,0,0.26,-0.49,-0.3,0.0,-0.25,0.15,0.1,-0.3,0.1,-0.1,0.15,-0.15,-0.1,0.2,-0.3,0.1,-0.1,0.1,0.4,0.2-0.1,0.05,-0.1,0.3,0.2,-0.2,0.1,-0.1,0.4,0.1,0.23,-0.3,0,0.3,0.31]

    N_Ax = [-4.6,-4.5,-4.45,-4.2,-4.1,-3.4,-3.35,-3.2,-3.15,-3.05,-3.0,-2.2,-2.1,-1.8,-1.7,-1.0,-0.9,-0.7,-0.6,0.3,0.4,0.5,0.55,0.8,0.9,1.5,1.6,1.8,1.9,2.1,2.8,2.9,3.2,3.3,3.5,3.6,4.4,4.5,4.6,4.7]
    N_Ay = [-0.28,0.28,-0.05,0.2,-0.24,-0.4,0,0.26,-0.48,-0.3,0,-0.25,0.15,0.1,-0.3,0.1,-0.1,0.15,-0.15,-0.1,0.2,-0.3,0.1,-0.1,0.1,0.3,-0.15,0.25,-0.15,0.25,0.2,-0.2,0.17,-0.19,0.2,-0.31,-0.33,0,0.3,-0.34]

    #N_Bx = [-45,-41,-35,-29,-23,-18,-11,-5,1,7,13,18,24,30,36,42,47,-37,-25,-15,-3,13,25,39,-49,-43,-38,-28,-21,-15,-8,-3,2,9,15,21,27,32,39,45]
    #N_By = [5.5,5.5,5.5,5.5,5.5,5.5,5.5,5.5,5.5,5.5,5.5,5.5,5.5,5.5,5.5,5.5,5.5,-0.1,0.1,-0.1,0.1,-0.1,0.1,-0.1,-5.5,-5.5,-5.5,-5.5,-5.5,-5.5,-5.5,-5.5,-5.5,-5.5,-5.5,-5.5,-5.5,-5.5,-5.5,-5.5]

    #N_Ax = [-47,-44,-45,-43,-41,-33,-32.5,-31.5-30,7,-30.5,-30,-22,-21,-18,-17,-10,-9,-7,-6,2,3,5,6,9,10,15,17,19,20,22,28,29,33,34,35,36,44,45,46,47]
    #N_Ay = [-3,2,-0.5,1.6,-2,-4,0,2.6,-6,-3,0,-2.5,1.5,1,-3,1,-1,1.5,-1.5,-1,2,-3,1,-1,1,4,2-1,0.5,-1,3,2,-2,1,-1,8,1,2.3,-3,0,3,3.1]

    N_msx = [-4.9]
    N_msy = [0]

    #print(len(N_Ay))

    for i in range(N_M):
        particle_array.append(new_particle.new_particle(N_msx[i],N_msy[i],3))

    for i in range(N_B):
        particle_array.append(new_particle.new_particle(N_Bx[i],N_By[i],2))

    for i in range(N_A):
        particle_array.append(new_particle.new_particle(N_Ax[i],N_Ay[i],1))

    #for i in range(N_B):
        #flag = True
        #while flag:
            #particle_array.append(new_particle.new_particle(random.uniform(-L_x/2-0.3,L_x/2-0.3),random.uniform(-L_y/4,L_y/4),2))
            #f1 = 0
            #for j in range(len(particle_array)):
                #if i != j:
                    #if np.sqrt((particle_array[i].getxcoords()-particle_array[j].getxcoords())**2+(particle_array[i].getycoords()-particle_array[j].getycoords())**2) < 1:
                        #f1 = 1
            #if f1 == 0:
                #flag = False
            #else:
                #particle_array.pop()

    #for i in range(N_A):
        #flag = True
        #i = N_B + i
        #while flag:
            #particle_array.append(new_particle.new_particle(random.uniform(-L_x/2-0.3,L_x/2-0.3),random.uniform(-L_y/4,L_y/4),1))
            #f1 = 0
            #for j in range(len(particle_array)):
                #if i != j:
                    #if np.sqrt((particle_array[i].getxcoords()-particle_array[j].getxcoords())**2+(particle_array[i].getycoords()-particle_array[j].getycoords())**2) < 1.5:
                        #f1 = 1
            #if f1 == 0:
                #flag = False
            #else:
                #particle_array.pop()

    #for i in range(N_M):
        #particle_array.append(new_particle.new_particle(random.uniform(-L_x/2-0.3,L_x/2-0.3),random.uniform(-L_y/2-0.3,L_y/2-0.3),3))

    for p in range(len(particle_array)):
        xplot.append(particle_array[p].getxcoords())
        yplot.append(particle_array[p].getycoords())
    
    fig, ax = plt.subplots()
    plt.cla()
    plt.xlim([-L_x/2, L_x/2])
    plt.ylim([-L_y/2, L_y/2])
    
    for l in range(len(xplot)):
        if particle_array[l].gettype() == 1:
            scat = plt.scatter(xplot[l],yplot[l],s = 0.6,c="black",marker='v')
        elif particle_array[l].gettype() == 2:
            scat = plt.scatter(xplot[l],yplot[l],s = 4,c="blue",marker=',')
        else:
            scat = plt.scatter(xplot[l],yplot[l],s = 6,c="red",marker='<')
        
    ani = animation.FuncAnimation(fig, update_plot, frames=range(int(t/tau)),
                                  fargs=(N_A,N_B,N,N_M,L_x,L_y,F_drx,r,t,tau,scat,particle_array,xplot,yplot,kb,T,Q_A,Q_B,Q_D,D_A,D_B,rho,mag_p,TotalAvA,TotalAvB),repeat=False)
    
    
    ax.set_aspect(1)
    
    #ani.save('Q1D_channel.gif', writer=writer)

    #ax.set_aspect(1/8)
    plt.show()


    print("Total Average A",Average(TotalAvA))
    print("Total Average B",Average(TotalAvB))
    print("Fdrx",F_drx)
    plt.close()

    xlol = []
    ylol = []



    for u in range(len(particle_array)):
        xlol.append(particle_array[u].getxcoords())
        ylol.append(particle_array[u].getycoords())

    #print(xlol)
    #print(ylol)

    return [Average(TotalAvB)*10,Average(TotalAvA)*10]

def update_plot(i,N_A,N_B,N,N_M,L_x,L_y,F_drx,r,t,tau,scat,particle_array,xplot,yplot,kb,T,Q_A,Q_B,Q_D,D_A,D_B,Rho,mag_p,TotalAvA,TotalAvB):
    print("time:",i,"of",int(t/tau))

    plt.cla()
    plt.xlim([-L_x/2, L_x/2])
    plt.ylim([-L_y/2, L_y/2])

    xplot = []
    yplot = []
    label = []
    both = []
    averagexA = []
    averagexB = []

    for j in range(1,N):
        both.append([particle_array[j].getxcoords(),particle_array[j].getycoords(),j])

    K = N-1
    R = [particle_array[0].getxcoords(),particle_array[0].getycoords()]
 
    triple = pClosest(both, K,R)

    for j in range(N_M):
        particle_array[j].setxcoords(particle_array[j].getnewxcoords())
        particle_array[j].setycoords(particle_array[j].getnewycoords())

        if particle_array[j].gettype() == 1:
            M = M_A
            Qi = Q_A
        elif particle_array[j].gettype() == 2:
            M = M_B
            Qi = Q_B
        elif particle_array[j].gettype() == 3:
            Qi = Q_D

        Vij = np.array([0,0])

        for p in range(N):

            xi = particle_array[j].getxcoords()%L_x
            xj = particle_array[p].getxcoords()%L_x
            yi = particle_array[j].getycoords()
            yj = particle_array[p].getycoords()

            if j != p:
                if particle_array[p].gettype() == 1:
                    M1 = M_A
                    Qj = Q_A
                    #Qj = Qi
                elif particle_array[p].gettype() == 2:
                    M1 = M_B
                    Qj = Q_B
                    #Qj = 8*Qi
                elif particle_array[p].gettype() == 3:
                    Qj = Q_D

                if np.sqrt((xi-xj)**2+(yi-yj)**2) < r: 
                    oVij = Repulsive_Dipole_interaction(Qi,Qj,xi,yi,xj,yj,r,1,L_x)                
                elif abs(np.sqrt((xi-xj)**2+(yi-yj)**2)-L_x) < r:
                    if xi < xj:
                        oVij = Repulsive_Dipole_interaction(Qi,Qj,xi,yi,xj,yj,r,2,L_x)
                    elif xi > xj:
                        oVij = Repulsive_Dipole_interaction(Qi,Qj,xi,yi,xj,yj,r,3,L_x)
                else:
                    oVij = np.array([0,0])
                #rint(oVij)

                #Vij = Vij + oVij
                Vij = Vij + np.array([0,0])
        
        #print("particle",j,Vij)

        if particle_array[j].gettype() == 1:
            D = D_A
        elif particle_array[j].gettype() == 2:
            D = D_B       
        elif particle_array[j].gettype() == 3:
            D = D_A

        if particle_array[j].gettype() == 3:
            F_d = F_drx
        else:
            F_d = 0  

        C = D/(k_b*T)

        particle_array[j].setnewxcoords(particle_array[j].getxcoords() + (1/C)*tau*(Vij[0])+F_d*tau+ np.sqrt(2*k_b*T*(1/D))*np.random.normal(0,0.1)) #lolololo
        particle_array[j].setnewycoords(particle_array[j].getycoords() + (1/C)*tau*(Vij[1]+0)+ np.sqrt(2*k_b*T*(1/D))*np.random.normal(0,0.1))

        if inside_boundary(particle_array[j].getnewycoords(),L_y):
            new_y_position(particle_array[j].getnewycoords(),L_y,j,particle_array)
        new_x_position(particle_array[j].getnewxcoords(),L_x,j,particle_array)

    for y in range(len(triple)):
        j = triple[y][2]
        label.append(j)
        particle_array[j].setxcoords(particle_array[j].getnewxcoords())
        particle_array[j].setycoords(particle_array[j].getnewycoords())

        if particle_array[j].gettype() == 1:
            M = M_A
            Qi = Q_A
        elif particle_array[j].gettype() == 2:
            M = M_B
            Qi = Q_B
        elif particle_array[j].gettype() == 3:
            Qi = Q_D

        Vij = np.array([0,0])

        for p in range(N):

            xi = particle_array[j].getxcoords()%L_x
            xj = particle_array[p].getxcoords()%L_x
            yi = particle_array[j].getycoords()
            yj = particle_array[p].getycoords()

            if j != p:
                if particle_array[p].gettype() == 1:
                    M1 = M_A
                    Qj = Q_A
                    #Qj = Qi
                elif particle_array[p].gettype() == 2:
                    M1 = M_B
                    Qj = Q_B
                    #Qj = 8*Qi
                elif particle_array[p].gettype() == 3:
                    Qj = Q_D

                if np.sqrt((xi-xj)**2+(yi-yj)**2) < r: 
                    oVij = Repulsive_Dipole_interaction(Qi,Qj,xi,yi,xj,yj,r,1,L_x)                
                elif abs(np.sqrt((xi-xj)**2+(yi-yj)**2)-L_x) < r:
                    if xi < xj:
                        oVij = Repulsive_Dipole_interaction(Qi,Qj,xi,yi,xj,yj,r,2,L_x)
                    elif xi > xj:
                        oVij = Repulsive_Dipole_interaction(Qi,Qj,xi,yi,xj,yj,r,3,L_x)
                else:
                    oVij = np.array([0,0])
                #print(oVij)
                Vij = Vij + oVij
        
        #print("particle",j,Vij)

        if particle_array[j].gettype() == 1:
            D = D_A
        elif particle_array[j].gettype() == 2:
            D = D_B       
        elif particle_array[j].gettype() == 3:
            D = D_A

        if particle_array[j].gettype() == 3:
            F_d = F_drx
        else:
            F_d = 0  

        C = D/(k_b*T)

        particle_array[j].setnewxcoords(particle_array[j].getxcoords() + (1/C)*tau*(Vij[0])+F_d*tau+ np.sqrt(2*k_b*T*(1/D)*tau)*np.random.normal(0,1))
        particle_array[j].setnewycoords(particle_array[j].getycoords() + (1/C)*tau*(Vij[1]+0)+ np.sqrt(2*k_b*T*(1/D)*tau)*np.random.normal(0,1))

        if inside_boundary(particle_array[j].getnewycoords(),L_y):
            new_y_position(particle_array[j].getnewycoords(),L_y,j,particle_array)
        new_x_position(particle_array[j].getnewxcoords(),L_x,j,particle_array)


        #print(np.sqrt((R[0]-particle_array[j].getnewxcoords())**2+(R[1]-particle_array[j].getnewycoords())**2))
        if np.sqrt((R[0]-particle_array[j].getnewxcoords())**2+(R[1]-particle_array[j].getnewycoords())**2) < 0.8:
            #print(y)
            if particle_array[j].gettype() == 1:
                averagexA.append((particle_array[j].getnewxcoords()%L_x-particle_array[j].getxcoords()%L_x)**2)
            elif particle_array[j].gettype() == 2:
                if (particle_array[j].getnewxcoords()%L_x-particle_array[j].getxcoords()%L_x)**2 < 100:
                    averagexB.append((particle_array[j].getnewxcoords()%L_x-particle_array[j].getxcoords()%L_x)**2)
        
        particle_array[j].setxcoords(particle_array[j].getnewxcoords())
        particle_array[j].setycoords(particle_array[j].getnewycoords())

    print("A",Average(averagexA))
    print("B",Average(averagexB))

    TotalAvA.append(Average(averagexA))
    TotalAvB.append(Average(averagexB))

    for p in range(len(particle_array)):
        xplot.append(particle_array[p].getxcoords())
        yplot.append(particle_array[p].getycoords())
        label.append(p)

    for l in range(len(xplot)):
        if particle_array[l].gettype() == 1:
            scat = plt.scatter(xplot[l],yplot[l],s = 10,c="red",marker='v')
        elif particle_array[l].gettype() == 2:
            scat = plt.scatter(xplot[l],yplot[l],s = 20,c="blue",marker=',')
            #scat = plt.text(xplot[l],yplot[l],label[l])
        else:
            scat = plt.scatter(xplot[l],yplot[l],s = 30,c="black",marker='<') 
            #circle1 = plt.Circle((xplot[l],yplot[l]),0.8,color=(0.1,0.2,0.3),fill=False)
            #scat = plt.gca().add_patch(circle1)

    #if i == 1:
    #    plt.show()

    return scat

N_A = 40
N_B = 40
N_M = 1

N = N_A + N_B + N_M

M_A = 0.5 #This value varies
M_B = 1
M_R = M_A/M_B

Q_A = 0.57#40
Q_B = 8*Q_A
Q_D = 8*Q_A

D_A = 3
D_B = 0.7

L_x = 10
L_y = 1.2

#F_drx = [0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3,3.2,3.4,3.6,3.8,4,4.2,4.4,4.6,4.8,5,5.2,5.4,5.6,5.8,6,6.2,6.4,6.6,6.8,7,7.2,7.4,7.6,7.8,8,8.2,8.4,8.6,8.8,9] #0.1, 0.4, 2, 10
#F_drx = [0.6]
F_drx = [5]

r = 6
t = 1
tau = 0.01    

T = 293
kb =  k_b = 1.38064852e-5
mag_p = 4*np.pi*10

rho = N/(L_x*L_y)
print("RHO:",rho)
print("fraction:",N_A,N_B)
print("M_R:",M_R)

TOTALAVERAGEB = []
TOTALAVERAGEA = []

for f_drx in F_drx:
    #print(f_drx)
    returned = simulation(N_A,N_B,N,N_M,L_x,L_y,f_drx,r,t,tau,kb,T,Q_A,Q_B,Q_D,D_A,D_B,rho,mag_p)
    TOTALAVERAGEB.append(returned[0])
    TOTALAVERAGEA.append(returned[1])
print(TOTALAVERAGEA)
print(TOTALAVERAGEB)

plt.plot(F_drx,TOTALAVERAGEB,'-o')
plt.show()