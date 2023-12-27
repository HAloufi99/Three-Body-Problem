import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.animation import FuncAnimation
from numpy import sqrt

def lorentz(X, t):
    r1=X[:3]
    r2=X[6:9]
    r3=X[12:15]
    
    v1=X[3:6]
    v2=X[9:12]
    v3=X[15:18]

    c1 = -G*m2*(r1-r2)/(np.linalg.norm((r1-r2), 3)) - G*m3* (r1-r3)/(np.linalg.norm((r1-r3), 3))
    c2 = -G*m3*(r2-r3)/(np.linalg.norm((r2-r3), 3)) - G*m1* (r2-r1)/(np.linalg.norm((r2-r1), 3))
    c3 = -G*m1*(r3-r1)/(np.linalg.norm((r3-r1), 3)) - G*m2* (r3-r2)/(np.linalg.norm((r3-r2), 3))
    dr1dt = v1
    dv1dt = c1

    dr2dt = v2
    dv2dt = c2
    dr3dt = v3
    dv3dt = c3

    return np.hstack((dr1dt, dv1dt, dr2dt, dv2dt, dr3dt, dv3dt))
G = 6.6743*10**(-11)
m1 = 1*10**10
m2 = 1*10**10
m3 = 1*10**10

### SOME COOL SETTINGS
# CIRCUILING a ball ##
# r1=[2,0,0]
# r2=[0,0,0]
# r3=[-2,0,0]

# v1=[1/2,sqrt(3)/2,0]
# v2=[0,0/4,0]
# v3=[-1/2,-sqrt(3)/2,0]
## 


##
r1=[3,0,0]
r2=[0,0,0]
r3=[2,2,0]

v1=[1/2,sqrt(3)/2,0]
v2=[-sqrt(3)/2,1/4,0]
v3=[-1/2,-sqrt(3)/2,0]

t = np.linspace(0, 30, 1800)  # 30 seconds with 60 FPS


Origin = np.hstack((r1, v1, r2, v2, r3, v3)) 

sol = odeint(lorentz, Origin, t, mxstep=1000)

# Animation parameters
r = 0.1  # Radius of the particle

fig, ax = plt.subplots(figsize=(16, 6))
particle1 = Circle((r1[0], r2[1]), radius=r, color='b')
particle2 = Circle((r2[0], r2[1]), radius=r, color='g')
particle3 = Circle((r3[0], r3[1]), radius=r, color='r')

ax.add_patch(particle1)
ax.add_patch(particle2)
ax.add_patch(particle3)

ax.set_xlim([-2*4, 4*4])
ax.set_ylim([-4, 4])

# Animation function
def animate(frame):
    r1 = sol[frame, :3]
    r2 = sol[frame, 6:9]
    r3 = sol[frame, 12:15]

    particle1.set_center((r1[0], r1[1]))
    particle2.set_center((r2[0], r2[1]))
    particle3.set_center((r3[0], r3[1]))

    return particle1, particle2, particle3


# Create the animation
ani = FuncAnimation(fig, animate, frames=len(t), interval=16.67, blit=True)  # 1000 ms / 60 FPS = 16.67 ms

# ani.save('gravity_simulation.gif', writer="ImageMagick" ,frames=60)

plt.show()
