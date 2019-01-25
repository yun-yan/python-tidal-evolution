import matplotlib
matplotlib.use('Qt4Agg')

import numpy as np
from matplotlib import pyplot
from matplotlib import animation
from math import sin

data = np.zeros((1000,4))

data[:,0] = [20*(1+sin(float(x)/200)) for x in range(1000)]
data[:,1] = [20*(1+sin(float(x)/100)) for x in range(1000)]
data[:,2] = [20*(1+sin(float(x)/50)) for x in range(1000)]
data[:,3] = [20*(1+sin(float(x)/25)) for x in range(1000)]

fig=pyplot.figure()
ax = pyplot.axes(xlim=(0, 40), ylim=(0, 40))

circle1=pyplot.Circle((data[0,0],1.0),0.2,fc='y')
circle2=pyplot.Circle((data[0,1],1.0),0.2,fc='g')
circle3=pyplot.Circle((data[0,2],1.0),0.2,fc='r')
circle4=pyplot.Circle((data[0,3],1.0),0.2,fc='b')

def init():
    circle1.center=(data[0,0],1)
    circle2.center=(data[0,1],1) 
    circle3.center=(data[0,2],1)
    circle4.center=(data[0,3],1)
    ax.add_patch(circle1)
    ax.add_patch(circle2)
    ax.add_patch(circle3)
    ax.add_patch(circle4)
    return circle1, circle2, circle3, circle4

def animate(i):
    # for state in data:
    circle1.center=(data[i,0],1)
    circle2.center=(data[i,1],1)
    circle3.center=(data[i,2],1)
    circle4.center=(data[i,3],1)
    return circle1, circle2, circle3, circle4


anim=animation.FuncAnimation(fig,animate,init_func=init,frames=1000,blit=True)

#pyplot.show()
anim.save('basic_animation2.mp4', fps=60)