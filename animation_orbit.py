import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import patches



plt.rcParams['figure.figsize'] = (8.0, 8.0)

pluto = plt.Circle((0, 0), 1, color='r')
pluto_ellipse =patches.Ellipse((0, 0), 2, 4,
                     angle=0, linewidth=2, fill=False, zorder=2)



tx=np.load('test.npy')

G=6.674e-11
Mp=870.3/G*1e9
Ms=101.4/G*1e9


tlist=tx[:,0]
xlist=tx[:,3]

def data_gen(t=0):
    cnt = 0
    while cnt < len(tlist)-1:
        cnt += 1
        #t += 0.1
        #yield t, np.sin(2*np.pi*t) * np.exp(-t/10.)
        t=tlist[cnt]
        n=(G*(Mp+Ms)/xlist[cnt]**3)**0.5/1000000
        x=xlist[cnt]*np.cos(n*t)
        y=xlist[cnt]*np.sin(n*t)
        yield x,y


def init():
    ax.set_ylim(-30, 30)
    ax.set_xlim(-30, 30)
    ax.add_artist(pluto)
    
    del xdata[:]
    del ydata[:]
    line.set_data(xdata, ydata)
    return line,

fig, ax = plt.subplots()
line, = ax.plot([], [], 'ro')
ax.grid()
xdata, ydata = [], []


def run(data):
    line.set_data(data)
    xd,yd=data
    pluto_ellipse.set_angle=xd
    return pluto_ellipse,line,

ani = animation.FuncAnimation(fig, run, data_gen, blit=False, interval=0.1,
                              repeat=False, init_func=init)
plt.show()
ani.save('basic_animation.mp4', fps=60)



