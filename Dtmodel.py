import numpy as np
from scipy.integrate import ode
import matplotlib.pyplot as plt

#Parameter Values
x0= (5.5,2,4,0)
G=6.674e-11
n=2*np.pi/(6.3867*24*60*60)
Mp=870.3/G*1e9
Ms=101.4/G*1e9
mu=(Mp*Ms)/(Mp+Ms)
k2p=0.058
k2s=0.006
Rp=1153e3
Rs=606.0e3
Dtp=600
Dts=600
Cp=0.383*Mp*Rp**2
Cs=0.4*Ms*Rs**2
ADt=10
AQ=1.15
Qp=100
Qs=Qp/AQ*(k2s/k2p)*(Mp/Ms)**2*(Rs/Rp)**5
t_end = 1e7
t_start = 1e-2
t_step_list =[10**(i/100000) for i in range(900001)]
#t_interval = np.arange(t_start, t_end, t_step)



def eq_system_e(t,x):
        

    Eqs= np.zeros((4))
    Eqs[0]= (-3*(n**2*(x[2]*Rp)**3/(Mp+Ms))/(Cp*x[2]**6)*k2p*Dtp*Ms**2*(Rp**-1)*((1+15/2*x[3]**2)*x[0]-(1+27/2*x[3]**2))
            +3/2*x[0]*(6*(n**2*(x[2]*Rp)**3/(Mp+Ms))/(mu*x[2]**8))*k2p*Dtp*Ms**2*(Rp**-3)*((1+27/2*x[3]**2)*(x[0]+ADt*x[1])-(1+23*x[3]**2)*(1+ADt)))*31536000
    Eqs[1]= (-3*(n**2*(x[2]*Rp)**3/(Mp+Ms))/(Cs*x[2]**6)*k2s*Dts*Mp**2*(Rp**5/Rs**6)*((1+15/2*x[3]**2)*x[1]-(1+27/2*x[3]**2))
            +3/2*x[0]*(6*(n**2*(x[2]*Rp)**3/(Mp+Ms))/(mu*x[2]**8))*k2p*Dtp*Ms**2*(Rp**-3)*((1+27/2*x[3]**2)*(x[0]+ADt*x[1])-(1+23*x[3]**2)*(1+ADt)))*31536000
    Eqs[2]= x[2]*(6*(n**2*(x[2]*Rp)**3/(Mp+Ms))/(mu*x[2]**8)*k2p*Dtp*Ms**2*(Rp**-3)*((1+27/2*x[3]**2)*(x[0]+ADt*x[1])-(1+23*x[3]**2)*(1+ADt)))*31536000
    Eqs[3]= (x[3]*(27*(n**2*(x[2]*Rp)**3/(Mp+Ms))/(mu*x[2]**8))*k2p*Dtp*Ms**2*(Rp**-3)*((11/18*(x[0]+ADt*x[1])-(1+ADt))))*31536000   
    
    

    
    return Eqs



ts = [t_start]
ys = [x0]

# BDF method suited to stiff systems of ODEs
#ode.set_integrator('vode',nsteps=500,method='bdf')


r_e =  ode(eq_system_e).set_integrator('lsoda', method='bdf')
r_e.set_initial_value(x0,t_start)


i=0
while ts[-1] < t_end:
#    t_step=t_step_list[i]
#    i+=1
    if ts[-1]<201:
        t_step=0.1
    elif ts[-1]<t_end-3000:
        t_step=1000
    else:    
        t_step=0.1       
#    t_step=100
    ts.append(r_e.t+t_step)
    ys.append(r_e.integrate(r_e.t+t_step))    



t = np.vstack(ts)
xs= np.vstack(ys)
tx=np.hstack((t,xs))
np.save('test',tx)



fig=plt.figure()
plt.subplot(2,2,3)
plt.plot(t,xs[:,0],label='x0')   
plt.xscale('log')
plt.subplot(2,2,4)
plt.plot(t,xs[:,1],label='x1')
plt.xscale('log')
plt.subplot(2,2,1)
plt.plot(t,xs[:,2],label='x2')
plt.xscale('log')
plt.subplot(2,2,2)
plt.plot(t,xs[:,3],label='x3')
plt.xscale('log')


