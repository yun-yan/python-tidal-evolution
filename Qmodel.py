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
t_end = 1e9
t_start = 1e-2
t_step = 1000
t_interval = np.arange(t_start, t_end, t_step)



def eq_system_e(t,x,Dp,Ep,Fp,Ds,Es,Fs,Sp,Ss):
        

    Eqs= np.zeros((4))

    Eqs[0]= (-3*(G)/((G*(Mp+Ms)/(x[2]*Rp)**3)**(0.5))/(2*Cp*x[2]**6)*k2p/Qp*Ms**2*(Rp**-1)*(Sp+x[3]**2*Dp)
            +3/2*x[0]*(3*(G*(Mp+Ms)/(x[2]*Rp)**3)**(0.5)*k2p/Qp*(Ms/Mp)/x[2]**5*(Sp+AQ*Ss+x[3]**2*(Ep+AQ*Es))))*31536000
    Eqs[1]=(-3*(G)/((G*(Mp+Ms)/(x[2]*Rp)**3)**(0.5))/(2*Cs*x[2]**6)*k2s/Qs*Mp**2*(Rs**5/Rp**6)*(Ss+x[3]**2*Ds)
            +3/2*x[0]*(3*(G*(Mp+Ms)/(x[2]*Rp)**3)**(0.5)*k2p/Qp*(Ms/Mp)/x[2]**5*(Sp+AQ*Ss+x[3]**2*(Ep+AQ*Es))))*31536000
    Eqs[2]= x[2]*(3*(G*(Mp+Ms)/(x[2]*Rp)**3)**(0.5)*k2p/Qp*(Ms/Mp)/x[2]**5*(Sp+AQ*Ss+x[3]**2*(Ep+AQ*Es)))*31536000
    Eqs[3]= x[3]*((G*(Mp+Ms)/(x[2]*Rp)**3)**(0.5)*k2p/Qp*(Ms/Mp)/x[2]**5*(Fp+AQ*Fs))*31536000


    
    return Eqs



ts = [t_start]
ys = [x0]

# BDF method suited to stiff systems of ODEs
#ode.set_integrator('vode',nsteps=500,method='bdf')


r_e =  ode(eq_system_e).set_integrator('lsoda', method='bdf')
r_e.set_initial_value(x0,t_start).set_f_params(15/2,51/4,57/8)


'''

while ts[-1] < t_end and ys[-1][0]>1:
    ts.append(r.t+t_step)
    ys.append(r.integrate(r.t+t_step))



'''
error=10**-4

while ts[-1] < t_end: 
    if ys[-1][0]>3/2: 
        if ys[-1][1]>3/2:   r_e.set_f_params(15/2,51/4,57/8,15/2,51/4,57/8,1,1)         
        elif ys[-1][1]==3/2:   r_e.set_f_params(15/2,51/4,57/8,-19/4,-45/8,-33/16,1,1)         
        elif ys[-1][1]<3/2 and ys[-1][1]>1+error: r_e.set_f_params(15/2,51/4,57/8,-17,-24,-45/4,1,1)
        elif abs(ys[-1][1]-1)<error: r_e.set_f_params(15/2,51/4,57/8,-12,-19,-21/2,1,0)
        else: r_e.set_f_params(15/2,51/4,57/8,-7,-14,-39/4,1,-1)
    elif ys[-1][0]==3/2: 
        if ys[-1][1]>3/2:   r_e.set_f_params(-19/4,-45/8,-33/16,15/2,51/4,57/8,1,1)         
        elif ys[-1][1]==3/2:   r_e.set_f_params(-19/4,-45/8,-33/16,-19/4,-45/8,-33/16,1,1)         
        elif ys[-1][1]<3/2 and ys[-1][1]>1+error: r_e.set_f_params(-19/4,-45/8,-33/16,-17,-24,-45/4,1,1)
        elif abs(ys[-1][1]-1)<error: r_e.set_f_params(-19/4,-45/8,-33/16,-12,-19,-21/2,1,0)
        else: r_e.set_f_params(-19/4,-45/8,-33/16,-7,-14,-39/4,1,-1)
    elif ys[-1][0]<3/2 and ys[-1][0]>1+error: 
        if ys[-1][1]>3/2:   r_e.set_f_params(-17,-24,-45/4,15/2,51/4,57/8,1,1)         
        elif ys[-1][1]==3/2:   r_e.set_f_params(-17,-24,-45/4,-19/4,-45/8,-33/16,1,1)         
        elif ys[-1][1]<3/2 and ys[-1][1]>1+error: r_e.set_f_params(-17,-24,-45/4,-17,-24,-45/4,1,1)
        elif abs(ys[-1][1]-1)<error: r_e.set_f_params(-17,-24,-45/4,-12,-19,-21/2,1,0)
        else: r_e.set_f_params(-17,-24,-45/4,-7,-14,-39/4,1,0)
    elif abs(ys[-1][0]-1)<=error: 
        if ys[-1][1]>3/2:   r_e.set_f_params(-12,-19,-21/2,15/2,51/4,57/8,0,1)         
        elif ys[-1][1]==3/2:   r_e.set_f_params(-12,-19,-21/2,-19/4,-45/8,-33/16,0,1)         
        elif ys[-1][1]<3/2 and ys[-1][1]>1+error: r_e.set_f_params(-12,-19,-21/2,-17,-24,-45/4,0,1)
        elif abs(ys[-1][1]-1)<error: r_e.set_f_params(-12,-19,-21/2,-12,-19,-21/2,0,0)
        else: r_e.set_f_params(-12,-19,-21/2,-7,-14,-39/4,0,0)
    else: 
        if ys[-1][1]>3/2:   r_e.set_f_params(-7,-14,-39/4,15/2,51/4,57/8,-1,1)         
        elif ys[-1][1]==3/2:   r_e.set_f_params(-7,-14,-39/4,-19/4,-45/8,-33/16,-1,1)         
        elif ys[-1][1]<3/2 and ys[-1][1]>1+error: r_e.set_f_params(-7,-14,-39/4,-17,-24,-45/4,-1,1)
        elif abs(ys[-1][1]-1)<error: r_e.set_f_params(-7,-14,-39/4,-12,-19,-21/2,-1,0)
        else: r_e.set_f_params(-7,-14,-39/4,-7,-14,-39/4,-1,0)
    if ts[-1]<10000:
        t_step=1
    else:
        t_step=100
    ts.append(r_e.t+t_step)
    ys.append(r_e.integrate(r_e.t+t_step))     
#    else:
#        ts.append(ts[-1]+t_step)
#        ys.append(ys[-1]) 
        
        
    





t = np.vstack(ts)
xs= np.vstack(ys)




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


