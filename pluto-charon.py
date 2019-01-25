import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt



def timescale(x,t):

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
    AQ=0.65
    Q=100
    
    
    
    return np.array([(-3*(G)/(Cp*x[2]**6)*k2p*Dtp*Ms**2*(Rp**-1)*((1+15/2*x[3]**2)*x[0]-(1+27/2*x[3]**2))
                     +3/2*x[0]*(6*(G)/(mu*x[2]**8))*k2p*Dtp*Ms**2*(Rp**-3)*((1+27/2*x[3]**2)*(x[0]+ADt*x[1])-(1+23*x[3]**2)*(1+ADt)))*31536000,
                     (-3*(G)/(Cs*x[2]**6)*k2s*Dts*Mp**2*(Rs**5/Rp**6)*((1+15/2*x[3]**2)*x[1]-(1+27/2*x[3]**2))
                     +3/2*x[0]*(6*(G)/(mu*x[2]**8))*k2p*Dtp*Ms**2*(Rp**-3)*((1+27/2*x[3]**2)*(x[0]+ADt*x[1])-(1+23*x[3]**2)*(1+ADt)))*31536000,
                     x[2]*(6*(G)/(mu*x[2]**8)*k2p*Dtp*Ms**2*(Rp**-3)*((1+27/2*x[3]**2)*(x[0]+ADt*x[1])-(1+23*x[3]**2)*(1+ADt)))*31536000,
                     (x[3]*(27*(G)/(mu*x[2]**8))*k2p*Dtp*Ms**2*(Rp**-3)*((11/18*(x[0]+ADt*x[1])-(1+ADt))))*31536000
                     ])
    

#    return np.array([(-3*(n*(x[2]*Rp)**3/(Mp+Ms))/(2*Cp*x[2]**6)*k2p/Q*Ms**2*(Rp**-1)*np.sign(x[0]-1)
#                     +3/2*x[0]*(3*n*k2p/Q*Ms/Mp/x[2]**5*(np.sign(x[0]-1)+AQ*np.sign(x[1]-1))))*31536000,
#                     (-3*(n*(x[2]*Rp)**3/(Mp+Ms))/(2*Cs*x[2]**6)*k2s/Q*Mp**2*(Rs**5/Rp**6)*(x[1]-1)
#                     +3/2*x[0]*(3*n*k2p/Q*Ms/Mp/x[2]**5*(np.sign(x[0]-1)+AQ*np.sign(x[1]-1))))*31536000,
#                     x[2]*(3*n*k2p/Q*Ms/Mp/x[2]**5*(np.sign(x[0]-1)+AQ*np.sign(x[1]-1)))*31536000,
#                     0
#                     ])

#    return np.array([(-3*(n*(x[2]*Rp)**3/(Mp+Ms))/(2*Cp*x[2]**6)*k2p/Q*Ms**2*(Rp**-1)*1
#                     +3/2*x[0]*(3*n*k2p/Q*Ms/Mp/x[2]**5*(1+AQ*1)))*31536000,
#                     (-3*(n*(x[2]*Rp)**3/(Mp+Ms))/(2*Cs*x[2]**6)*k2s/Q*Mp**2*(Rp**5/Rs**-6)*(x[1]-1)
#                     +3/2*x[0]*(3*n*k2p/Q*Ms/Mp/x[2]**5*(1+AQ*1)))*31536000,
#                     x[2]*(3*n*k2p/Q*Ms/Mp/x[2]**5*(1+AQ*1))*31536000,
#                     0
#                     ])
#    return np.array([-np.sign(x[0]-1),
#                     0,
#                     0,
#                     0
#                     ])
                     
                     
                     
def eq_system_e(t,x):
     
     
     
    G=6.674e-11
    n=2*np.pi/(6.3867*24*60*60)
    Mp=870.3/G
    Ms=101.4/G
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
    AQ=0.65
    Q=100
    
    x=xinit
    
    if x[0]>3/2:
        Dp=15/2
        Ep=51/4
        Fp=57/8
    elif x[0]==3/2:
        Dp=-19/4
        Ep=-45/8
        Fp=-33/16
    elif x[0]<3/2 and x[0]>1:
        Dp=-17
        Ep=-24
        Fp=-45/4
    elif x[0]==1:
        Dp=-12
        Ep=-19
        Fp=-21/2
    elif x[0]<1 and x[0]>1/2:
        Dp=-7
        Ep=-14
        Fp=-39/4
    else:
        Dp=-7.5
        Ep=-14.25
        Fp=-75/8
    
    if x[1]>3/2:
        Ds=15/2
        Es=51/4
        Fs=57/8
    elif x[1]==3/2:
        Ds=-19/4
        Es=-45/8
        Fs=-33/16
    elif x[1]<3/2 and x[1]>1:
        Ds=-17
        Es=-24
        Fs=-45/4
    elif x[1]==1:
        Ds=-12
        Es=-19
        Fs=-21/2
    elif x[1]<1 and x[1]>1/2:
        Ds=-7
        Es=-14
        Fs=-39/4
    else:
        Ds=-7.5
        Es=-14.25
        Fs=-75/8

    Eqs= np.zeros((4))
    Eqs[0]=(-3*(n*(x[2]*Rp)**3/(Mp+Ms))/(2*Cp*x[2]**6)*k2p/Q*Ms**2*(Rp**-1)*(np.sign(x[0]-1)+x[3]**2*Dp)
            +3/2*x[0]*(3*n*k2p/Q*(Ms/Mp)/x[2]**5*(np.sign(x[0]-1)+AQ*np.sign(x[1]-1)+x[3]**2*(Ep+AQ*Es))))*31536000
#    Eqs[1]=(-3*(n*(x[2]*Rp)^3/(Mp+Ms))/(2*Cs*x[2]^6)*k2s/Q*Mp^2*(Rp^5/Rs^6)*(np.sign(x[1]-1)+x[3]^2*Ds)
#            +3/2*x[0]*(3*n*k2p/Q*(Ms/Mp)/x[2]^5*(np.sign(x[0]-1)+AQ*np.sign(x[1]-1)+x[3]^2*(Ep+AQ*Es))))*31536000
    Eqs[1]= (-3*(n**2*(x[2]*Rp)**3/(Mp+Ms))/(Cs*x[2]**6)*k2s*Dts*Mp**2*(Rp**5/Rs**6)*((1+15/2*x[3]**2)*x[1]-(1+27/2*x[3]**2))
            +3/2*x[0]*(6*(n**2*(x[2]*Rp)**3/(Mp+Ms))/(mu*x[2]**8))*k2p*Dtp*Ms**2*(Rp**-3)*((1+27/2*x[3]**2)*(x[0]+ADt*x[1])-(1+23*x[3]**2)*(1+ADt)))*31536000
    Eqs[2]=x[2]*(3*n*k2p/Q*(Ms/Mp)/x[2]**5*(np.sign(x[0]-1)+AQ*np.sign(x[1]-1)+x[3]**2*(Ep+AQ*Es)))*31536000
    Eqs[3]=x[3]*(n*k2p/Q*(Ms/Mp)/x[2]**5*(Fp+AQ*Fs))*31536000
                     
    
    return Eqs       

time=np.linspace(1e-2,10e10,num=10000000)
xinit=np.array([5.5,2,4,0])
x = odeint(timescale,xinit,time)
#x = odeint(eq_system_e,xinit,time)


fig=plt.figure()
plt.subplot(2,2,3)
plt.plot(time,x[:,0],label='x0')   
plt.xscale('log')
plt.subplot(2,2,4)
plt.plot(time,x[:,1],label='x1')
plt.xscale('log')
plt.subplot(2,2,1)
plt.plot(time,x[:,2],label='x2')
plt.xscale('log')
plt.subplot(2,2,2)
plt.plot(time,x[:,3],label='x3')
plt.xscale('log')
#xlabel('t')
#ylabel('y')
#legend()

plt.show()