import matplotlib.pyplot as plt
import numpy as np
def ode_Euler(f,t0,tf,y0=0,n=100):
    """
    First order ODE (y' = f(t,y)) Solver using Euler method
     
    t0: initial value of independent variable
    tf: final value of independent variable
    y0: initial value of dependent variable
    n : number of steps
    f : function of f(t,y)
 
    The return is list of y(t)
    """
 
    t = np.linspace(t0,tf,n)
    y = list([y0])
    for i in range(n-1):
        h = t[i+1]-t[i]
        y.append(y[-1]+h*f(t[i],y[-1]))
 
    y = np.array(y)
     
    return y

def ode_RK4(f,t0,tf,y0,n):
    """
    First order ODE (y' = f(t,y)) Solver using RK4 method
     
    t0: initial value of independent variable
    tf: final value of independent variable
    y0: initial value of dependent variable
    n : number of steps
    f : function of f(t,y)
 
    The return is list of y(t)
    """
 
    t = np.linspace(t0,tf,n)
    y = list([y0])
    for i in range(n-1):
        h = t[i+1]-t[i]
        k1 = h*f(t[i],y[-1])
        k2 = h*f(t[i]+h/2.0,y[-1]+k1/2.0)
        k3 = h*f(t[i]+h/2.0,y[-1]+k2/2.0)
        k4 = h*f(t[i]+h,y[-1]+k3)
        y.append(y[-1]+(k1+2*k2+2*k3+k4)/6)
 
    y = np.array(y)
     
    return y

def fun(t,x):
    
    
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

#    Eqs= np.zeros((4))
#    Eqs[0]= (-3*(n*(x[2]*Rp)**3/(Mp+Ms))/(2*Cp*x[2]**6)*k2p/Q*Ms**2*(Rp**-1)*np.sign(x[0]-1)
#                     +3/2*x[0]*(3*n*k2p/Q*Ms/Mp/x[2]**5*(np.sign(x[0]-1)+AQ*np.sign(x[1]-1))))*31536000
#    Eqs[1]= (-3*(n**2*(x[2]*Rp)**3/(Mp+Ms))/(Cs*x[2]**6)*k2s*Dts*Mp**2*(Rp**5/Rs**6)*((1+15/2*x[3]**2)*x[1]-(1+27/2*x[3]**2))
#            +3/2*x[0]*(6*(n**2*(x[2]*Rp)**3/(Mp+Ms))/(mu*x[2]**8))*k2p*Dtp*Ms**2*(Rp**-3)*((1+27/2*x[3]**2)*(x[0]+ADt*x[1])-(1+23*x[3]**2)*(1+ADt)))*31536000
#    Eqs[2]= x[2]*(3*n*k2p/Q*Ms/Mp/x[2]**5*(np.sign(x[0]-1)+AQ*np.sign(x[1]-1)))*31536000
#    Eqs[3]= 0
    

    Eqs= np.zeros((4))
    Eqs[0]= (-3*(G)/(Cp*x[2]**6)*k2p*Dtp*Ms**2*(Rp**-1)*((1+15/2*x[3]**2)*x[0]-(1+27/2*x[3]**2))
            +3/2*x[0]*(6*(G)/(mu*x[2]**8))*k2p*Dtp*Ms**2*(Rp**-3)*((1+27/2*x[3]**2)*(x[0]+ADt*x[1])-(1+23*x[3]**2)*(1+ADt)))*31536000
    Eqs[1]= (-3*(G)/(Cs*x[2]**6)*k2s*Dts*Mp**2*(Rp**5/Rs**6)*((1+15/2*x[3]**2)*x[1]-(1+27/2*x[3]**2))
            +3/2*x[0]*(6*(G)/(mu*x[2]**8))*k2p*Dtp*Ms**2*(Rp**-3)*((1+27/2*x[3]**2)*(x[0]+ADt*x[1])-(1+23*x[3]**2)*(1+ADt)))*31536000
    Eqs[2]= x[2]*(6*(G)/(mu*x[2]**8)*k2p*Dtp*Ms**2*(Rp**-3)*((1+27/2*x[3]**2)*(x[0]+ADt*x[1])-(1+23*x[3]**2)*(1+ADt)))*31536000
    Eqs[3]= (x[3]*(27*(G)/(mu*x[2]**8))*k2p*Dtp*Ms**2*(Rp**-3)*((11/18*(x[0]+ADt*x[1])-(1+ADt))))*31536000   
     
    
    return  Eqs

y = ode_RK4(fun,1e-2,10e7,[5.5,2,4,0],100000)
#y = ode_Euler(fun,1,10e7,[5.5,2,4,0],10)
#x=np.array(x)

time =np.linspace(1e-2,1e7,100000)

fig=plt.figure()
plt.subplot(2,2,3)
plt.plot(time,y[:,0],label='x0')   
plt.xscale('log')
plt.subplot(2,2,4)
plt.plot(time,y[:,1],label='x1')
plt.xscale('log')
plt.subplot(2,2,1)
plt.plot(time,y[:,2],label='x2')
plt.xscale('log')
plt.subplot(2,2,2)
plt.plot(time,y[:,3],label='x3')
plt.xscale('log')
#xlabel('t')
#ylabel('y')
#legend()

plt.show()

'''


test


def f(t,y):
    return y-t**2+1
y = ode_RK4(f,0,2,0.5,11)

'''