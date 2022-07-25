import numpy as np
from numpy import fft as fft

## Numerical differentiation
def NDiffFd1(X,Y):
    Dx  = X[1]-X[0]
    Yd  = (np.roll(Y,-1)-Y)/Dx
    return Yd

def NDiffFd2(X,Y):
    Dx  = X[1]-X[0]
    Yd  = (np.roll(Y,-1)-np.roll(Y,1))/(2*Dx)
    return Yd

def NDiffFd4(X,Y):
    Dx  = X[1]-X[0]
    Yd  = (np.roll(Y,-1)-np.roll(Y,1))*(2/3/Dx)-(np.roll(Y,-2)-np.roll(Y,2))*(1/12/Dx)
    return Yd

def NDiffFd6(X,Y):
    Dx  = X[1]-X[0]
    Yd  = (np.roll(Y,-1)-np.roll(Y,1))*(3/4/Dx)-(np.roll(Y,-2)-np.roll(Y,2))*(3/20/Dx)+(np.roll(Y,-3)-np.roll(Y,3))*(1/60/Dx)
    return Yd

def NDiffFft(X,Y,od=1):
    N = len(X)
    L = X[-1]-X[0]+X[1]-X[0]
    if np.mod(len(Y),2) == 0:
        omega = fft.fftshift(1j*(2*np.pi/L)*np.arange(-(N/2),(N/2),1))
    else:
        omega = fft.fftshift(1j*(2*np.pi/L)*np.arange(-np.floor(N/2),np.floor(N/2)+1,1))
    YHat  = fft.fft(Y)
    YdHat = omega**od*YHat
    Yd    = np.real(fft.ifft(YdHat))
    return Yd

## Numerical integration
def NIntgRk4(t,Q,Dt,f):
    k1 = f(t,Q)
    k2 = f(t+Dt/2,Q+(1/2)*k1*Dt)
    k3 = f(t+Dt/2,Q+(1/2)*k2*Dt)
    k4 = f(t+Dt,Q+k3*Dt)
    Q_next = Q+(k1+2*k2+2*k4+k4)/6*Dt
    return Q_next

def NIntgFd1(t,Q,Dt,f):
    k1 = f(t,Q)
    Q_next = Q+k1*Dt
    return Q_next

# Error Calculation
def rms(A):
    RMS = np.sqrt(np.sum(A**2)/np.prod((A).shape))
    return RMS

def CalErrA(F,G):
    Err = (np.mean(F)-np.mean(G))**2+(np.std(F)-np.std(G))**2;
    return Err

def CalErrF(F,G):
    Err = 2*(1-np.corrcoef(F,G)[1])*np.std(F)*np.std(G);
    return Err