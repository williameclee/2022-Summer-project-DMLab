import numpy as np
from numpy import fft as fft

def NDiffFd1(X,Y):
    Dx  = X[1]-X[0]
    Yd  = (np.roll(Y,1)-Y)/Dx
    return Yd

def NDiffFd2(X,Y):
    Dx  = X[1]-X[0]
    Yd  = (np.roll(Y,1)-np.roll(Y,-1))/(2*Dx)
    return Yd

def NDiffFd4(X,Y):
    Dx  = X[1]-X[0]
    Yd  = (np.roll(Y,1)-np.roll(Y,-1))*(2/3/Dx)-(np.roll(Y,2)-np.roll(Y,-2))*(1/12/Dx)
    return Yd

def NDiffFd6(X,Y):
    Dx  = X[1]-X[0]
    Yd  = (np.roll(Y,1)-np.roll(Y,-1))*(3/4/Dx)-(np.roll(Y,2)-np.roll(Y,-2))*(3/20/Dx)+(np.roll(Y,3)-np.roll(Y,-3))*(1/60/Dx)
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

def NIntgRk4(t,Q,dT,f):
    k1 = f(t,Q)
    k2 = f(t+dT/2,Q+(1/2)*k1*dT)
    k3 = f(t+dT/2,Q+(1/2)*k2*dT)
    k4 = f(t+dT,Q+k3*dT)
    Q_next = Q+(k1+2*k2+2*k4+k4)/6*dT
    return Q_next
