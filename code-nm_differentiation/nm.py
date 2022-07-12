import numpy as np
from numpy import fft as fft

def NDiffFd1(X,Y):
    Dx  = X[1]-X[0]
    Yp1 = np.append(Y[1:],Y[0])
    Yd  = (Yp1-Y)/Dx
    return Yd

def NDiffFd2(X,Y):
    Dx  = X[1]-X[0]
    Yp1 = np.append(Y[1:],Y[0])
    Yn1 = np.append(Y[-1],Y[:-1])
    Yd  = (Yp1-Yn1)/(2*Dx)
    return Yd

def NDiffFd4(X,Y):
    Dx  = X[1]-X[0]
    Yp1 = np.append(Y[1:],Y[0])
    Yn1 = np.append(Y[-1],Y[:-1])
    Yp2 = np.append(Y[2:],Y[0:2])
    Yn2 = np.append(Y[-2:],Y[:-2])
    Yd  = (Yp1-Yn1)*(2/3/Dx)-(Yp2-Yn2)*(1/12/Dx)
    return Yd

def NDiffFd6(X,Y):
    Dx  = X[1]-X[0]
    Yp1 = np.append(Y[1:],Y[0])
    Yn1 = np.append(Y[-1],Y[:-1])
    Yp2 = np.append(Y[2:],Y[0:2])
    Yn2 = np.append(Y[-2:],Y[:-2])
    Yp3 = np.append(Y[3:],Y[0:3])
    Yn3 = np.append(Y[-3:],Y[:-3])
    Yd  = (Yp1-Yn1)*(3/4/Dx)-(Yp2-Yn2)*(3/20/Dx)+(Yp3-Yn3)*(1/60/Dx)
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