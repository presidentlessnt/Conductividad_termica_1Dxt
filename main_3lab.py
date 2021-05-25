#! /usr/bin/python3
import sys
import math
import numpy as np
import scipy.stats
import scipy.special as sp
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import warnings
warnings.filterwarnings("ignore")
import modulo as mtd

#def f1(x):
#    return x

def f1(x):
    return -.0005*x**3+.0068*x**2+18.366*x+1.968

def limp(A,col,t):
    l=[]
    for i in range(len(A)):
        if A[i,col]==0 and i==0:
            if t==0:
                l.append(int(A[i,col]))
            if t==1:
                l.append(A[i,col])
        elif A[i,col]!=0:
            if t==0:
                l.append(int(A[i,col]))
            if t==1:
                l.append(A[i,col])
    return l


#obtener picos
def findp(mM,radio,tt,fun): # M:data, radio, tt:minimo deltaT entre picos subsiguientes
    M=np.copy(mM)
    M[:,1]=fun(M[:,1])
    M[:,2]=fun(M[:,2])
    m=len(M[0])
    A=np.zeros((20,m-1))  # valor T de pico
    B=np.zeros((20,m-1))  # ordinal (posicion en matriz) de voltaje pico
    for j in range(1,m):
        cont=0
        val=0
        for i in range(radio,len(M)-radio):
            if M[i,j]==max(M[i-radio:i+radio,j]):
                if cont==0: # busca primer pico
                    A[cont,j-1]=M[i,0]
                    B[cont,j-1]=i
                    val=M[i,j]
                    cont+=1
                elif cont>0 and abs(M[i,0]-A[cont-1,j-1])>tt: # busca demas picos
                    A[cont,j-1]=M[i,0]
                    B[cont,j-1]=i
                    val=M[i,j]
                    cont+=1    
        for i in range(radio):
            if M[i,j]==max(M[:i+radio,j]) and abs(M[i,0]-A[0,j-1])>tt:
                    for k in range(cont,0,-1):
                        A[k,j-1]=A[k-1,j-1]
                        B[k,j-1]=B[k-1,j-1]
                    A[0,j-1]=M[i,0]
                    B[0,j-1]=i
    lA1=limp(A,0,1)
    lA2=limp(A,1,1)
    lB1=limp(B,0,0)
    lB2=limp(B,1,0)
    DT=np.array(lA2)-np.array(lA1)
#    print(DT)
#    print(np.mean(DT))
#    print(np.std(DT),'\n')
    P1=np.zeros((len(lB1),2))
    P2=np.zeros((len(lB2),2))
    for i in range(len(P1)):
        P1[i,0]=lA1[i]
        P1[i,1]=M[lB1[i],1]
    for i in range(len(P2)):
        P2[i,0]=lA2[i]
        P2[i,1]=M[lB2[i],2]
#    print(P2[:,0]-P1[:,0])
#    print((P2[:,1]+P1[:,1])*.5,'\n')
    return M,P1,P2


def deltapm(M,fun):
    Mm=np.copy(M)
    Mm[:,1]=M[:,1]*.8/100+.03
    Mm[:,2]=M[:,2]*.8/100+.03
    Mm[:,1]=fun(Mm[:,1])
    Mm[:,2]=fun(Mm[:,2])
    return Mm


## Obtener datos del .txt ##
D1=np.loadtxt('D1.txt')
D2=np.loadtxt('D2.txt')
D3=np.loadtxt('D3.txt')
D4=np.loadtxt('D4.txt')
D5=np.loadtxt('D5.txt')
D6=np.loadtxt('D6.txt')


M1=findp(D1,4,1,f1)
M2=findp(D2,4,2,f1)
M3=findp(D3,4,2,f1)
M4=findp(D4,5,4,f1)
M5=findp(D5,4,2,f1)
M6=findp(D6,3,4,f1)


K1=deltapm(D1,f1)
K2=deltapm(D2,f1)
K3=deltapm(D3,f1)
K4=deltapm(D4,f1)
K5=deltapm(D5,f1)
K6=deltapm(D6,f1)



## Graficas ##

def custom_plot(M1,K1,i,p,mm,ll,ax=None,**plt_kwargs):
    if ax is None:
        ax=plt.figure(i)
    plt.plot(M1[0][:,0],M1[0][:,1],'o-',label=r'$A$',markersize=2,lw=.5)
    plt.plot(M1[0][:,0],M1[0][:,2],'o-',label=r'$B$',markersize=2,lw=.3)
    plt.plot(M1[1][:,0],M1[1][:,1],'ks',markersize=5,label='Picos')
    plt.plot(M1[2][:,0],M1[2][:,1],'ks',markersize=5)
    if p==1:
        plt.plot(M1[0][:,0],M1[0][:,1]+K1[:,1],'--',lw=.5,color='C0')
        plt.plot(M1[0][:,0],M1[0][:,1]-K1[:,1],'--',lw=.5,color='C0')
        plt.plot(M1[0][:,0],M1[0][:,2]+K1[:,2],'--',lw=.5,color='C1')
        plt.plot(M1[0][:,0],M1[0][:,2]-K1[:,2],'--',lw=.5,color='C1')
        lines=plt.gca().get_lines()
        legend1=plt.legend([lines[i] for i in [5,7]],[
'A - $\Delta T_{min}:\pm$%1.2f°C, $\Delta T_{max}:\pm$%1.2f°C'%(min(K1[:,1]),max(K1[:,1])),
'B - $\Delta T_{min}:\pm$%1.2f°C, $\Delta T_{max}:\pm$%1.2f°C'%(min(K1[:,2]),max(K1[:,2]))],loc=ll)
    if p==0:
        lines=plt.gca().get_lines()
        legend1=plt.legend([lines[i] for i in [0,1]],[
'A - $\Delta T_{min}:\pm$%1.2f°C, $\Delta T_{max}:\pm$%1.2f°C'%(min(K1[:,1]),max(K1[:,1])), 
'B - $\Delta T_{min}:\pm$%1.2f°C, $\Delta T_{max}:\pm$%1.2f°C'%(min(K1[:,2]),max(K1[:,2]))],loc=ll)
    plt.xlabel('$t\;[min]$')
    plt.ylabel('$T\;[°C]$')
    plt.ylim(min(M1[0][:,2])-15,max(M1[0][:,1])+15)
    plt.grid(ls=':',color='grey',alpha=.5)
    plt.legend(loc=mm)
    plt.gca().add_artist(legend1)
    return ax


custom_plot(M1,K1,1,1,1,7)
custom_plot(M2,K2,2,1,1,7)
custom_plot(M3,K3,3,1,1,4)
custom_plot(M4,K4,4,1,1,4)
custom_plot(M5,K5,5,1,1,4)
custom_plot(M6,K6,6,1,1,4)


nn=4
X=np.linspace(0,20,nn*20)
T=np.linspace(0,40,nn*40)
X,T=np.meshgrid(X,T)
a,b,c,d=[25,-.5,120,1]
w,D=[2*np.pi/10,12.73]
Z=a*np.exp(-np.sqrt(w/2/D)*X)*np.cos(w*T-np.sqrt(w/2/D)*X)+c+b*X

#def teo(x):
#    pp=10
#    return 25*np.exp(-np.sqrt(w/2/D)*pp)*np.cos(w*x-np.sqrt(w/2/D)*pp)*np.exp(-.06*x)+200


# Grafica 3D
fig1=plt.figure()
axs=fig1.gca(projection='3d')
surf1=axs.plot_surface(X,T,Z,cmap=cm.viridis,linewidth=0,antialiased=False,alpha=.5)
axs.contour(X,T,Z,zdir='x',levels=[9],lw=1,colors='red',ls='solid',offset=22)
axs.contour(X,T,Z,zdir='x',levels=[14],lw=1,ls='solid',offset=22)
axs.contour(X,T,Z,zdir='x',levels=[9],linewidths=3,colors='red',linestyles='dashed')
axs.contour(X,T,Z,zdir='x',levels=[14],linewidths=3,linestyles='dashed')
cbar=fig1.colorbar(surf1,shrink=.5,aspect=5)
cbar.ax.set_title('$T\;[°C]$')
axs.set_xlabel(r'$x$')
axs.set_ylabel(r'$t$')
axs.set_zlabel(r'$T[°C]$')


plt.show()



