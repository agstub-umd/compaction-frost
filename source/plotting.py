import matplotlib.pyplot as plt 
import matplotlib.pylab as pl
import matplotlib as mpl
import numpy as np 
from constitutive import get_fields, Phi 
from params import z_f

def plot(N,z,timesteps):
    plt.figure(figsize=(12,6))
    phi = Phi(N,log=np.log)
    colors = pl.cm.plasma_r(timesteps/timesteps.max())
    T,S,k = get_fields(z)

    j = np.argmax(z[:,-1])

    N_f = N[0,0]

    plt.subplot(121)
    for i in range(timesteps.size):
        plt.plot(phi[i,:]/phi[0,0].mean(),z[i,:]-z_f,color=colors[i],linewidth=1)
    plt.axhline(0,linestyle='--',color='k',linewidth=1)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.ylabel(r'$z-z_\mathrm{f}$',fontsize=16)
    plt.xlabel(r'$\phi\,/\,\phi_\mathrm{f}$',fontsize=20)

    plt.subplot(122)
    for i in range(timesteps.size):
        plt.plot(N[i,:]/N_f,z[i,:]-z_f,color=colors[i])
    plt.axhline(0,linestyle='--',color='k',linewidth=1)
    plt.axvline(0,linestyle=':',color='k',linewidth=1)
    plt.xticks(fontsize=16)
    plt.gca().yaxis.set_label_position("right")
    plt.gca().yaxis.tick_right()
    plt.yticks(fontsize=16)
    plt.ylabel(r'$z-z_\mathrm{f}$',fontsize=16)
    plt.xlabel(r'$N\,/\,N_\mathrm{f}$',fontsize=20)
    plt.xlim(-1,N.max()/N_f+1)
 

    plt.show()
    plt.close()


def plot_steady(N,z):
    T,S,k = get_fields(z)
 
    N_l = ((1-Phi(N,log=np.log)*S)*(1+T))[-1]
    phi_l = Phi(N_l,log=np.log)

    N_f = N[np.argmin(np.abs(z-z_f))]

    z_n = z[np.argmin(N)]-z_f
    phi = Phi(N,log=np.log)
 
    plt.figure(figsize=(8,6))
    plt.subplot(121)
    plt.plot(phi/phi[0],z-z_f,color='indigo',label=r'$\phi$')
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.ylabel(r'$z-z_\mathrm{f}$',fontsize=16)
    plt.xlabel(r'$\phi\,/\,\phi_\mathrm{f}$',fontsize=20)
    plt.legend(fontsize=12,loc='lower right')
 
    plt.subplot(122)
    plt.plot(N/N_f,z-z_f,color='indigo',label=r'$N$')
    plt.axhline(z_n,linestyle='--',color='springgreen',linewidth=2,label=r'$z_n$')
    plt.axvline(0,linestyle=':',color='k',linewidth=1)
    plt.xlim(-1,N.max()/N_f+1)
    plt.xticks(fontsize=16)
    plt.gca().yaxis.set_label_position("right")
    plt.gca().yaxis.tick_right()
    plt.yticks(fontsize=16)
    plt.ylabel(r'$z-z_\mathrm{f}$',fontsize=16)
    plt.xlabel(r'$N\,/\,N_\mathrm{f}$',fontsize=20)
    plt.legend(fontsize=12,loc='lower right')
    plt.show()
    plt.close()

