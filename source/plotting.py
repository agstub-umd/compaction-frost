import matplotlib.pyplot as plt 
import matplotlib.pylab as pl
import matplotlib as mpl
import numpy as np 
from constitutive import get_fields, Phi

def plot(N,heave,z,timesteps):
    plt.figure(figsize=(12,6))
    phi = Phi(N,log=np.log)

    ind = np.arange(0,timesteps.size,25)

    colors = pl.cm.winter_r(timesteps/timesteps.max())
    T,S,k = get_fields(z)

    j = np.argmax(z[:,-1])

    N_f = N[0,0]

    z_n = z[ind[-1],np.argmin(N[ind[-1],:])]
    z_l = z[ind[-1],-1]


    plt.subplot(132)
    plt.title(r'(b)',fontsize=20,loc='left')


    for i in ind:
        plt.plot(phi[i,:]/phi[0,0].mean(),z[i,:],color=colors[i],linewidth=2)
    # plt.axhline(0,linestyle='--',color='k',linewidth=1)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlabel(r'$\phi\,/\,\phi_\mathrm{f}$',fontsize=20)
    plt.gca().yaxis.set_ticklabels([])
    plt.ylim(0,10)

    plt.subplot(131)
    plt.title(r'(a)',fontsize=20,loc='left')
    plt.ylabel(r'$z$',fontsize=20)

    plt.axhline(z_n,linestyle='--',color='crimson',linewidth=1.5,label=r'$z_\mathrm{n}$')
    plt.axhline(z_l,linestyle='--',color='darkblue',linewidth=1.5,label=r'$z_\ell^\star$')

    for i in ind:
        plt.plot(N[i,:]/N_f,z[i,:],color=colors[i],linewidth=2)


    # plt.axhline(0,linestyle='--',color='k',linewidth=1)
    plt.axvline(0,linestyle=':',color='k',linewidth=1)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)

    plt.xlabel(r'$N\,/\,N_\mathrm{f}$',fontsize=20)
    plt.xlim(-1,N.max()/N_f+1)
    plt.legend(fontsize=20,loc='lower right')
    plt.ylim(0,10)

 
    plt.subplot(133)
    plt.title(r'(c)',fontsize=20,loc='left')
    for i in ind:
        plt.plot(heave[i,:],z[i,:],color=colors[i],linewidth=2)
    # plt.axhline(0,linestyle='--',color='k',linewidth=1)
    plt.axvline(0,linestyle=':',color='k',linewidth=1)
    plt.xticks(fontsize=16)
    plt.gca().yaxis.set_label_position("right")
    plt.gca().yaxis.tick_right()
    plt.yticks(fontsize=16)
    plt.ylabel(r'$z$',fontsize=20)
    plt.xlabel(r'$v_{\ast}\,/\,v_\mathrm{i}$',fontsize=20)
    plt.xlim(0,0.6)
    plt.ylim(0,10)
    plt.savefig('profiles')
    plt.show()
    plt.close()


def plot_steady(N,z):
    T,S,k = get_fields(z)
 
    N_l = ((1-Phi(N,log=np.log)*S)*(1+T))[-1]

    N_f = N[np.argmin(np.abs(z))]

    z_n = z[np.argmin(N)]
    phi = Phi(N,log=np.log)
 
    plt.figure(figsize=(8,6))
    plt.subplot(121)
    plt.title(r'(a)',loc='left',fontsize=20)
    plt.plot(phi/phi[0],z,color='indigo',label=r'$\phi$')
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.ylabel(r'$z$',fontsize=20)
    plt.xlabel(r'$\phi\,/\,\phi_\mathrm{f}$',fontsize=20)
    # plt.legend(fontsize=12,loc='lower right')
 
    plt.subplot(122)
    plt.title(r'(b)',loc='left',fontsize=20)
    plt.plot(N/N_f,z,color='indigo')
    plt.axhline(z_n,linestyle='--',color='indianred',linewidth=2,label=r'$z_n$')
    plt.axvline(0,linestyle=':',color='k',linewidth=1)
    # plt.axvline(N_l/N_f,linestyle=':',color='goldenrod',linewidth=3)
    plt.xlim(-1,N.max()/N_f+1)
    plt.xticks(fontsize=16)
    plt.gca().yaxis.set_label_position("right")
    plt.gca().yaxis.tick_right()
    plt.yticks(fontsize=16)
    plt.ylabel(r'$z$',fontsize=20)
    plt.xlabel(r'$N\,/\,N_\mathrm{f}$',fontsize=20)
    plt.legend(fontsize=12,loc='lower right')
    plt.savefig('initial')
    plt.show()
    plt.close()

