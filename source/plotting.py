import matplotlib.pyplot as plt 
import matplotlib.pylab as pl
import numpy as np 
from constitutive import sat,temp 
from params import zf

def plot(ws,phi,N,Nr,z,timesteps):
    colors = pl.cm.plasma_r(timesteps/timesteps.max())
    plt.figure(figsize=(8,6))
  
    plt.subplot(131)
    plt.axhline(zf,color='k',linestyle='--')
    for i in range(timesteps.size):
        plt.plot(ws[i,:],z[i,:],color=colors[i])  
    plt.xlabel(r'$v_\mathrm{s}^*$',fontsize=16)
    plt.ylabel(r'$z$',fontsize=16)

    plt.subplot(132)
    plt.axhline(zf,color='k',linestyle='--')
    plt.axvline(0,color='indigo',linestyle='--')
    for i in range(timesteps.size):
        plt.plot(N[i,:],z[i,:],color=colors[i])
        # plt.plot(Nr[i,:],z[i,:],color='k')  
    plt.xlabel(r'$N$',fontsize=16)
    plt.gca().set_yticks([])


    plt.subplot(133)
    plt.axhline(zf,color='k',linestyle='--')
    for i in range(timesteps.size):
        plt.plot(phi[i,:],z[i,:],color=colors[i])
 
    plt.xlim(-0.1,1.1)
    plt.xlabel(r'$\phi$',fontsize=16)
    plt.gca().yaxis.set_label_position("right")
    plt.gca().yaxis.tick_right()
    plt.ylabel(r'$z$',fontsize=16)
    plt.show()
    plt.close()



def plot_movie(ws,phi,N,Nr,z,timesteps):    
    ind = np.arange(0,timesteps.size,2)
    for i in range(ind.size):
        plt.figure(figsize=(8,6))
        plt.suptitle(r'$t/t_\mathrm{max}=$'+'{:.2f}'.format(timesteps[ind[i]]/timesteps.max()),fontsize=24)
        plt.subplot(131)
        plt.plot(ws[ind[i],:],z[ind[i],:],color='indigo')  
        plt.xlabel(r'$v_\mathrm{s}^*$',fontsize=16)
        plt.ylabel(r'$z$',fontsize=16)
        plt.xlim(ws.min()-0.1,ws.max()+0.1)
        plt.ylim(z[ind,:].min(),z[ind,:].max())
        plt.axhline(y=zf,color='k',linewidth=1,linestyle='-')

        plt.subplot(132)
        plt.plot(N[ind[i],:],z[ind[i],:],color='indigo')  
        plt.xlabel(r'$N$',fontsize=16)
        plt.gca().set_yticks([])
        plt.xlim(-0.1,+N.max())
        plt.ylim(z[ind,:].min(),z[ind,:].max())
        plt.axvline(x=0,color='indigo',linewidth=1,linestyle='--')
        plt.axhline(y=zf,color='k',linewidth=1,linestyle='-')


        plt.subplot(133)
        plt.plot(phi[ind[i],:],z[ind[i],:],color='indigo')
        plt.axhline(y=zf,color='k',linewidth=1,linestyle='-')
        plt.xlim(0,1)
        plt.ylim(z[ind,:].min(),z[ind,:].max())
        plt.xlabel(r'$\phi$',fontsize=16)
        plt.gca().yaxis.set_label_position("right")
        plt.gca().yaxis.tick_right()
        plt.ylabel(r'$z$',fontsize=16)
        plt.savefig('./movie/'+str(i)+'.png')
        plt.close()



# t = np.outer(np.linspace(0,t_f,nt),np.ones(nz+1))
# ind = np.arange(0,nt,1)
# colors = pl.cm.plasma_r(np.linspace(0,1,ind.size))
# plt.figure(figsize=(14,4))
# plt.suptitle(r'time dependent solutions',fontsize=22)
# plt.subplot(151)
# plt.contourf(t,z,ws/np.max(np.abs(ws)),cmap='RdBu_r',levels=np.linspace(-1,1,99),extend='both')
# plt.xlabel(r'$t$',fontsize=16)
# # plt.xlabel(r'$w_\mathrm{s}$',fontsize=16)
# plt.ylabel(r'$z$',fontsize=16)
# cbar = plt.colorbar(orientation='horizontal',pad=0.2,ticks=np.linspace(-1,1,5))
# cbar.set_label(r'$w_\mathrm{s}$',fontsize=20)
# cbar.ax.tick_params(labelsize=12)

# plt.subplot(152)
# plt.contourf(t,z,wi/np.max(np.abs(wi)),cmap='RdBu_r',levels=np.linspace(-1,1,99),extend='both')    
# plt.xlabel(r'$t$',fontsize=16)
# plt.axhline(zf,linestyle='--',color='k',linewidth=1)
# plt.gca().set_yticks([])
# cbar = plt.colorbar(orientation='horizontal',pad=0.2,ticks=np.linspace(-1,1,5))
# cbar.set_label(r'$w_\mathrm{i}$',fontsize=20)
# cbar.ax.tick_params(labelsize=12)

# plt.subplot(153)
# S = sat(temp(z),phi)
# q = (1-phi*S)*(wi-ws)    
# plt.contourf(t,z,q/np.max(np.abs(q)),cmap='RdBu_r',levels=np.linspace(-1,1,99),extend='both')        
# plt.xlabel(r'$t$',fontsize=16)
# plt.axhline(zf,linestyle='--',color='k',linewidth=1)
# plt.gca().set_yticks([])
# cbar = plt.colorbar(orientation='horizontal',pad=0.2,ticks=np.linspace(-1,1,5))
# cbar.set_label(r'$q$',fontsize=20)
# cbar.ax.tick_params(labelsize=12)

# plt.subplot(154)
# plt.contourf(t,z,N_/np.max(np.abs(N_)),cmap='RdBu_r',levels=np.linspace(-1,1,99),extend='both')        
# plt.xlabel(r'$t$',fontsize=16)
# plt.axhline(zf,linestyle='--',color='k',linewidth=1)
# plt.gca().set_yticks([])
# cbar = plt.colorbar(orientation='horizontal',pad=0.2,ticks=np.linspace(-1,1,5))
# cbar.set_label(r'$N$',fontsize=20)
# cbar.ax.tick_params(labelsize=12)

# plt.subplot(155)
# plt.contourf(t,z,phi,cmap='Blues_r',levels=np.linspace(0,1,99),extend='both')
# plt.xlabel(r'$t$',fontsize=16)
# plt.gca().yaxis.set_label_position("right")
# plt.gca().yaxis.tick_right()
# plt.ylabel(r'$z$',fontsize=16)
# cbar = plt.colorbar(orientation='horizontal',pad=0.2,ticks=np.linspace(0,1,5))
# cbar.set_label(r'$\phi$',fontsize=20)
# cbar.ax.tick_params(labelsize=12)
# plt.show()
# plt.close()