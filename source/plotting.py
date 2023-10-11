import matplotlib.pyplot as plt 
import matplotlib.pylab as pl
import numpy as np 
from constitutive import sat,temp
from params import zf

def plot(ws,wi,phi,z,N_,timesteps):
    ind = np.arange(0,timesteps.size,1)
    colors = pl.cm.plasma_r(timesteps[ind]/timesteps.max())
    plt.figure(figsize=(12,6))
    disp = 0*ind
        
    plt.subplot(151)
    for i in range(ind.size):
        plt.plot(ws[ind[i],:]+disp[i],z[ind[i],:],color=colors[i],alpha=0.5)  
    plt.xlabel(r'$w_\mathrm{s}$',fontsize=16)
    plt.ylabel(r'$z$',fontsize=16)

    plt.subplot(152)
    for i in range(ind.size):
        plt.plot(wi[ind[i],:]+disp[i],z[ind[i],:],color=colors[i],alpha=0.5)  
    plt.xlabel(r'$w_\mathrm{i}$',fontsize=16)
    plt.gca().set_yticks([])

    plt.subplot(153)
    for i in range(ind.size):
        S = sat(temp(z[ind[i],:]))
        q = (1-phi[ind[i],:]*S)*(wi[ind[i],:]-ws[ind[i],:])
        plt.plot(q+disp[i],z[ind[i],:],color=colors[i],alpha=0.5)  
    plt.xlabel(r'$q$',fontsize=16)
    plt.gca().set_yticks([])

    plt.subplot(154)
    for i in range(ind.size):
        plt.plot(N_[ind[i],:]+disp[i],z[ind[i],:],color=colors[i],alpha=0.5)  
    plt.xlabel(r'$N$',fontsize=16)
    plt.gca().set_yticks([])

    plt.subplot(155)
    for i in range(ind.size):
        S = sat(temp(z[ind[i],:]))
        plt.plot(phi[ind[i],:]+disp[i],z[ind[i],:],color=colors[i],alpha=0.5)
        # plt.plot(S,z[ind[i],:],'crimson')

    plt.xlim(-0.1,1.1)
    # plt.axvline(x=phi_c,color='k',linestyle='--',linewidth=1)
    plt.xlabel(r'$\phi$',fontsize=16)
    plt.gca().yaxis.set_label_position("right")
    plt.gca().yaxis.tick_right()
    plt.ylabel(r'$z$',fontsize=16)
    plt.show()
    plt.close()



def plot_movie(ws,wi,phi,z,N_,timesteps):    
    q = (1-phi*sat(temp(z)))*(wi-ws)
    ind = np.arange(0,timesteps.size,6)
    for i in range(ind.size):
        plt.figure(figsize=(12,6))
        plt.suptitle(r'$t/t_\mathrm{max}=$'+'{:.2f}'.format(timesteps[ind[i]]/timesteps.max()),fontsize=24)
        plt.subplot(151)
        plt.plot(ws[ind[i],:],z[ind[i],:],color='indigo')  
        plt.xlabel(r'$w_\mathrm{s}$',fontsize=16)
        plt.ylabel(r'$z$',fontsize=16)
        plt.xlim(ws.min()-0.1,ws.max()+0.1)
        plt.ylim(z[ind,:].min(),z[ind,:].max())
        # plt.axhline(y=z[ind[i],:].max(),color='lightsteelblue',linewidth=4,linestyle='-')
        # plt.axhline(y=zf,color='burlywood',linewidth=4,linestyle='-')

        plt.subplot(152)
        plt.plot(wi[ind[i],:],z[ind[i],:],color='indigo')  
        plt.xlabel(r'$w_\mathrm{i}$',fontsize=16)
        plt.gca().set_yticks([])
        plt.xlim(wi.min()-0.1,wi.max()+0.1)
        plt.ylim(z[ind,:].min(),z[ind,:].max())
        # plt.axhline(y=z[ind[i],:].max(),color='lightsteelblue',linewidth=4,linestyle='-')
        # plt.axhline(y=zf,color='burlywood',linewidth=4,linestyle='-')


        plt.subplot(153)
        S = sat(temp(z[ind[i],:]))
        plt.plot(q[ind[i],:],z[ind[i],:],color='indigo')  
        plt.xlabel(r'$q$',fontsize=16)
        plt.gca().set_yticks([])
        plt.xlim(q.min()-0.1,q.max()+0.1)
        plt.ylim(z[ind,:].min(),z[ind,:].max())
        # plt.axhline(y=z[ind[i],:].max(),color='lightsteelblue',linewidth=4,linestyle='-')
        # plt.axhline(y=zf,color='burlywood',linewidth=4,linestyle='-')



        plt.subplot(154)
        plt.plot(N_[ind[i],:],z[ind[i],:],color='indigo')
        plt.xlabel(r'$N$',fontsize=16)
        plt.gca().set_yticks([])
        plt.xlim(N_.min()-0.1,N_.max()+0.1)
        plt.ylim(z[ind,:].min(),z[ind,:].max())
        plt.axvline(x=0,linestyle='--',color='k')
        # plt.axhline(y=z[ind[i],:].max(),color='lightsteelblue',linewidth=4,linestyle='-')
        # plt.axhline(y=zf,color='burlywood',linewidth=4,linestyle='-')


        plt.subplot(155)
        plt.plot(phi[ind[i],:],z[ind[i],:],color='indigo')
        # plt.axhline(y=z[ind[i],:].max(),color='lightsteelblue',linewidth=4,linestyle='-')
        # plt.axhline(y=zf,color='burlywood',linewidth=4,linestyle='-')
        plt.xlim(0,1)
        plt.ylim(z[ind,:].min(),z[ind,:].max())
        plt.xlabel(r'$\phi$',fontsize=16)
        plt.gca().yaxis.set_label_position("right")
        plt.gca().yaxis.tick_right()
        plt.ylabel(r'$z$',fontsize=16)
        plt.savefig('./movie/'+str(i))
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