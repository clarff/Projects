import matplotlib.pyplot as plt
import numpy as np
from pylab import *

#constants
T = 303 #temperature in K
C_elyte = 5e-3 #5 milimolar solution, from Battistel et al 2012, natural brine conc
dY = 30e-6 #meters, distance between nodes
R = 8.3145 # J/mol-K, gas constant
F = 96485  # C/mol equiv

#==========================================================================================
#==========Select (comment out unwanted code) Case1/2/3 BELOW===============================
#mole fractions, respectively Li+, Cl-, Na+, Ca2+, solvent

#Case 1, Li+, Cl-, solvent
X_k_1 = np.array([0.03, 0.03, 0.94])
X_k_2 = np.array([0.06,0.06,0.88])

#Case 2, Li+, Cl-, Na+, solvent
# X_k_1 = np.array([0.0003, 0.03, 0.03, 0.9397])
# X_k_2 = np.array([0.0006,0.06, 0.06, 0.8794])

# #Case 3, Li+, Cl-, Na+, Ca2+, solvent
# X_k_1 = np.array([0.0003, 0.03, 0.03, 0.03, 0.9097])
# X_k_2 = np.array([0.0006,0.06, 0.06, 0.06, 0.8194])

#charge equivalents, same order as above for each case

#Case 1
z_k = np.array([1., -1., 0.])
#
# # #Case 2
#z_k = np.array([1., -1., 1., 0.])
#
# Case 3
# z_k = np.array([1., -1., 1., 2., 0.])

#diffusion coefficients
#Case 1
D_k = np.array([1.52e-10, 0.25e-10, 1e-12])

#Case 2
#D_k = np.array([1.52e-10, 0.25e-10, 0.5e-10, 1e-12])

#Case 3
# D_k = np.array([1.52e-10, 0.25e-10, 0.5e-10, 0.5e-10, 1e-12])

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#===========VARIABLES BELOW DO NOT CHANGE FOR EACH CASE===================

#electric potentials at node 1 and node 2
phi_1 = 0.9
phi_2 = 0.5

#active material particle diameter
d_part = 5e-6

#porosity
eps_elyte = .4

#tortuosity
n_Brugg = -0.5


#====DICTIONARIES====================================================


# State variables for node 1:
s1 = {'X_k':X_k_1, 'T':T, 'phi':phi_1, 'C_k':X_k_1*C_elyte}
# State variables for node 2:
s2 = {'X_k':X_k_2, 'T':T, 'phi':phi_2, 'C_k':X_k_2*C_elyte}

#interface concentrations
C_k_int = (s1['C_k']+s2['C_k'])/2

# Geometric and microstructure parameters:
geom = {'eps_elyte':eps_elyte, 'n_Brugg':n_Brugg, 'd_part':d_part, 'dY':dY}

# Electrolyte properties
elyte_pars = {'D_k':D_k, 'z_k':z_k, 'C_elyte':C_elyte}


#====TRANSPORT FUNCTION=================================================
subplot(1,2,1)
def electrolyte_transport(state1, state2, geom, elyte_pars):
    N_k = np.zeros_like(s1['X_k'])
    #gradC = (s2['C_k']-s1['C_k'])/geom['dY']
    gradX_k = (s2['X_k']-s1['X_k'])/geom['dY']
    gradPhi = (s2['phi'] - s1['phi'])/geom['dY']
    D_k_eff = geom['eps_elyte']**1.5*elyte_pars['D_k']
    D_k_mig = D_k_eff*C_k_int*elyte_pars['z_k']*F/R/T
    X_k_int = (s1['X_k']+s2['X_k'])/2
    #per second
    N_k = -D_k_eff*gradX_k*C_k_int/X_k_int-D_k_mig*gradPhi
    #per day
    #N_k = -D_k_eff*gradX_k*C_k_int/X_k_int-D_k_mig*gradPhi*86400

    return N_k

#vary cell potential from 0-1.1 Volts
dPhi = np.linspace(0,1.1,25)
currents = np.zeros_like(dPhi)
N_k = np.zeros((len(dPhi), len(z_k)))

#===========loop through different Node 2 potentials between 0-1.1 ============
for j, phi in enumerate(dPhi):
    s2['phi'] = phi
    N_k[j,:] = electrolyte_transport(s1,s2, geom, elyte_pars)
    currents[j] = np.dot(z_k,N_k[j,:])*F

#=============PLOT: where current will be equal to zero============
plt.plot(dPhi, currents, 'k', color='green')
# plt.plot(dPhi, current_check, 'ro', markerfacecolor=None)
plt.plot(dPhi, np.zeros_like(dPhi),'--',color='0.5')
plt.xlabel('Electric potential difference (V)',fontsize=10)
plt.ylabel('Current density (A/m$^2$)',fontsize=12)
plt.title("Charge transfer in brine")
zero=np.interp(0, np.flip(currents), np.flip(dPhi))
print('Zero current at dPhi = ',zero)
# plt.savefig("IV_curve.png")


#=================MOLAR FLUX PLOT===================
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

subplot(1,2,2)
plt.plot(dPhi, N_k[:,0]*86400,linewidth=2.5)
plt.plot(dPhi, N_k[:,1]*86400,linewidth=2.5)
# plt.plot(dPhi, N_k[:,2]*86400,linewidth=2.5)
# plt.plot(dPhi, N_k[:,3]*86400,linewidth=2.5)


plt.xlabel('Electric potential difference (V)',fontsize=10)
plt.ylabel('Molar flux (mol/m$^2$-day)',fontsize=12)
plt.title('Molar flux of ions in brine')
plt.legend(['Li$^+$','Cl$^-$'],frameon=False,fontsize=12)
#Case3 plt.legend(['Li$^+$','Cl$^-$','Na$^+$','Ca$^{2+}$'],frameon=False,fontsize=12)


plt.plot([zero,zero],[N_k[-1,0],N_k[0,0]],'--',color='0.5')
# plt.axvline(0.8868918573910756, color = 'cyan')
plt.plot([0,1],[0,0],'--',color='0.5')

#give plot a meta title
# plt.suptitle('Case 3', fontsize=14)

plt.subplots_adjust(left=None, bottom= None, right=None, top=None, wspace=.75, hspace=None)
pngtitle = input('type in your png title: ')
plt.savefig(str(pngtitle), bbox_inches = 'tight',dpi=300)

plt.show()

