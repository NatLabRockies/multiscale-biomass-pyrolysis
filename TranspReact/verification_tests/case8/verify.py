import yt
from sys import argv
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve

root_dirctory = os.path.abspath(os.path.join(argv[1], os.pardir))+"/"
ds=yt.load(argv[1])
axialdir=np.argmax(ds.domain_dimensions)
prob_lo=ds.domain_left_edge.d
prob_hi=ds.domain_right_edge.d
lengths=prob_hi-prob_lo
clength   = lengths[axialdir]
print("axialdir:",axialdir)

ncells=ds.domain_dimensions
dx=lengths[:]/ncells[:]
res=ncells[axialdir]
eps=1e-10

loend_1d = prob_lo[axialdir]+0.5*dx[axialdir]+eps
hiend_1d = prob_hi[axialdir]-0.5*dx[axialdir]-eps
xarr = np.linspace(loend_1d,hiend_1d,res)

loend_3d=(0.5+eps)*(prob_lo+prob_hi)
hiend_3d=(0.5+eps)*(prob_lo+prob_hi)
loend_3d[axialdir]=loend_1d
hiend_3d[axialdir]=hiend_1d

lb = yt.LineBuffer(ds, (loend_3d[0],loend_3d[1],loend_3d[2]), \
        (hiend_3d[0],hiend_3d[1],hiend_3d[2]), res)
fld_pot=lb["Potential"]

#=======================================
#exact solution
#=======================================
E0=0.2
jco=1.0
alpha=50.0
sigma_c=0.0002
sigma_e=0.0001

J0=4.0
x_interface=0.5*clength

A=J0/sigma_c
B=J0/sigma_e
C=-1/alpha*np.log(jco/J0)-E0-J0*(1/sigma_e-1/sigma_c)*x_interface

phi_exact=np.zeros(res)

for i in range(res):
    if(xarr[i]<x_interface):
        phi_exact[i]=A*xarr[i]
    else:
        phi_exact[i]=B*xarr[i]+C

np.savetxt('exactsoln',np.transpose(np.vstack((xarr,phi_exact))),delimiter=' ')
#=======================================
#Plot solutions
#=======================================
fig,ax=plt.subplots(1,1,figsize=(4,4))
ax.plot(xarr,fld_pot,'k-',label="echemAMR",markersize=2)
ax.plot(xarr,phi_exact,'r^',label="exact",markersize=2)
ax.legend(loc="best")

dir_char=chr(ord('x')+int(axialdir))
fig.suptitle("Potential solution along "+dir_char+" direction (AE) ")
plt.savefig(root_dirctory+"pot_"+dir_char+"_AE.png")
plt.show()
#=======================================
