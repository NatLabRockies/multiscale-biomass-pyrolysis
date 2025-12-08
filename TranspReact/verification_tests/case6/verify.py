import yt
from sys import argv
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve

def nonlinsolve_func_left(a2,data):
    phiL=data[0]
    phiR=data[1]
    sigma_a=data[2]
    sigma_e=data[3]
    ocp=data[4]
    j0=data[5]
    phi0=data[6]
    xi=data[7]
    L=data[8]
    lhs=np.arcsinh(-sigma_e*a2/j0)*phi0
    rhs=sigma_e/sigma_a*a2*xi+phiL-a2*(xi-L)-phiR-ocp
    return(lhs-rhs)

def nonlinsolve_func_right(a2,data):
    phiL=data[0]
    phiR=data[1]
    sigma_a=data[2]
    sigma_e=data[3]
    ocp=data[4]
    j0=data[5]
    phi0=data[6]
    xi=data[7]
    L=data[8]
    lhs=np.arcsinh(sigma_a*a2/j0)*phi0
    rhs=a2*(xi-L)+phiR-sigma_a/sigma_e*a2*xi-phiL-ocp
    return(lhs-rhs)

root_dirctory = os.path.abspath(os.path.join(argv[1], os.pardir))+"/"
ds=yt.load(argv[1])
leftcase=int(argv[2])
axialdir=np.argmax(ds.domain_dimensions)
prob_lo=ds.domain_left_edge.d
prob_hi=ds.domain_right_edge.d
lengths=prob_hi-prob_lo
clength   = lengths[axialdir]
cwidth    = lengths[(axialdir+1)%3]
cdepth    = lengths[(axialdir+2)%3]
ar     = (clength/cwidth)

axialdir_char=chr(ord('x')+axialdir)
res=ds.domain_dimensions[axialdir]
slicedepth = cdepth/2
slc = ds.slice((axialdir+2)%3,slicedepth)
frb = slc.to_frb((1.0,'cm'),res)
x = np.linspace(0,clength,res)
fld_pot = np.array(frb["Potential"])[res//2,:]

#=======================================
#exact solution
#=======================================
ocp=0.2
sigma_a=3.0
sigma_e=1.0
j0=3.0
phi0=1.0

phiL=1.0
phiR=0.0
x_interface=0.2

if(not leftcase):
    phiL=0.0
    phiR=1.0
    x_interface=0.8

data=[phiL,phiR,sigma_a,sigma_e,ocp,j0,phi0,x_interface,clength]

a2=0.0
a1=0.0

if(leftcase):
    a2=fsolve(nonlinsolve_func_left,0.0,args=data)[0]
    a1=sigma_e/sigma_a*a2
else:
    a2=fsolve(nonlinsolve_func_right,0.0,args=data)[0]
    a1=sigma_a/sigma_e*a2

phi_exact=np.zeros(res)

for i in range(res):
    if(x[i]<x_interface):
        phi_exact[i]=a1*x[i]+phiL
    else:
        phi_exact[i]=a2*(x[i]-clength)+phiR


#=======================================
#Plot solutions
#=======================================
fig,ax=plt.subplots(1,1,figsize=(4,4))
ax.plot(x,fld_pot,'k-',label="echemAMR",markersize=2)
ax.plot(x,phi_exact,'r^',label="exact",markersize=2)
ax.legend(loc="best")

dir_char=chr(ord('x')+int(axialdir))
fig.suptitle("Potential solution along "+dir_char+" direction (AE) ")
plt.savefig(root_dirctory+"pot_"+dir_char+"_AE.png")
plt.show()
#=======================================

