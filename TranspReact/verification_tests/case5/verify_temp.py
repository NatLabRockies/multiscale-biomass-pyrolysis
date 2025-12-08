import yt
from sys import argv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
    

clength   = 1.0
cwidth    = 0.0625
cdepth    = 0.0625
ar     = (clength/cwidth)


ds=yt.load(argv[1])
axialdir=np.argmax(ds.domain_dimensions)
axialdir_char=chr(ord('x')+axialdir)
res=100
slicedepth = cdepth/2
slc = ds.slice((axialdir+2)%3,slicedepth)
frb = slc.to_frb((1.0,'cm'),res)
x = np.linspace(0,clength,res)
fld_temp = np.array(frb["Temperature"])[res//2,:]
x1=0.2
x2=0.6
k1=0.1
k2=1.0
k3=0.2

R1=x1/k1
R2=(x2-x1)/k2
R3=(1.0-x2)/k3

Tleft=1.0
Tright=0.0

#Ohm's law for total current
I_tot=(Tleft-Tright)/(R1+R2+R3)

T0=1.0
T12=Tleft-I_tot*R1
T23=T12-I_tot*R2

exactsoln=np.zeros(res)
for i in range(res):
    if(x[i]<x1):
        exactsoln[i]=T0-(T0-T12)/x1*x[i]
    elif(x[i]>=x1 and x[i]<x2):
        exactsoln[i]=T12-(T12-T23)/(x2-x1)*(x[i]-x1)
    else:
        exactsoln[i]=T23-(T23-Tright)/(1.0-x2)*(x[i]-x2)

#=======================================

#=======================================
#Plot solutions
#=======================================
fig,ax=plt.subplots(1,1,figsize=(4,4))
ax.plot(x,exactsoln,'r-',label="Exact solution")
ax.plot(x,fld_temp,'k*',label="TR",markersize=2)
ax.legend(loc="best")

dir_char=axialdir_char
fig.suptitle("T solution along "+dir_char+" direction ")
plt.savefig("temperature_"+dir_char+".png")
plt.show()
#=======================================

