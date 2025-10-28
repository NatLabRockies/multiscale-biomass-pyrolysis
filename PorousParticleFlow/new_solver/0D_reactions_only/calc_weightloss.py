# trace generated using paraview version 5.8.1
#
# To ensure correct image size when batch processing, please search 
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
import numpy as np
from paraview import simple as pv
import vtk.numpy_interface.dataset_adapter as dsa 
import sys
from sys import argv

MW_CELL = 123.7/1000
MW_HCELL = 119.2/1000
MW_LIG = 116.5/1000
MW_ACELL = 123.7/1000
MW_AHCELL = 119.2/1000
MW_ALIG = 116.5/1000
MW_TRC = 100/1000
MW_TRH = 100/1000
MW_TRL = 100/1000
MW_SGC = 30/1000
MW_SGH = 30/1000
MW_SGL = 30/1000
MW_CHAR = 12.0/1000


speclist = ['CELL.solid', 'HCELL.solid', 'LIG.solid', 'ACELL.solid', 'AHCELL.solid', 'ALIG.solid', 'CHAR.solid']

MWlist = [MW_CELL, MW_HCELL, MW_LIG, MW_ACELL, MW_AHCELL, MW_LIG, MW_CHAR]

calc_fn = ''
for s, mw in zip(speclist, MWlist):
    calc_fn+=f'"{s}"*{mw}+'

print(calc_fn[:-1])

poro_part = 0.1
solfoam = pv.OpenFOAMReader(FileName = './soln.foam') # just need to provide folder
solfoam.CaseType = 'Reconstructed Case'
solfoam.MeshRegions = ['internalMesh']
solfoam.CellArrays = ['porosity', 'CELL.solid', 'HCELL.solid', 'LIG.solid', 'ACELL.solid', 'AHCELL.solid', 'ALIG.solid', 'CHAR.solid', 'T', 'T.solid']
t = np.array(solfoam.TimestepValues)
N=t.size
print(N)
# create a new 'Threshold'
threshold1 = pv.Threshold(registrationName='Threshold1', Input=solfoam)
threshold1.Scalars = ['POINTS', 'porosity']
threshold1.LowerThreshold = 0.0
#threshold1.UpperThreshold = poro_part+0.01
threshold1.UpperThreshold = 0.9

calculator1 = pv.Calculator(Input=threshold1)
calculator1.AttributeType = 'Point Data'
calculator1.ResultArrayName = 'rho_solid'
calculator1.Function = calc_fn[:-1]
# create a new 'Integrate Variables'
int1 = pv.IntegrateVariables(Input=calculator1)
int2 = pv.IntegrateVariables(Input=solfoam)

outfile=open("weightloss.dat","w")
outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%('t[i]','vol','poro','Tgas','Tsolid','solid_m(0)','solid_m(t)','wloss'))
for i in range(N):
    try:
        idat    = dsa.WrapDataObject(pv.servermanager.Fetch(int1) )
        idat2    = dsa.WrapDataObject(pv.servermanager.Fetch(int2) )
        #print(idat.CellData.keys())
        pv.UpdatePipeline(time=t[i], proxy=int1)
        pv.UpdatePipeline(time=t[i], proxy=int2)
        if i==0:
            m_init = idat.PointData['rho_solid'].item()
        vol       = idat.CellData['Volume'].item()
        vol2       = idat2.CellData['Volume'].item()
        poro   = idat.CellData['porosity'].item()/vol
        gasT   = idat2.CellData['T'].item()/vol2
        solidT   = idat.CellData['T.solid'].item()/vol
        # division by volume cancels out with conversion to kg
        m_t   = idat.PointData['rho_solid'].item()
        wloss   = m_t / m_init
        #print(idat.CellData.keys())
        print("processing time =",t[i],vol,poro,gasT,solidT,m_init,m_t,wloss)
        outfile.write("%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n"%(t[i],vol,poro,gasT,solidT,m_init,m_t,wloss))
    except Exception as e:
        print(f'Error at time {t[i]}: {e}')
        wloss=0
        print('setting weight loss to 0')
        outfile.write("%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n"%(t[i],vol,poro,gasT,solidT,m_init,m_t,wloss))

outfile.close()
