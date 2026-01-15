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

MW_CELL = 162.13940/1000
MW_GMSW = 132.11372/1000
MW_XYHW = MW_GMSW
MW_XYGR = MW_GMSW
MW_HCE1 = MW_GMSW
MW_HCE2 = MW_GMSW
MW_LIG = 258.27076/1000
MW_LIGC = MW_LIG
MW_LIGCC = MW_LIG
MW_LIGH = 436.45252/1000
MW_LIGO = 422.38248/1000
MW_LIGOH = 378.37348/1000
MW_TGL = 897.40400/1000
MW_TANN = 304.25208/1000
MW_ITANN = 164.11536/1000
MW_ACQUA = 18.01468/1000
MW_CELLA = MW_CELL
MW_H2S = 2.01568/1000
MW_CH4S =  16.04236/1000
MW_C2H4S = 28.05336/1000
MW_C2H6S = 30.06904/1000
MW_COS = 28.01000/1000
MW_CH2OS = 30.02568/1000
MW_COH2S = MW_CH2OS
MW_CO2S = 44.00900/1000
MW_CH3OHS = 32.04136/1000
MW_C6H5OHS = 94.11204/1000
MW_LVGS = 162.13940/1000
MW_COSTIFFS = 28.01000/1000
MW_CHARO = 28.01000/1000
MW_CHAR = 12.01100/1000
MW_ASH = 1/1000
# debiagi 2018 solid species
#speclist = ['CELL.solid', 'GMSW.solid', 'XYHW.solid', 'XYGR.solid', 'HCE1.solid', 'HCE2.solid','LIG.solid', 'LIGC.solid', 'LIGCC.solid', 'LIGH.solid', 
#                      'LIGO.solid', 'LIGOH.solid', 'TGL.solid', 'TANN.solid', 'ITANN.solid', 'ACQUA.solid', 'CELLA.solid', 'H2S.solid', 'CH4S.solid', 'C2H4S.solid',
#                      'C2H6S.solid', 'COS.solid', 'CH2OS.solid', 'COH2S.solid', 'CO2S.solid', 'CH3OHS.solid', 'C6H5OHS.solid', 'ASH.solid', 'CHAR.solid']
#MWlist = [MW_CELL, MW_GMSW, MW_XYHW, MW_XYGR, MW_HCE1, MW_HCE2, MW_LIG, MW_LIGC, MW_LIGCC, MW_LIGH, 
#           MW_LIGO, MW_LIGOH, MW_TGL, MW_TANN, MW_ITANN, MW_ACQUA, MW_CELLA, MW_H2S, MW_CH4S, MW_C2H4S,
#           MW_C2H6S, MW_COS, MW_CH2OS, MW_COH2S, MW_CO2S, MW_CH3OHS, MW_C6H5OHS, MW_ASH, MW_CHAR]

# new debiagi solid specs
#speclist = ['CELL.solid', 'GMSW.solid', 'XYHW.solid', 'XYGR.solid', 'HCE1.solid', 'HCE2.solid','LIG.solid', 'LIGC.solid', 'LIGCC.solid', 'LIGH.solid', 
#                      'LIGO.solid', 'LIGOH.solid', 'TGL.solid', 'TANN.solid', 'ITANN.solid', 'ACQUA.solid', 'CELLA.solid', 'H2S.solid', 'CH4S.solid', 'C2H4S.solid',
#                      'C2H6S.solid', 'COS.solid', 'CH2OS.solid', 'COH2S.solid', 'CO2S.solid', 'CH3OHS.solid', 'C6H5OHS.solid', 'LVGS.solid', 'COSTIFFS.solid', 'CHARO.solid', 'CHAR.solid']
#MWlist = [MW_CELL, MW_GMSW, MW_XYHW, MW_XYGR, MW_HCE1, MW_HCE2, MW_LIG, MW_LIGC, MW_LIGCC, MW_LIGH, 
#           MW_LIGO, MW_LIGOH, MW_TGL, MW_TANN, MW_ITANN, MW_ACQUA, MW_CELLA, MW_H2S, MW_CH4S, MW_C2H4S,
#           MW_C2H6S, MW_COS, MW_CH2OS, MW_COH2S, MW_CO2S, MW_CH3OHS, MW_C6H5OHS, MW_LVGS, MW_COSTIFFS, MW_CHARO, MW_CHAR]

MWlist = [MW_CELL, MW_GMSW, MW_XYHW, MW_XYGR, MW_HCE1, MW_HCE2, MW_LIG, MW_LIGC, MW_LIGCC, MW_LIGH, MW_LIGO, MW_LIGOH, MW_TGL, MW_TANN, MW_ITANN, MW_ACQUA, MW_CELLA, MW_H2S, MW_CH4S, MW_C2H4S, MW_C2H6S, MW_COS, MW_CH2OS, MW_COH2S, MW_CO2S, MW_CH3OHS, MW_C6H5OHS, MW_CHAR]
speclist = ['CELL.solid', 'GMSW.solid', 'XYHW.solid', 'XYGR.solid', 'HCE1.solid', 'HCE2.solid', 'LIG.solid', 'LIGC.solid', 'LIGCC.solid', 'LIGH.solid','LIGO.solid', 'LIGOH.solid', 'TGL.solid', 'TANN.solid', 'ITANN.solid', 'ACQUA.solid', 'CELLA.solid', 'H2S.solid', 'CH4S.solid', 'C2H4S.solid', 'C2H6S.solid', 'COS.solid', 'CH2OS.solid', 'COH2S.solid', 'CO2S.solid', 'CH3OHS.solid', 'C6H5OHS.solid', 'CHAR.solid']
calc_fn = ''
for s, mw in zip(speclist, MWlist):
    calc_fn+=f'"{s}"+'

print(calc_fn[:-1])

poro_part = 0.16
solfoam = pv.OpenFOAMReader(FileName = './soln.foam') # just need to provide folder
solfoam.CaseType = 'Reconstructed Case'
solfoam.MeshRegions = ['internalMesh']
# new debiagi scheme solid species
#solfoam.CellArrays = ['porosity', 'CELL.solid', 'GMSW.solid', 'XYHW.solid', 'XYGR.solid', 'HCE1.solid', 'HCE2.solid','LIG.solid', 'LIGC.solid', 'LIGCC.solid', 'LIGH.solid', 
#                      'LIGO.solid', 'LIGOH.solid', 'TGL.solid', 'TANN.solid', 'ITANN.solid', 'ACQUA.solid', 'CELLA.solid', 'H2S.solid', 'CH4S.solid', 'C2H4S.solid',
#                      'C2H6S.solid', 'COS.solid', 'CH2OS.solid', 'COH2S.solid', 'CO2S.solid', 'CH3OHS.solid', 'C6H5OHS.solid', 'LVGS.solid', 'COSTIFFS.solid', 'CHARO.solid', 'CHAR.solid']

# debiagi 2018 solid species
#solfoam.CellArrays = ['porosity', 'CELL.solid', 'GMSW.solid', 'XYHW.solid', 'XYGR.solid', 'HCE1.solid', 'HCE2.solid','LIG.solid', 'LIGC.solid', 'LIGCC.solid', 'LIGH.solid', 
#                      'LIGO.solid', 'LIGOH.solid', 'TGL.solid', 'TANN.solid', 'ITANN.solid', 'ACQUA.solid', 'CELLA.solid', 'H2S.solid', 'CH4S.solid', 'C2H4S.solid',
#                      'C2H6S.solid', 'COS.solid', 'CH2OS.solid', 'COH2S.solid', 'CO2S.solid', 'CH3OHS.solid', 'C6H5OHS.solid', 'ASH.solid', 'CHAR.solid']

solfoam.CellArrays = ['porosity', 'CELL.solid', 'GMSW.solid', 'XYHW.solid', 'XYGR.solid', 'HCE1.solid', 'HCE2.solid', 'LIG.solid', 'LIGC.solid', 'LIGCC.solid', 'LIGH.solid', 'LIGO.solid', 'LIGOH.solid', 'TGL.solid', 'TANN.solid', 'ITANN.solid', 'ACQUA.solid', 'CELLA.solid', 'H2S.solid', 'CH4S.solid', 'C2H4S.solid', 'C2H6S.solid', 'COS.solid', 'CH2OS.solid', 'COH2S.solid', 'CO2S.solid', 'CH3OHS.solid', 'C6H5OHS.solid', 'CHAR.solid', 'T', 'T.solid', 'rho.solid']
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
calculator1.ResultArrayName = 'm_solid'
calculator1.Function = '"rho.solid"'

# threshold to select gas region
threshold2 = pv.Threshold(registrationName='Threshold2', Input=solfoam)
threshold2.Scalars = ['POINTS', 'porosity']
threshold2.LowerThreshold = 1
threshold2.UpperThreshold = 1
calculator2 = pv.Calculator(Input=threshold2)
calculator2.AttributeType = 'Cell Data'
calculator2.ResultArrayName = 'Tgas'
calculator2.Function = '"T"'

#calculator1.Function = calc_fn[:-1]
# create a new 'Integrate Variables'
int1 = pv.IntegrateVariables(Input=calculator1)
int2 = pv.IntegrateVariables(Input=calculator2)

outfile=open("weightloss_rho.dat","w")
outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%('t[i]','vol','poro','Tgas','Tsolid','rhosolid','solid_m(0)','solid_m(t)','wloss'))
for i in range(N):
    try:
        idat    = dsa.WrapDataObject(pv.servermanager.Fetch(int1) )
        idat2    = dsa.WrapDataObject(pv.servermanager.Fetch(int2) )
        #print(idat2.CellData.keys())
        pv.UpdatePipeline(time=t[i], proxy=int1)
        vol       = idat.CellData['Volume'].item()
        rhosolid = idat.CellData['rho.solid'].item()/vol
        if i==0:
            m_init = rhosolid
        pv.UpdatePipeline(time=t[i], proxy=int2)
        vol2       = idat2.CellData['Volume'].item()
        poro   = idat.CellData['porosity'].item()/vol
        gasT   = idat2.CellData['Tgas'].item()/vol2
        solidT   = idat.CellData['T.solid'].item()/vol
        # division by volume cancels out with conversion to kg
        #m_t   = idat.PointData['m_solid'].item()*vol
        m_t   = rhosolid
        wloss   = m_t / m_init
        #print(idat.CellData.keys())
        print("processing time =",t[i],vol,poro,gasT,solidT,rhosolid,m_init,m_t,wloss)
        outfile.write("%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n"%(t[i],vol,poro,gasT,solidT,rhosolid,m_init,m_t,wloss))
    except Exception as e:
        print(f'Error at time {t[i]}: {e}')
        wloss=0
        print('setting weight loss to 0')
        outfile.write("%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n"%(t[i],vol,poro,gasT,solidT,rhosolid,m_init,m_t,wloss))

outfile.close()
