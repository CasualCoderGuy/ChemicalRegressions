import symfit
from lmfit import Parameters, fit_report, minimize
import scipy
from numpy import exp, sign, sin, pi, array, concatenate
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import uncertainties
import asteval
import csv

def residual(pars, CL, CT, data1=None, data2=None, data3=None, data4=None):#, eps=None):#t=x, 
    # unpack parameters: extract .value attribute for each parameter
    parvals = pars.valuesdict()

    Ka = parvals['Ka']
    N = parvals['N']
    RF1 = parvals['RF1']
    RF2 = parvals['RF2']
    RF3 = parvals['RF3']
    RF4 = parvals['RF4']
    RD1 = parvals['RD1']
    RD2 = parvals['RD2']
    RD3 = parvals['RD3']
    RD4 = parvals['RD4']
    RB1 = parvals['RB1']
    RB2 = parvals['RB2']
    RB3 = parvals['RB3']
    RB4 = parvals['RB4']

    ROBS1 = (((Ka*CT+N*CL*Ka+1)-((Ka*CT+N*CL*Ka+1)**2-4*Ka**2*CT*N*CL)**(1/2))/(2*Ka))*1000*(RB1-RF1)+(RF1*CT*1000)+RD1
    ROBS2 = (((Ka*CT+N*CL*Ka+1)-((Ka*CT+N*CL*Ka+1)**2-4*Ka**2*CT*N*CL)**(1/2))/(2*Ka))*1000*(RB2-RF2)+(RF2*CT*1000)+RD2
    ROBS3 = (((Ka*CT+N*CL*Ka+1)-((Ka*CT+N*CL*Ka+1)**2-4*Ka**2*CT*N*CL)**(1/2))/(2*Ka))*1000*(RB3-RF3)+(RF3*CT*1000)+RD3
    ROBS4 = (((Ka*CT+N*CL*Ka+1)-((Ka*CT+N*CL*Ka+1)**2-4*Ka**2*CT*N*CL)**(1/2))/(2*Ka))*1000*(RB4-RF4)+(RF4*CT*1000)+RD4

    if data1 is None:
        resid1 = ROBS1
    else:
        resid1 = ROBS1-data1
    if data2 is None:
        resid2 = ROBS2
    else:
        resid2 = ROBS2-data2
    if data3 is None:
        resid3 = ROBS3
    else:
        resid3 = ROBS3-data3
    if data4 is None:
        resid4 = ROBS4
    else:
        resid4 = ROBS4-data4
    return concatenate((resid1, resid2, resid3, resid4))

csvfile = open("D://Gergo//Research//Vegyületek//BPPA//relax//HSA//affi2.csv")
csvreader = csv.reader(csvfile)
cl_dat = []
ct_dat = []
r1_dat = []
r2_dat = []
r3_dat = []
r4_dat = []
next(csvreader)
for row in csvreader:
    cl_dat.append(float(row[0]))
    ct_dat.append(float(row[1]))
    r1_dat.append(float(row[2]))
    r2_dat.append(float(row[3]))
    r3_dat.append(float(row[4]))
    r4_dat.append(float(row[5]))

cldata = array(cl_dat)
ctdata = array(ct_dat)
r1data = array(r1_dat)
r2data = array(r2_dat)
r3data = array(r3_dat)
r4data = array(r4_dat)

fit_params = Parameters()
fit_params.add('N', value=1.0, min=0, max=6, vary=False)
fit_params.add('RF1', value=2.74, vary=False) # 2.74 | 2.05 @60MHz T1
fit_params.add('RF2', value=7.28, vary=False) # 7.28 | 5.69 @60MHz T2
fit_params.add('RF3', value=3.11, vary=False) # 3.11 | 2.48 @20MHz T1
fit_params.add('RF4', value=4.30, vary=False) # 4.30 | 3.52 @20MHz T2
fit_params.add('RD1', value=0.38, vary=False) # 0.38 | 0.29
fit_params.add('RD2', value=0.50, vary=False) # 0.50 | 0.41
fit_params.add('RD3', value=0.34, vary=False) # 0.34 | 0.26
fit_params.add('RD4', value=0.49, vary=False) # 0.49 | 0.43
fit_params.add('RB1', value=100, min=1.0, max=10000)
fit_params.add('RB2', value=100, min=1.0, max=10000)
fit_params.add('RB3', value=100, min=1.0, max=10000)
fit_params.add('RB4', value=100, min=1.0, max=10000)
fit_params.add('Ka', value=68, min=1.0, max=10000) #, max=6

out = minimize(residual, fit_params, args=(cldata,ctdata,), kws={'data1': r1data, 'data2': r2data, 'data3': r3data, 'data4': r4data}, method='differential_evolution')

print(fit_report(out))