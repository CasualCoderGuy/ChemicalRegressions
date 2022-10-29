
import symfit
import lmfit
import scipy
import numpy
import matplotlib
import matplotlib.pyplot as plt
import uncertainties
import asteval
import csv

#sklearn

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

cldata = numpy.array(cl_dat)
ctdata = numpy.array(ct_dat)
r1data = numpy.array(r1_dat)
r2data = numpy.array(r2_dat)
r3data = numpy.array(r3_dat)
r4data = numpy.array(r4_dat)


ROBS1, ROBS2, ROBS3, ROBS4, CL, CT = symfit.variables('ROBS1, ROBS2, ROBS3, ROBS4, CL, CT')
Ka = symfit.Parameter('Ka', value=68.0, min=1.0, fixed=False)#, max=6
N = symfit.Parameter('N', value=1.0, fixed=True)
RF1 = symfit.Parameter('RF1', value=2.74, fixed=True) # 2.74 | 2.05 @60MHz T1
RF2 = symfit.Parameter('RF2', value=7.28, fixed=True) # 7.28 | 5.69 @60MHz T2
RF3 = symfit.Parameter('RF3', value=3.11, fixed=True) # 3.11 | 2.48 @20MHz T1
RF4 = symfit.Parameter('RF4', value=4.30, fixed=True) # 4.30 | 3.52 @20MHz T2
RD1 = symfit.Parameter('RD1', value=0.38, fixed=True) # 0.38 | 0.29
RD2 = symfit.Parameter('RD2', value=0.50, fixed=True) # 0.50 | 0.41
RD3 = symfit.Parameter('RD3', value=0.34, fixed=True) # 0.34 | 0.26
RD4 = symfit.Parameter('RD4', value=0.49, fixed=True) # 0.49 | 0.43
RB1, RB2, RB3, RB4 = symfit.parameters('RB1, RB2, RB3, RB4', value=50.0, min=0.1, max=10000.0, fixed=False)

model = symfit.Model({
    ROBS1: (((Ka*CT+N*CL*Ka+1)-((Ka*CT+N*CL*Ka+1)**2-4*Ka**2*CT*N*CL)**(1/2))/(2*Ka))*1000*(RB1-RF1)+(RF1*CT*1000)+RD1,
    ROBS2: (((Ka*CT+N*CL*Ka+1)-((Ka*CT+N*CL*Ka+1)**2-4*Ka**2*CT*N*CL)**(1/2))/(2*Ka))*1000*(RB2-RF2)+(RF2*CT*1000)+RD2,
    ROBS3: (((Ka*CT+N*CL*Ka+1)-((Ka*CT+N*CL*Ka+1)**2-4*Ka**2*CT*N*CL)**(1/2))/(2*Ka))*1000*(RB3-RF3)+(RF3*CT*1000)+RD3,
    ROBS4: (((Ka*CT+N*CL*Ka+1)-((Ka*CT+N*CL*Ka+1)**2-4*Ka**2*CT*N*CL)**(1/2))/(2*Ka))*1000*(RB4-RF4)+(RF4*CT*1000)+RD4
})

fit = symfit.Fit(model, ROBS1=r1data, ROBS2=r2data, ROBS3=r3data, ROBS4=r4data, CL=cldata, CT=ctdata)
fit_result = fit.execute()

print(fit_result)
