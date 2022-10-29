import symfit
from lmfit import Parameters, fit_report, minimize
import scipy
from numpy import exp, sign, sin, pi, array
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import uncertainties
import asteval
import csv

def residual(pars, cH, cCu, data=None, eps=None):#t=x, 
    # unpack parameters: extract .value attribute for each parameter
    parvals = pars.valuesdict()

    k0 = parvals['k0']
    k1 = parvals['k1']
    k2 = parvals['k2']
    #k3 = parvals['k3']

    model = k0 + k1 * cH + k2 * cH**2 #+ k3 * cCu

    if data is None:
        return model
    if eps is None:
        return model - data
    return (model-data) / eps

def recalc(pars, t):#t=x, 
    # unpack parameters: extract .value attribute for each parameter
    parvals = pars.valuesdict()
    k = parvals['k']
    A0 = parvals['A0']
    Ainf = parvals['Ainf']
    result = (A0 - Ainf) * exp(- k * t) + Ainf
    return result

csvfile = open("D://Gergo//Research//Vegyületek//BPPA//spektofotometria//mn-cu_kin//kintest3d.csv")
csvreader = csv.reader(csvfile)
xdataL = []
ydataL = []
zdataL = []
next(csvreader)
for row in csvreader:
    xdataL.append(float(row[0]))
    ydataL.append(float(row[1]))
    zdataL.append(float(row[2]))

xdata = array(xdataL)
ydata = array(ydataL)
data = array(zdataL)

fit_params = Parameters()
fit_params.add('k0', value=1)
fit_params.add('k1', value=1)
fit_params.add('k2', value=1000)
#fit_params.add('k3', value=0.1)

out = minimize(residual, fit_params, args=(xdata,ydata,), kws={'data': data})

print(fit_report(out))

#plt.plot(xdata, recalc(out.params, xdata), c='orange')
#plt.scatter(xdata, data, s=2)
#plt.show()