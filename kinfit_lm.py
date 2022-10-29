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

def residual(pars, t, data=None, eps=None):#t=x, 
    # unpack parameters: extract .value attribute for each parameter
    parvals = pars.valuesdict()

    k = parvals['k']
    A0 = parvals['A0']
    Ainf = parvals['Ainf']

    model = (A0 - Ainf) * exp(- k * t) + Ainf

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

csvfile = open("D://Gergo//Research//Vegyületek//BPPA//spektofotometria//mn-cu_kin//kintest2.csv")
csvreader = csv.reader(csvfile)
xdataL = []
ydataL = []
next(csvreader)
for row in csvreader:
    xdataL.append(float(row[0]))
    ydataL.append(float(row[1]))

xdata = array(xdataL)
data = array(ydataL)

fit_params = Parameters()
fit_params.add('k', value=0.001)
fit_params.add('A0', value=0.1)
fit_params.add('Ainf', value=0.8)

out = minimize(residual, fit_params, args=(xdata,), kws={'data': data})

print(fit_report(out))

plt.plot(xdata, recalc(out.params, xdata), c='orange')
plt.scatter(xdata, data, s=2)
plt.show()