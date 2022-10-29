
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

csvfile = open("D://Gergo//Research//Vegyületek//BPPA//spektofotometria//mn-cu_kin//kintest.csv")
csvreader = csv.reader(csvfile)
xdataL = []
ydataL = []
next(csvreader)
for row in csvreader:
    xdataL.append(float(row[0]))
    ydataL.append(float(row[1]))

xdata = numpy.array(xdataL)
ydata = numpy.array(ydataL)
    

"""next(csvreader)
row_count= len(list(csvreader))
print(row_count)
xdata = numpy.zeros(row_count)
ydata = numpy.zeros(row_count)
index = 0
for row in csvreader:
    xdata[index] = (float(row[0]))
    ydata[index] = (float(row[1]))
    index = index+1"""


t, A = symfit.variables('t, A')
k = symfit.Parameter('k', value=0.001, fixed=False)#, min=3, max=6
A0 = symfit.Parameter('A0', value=0.5, fixed=False)
Ainf = symfit.Parameter('Ainf', value=0.5, fixed=False)



model = symfit.Model({
    A: (A0 - Ainf) * symfit.exp(- k * t) + Ainf#,
#    y_2: y0 + a_2 * exp(- b_2 * x_2),
})

fit = symfit.Fit(model, t=xdata, A=ydata)
fit_result = fit.execute()

print(fit_result)

y = model(t=xdata, **fit_result.params)#A0=fit_result.value(A0), Ainf=fit_result.value(Ainf), k=fit_result.value(k))
plt.plot(xdata, y.output[0], c='orange')
plt.scatter(xdata, ydata, s=2)
plt.show()
