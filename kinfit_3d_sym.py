
import symfit
import lmfit
import scipy
import numpy
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import uncertainties
import asteval
import csv

#sklearn

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

xdata = numpy.array(xdataL)
ydata = numpy.array(ydataL)
zdata = numpy.array(zdataL)
    

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


cH, cCu, kobs = symfit.variables('cH, cCu, kobs')
k0 = symfit.Parameter('k0', value=10000, fixed=False)#, min=3, max=6
k1 = symfit.Parameter('k1', value=10000, fixed=False)#, min=3, max=6
k2 = symfit.Parameter('k2', value=10000, fixed=False)#, min=3, max=6
k3 = symfit.Parameter('k3', value=10000, fixed=False)#, min=3, max=6



model = symfit.Model({
    kobs: k0 + k1 * cH + k2 * cH**2 + k3 * cCu  #,
#    y_2: y0 + a_2 * exp(- b_2 * x_2),
})

fit = symfit.Fit(model, cH=xdata, cCu=ydata, kobs=zdata)
fit_result = fit.execute()

print(fit_result)

z = model(cH=xdata, cCu=ydata, **fit_result.params)#A0=fit_result.value(A0), Ainf=fit_result.value(Ainf), k=fit_result.value(k))

#fig = plt.figure()
#ax = fig.add_subplot(projection='3d')
#Axes3D.plot_surface(xdata, ydata, z.output[0])
#Axes3D.plot_surface(xdata, ydata, zdata)
#Axes3D.show()
