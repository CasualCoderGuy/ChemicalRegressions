import symfit
import lmfit
import scipy
import numpy
import matplotlib
import uncertainties
import asteval
import csv
from sympy.concrete.summations import Sum

csvfile = open("D://Gergo//Research//Vegyületek//BPPA//bppa.csv")
csvreader = csv.reader(csvfile)
xdataL = []
ydataL = []
next(csvreader)
for row in csvreader:
    xdataL.append(float(row[0]))
    ydataL.append(float(row[1]))

xdata = numpy.array(xdataL)
ydata = numpy.array(ydataL)

comp_no = 2
spec_no = 6
assoc_no = comp_no - spec_no
no_data_points = len(xdata)

measured_comp = 1
A_factor = 0.0 #Irving factor
V0 = 0.0 #starting volume
c0_data = numpy.zeros(comp_no) #starting component conc.
V_data = numpy.zeros(no_data_points) #volume data
Lct_data = numpy.zeros((no_data_points, comp_no)) #Lct = ln(ct) where ct is the total component conc.
b_data = numpy.zeros(assoc_no) #b=ln(beta) where beta is the form. const. of associates
fixed_b = [False]*assoc_no
fixed_b[0] = True
a_data = numpy.zeros((assoc_no, comp_no)) #alfa(j,t) :: component matrix, assign as consts
L_data = numpy.zeros((no_data_points, comp_no)) #L = ln(c) where c is the eq. conc. of components
x_data = numpy.zeros(no_data_points) # conc. of comp. without correction (eg.: 10^-pH)
Lct = []
b = []
L = [] #initial guess: ct
Lcts = []
Ls = []
x = []
#Lcts = [Lct for l in range(no_data_points)]
#Ls = [L for l in range(no_data_points)]
for j in range(assoc_no):
    b.append(symfit.Parameter('b'+str(j), value=b_data[j], fixed=fixed_b[j]))

for l in range(no_data_points):
    Lcts.append([])
    x.append(symfit.Variable('ct'+str(l)+'x'+str(t), value=Lct_data[l, t], fixed=True))
    for t in range(comp_no):
        Lcts[l].append(symfit.Parameter('ct'+str(l)+'x'+str(t), value=Lct_data[l, t], fixed=True))
        Ls[l].append(symfit.Parameter('L'+str(l)+'x'+str(t), value=Lct_data[l, t], fixed=(comp_no == measured_comp)))


model_dict = {
    z: a/(y * b) *  exp(- a * x)
        for x, y, z in zip(xs, ys, zs)
}