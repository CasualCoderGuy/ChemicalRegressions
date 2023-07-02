import symfit
import lmfit
import scipy
import numpy
import matplotlib
import uncertainties
import asteval
import csv
from sympy.concrete.summations import Sum
import pandas
import itertools



def ResidualConcentration(LnCeqRow, LnBetaColumn, AlphaMat, CtTotal, Weights=None, nonconverging=False):#pars:LnCeqRow
    #parvals = pars.valuesdict()
    #k = parvals['k']
    #LnCeq:row, kibővíteni, minden row ua
    LnCeqMat = LnCeqRow
    for j in range(len(LnBetaColumn)-1):
        LnCeqMat = numpy.vstack((LnCeqMat, LnCeqRow))
    SjColumn = numpy.exp(LnBetaColumn)*(numpy.exp(numpy.sum((LnCeqMat*AlphaMat), axis=1).reshape(-1,1)))
    #Sj:column, kibővíteni, minden column ua
    SjMat = SjColumn
    for j in range(len(LnCeqRow)-1):
        SjMat = numpy.hstack((SjMat, SjColumn))
    CtRow = numpy.sum(AlphaMat*SjMat, axis=0)
    if nonconverging:
        #if isH:
        #    return 
        return numpy.log(CtRow/CtTotal)
    return numpy.abs((CtTotal-CtRow)/CtTotal) #CtRow-CtTotal #
    #(model-data) / Weights

def ConcentrationJacobian(LnCeqRow, LnBetaColumn, AlphaMat, CtTotal, Weights=None, nonconverging=False):#pars:LnCeqRow
    #parvals = pars.valuesdict()
    #k = parvals['k']
    #LnCeq:row, kibővíteni, minden row ua
    LnCeqMat = LnCeqRow
    for j in range(len(LnBetaColumn)-1):
        LnCeqMat = numpy.vstack((LnCeqMat, LnCeqRow))
    SjColumn = numpy.exp(LnBetaColumn)*(numpy.exp(numpy.sum((LnCeqMat*AlphaMat), axis=1).reshape(-1,1)))
    #Sj:column, kibővíteni, minden column ua
    SjMat = SjColumn
    matsize = len(LnCeqRow)
    diag = numpy.zeros((matsize,matsize))
    diag[0,0] = numpy.exp(LnCeqRow[0])
    for j in range(matsize-1):
        SjMat = numpy.hstack((SjMat, SjColumn))
        diag[j+1,j+1] = numpy.exp(LnCeqRow[j+1])
    CtRow = numpy.sum(AlphaMat*SjMat, axis=0)
    jac = numpy.dot(numpy.transpose(AlphaMat*SjMat),AlphaMat)+diag#+diag(exp(compEqConcLN));
    return jac

def CalculateConcentration(LnCeqRow, LnBetaColumn, AlphaMat, CtTotal, Weights=None, nonconverging=False): #ret: conc, jac, res
    #parvals = pars.valuesdict()
    #k = parvals['k']
    #LnCeq:row, kibővíteni, minden row ua
    LnCeqMat = LnCeqRow
    for j in range(len(LnBetaColumn)-1):
        LnCeqMat = numpy.vstack((LnCeqMat, LnCeqRow))
    SjColumn = numpy.exp(LnBetaColumn)*(numpy.exp(numpy.sum((LnCeqMat*AlphaMat), axis=1).reshape(-1,1)))
    #Sj:column, kibővíteni, minden column ua
    SjMat = SjColumn
    matsize = len(LnCeqRow)
    diag = numpy.zeros((matsize,matsize))
    diag[0,0] = numpy.exp(LnCeqRow[0])
    for j in range(matsize-1):
        SjMat = numpy.hstack((SjMat, SjColumn))
        diag[j+1,j+1] = numpy.exp(LnCeqRow[j+1])
    CtRow = numpy.sum(AlphaMat*SjMat, axis=0)
    jac = numpy.dot(numpy.transpose(AlphaMat*SjMat),AlphaMat)+diag#+diag(exp(compEqConcLN));
    res = numpy.abs((CtTotal-CtRow)/CtTotal) #CtRow-CtTotal #
    if nonconverging:
        #if isH:
        #    return 
        res = numpy.log(CtRow/CtTotal)
    return res, jac, CtRow

def CalcTotalConcentrationMatrix(CinitialRow, CcompInTitrantRow, AddedVolumeCol, Vstart):
    return (numpy.dot(AddedVolumeCol.reshape(-1,1),CcompInTitrantRow.reshape(1,-1))+CinitialRow*Vstart)/((AddedVolumeCol+Vstart)[:,None]) #create mat from CinitialRow

def StartCeqCalc(CtotalMat, C0min):
    mat = numpy.where(CtotalMat < C0min, C0min, CtotalMat)
    return numpy.log(mat/2)

def NewtonRaphson(LnCeqRow, LnBetaColumn, AlphaMat, CtTotal, maxit, maxresidual):
    print("NR")
    #CtCalc = Concentration(LnCeqRow, LnBetaColumn, AlphaMat, CtTotal)
    CtCalc, jac, res = CalculateConcentration(LnCeqRow, LnBetaColumn, AlphaMat, CtTotal)
    for it in range(maxit):
        dX=numpy.dot(jac.T,(CtCalc-CtTotal))
        dX = numpy.where(dX>4.6, 4.6, dX)
        dX = numpy.where(dX<4.6, -4.6, dX)
        dX = numpy.where((dX<0).any() and (dX>-1e-6).any(), -1e-2, dX)
        dX = numpy.where((1e-6>dX).any() and (dX>0).any(), 1e-2, dX)
        LnCeqRow = LnCeqRow-dX
        #CtCalc = Concentration(LnCeqRow, LnBetaColumn, AlphaMat, CtTotal)
        res, jac, CtCalc = CalculateConcentration(LnCeqRow, LnBetaColumn, AlphaMat, CtTotal)
        if all(i < maxresidual for i in res):
            break
    return LnCeqRow


excel_data = pandas.read_excel("D://Gergo//psequad//PROGRAM//MATLAB//psequad_test_import2.xlsx", sheet_name=None)
AlphaMat = numpy.array(excel_data['CompMat'])[:,1:].astype(float)
print(AlphaMat)
CinitialRow = numpy.array(excel_data['Components'])[:,0].astype(float)
print(CinitialRow)
Beta = numpy.array(excel_data['Species'])[:,0].astype(float)
print(Beta)
LnBetaColumn = (numpy.log(10**Beta)).reshape(-1,1)
print(LnBetaColumn)
DisabledCol = numpy.array(excel_data['Species'])[:,1].astype(float)
print(DisabledCol)
TitrData = numpy.array(excel_data['Titration']).astype(float)
CcompInTitrantRow = numpy.array([-0.149987, 0,0]).astype(float)
print(CcompInTitrantRow)
totC = CalcTotalConcentrationMatrix(CinitialRow,CcompInTitrantRow,TitrData[:,0],6.0117)
print(totC[0,:])
ceq = StartCeqCalc(totC, 0.000001)
#ceq = numpy.array([2.06155E-02, 1.23836E-17, 1.27912E-06])
#ceq = numpy.log(ceq)
print(ceq[0,:]) #[0,:]


#result = scipy.optimize.root(CalculateConcentration, ceq[0,:], method='lm', jac=True, args=(LnBetaColumn, AlphaMat, totC[0,:]), options={})
#result = scipy.optimize.root(ResidualConcentration, ceq[0,:], jac=ConcentrationJacobian, args=(LnBetaColumn, AlphaMat, totC[0,:]))
#result, res2, res3 = CalculateConcentration(ceq, LnBetaColumn, AlphaMat, totC[0,:])
result = NewtonRaphson(ceq[0,:], LnBetaColumn, AlphaMat, totC[0,:], 300000, 0.3)
print(result)
#print(numpy.exp(result.x))


xdata=[1, 1]
comp_no = 2
spec_no = 6
assoc_no = spec_no-comp_no
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
"""for j in range(assoc_no):
    b.append(symfit.Parameter('b'+str(j), value=b_data[j], fixed=fixed_b[j]))

for l in range(no_data_points):
    Lcts.append([])
    x.append(symfit.Variable('ct'+str(l)+'x'+str(t), value=Lct_data[l, t], fixed=True))
    for t in range(comp_no):
        Lcts[l].append(symfit.Parameter('ct'+str(l)+'x'+str(t), value=Lct_data[l, t], fixed=True))
        Ls[l].append(symfit.Parameter('L'+str(l)+'x'+str(t), value=Lct_data[l, t], fixed=(comp_no == measured_comp)))"""


