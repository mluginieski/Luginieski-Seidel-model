## IdFLSm: Id current Fit of Luginieski-Seidel model
#                   Main code
#
#  This code uses the code IDcLSm, with the drain current equations from
#      Luginieski-Seidel model.
#
#  IdcLSm uses two auxiliary codes: constants, with some important physical constants
#      and the unity_converter, with is responsible for converting the input values to SI.
#
#  Input Experimental data: the code is set to read a csv file with ; as column separator;
#         By default, the code considers a data file without headers, and only two columns: Vds and Ids.
#         If the data file have headers, the number of header rows shoud be informed in the following code line:
#         data = pd.read_csv(fileName, sep=';', header=0, names=['Vd','Id'])
#
#         If the data file have more than 2 columns, the array names must be corrected.
#         The user is recommended to read the Pandas.read_csv docummentation:
#         https://pandas.pydata.org/docs/reference/api/pandas.read_csv.html
#
#  Guess: user shoud enter the guesses of three parameters in the following order:
#       >> Dc
#       >> alfa
#       >> gamma
#
#  Bounds: if user want to use a bound, the min and max values shoud be entered as:
#       >> min max
#          with a space between min and max values. The order is the same as for the guesses;
#          if user does not want to set a bound, it is only necessary to type:
#       >> inf
#
#  The authors acknowledge for the citation:
#  Adv. Theory Simul. 2023, 2200852
#  DOI: 10.1002/adts.202200852 
#
#  code created by: Msc. Marcos Luginieski
#          contact: mluginieski@ifsc.usp.com
#             date: 25/11/2022

import os
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
import IdcLSm as LSm
import warnings
from scipy.optimize import curve_fit
from scipy.optimize import differential_evolution

# Reads the file name and check if it exists
fileName = input('Enter the file name: ')
if os.path.isfile(fileName) != 1:
    print('File does not exist!')
    while os.path.isfile(fileName) != 1:
        fileName = input('Try a valid file name: ')
    
data = pd.read_csv(fileName, sep=';', header=0, names=['Vd','Id']) # Reads the experimental data file

expVd      = abs(data.Vd) # Experimental drain voltage
expId      = abs(data.Id) # Experimental drain current
lenExpData = len(expVd)   # Lenght of list expVd (needed for the IdcLSm code)
maxExpVd   = max(expVd)   # Maximum experimental Vd (needed for the IdcLSm code)

# User enter the initial guess of the fitted parameters
Dc_guess   = float(input('Initial guess of Dc: '))
a_guess    = float(input('Initial guess of a : '))
y_guess    = float(input('Initial guess of y : '))

# Initial guess array
p0 = []
p0.append(Dc_guess)
p0.append(a_guess)
p0.append(y_guess)

# User enter the bounds of the three parameters
print('\nBounds: enter the min and max bounds space-separated like: min  max')
Dc_bounds = tuple(input('Enter the Dc bounds: ').split())
a_bounds  = tuple(input('Enter the a bounds : ').split())
y_bounds  = tuple(input('Enter the y bounds : ').split())

if Dc_bounds[0] == 'inf':
    Dc_bounds = (-np.inf, np.inf)
if a_bounds[0] == 'inf':
    a_bounds = (-np.inf, np.inf)
if y_bounds[0] == 'inf':
    y_bounds = (-np.inf, np.inf)

bds       = ((Dc_bounds[0], a_bounds[0], y_bounds[0]), (Dc_bounds[1], a_bounds[1], y_bounds[1])) # Bounds array

# User inform if want to see the test parameters. ONLY LOWERCASE!
testPar   = input('\nWant to see the test parameters? (y/n): ')

# Function used for the fitting
def func(x, Dc, a, y):
    return LSm.Id_current(Dc, a, y, maxExpVd, lenExpData)

# Sum of Squared Error used for the minimization in the differential_evolution module
def sumOfSquaredError(parameterTuple):
    warnings.filterwarnings('ignore')
    val = func(expVd, *parameterTuple)
    return np.sum((expId - val) ** 2.0)

# Initial parameters generator for the differential_evolution module
def generate_Initial_Parameters():
    if Dc_bounds[1] == np.inf:
        parametersBounds = [(0.0, 1e-9), (1.0, 2.0), (0.1, 5.0)]
    else:
        parametersBounds = [Dc_bounds, a_bounds, y_bounds]
        
    guess  = [p0, p0, p0, p0, p0, p0]
    result = differential_evolution(sumOfSquaredError, parametersBounds, tol=1e-10, seed=3, polish=True, init=guess, updating='deferred')
  
    return result.x

# Call for the initial parameters generator
geneticParameters = generate_Initial_Parameters()

# Test parameters section
if testPar == 'y':

    # Test fitted parameters, standard error and lc
    testParameters, test_pcov = curve_fit(func, expVd, expId, geneticParameters)
    test_standardError        = np.sqrt(np.diag(test_pcov))
    lc                        = LSm.criticalThickness(testParameters[2])

    print('\nTest Parameters:')
    print('Dc = %.4f +- %.5f nm' % (testParameters[0] / 1e-9, test_standardError[0] / 1e-9))
    print('a  = %.3f +- %.4f' % (testParameters[1], test_standardError[1]))
    print('y  = %.3f +- %.4f' % (testParameters[2], test_standardError[2]))
    print('lc = %.4f nm' % (lc / 1e-9))

    # Calculation of Id curve for the graph
    xModel_test = np.linspace(min(expVd), max(expVd), lenExpData)
    yModel_test = func(xModel_test, *testParameters)

    # Graph of test and experimental data
    plt.plot(expVd, expId/1e-6, 'D', label='Exp')                       # Test data
    plt.plot(xModel_test, yModel_test/1e-6, linewidth='3', label='Fit') # Exp data
    plt.show()

# Fitted parameters section
fittedParameters, pcov = curve_fit(func, expVd, expId, p0, bounds=bds) # Calculation of fitted parameters
standardError          = np.sqrt(np.diag(pcov))                        # Standard error
lc                     = LSm.criticalThickness(fittedParameters[2])    # lc
modelPredictions       = func(expVd, *fittedParameters)                # Predicted Id current
absError               = modelPredictions - expId                      # Absolute Error
SE                     = np.sqrt(absError)                             # Squared Error
MSE                    = np.mean(SE)                                   # Mean squared error
RMSE                   = np.sqrt(MSE)                                  # Root mean squared error
Rsquared               = 1.0 - (np.var(absError) / np.var(expId))      # R squared

print('\nFitted parameters:')
print('Dc   = %.4f +- %.5f nm' % (fittedParameters[0] / 1e-9, standardError[0] / 1e-9))
print('a    = %.3f +- %.4f' % (fittedParameters[1], standardError[1]))
print('y    = %.3f +- %.4f' % (fittedParameters[2], standardError[2]))
print('lc   = %.4f nm' % (lc / 1e-9))
print('RMSE = ', RMSE)
print('r²   = ',Rsquared)

# Calculation of Id curve with the fitted parameters for the graph
xModel = np.linspace(min(expVd), max(expVd), lenExpData)
yModel = func(xModel, *fittedParameters)

# Graph of fit and exp data
plt.plot(expVd, expId/1e-6, 'D', label='Exp')             # Exp data
plt.plot(xModel, yModel/1e-6, linewidth='3', label='Fit') # Fit data
plt.show()
