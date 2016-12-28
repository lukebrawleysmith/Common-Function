# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 08:25:33 2016

@author: luksmi
"""
import sys
#sys.path.append('/local/projects/luksmi/PythonScripts')
sys.path.append('/home/local/ANT/luksmi/Documents/PythonScripts')
from CommonFunctions import *

pd.set_option('html', False)
pd.set_option('max_columns', 30)
pd.set_option('max_rows', 20)
np.set_printoptions(linewidth= 200)


# checkInteger
negTwoInt = int(-2)
zeroInt = int(0)
twoInt = int(2)
intArrayNeg = np.array(range(0, -10, -2), dtype = np.int64)
intArrayZero = np.zeros(shape = (2, 2,), dtype = np.int64)
intArray = np.array(3, dtype = np.int64)

checkInteger(negTwoInt, "negTwoInt")
checkInteger(negTwoInt, "negTwoInt", nonNegative = True)
checkInteger(zeroInt, "zeroInt", nonNegative = True)
checkInteger(zeroInt, "zeroInt", nonNegative = True, positive = True)
checkInteger(twoInt, "twoInt", nonNegative = True, positive = True)
checkInteger(twoInt, "twoInt", 
  nonNegative = True, positive = True, isArray = True)

checkInteger(intArrayNeg, "intArrayNeg", isArray = False)
checkInteger(intArrayNeg, "intArrayNeg", isArray = True)
checkInteger(intArrayNeg, "intArrayNeg", 
   isArray = True, lengthOne = False)
checkInteger(intArrayNeg, "intArrayNeg", 
   isArray = True, lengthOne = False, nonNegative = True)
checkInteger(intArrayZero, "intArrayZero", 
   isArray = True, lengthOne = False, nonNegative = True)
checkInteger(intArrayZero, "intArrayZero", 
   isArray = True, lengthOne = False, nonNegative = True, positive = True)
checkInteger(intArray, "intArray", 
   isArray = True, lengthOne = False, nonNegative = True, positive = True)

# checkFloat
negTwoFloat = float(-2)
zeroFloat = float(0)
twoFloat = float(2)
nanFloat = np.nan
floatArrayNeg = np.array(range(0, -10, -2), dtype = float)
floatArrayNeg[4] = np.nan
floatArrayZero = np.zeros(shape = (2, 2), dtype = float)
floatArray = np.array(3.0, dtype = float)
floatArray2 = np.array([3.0, np.nan], dtype = float)

checkFloat(negTwoFloat, "negTwoFloat", isArray = False)
checkFloat(negTwoFloat, "negTwoFloat", nonNegative = True)
checkFloat(zeroFloat, "zeroFloat", nonNegative = True)
checkFloat(zeroFloat, "zeroFloat", positive = True)
checkFloat(twoFloat, "twoFloat", positive = True)
checkFloat(nanFloat, "nanFloat")
checkFloat(nanFloat, "nanFloat", noMissing = False)
checkFloat(nanFloat, "nanFloat", noMissing = False, nonNegative = True)
checkFloat(nanFloat, "nanFloat", noMissing = False, positive = True)
checkFloat(nanFloat, "nanFloat", isArray = True)

checkFloat(floatArrayNeg, "floatArrayNeg", isArray = False)
checkFloat(floatArrayNeg, "floatArrayNeg", isArray = True)
checkFloat(floatArrayNeg, "floatArrayNeg", 
 isArray = True, noMissing = False)
checkFloat(floatArrayNeg, "floatArrayNeg", 
 isArray = True, noMissing = False, lengthOne = False)
checkFloat(floatArrayNeg, "floatArrayNeg", 
 isArray = True, noMissing = False,  
 lengthOne = False, nonNegative = True)
checkFloat(floatArrayZero, "floatArrayZero", 
 isArray = True, noMissing = False,  
 lengthOne = False, nonNegative = True)
checkFloat(floatArrayZero, "floatArrayZero", 
 isArray = True, noMissing = False,  
 lengthOne = False, positive = True)
checkFloat(floatArray, "floatArray", 
 isArray = True, noMissing = False,  
 lengthOne = False, positive = True)
checkFloat(floatArray, "floatArray", 
 isArray = True, noMissing = False,  
 lengthOne = False, positive = True)
checkFloat(floatArray2, "floatArray2", 
 isArray = True, noMissing = False,  
 lengthOne = False, positive = True)
checkFloat(floatArray2, "floatArray2", 
 isArray = True, noMissing = True,  
 lengthOne = False, positive = True)

# checkUnitInterval
zeroFloat = 0.0
oneFloat = 1.0
midFloat = np.random.uniform(size = 1)
v = np.random.uniform(size = 100)
vName = "midFloat"

v = 0.0
vName = "zeroFloat"
isArray = False
noMissing = True
lengthOne = True
equal = False

checkUnitInterval(zeroFloat, "zeroFloat", isArray = False, 
  noMissing = True, lengthOne = True, equal = False)
checkUnitInterval(zeroFloat, "zeroFloat", isArray = False, 
  noMissing = True, lengthOne = True, equal = True)
checkUnitInterval(oneFloat, "oneFloat", isArray = False, 
  noMissing = True, lengthOne = True, equal = False)
checkUnitInterval(zeroFloat, "oneFloat", isArray = False, 
  noMissing = True, lengthOne = True, equal = True)
checkUnitInterval(midFloat, "midFloat", isArray = True, 
  noMissing = True, lengthOne = True, equal = False)

# invert symmetric matrix
n = 5
v = scipy.stats.invwishart.rvs(df = 21, scale = np.identity(n))
out = invertSymmetricMatrix(v)
vInv = out[0]
vInvLogDet = out[1]
np.ndarray.max(abs(v.dot(vInv) - np.identity(n)))
vInvLogDet - math.log(np.linalg.det(vInv))

#bSpline - R code below for the first check
#degree <- 3
#interiorKnots <- c(0.25, 0.5, 0.75, 0.80)
#x <- c(0, 0.01, 0.25, 0.26, 0.5, 0.51, 0.75, 0.76, 0.99, 1.0)
#boundary.knots <- c(0, 1)
#bs(x, df = degree, knots = interiorKnots)

degree = 3
interiorKnots = np.array([0.25, 0.5, 0.75, 0.80])
x = np.array([0, 0.01, 0.25, 0.26, 0.5, 0.51, 0.75, 0.76, 0.99, 1.0])
boundaryKnots = np.array([0, 1])
bSpline(x, degree, boundaryKnots, interiorKnots, intercept = False)
bSpline2(x, degree, boundaryKnots, interiorKnots, intercept = False)


x = np.random.uniform(size = int(1E5))
start = time.time() 
out1 = bSpline(x, degree, boundaryKnots, interiorKnots, intercept = False)
end = time.time() 
print(end - start)
start = time.time() 
out2 = bSpline2(x, degree, boundaryKnots, interiorKnots, intercept = False)
end = time.time() 
print(end - start)
start = time.time() 
out3 = bSpline3(x, degree, boundaryKnots, interiorKnots, intercept = False)
end = time.time() 
print(end - start)

start = time.time() 
out4 = bSpline4(x, degree, boundaryKnots, interiorKnots, intercept = False)
end = time.time() 
print(end - start)




np.max(np.abs(out3 - out2))












