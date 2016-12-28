# -*- coding: utf-8 -*-
"""
Created on Wed Oct  5 17:37:30 2016

@author: luksmi

"""
import pandas as pd
import numpy as np
import scipy
import math
import time
from functools import partial

def checkInteger(v, vName, isArray = False, lengthOne = True,
  nonNegative = False, positive = False):
  """Checks if argument v is an integer.  
  Note that np.nan is considered a float, so nan cannot be in an integer array.
  Args:
    v (int):
    vName (string): the name of v
    isArray (boolean): if v should be an array
    lengthOne (boolean): should v be of length one?
    nonNegative (boolean): should v be nonNegative?
    positive (boolean): is v restricted to be positive?
  Returns nothing.
  """
  # Error handling:  
  errorString = 'in checkInteger vName must be a float.'
  assert isinstance(vName, str), errorString
  errorString = 'in checkInteger isArray must be boolean.'
  assert isinstance(isArray, bool), errorString
  errorString = 'in checkInteger lengthOne must be boolean.'
  assert isinstance(lengthOne, bool), errorString
  errorString = 'in checkInteger nonNegative must be boolean.'
  assert isinstance(nonNegative, bool), errorString
  errorString = 'in checkInteger positive must be boolean.'
  assert isinstance(positive, bool), errorString
  if not isArray:
    #assumes v is of length 1    
    errorString = vName + ' must contain integers.'  
    assert isinstance(v, int), errorString
    if nonNegative:
      errorString = vName + ' cannot be negative.'
      assert v >= 0, errorString
    if positive:
      errorString = vName + ' cannot be nonNegative.'
      assert v > 0, errorString
  if isArray:
    errorString = vName + ' must be a numpy array.'
    assert isinstance(v, np.ndarray), errorString
    errorString = vName + ' must contain integers.'
    assert v.dtype == 'int64' or v.dtype == 'int32', errorString
    if lengthOne:
      errorString = vName + ' must be of length 1'
      assert v.size == 1, errorString
    if nonNegative:
      errorString = vName + ' cannot contain negative values.'
      assert v[v < 0].size == 0, errorString
    if positive:
      errorString = vName + ' cannot contain nonNegative values.'
      assert v[v <= 0].size == 0, errorString

def checkFloat(v, vName, isArray = False, noMissing = True,
  lengthOne = True,  nonNegative = False, positive = False):
  """Checks if argument v is a float.  
  Note that np.nan is considered a float.
  Args:
    v (int):
    vName (string): the name of v
    isArray (boolean): if v should be an array
    noMissing (boolean): if v should not contain missing
    lengthOne (boolean): should v be of length one?
    nonNegative (boolean): should v be nonNegative?
    positive (boolean): is v restricted to be positive?
  Returns nothing.
  """
  # Error handling:  
  errorString = 'in checkFloat vName must be a float.'
  assert isinstance(vName, str), errorString
  errorString = 'in checkFloat isArray must be boolean.'
  assert isinstance(isArray, bool), errorString
  errorString = 'in checkFloat noMissing must be boolean.'
  assert isinstance(noMissing, bool), errorString
  errorString = 'in checkFloat lengthOne must be boolean.'
  assert isinstance(lengthOne, bool), errorString
  errorString = 'in checkFloat nonNegative must be boolean.'
  assert isinstance(nonNegative, bool), errorString
  errorString = 'in checkFloat positive must be boolean.'
  assert isinstance(positive, bool), errorString
  if not isArray:
    #assumes v is of length 1    
    errorString = vName + ' must contain floats.'  
    assert isinstance(v, float), errorString
    if noMissing:
      errorString = vName + ' cannot be missing.'
      assert not np.isnan(v), errorString
    if nonNegative:
      errorString = vName + ' cannot be negative.'
      assert v >= 0.0, errorString
    if positive:
      errorString = vName + ' cannot be nonNegative.'
      assert v > 0.0, errorString
  if isArray:
    errorString = vName + ' must be a numpy array.'
    assert isinstance(v, np.ndarray), errorString
    errorString = vName + ' must contain floats.'
    assert v.dtype == 'float64' or v.dtype == 'float32', errorString
    if noMissing:
      errorString = vName + ' cannot contain missing values.'
      assert not np.isnan(v).any(), errorString
    if lengthOne:
      errorString = vName + ' must be of length 1'
      assert v.size == 1, errorString
    if nonNegative:
      errorString = vName + ' cannot contain negative values.'
      assert v[v < 0.0].size == 0, errorString
    if positive:
      errorString = vName + ' cannot contain nonNegative values.'
      assert v[v <= 0.0].size == 0, errorString

def checkUnitInterval(v, vName, isArray = False, noMissing = True,
  lengthOne = True, equal = False):
  """Checks if argument v is in (0, 1) or [0, 1].  
  Args:
    v (int):
    vName (string): the name of v
    isArray (boolean): if v should be an array
    noMissing (boolean): if v should not contain missing
    lengthOne (boolean): should v be of length one?
    equal (boolean): [0, 1] is interval of interest
  Returns nothing.
  """
  # Error handling:  
  errorString = 'in checkFloat vName must be a float.'
  assert isinstance(vName, str), errorString
  errorString = 'in checkFloat isArray must be boolean.'
  assert isinstance(isArray, bool), errorString
  errorString = 'in checkFloat noMissing must be boolean.'
  assert isinstance(noMissing, bool), errorString
  errorString = 'in checkFloat lengthOne must be boolean.'
  assert isinstance(lengthOne, bool), errorString
  errorString = 'in checkFloat equal must be boolean.'
  assert isinstance(equal, bool), errorString
  if equal: positive = False
  else: positive = True
  checkFloat(v, vName, isArray, noMissing,
    lengthOne, True, positive)
  if equal:
    errorString = vName + ' must contain elements less than or equal to 1.0.'
    if isArray:
      assert v[v > 1.0].size == 0, errorString
    else:
      errorString = vName + ' must be less than or equal to 1.0.'
      assert v <= 1.0, errorString
  else:
    if isArray:
      errorString = vName + ' must contain elements less than 1.0.'
      assert v[v >= 1.0].size == 0, errorString
    else:
      errorString = vName + ' must be less than 1.0.'
      assert v < 1.0, errorString

def invertSymmetricMatrix(v):
  """Inverts symmetric matrix v.  
  Args:
    v : symmetric invertible matrix. 
  Returns:
    vInv: inverse of v.
    vInvLogDet: log determinant of vInv.
  """
  # Error handling:  
  assert isinstance(v, np.ndarray), 'v must be a numpy array.'
  vChol = scipy.linalg.cholesky(v).T
  n = v.shape[0]
  vInvLogDet = sum(math.log(vChol[i][i]) for i in range(n))
  vInvLogDet *= -2.0
  vCholInv = scipy.linalg.solve_triangular(vChol.T, b = np.identity(n))
  vInv = vCholInv.dot(vCholInv.T)
  return([vInv, vInvLogDet])


def makeAbsoluteDistance(timepoints):
  """Creates a matrix of absolute distances.  
  Args:
    timepoints (np.ndarray): vector of length n.
  Returns d, the n x n distance matrix.
  """    
  assert isinstance(timepoints, np.ndarray), 'timepoints must be a numpy array.'  
  assert len(timepoints.shape) == 1, 'timepoints must be a 1-dimensional array.'
  n = timepoints.shape[0]
  k = len(timepoints)
  d1 = np.zeros(shape = (k, k))
  d2 = np.zeros(shape = (k, k))
  for i in range(k):
    d1[:, i] = timepoints
    d2[i, :] = timepoints
  d = abs(d1 - d2)
  return d 




def bSpline(x, degree = 3, boundaryKnots = [0, 1], interiorKnots = None, intercept = False):
  """bSpline in loop.  
  Args:
    x: points at which to evaluate the spline basis.
    degree: spline degree (default cubic).
    boundaryKnots: 
    noMissing (boolean): if v should not contain missing
    lengthOne (boolean): should v be of length one?
    equal (boolean): [0, 1] is interval of interest
  Returns nothing.
  """
  # Error handling:  
  assert isinstance(x, np.ndarray), 'x must be a numpy array.'
  assert len(x.shape) == 1, 'if x has length greater than 1, x be a vector.'
  checkInteger(degree, 'degree', positive = True)  
  assert isinstance(boundaryKnots, np.ndarray), 'boundaryKnots must be a numpy array.'
  assert len(boundaryKnots) == 2, 'boundaryKnots must be of length 2.'  
  boundaryKnots = sorted(boundaryKnots)
  if interiorKnots is None:
    lenInteriorKnots = 0
    knots = boundaryKnots
  else:
    assert isinstance(interiorKnots, np.ndarray), 'if not empty, interiorKnots must be a numpy array.'
    interiorKnots = sorted(interiorKnots)
    lenInteriorKnots = len(np.atleast_1d(interiorKnots))
    assert np.min(interiorKnots) > boundaryKnots[0], 'boundaryKnots[0] must be smaller than all interiorKnots.'    
    assert np.max(interiorKnots) < boundaryKnots[1], 'boundaryKnots[1] must be larger than all interiorKnots.'    
    knots = np.concatenate([np.repeat(boundaryKnots[0], degree + 1), interiorKnots, np.repeat(boundaryKnots[1], degree + 1)], axis = 0)
  errorString = 'in bSplineMatrix intercept must be boolean.'
  assert isinstance(intercept, bool), errorString 
  #create function that evaluates bSpline for a scalar x  
  def bSplineScalar(x, degree, i, knots):
    if(degree == 0):
      if (knots[i] <= x) & (x < knots[i + 1]):
        v = 1.0
      else:
        v = 0.0
    else:
      if((knots[degree + i] - knots[i]) == 0.0):
        alpha1 = 0.0
      else:
        alpha1 = (x - knots[i]) / (knots[degree + i] - knots[i])
      if((knots[degree + i + 1] - knots[i + 1]) == 0):
        alpha2 = 0.0
      else:
        alpha2 = (knots[degree + i + 1] - x) / (knots[i + degree + 1] - knots[i + 1])      
      v = alpha1 * bSplineScalar(x, degree - 1, i, knots) + alpha2 * bSplineScalar(x, degree - 1, i + 1, knots)
    return v
  nColumnsB = lenInteriorKnots + degree + 1
  xLength = len(np.atleast_1d(x))
  vMat = np.zeros(shape = (xLength, nColumnsB))
  for i in range(nColumnsB):
    for j in range(xLength):    
      vMat[j, i] = bSplineScalar(x[j], degree, i, knots)      
  vMat[x == boundaryKnots[1], nColumnsB - 1] = 1.0
  if not intercept:
    vMat = vMat[:, 1:nColumnsB]
  return vMat

def bSpline2(x, degree = 3, boundaryKnots = [0, 1], interiorKnots = None, intercept = False):
  """#vectorized bSpline.  
  Args:
    x: points at which to evaluate the spline basis.
    degree: spline degree (default cubic).
    boundaryKnots: 
    noMissing (boolean): if v should not contain missing
    lengthOne (boolean): should v be of length one?
    equal (boolean): [0, 1] is interval of interest
  Returns nothing.
  """
  # Error handling:  
  assert isinstance(x, np.ndarray), 'x must be a numpy array.'
  assert len(x.shape) == 1, 'if x has length greater than 1, x be a vector.'
  checkInteger(degree, 'degree', positive = True)  
  assert isinstance(boundaryKnots, np.ndarray), 'boundaryKnots must be a numpy array.'
  assert len(boundaryKnots) == 2, 'boundaryKnots must be of length 2.'  
  boundaryKnots = sorted(boundaryKnots)
  if interiorKnots is None:
    lenInteriorKnots = 0
    knots = boundaryKnots
  else:
    assert isinstance(interiorKnots, np.ndarray), 'if not empty, interiorKnots must be a numpy array.'
    interiorKnots = sorted(interiorKnots)
    lenInteriorKnots = len(np.atleast_1d(interiorKnots))
    assert np.min(interiorKnots) > boundaryKnots[0], 'boundaryKnots[0] must be smaller than all interiorKnots.'    
    assert np.max(interiorKnots) < boundaryKnots[1], 'boundaryKnots[1] must be larger than all interiorKnots.'    
    knots = np.concatenate([np.repeat(boundaryKnots[0], degree + 1), interiorKnots, np.repeat(boundaryKnots[1], degree + 1)], axis = 0)
  errorString = 'in bSplineMatrix intercept must be boolean.'
  assert isinstance(intercept, bool), errorString 
  #create function that evaluates bSpline for a scalar x  
  def bSplineScalar(x, degree, i, knots):
    if(degree == 0):
      if (knots[i] <= x) & (x < knots[i + 1]):
        v = 1.0
      else:
        v = 0.0
    else:
      if((knots[degree + i] - knots[i]) == 0.0):
        alpha1 = 0.0
      else:
        alpha1 = (x - knots[i]) / (knots[degree + i] - knots[i])
      if((knots[degree + i + 1] - knots[i + 1]) == 0):
        alpha2 = 0.0
      else:
        alpha2 = (knots[degree + i + 1] - x) / (knots[i + degree + 1] - knots[i + 1])      
      v = alpha1 * bSplineScalar(x, degree - 1, i, knots) + alpha2 * bSplineScalar(x, degree - 1, i + 1, knots)
    return v
  nColumnsB = lenInteriorKnots + degree + 1
  vMat = np.zeros(shape = (len(np.atleast_1d(x)), nColumnsB))
  for i in range(nColumnsB):
    bSplineVector = partial(bSplineScalar, degree = degree, i = i, knots = knots)  
    bSplineVector2 = np.vectorize(bSplineVector)
    vMat[:, i] = bSplineVector2(x)
  vMat[x == boundaryKnots[1], nColumnsB - 1] = 1.0
  if not intercept:
    vMat = vMat[:, 1:nColumnsB]
  return vMat


def bSpline3(x, degree = 3, boundaryKnots = [0, 1], interiorKnots = None, intercept = False):
  """List comprehension version.  
  Args:
    x: points at which to evaluate the spline basis.
    degree: spline degree (default cubic).
    boundaryKnots: 
    noMissing (boolean): if v should not contain missing
    lengthOne (boolean): should v be of length one?
    equal (boolean): [0, 1] is interval of interest
  Returns nothing.
  """
  # Error handling:  
  assert isinstance(x, np.ndarray), 'x must be a numpy array.'
  assert len(x.shape) == 1, 'if x has length greater than 1, x be a vector.'
  checkInteger(degree, 'degree', positive = True)  
  assert isinstance(boundaryKnots, np.ndarray), 'boundaryKnots must be a numpy array.'
  assert len(boundaryKnots) == 2, 'boundaryKnots must be of length 2.'  
  boundaryKnots = sorted(boundaryKnots)
  if interiorKnots is None:
    lenInteriorKnots = 0
    knots = boundaryKnots
  else:
    assert isinstance(interiorKnots, np.ndarray), 'if not empty, interiorKnots must be a numpy array.'
    interiorKnots = sorted(interiorKnots)
    lenInteriorKnots = len(np.atleast_1d(interiorKnots))
    assert np.min(interiorKnots) > boundaryKnots[0], 'boundaryKnots[0] must be smaller than all interiorKnots.'    
    assert np.max(interiorKnots) < boundaryKnots[1], 'boundaryKnots[1] must be larger than all interiorKnots.'    
    knots = np.concatenate([np.repeat(boundaryKnots[0], degree + 1), interiorKnots, np.repeat(boundaryKnots[1], degree + 1)], axis = 0)
  errorString = 'in bSplineMatrix intercept must be boolean.'
  assert isinstance(intercept, bool), errorString 
  #create function that evaluates bSpline for a scalar x  
  def bSplineScalar(x, degree, i, knots):
    if(degree == 0):
      if (knots[i] <= x) & (x < knots[i + 1]):
        v = 1.0
      else:
        v = 0.0
    else:
      if((knots[degree + i] - knots[i]) == 0.0):
        alpha1 = 0.0
      else:
        alpha1 = (x - knots[i]) / (knots[degree + i] - knots[i])
      if((knots[degree + i + 1] - knots[i + 1]) == 0):
        alpha2 = 0.0
      else:
        alpha2 = (knots[degree + i + 1] - x) / (knots[i + degree + 1] - knots[i + 1])      
      v = alpha1 * bSplineScalar(x, degree - 1, i, knots) + alpha2 * bSplineScalar(x, degree - 1, i + 1, knots)
    return v
  nColumnsB = lenInteriorKnots + degree + 1
  xLength = len(np.atleast_1d(x))
  vMat = np.zeros(shape = (xLength, nColumnsB))
  for i in range(nColumnsB):
    vMat[:, i] = [bSplineScalar(y, degree, i, knots) for y in x]      
  vMat[x == boundaryKnots[1], nColumnsB - 1] = 1.0
  if not intercept:
    vMat = vMat[:, 1:nColumnsB]
  return vMat

