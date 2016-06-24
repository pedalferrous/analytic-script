from numpy import *
from collections import *
from itertools import chain
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import os
from scipy.interpolate import interp1d

"""Script to model E-field, transmission and reflection through a number 
of layers of dielectric.  Matrix method used as described in Pettersson 
et al 1999 or Guenther appendix 4-A (beware of Guenther - some mistakes 
in the text. Pettersson is a more reliable reference.  Using P instead 
of T for the propagation matrix to avoid confusion).

User specifies stack of layers with thicknesses, indices, and labels.  
The propagation (P) and transmission (I) matrices are calculated first 
for each layer for each wavelength, layer thickness (average over 
roughness), and incidence angle.  Then quantities including total 
transmission (T), reflection (R), and E^2 (across thickness and total 
integrated) are calculated.

The parameters that are allowed to vary in this script represent external
parameters (incidence angle, wavelength) for one device.  The paramVary.py
script, on the other hand, is intended to sweep devices (eg, consider many
devices with different layer thicknesses).
"""

#Create a named tuple to hold layer information
#depth (along with wavelength) is in nm
#Index is either index constant, name of lookup file, or information
#for the Drude/Lorentz models
#active indicates whether the layer is active (whether integrated E^2
#will be calculated)
film = namedtuple('layer',['depth','index','name','active'])

def main():

  ###################################################################
  # ENTER SAVE AND OUTPUT OPTIONS HERE
  ###################################################################  

  save = False
  saveFileName = 'test'

  #options to plot T/R as function of angle or incident wavelength
  plotTRangle = False
  plotTRspectrum = False

  #option to evaluate E^2 integral in active layer and plot spectrum
  evalESQint = False
  #option to plot E^2 throughout structure
  plotESQ = True

  ###################################################################
  # ENTER INCIDENT LIGHT PARAMS AND ENVIRONMENT INDICES HERE
  ###################################################################

  #wls = linspace(1000, 7000, 1000) #wavelengths (nm)
  wls = [4000] #nm
  angles = (pi/180)*array([0.0])
  #angles = (pi/180)*linspace(0,40,10)
  #pols: 0 is s- [normal to plane of incidence], 1 is p- [parallel]
  pols = [0,1]
  #incident and final indices
  n_i = 1.0
  n_f = 1.0

  ###################################################################
  # ENTER STACK AND MATERIALS PROPERTIES HERE
  ###################################################################

  #list of the parameters for the ITO Drude model
  #these will be used by indexLookup
  #these parameters are: restricted fit, both constraints (see notes 3/05)
  ito_omega_p = 2.73649315189e+14 #s^-1
  ito_epsilon_inf = 3.24
  ito_gamma = 5.96450643459e+12 #s^-1
  itoDrudeParams = [ito_omega_p, ito_epsilon_inf, ito_gamma]

  #EXAMPLE STACKS
  """
  stack = [film(300.0, 2.4, 'HgTe dots',True),
        film(10.0, 'ni.csv','Ni',False)]
  """

  """
  stack = [
      film(840, 1.399, 'SiO2', False),
      film(50.0, itoDrudeParams, 'ITO', False),
      film(300,2.4+0.106j,'HgTe CQDs',True),
      film(50.0,'au','Au', False),
      ]
  """

  stack = [
    film(3500, 1.0, '1.0', False),
    film(10000, 0.5+0.01j, '0.5', False),
    film(3000, 0.8, '1', False)
    ]

  """
  stack = [
    film(6.0, 4+13.5j, 'NiCr approximation', False),
    film(300.0, 2.4+0.046j, 'HgTe CQD film', True)
    ]
  """

  ################################################################### 
  ###################################################################

  if save:
    os.makedirs('results/{0}'.format(saveFileName))
    saveFile = open('results/{0}/{0}.txt'.format(saveFileName), 'w')
    saveFile.write('stack:\n{0}\n'.format(stack))
    saveFile.write('wavelengths (nm): {0}\n'.format(wls))
    saveFile.write('angles (rad): {0}\n'.format(angles))
    saveFile.write('n_i (initial index): {0}\n'.format(n_i))
    saveFile.write('n_f (final index): {0}\n'.format(n_f))
    saveFile.write('polarizations (0=normal [s], 1=parallel [p]): {0}\n'.format(pols))

  #indices is n_layers x n_wavelengths
  indices = [[indexLookup(layer.index, wl) for wl in wls] for layer in stack]
  #t_angles is an (n_layers+1 x n_wavelengths x n_angles) matrix
  t_angles = snell(indices, angles, n_i, n_f)
  (I,P) = genMatrices(stack, wls, angles, n_i, n_f, indices, t_angles)

  (T,R,TAvg,RAvg,E_i) = evalMatrices(stack, wls, angles, pols,
    n_i, n_f, indices, t_angles, P, I, addBulkT=True)
  print 'TAvg = {0}'.format(TAvg)
  print 'RAvg = {0}'.format(RAvg)
  print 'RAvg + TAvg = {0}'.format(RAvg+TAvg)

  if save:
    saveFile.write('TAvg = {0}\n'.format(TAvg))
    saveFile.write('RAvg = {0}\n'.format(RAvg))
    saveFile.write('RAvg + TAvg = {0}\n'.format(RAvg+TAvg))

  if plotTRangle:
    TRAnglePlot(T,R,wls,angles,save,saveFileName)

  if plotTRspectrum:
    TRSpectrumPlot(T,R,wls,angles,save,saveFileName)

  if evalESQint:
    (ESqInt,ESqIntAvg) = ESqIntEval(stack,wls,angles,pols,indices,E_i)

    for i in range(len(stack)):
      if stack[i].active:
        print 'Average E^2 integral in active layer ({0}): {1}'.format(
          stack[i].name,ESqIntAvg[i])
        if save:
          saveFile.write('Average E^2 integral in active layer ({0}): {1}\n'.format(
          stack[i].name,ESqIntAvg[i]))

    ESqIntSpectrumPlot(ESqInt, stack, wls, angles, pols, indices, save, saveFileName)

  if plotESQ:  
    ESqPlot(E_i, stack, wls, angles, pols, indices, save, saveFileName)

"""Look up the index of refraction for a given wavelength and index look-up 
file.  If the index is given as a float, return the value.  If the index
is a tuple (of three lists, a wavelength table, one for real, one for imaginary)
then interpolate the index for the desired wavelength.  If the index is a 
list of Drude model parameters, feed them to drude() to get the index

Note on units: the save files store wavelengths in microns and the rest of the
program uses wavelengths in nm so there is a conversion factor"""
def indexLookup(index, wl):
  if isinstance(index, float) or isinstance(index, complex):
    return index
  elif isinstance(index, tuple):
    wlTable = index[0]
    nTable = index[1]
    kTable = index[2]
    nInterp = interp1d(wlTable, nTable)
    kInterp = interp1d(wlTable, kTable)
    return nInterp(wl) + (1j)*kInterp(wl)
  elif isinstance(index, list):
    #in this case the argument is parameters for the Drude model
    #in this case expect the parameters to be in order:
    #omega_p, epsilon_inf, gamma
    omega_p = index[0]
    epsilon_inf = index[1]
    gamma = index[2]
    k = (1.0e9/wl) #wavenumber in m^-1
    (_, n_r, n_i) = drude(k, wl, omega_p, epsilon_inf, gamma)
    return n_r + (1j)*n_i
  else:
    ntable = loadtxt('{0}-n.csv'.format(index), skiprows=1, unpack=True, delimiter=',')
    ktable = loadtxt('{0}-k.csv'.format(index), skiprows=1, unpack=True, delimiter=',')
    layerIndices = []
    nInterp = interp1d(1.0e3*ntable[0], ntable[1])
    kInterp = interp1d(1.0e3*ktable[0], ktable[1])
    return nInterp(wl)+(1j)*kInterp(wl)

"""Function to give the index of refraction from the Drude model
for a given omega_p (plasma frequency) and gamma (=1/tau the relaxation
time).  This function takes redundant information (wavelengths as well
as wavenumbers) since k is used to calculate while wls is used by the
indexLookup function to interpolate.  There is a term for the high-frequency
dielectric constant."""
def drude(k, wls, omega_p, epsilon_inf, gamma):
  c = 3.0e8 #m/s
  omega = c*k #expect k in m^-1 -> this will be s^-1
  n = sqrt(epsilon_inf-((omega_p**2)/(omega**2 + (1j)*omega*gamma)))
  n_r = real(n)
  n_i = imag(n)
  return (wls, n_r, n_i)

"""Given a 2D list of indices (n_layers x n_wavelengths), initial and final 
indices, and a list of angles incident on the structure, give a 3D list of 
the angles transmitted in each layer (n_layers+1 x n_wavelengths x n_angles).  
In this case the layer index refers to the stack layer with that index (eg 
t_angles[0] gives transmission angles in the first layer).  The final layer
index gives the transmission angle in the final medium (n_f).

This could also be done implicitly in the calculation of the matrices 
(see Pettersson) but I don't think there is a computational advantage, and
this way the list of angles is accessible in case it is ever useful
to look at it directly"""
def snell(indices, angles, n_i, n_f):
  t_angles = []
  #index over number of layers
  for i in range(len(indices)+1):
    t_angles.append([])
    #index over number of wavelengths
    for j in range(len(indices[0])):
      t_angles[i].append([])
      #index over number of angles
      for k in range(len(angles)):
        t_angles[i][j].append([])
        if (i == 0):
          t_angles[i][j][k] = arcsin(n_i*sin(angles[k])/indices[i][j])
        elif (i == len(indices)):
          t_angles[i][j][k] = arcsin(indices[i-1][j]*sin(
            t_angles[i-1][j][k])/n_f)
        else:
          t_angles[i][j][k] = arcsin(indices[i-1][j]*sin(
            t_angles[i-1][j][k])/indices[i][j])
  return t_angles

"""Returns 4D lists of P and I matrices indexed by (layer number, wavelength, 
incidence angle, polarization [for I matrices]).  

Because there are (n_layers + 1) I-matrices, the first index of the list of 
I matrices runs from 0 to n_layers and refers to the matrix before (entering) 
the layer to which the index corresponds, eg I[0] is for the interface from the
inital medium into the first layer.

P matrix index corresponds to the layer in which propagation occurs"""
def genMatrices(stack, wls, angles, n_i, n_f, indices, t_angles):
  I = []  
  P = []
  """
  #indices is an (n_layers x n_wavelengths) matrix
  indices = [[indexLookup(layer.index, wl) for wl in wls] for layer in stack]
  #t_angles is an (n_layers+1 x n_wavelengths x n_angles) matrix
  t_angles = snell(indices, angles, n_i, n_f)
  """
  # generate the I matrices (n_layers+1 x n_wls x n_angles x 2)
  for i in range(len(stack)+1):
    I.append([])
    for j in range(len(wls)):
      I[i].append([])
      for k in range(len(angles)):
        I[i][j].append([])
        #0 is polarization normal to the plane of incidence (s)
        #1 is polarization parallel to the plane of incidence (p)
        for l in [0,1]:
          if (i == 0):
            I[i][j][k].append(IMatrix(n_i, indices[i][j], angles[k],
              t_angles[i][j][k], l))
          elif (i == len(stack)):
            I[i][j][k].append(IMatrix(indices[i-1][j], n_f, t_angles[i-1][j][k], 
              t_angles[i][j][k], l))
          else:
            I[i][j][k].append(IMatrix(indices[i-1][j], indices[i][j], 
              t_angles[i-1][j][k], t_angles[i][j][k], l))
  
  # I = [[[[
  #   IMatrix(n_i, indices[i][j], angles[k], t_angles[i][j][k], l) if i == 0
  #   else IMatrix(indices[i-1][j], n_f, t_angles[i-1][j][k], t_angles[i][j][k], l) if i == len(stack)
  #   else IMatrix(indices[i-1][j], indices[i][j], t_angles[i-1][j][k], t_angles[i][j][k], l)
  #   for l in [0,1]] 
  #   for k in range(len(angles))] 
  #   for j in range(len(wls))] 
  #   for i in range(len(stack)+1)]


  #generate the P matrices (n_layers x n_wls x n_angles)
  P = [[[PMatrix(indices[i][j], wls[j], stack[i].depth, t_angles[i][j][k]) 
    for k in range(len(angles))] 
    for j in range(len(wls))] 
    for i in range(len(stack))]
  return (I,P)

"""method to compute matrices for transmission across interfaces"""
def IMatrix(n_1, n_2, theta_1, theta_2, pol):
  #Uses t,r consistent with Petersson et al 1999
  #The indices n_1 and n_1 input will be complex
  if (pol==0): #perpendicular (s-) polarization
    r = (((n_1*cos(theta_1))-(n_2*cos(theta_2)))/
      ((n_1*cos(theta_1))+(n_2*cos(theta_2))))
    t = ((2*n_1*cos(theta_1))/
      ((n_1*cos(theta_1))+(n_2*cos(theta_2))))
  elif (pol==1): #parallel (p-) polarization
    r = (((n_1*cos(theta_2))-(n_2*cos(theta_1)))/
      ((n_2*cos(theta_1))+(n_1*cos(theta_2))))
    t = ((2*n_1*cos(theta_1))/
      ((n_2*cos(theta_1))+(n_1*cos(theta_2))))
  return array([[1/t, r/t],[r/t, 1/t]])

"""method to compute matrices for transmission through layers
(addition of a phase from propagation)"""
def PMatrix(n, wl, d, theta): 
  delta = (2*pi*n*d*cos(theta))/(wl)
  return array([[exp((-1j)*delta),0],[0,exp(1j*delta)]])

"""Given the modeling parameters and the transmission matrices generate
the E-field distribution throughout the device and return the integral
of E^2 in each layer.  Average over all given parameters (wavelengths,
angles, polarizations [expected to be a list of the form [0] or [0,1] where
0 is normal to the plane of incidence and 1 is parallel to the plane of 
incidence])."""

def matrixCreate(startLayer, endLayer, pol, angle, wLength, P, I):
    M_partial = identity(2)
    for index in range(startLayer, endLayer):
        M_partial = dot(M_partial,dot(P[index][wLength][angle], I[index+1][wLength][angle][pol]))
    return M_partial

def evalMatrices(stack, wls, angles, pol, n_i, n_f, indices, t_angles, P, I,
  addBulkT=False):
  #Evaluate the M matrix for each layer - eg E_i[i] = M[i]E_f where
  #E_i[i] is the E-field vector after transmission into layer i but
  #before propagation.
  #M indexed by layer, wavelength, angle, polarization.



  # M = []
  # for i in range(len(stack)):
  #   M.append([])
  #   for j in range(len(wls)):
  #     M[i].append([])
  #     for k in range(len(angles)):
  #       M[i][j].append([])
  #       for l in range(len(pol)):
  #         M_partial = array([[1,0],[0,1]])
  #         # indicates propagation from i_th layer to terminal layer 
  #         # (excluding initial)
  #         for index in range(i, len(stack)):
  #           M_partial = dot(M_partial,dot(P[index][j][k], I[index+1][j][k][l]))
  #         M[i][j][k].append(M_partial)
  # #Calculate the total E_f and E_r for the stack using M_0
  # M_0 = [[[dot(I[0][j][k][l], M[0][j][k][l]) 
  #   for l in range(len(pol))] 
  #   for k in range(len(angles))] 
  #   for j in range(len(wls))]
  # #E_0 is the E-field vector before transmission into the first layer
  # # indicating no reverse incident wave (i.e. from the right)
  # # this operation yields non normalized pre-boundary field
  # E_0_noNorm = [[[dot(M_0[j][k][l],array([1,0])) 
  #   for l in range(len(pol))] 
  #   for k in range(len(angles))] 
  #   for j in range(len(wls))]

  # E_0 = [[[array([1.0, E_0_noNorm[j][k][l][1]/E_0_noNorm[j][k][l][0]]) 
  #   for l in range(len(pol))] 
  #   for k in range(len(angles))] 
  #   for j in range(len(wls))]
  # # representative of final E field, now normalized akin to incident E field
  # # solely forward propagating
  # E_f = [[[array([1.0/E_0_noNorm[j][k][l][0],0.0])
  #   for l in range(len(pol))]
  #   for k in range(len(angles))] 
  #   for j in range(len(wls))]
  # E_i = [[[[dot(M[i][j][k][l], E_f[j][k][l])
  #   for l in range(len(pol))] 
  #   for k in range(len(angles))]
  #   for j in range(len(wls))]
  #   for i in range(len(stack))]

  M = [[[[
      matrixCreate(i, len(stack),l, k, j, P, I)
      for l in range(len(pol))] 
      for k in range(len(angles))] 
      for j in range(len(wls))] 
      for i in range(len(stack))]
  
  #Calculate the total E_f and E_r for the stack using M_0
  M_0 = [[[dot(I[0][j][k][l], M[0][j][k][l]) 
      for l in range(len(pol))] 
      for k in range(len(angles))] 
      for j in range(len(wls))]

  # incident and reflected
  E_0 = [[[(1/(dot(M_0[j][k][l],array([1,0]))[0]))*dot(M_0[j][k][l],array([1,0])) 
      for l in range(len(pol))] 
      for k in range(len(angles))] 
      for j in range(len(wls))]

  # exiting field
  E_f = [[[1/(dot(M_0[j][k][l],array([1,0]))[0])*array([1.0,0.0])
      for l in range(len(pol))]
      for k in range(len(angles))] 
      for j in range(len(wls))]

  # field as indexed by layer, back-extrapolating from exiting field
  E_i = [[[[dot(M[i][j][k][l], E_f[j][k][l])
      for l in range(len(pol))] 
      for k in range(len(angles))]
      for j in range(len(wls))]
      for i in range(len(stack))]

  #overall transmittance, reflectance coefficients (following Hecht)
  #there is also the option to add a bulk transmission coefficient
  #for a thick substrate (if addBulkT = True).  In this case the bulk
  #coefficient added will be for the final index (n_f) to air (n=1).
  #this expression is again Hecht 4.68 for normal incidence
  bulkT = (4.0*n_f)/((n_f+1.0)**2)
  if addBulkT:
    T = [[[real(bulkT*(n_f/n_i)*(cos(t_angles[len(stack)][j][k])/cos(angles[k]))*abs(
      E_f[j][k][l][0])**2)
      for l in range(len(pol))]
      for k in range(len(angles))]
      for j in range(len(wls))]
  else:
    T = [[[real((n_f/n_i)*(cos(t_angles[len(stack)][j][k])/cos(angles[k]))*abs(
      E_f[j][k][l][0])**2)
      for l in range(len(pol))]
      for k in range(len(angles))]
      for j in range(len(wls))]
  R = [[[(abs(E_0[j][k][l][1])/abs(E_0[j][k][l][0]))**2
    for l in range(len(pol))]
    for k in range(len(angles))]
    for j in range(len(wls))]

  TAvg = mean(list(chain.from_iterable(chain.from_iterable(T))))
  # TAvg = (sum(T[j][k][l]
  #   for l in range(len(pol))
  #   for k in range(len(angles))
  #   for j in range(len(wls)))
  #   /(len(pol)*len(angles)*len(wls)))
  RAvg = (sum(R[j][k][l]
    for l in range(len(pol))
    for k in range(len(angles))
    for j in range(len(wls)))
    /(len(pol)*len(angles)*len(wls)))
  return (T,R,TAvg,RAvg,E_i)

"""Return the integral of Re(E)^2 in any layer given the initial E-field
vector.  This function accounts for absorption and evaluates the integral
numerically."""
def ESqIntegral(E_0, index, wl, d, theta, pol):
  return integrate.quad(lambda x: ESqEval(E_0, index, wl, x, theta, pol),
    0,d)[0]

"""Evaluate E^2 at a point (given by x) within a desired layer using the
initial E-field within that layer"""
def ESqEval(E_0, index, wl, x, theta, pol):
  E_f = E_0[0]
  E_r = E_0[1]
  delta = (2*pi*index*x*cos(theta))/wl

  #s-polarization
  if (pol == 0):
    E = E_f*exp(1j*delta) + E_r*exp((-1j)*delta)
    return abs(E)**2
  #p-polarization
  elif (pol == 1): 
    return ((sin(theta)*abs((E_r*exp((-1j)*delta))-(E_f*exp(1j*delta))))**2 
      + (cos(theta)*abs((E_r*exp((-1j)*delta))+(E_f*exp((1j)*delta))))**2)

"""Return the integral of Re(E)^2 in a layer with no absorption 
given the thickness and the initial E-field vector.  Integral is 
evaluated analytically.  This method is not currently used."""
def zeroKESqIntegral(E_0, n, wl, d, theta):
  phi = (2*pi*n*d)/(wl*cos(theta))
  #E_0 for forward (F) and reverse (R) components decomposed into real (R)
  #and imaginary (I) components to make expression easier to read
  E_0FR = real(E_0[0])
  E_0FI = imag(E_0[0])
  E_0RR = real(E_0[1])
  E_0RI = imag(E_0[1])
  A = E_0RR + E_0FR
  B = E_0RI - E_0FI
  return ((0.5*(phi**2)*(A**2 + B**2)) + (0.25*sin(2*phi)*(A**2 - B**2))
   - (0.5*A*B*cos(2*phi)))

"""Return an (n_layers x n_wls x n_angles x n_polarizations) list of the
integral of E^2 in each layer evaluated analytically using the ESqIntegral
method along with the average of the integral over wavelengths, angles,
and polarizations.  Use the initial E-field vectors E_i from evalMatrices()"""
def ESqIntEval(stack, wls, angles, pols, indices, E_i):
  ESqInt = [[[[ESqIntegral(E_i[i][j][k][l], indices[i][j], 
      wls[j], stack[i].depth, angles[k], pols[l])
    for l in range(len(pols))] 
    for k in range(len(angles))]
    for j in range(len(wls))]
    for i in range(len(stack))]
  ESqIntAvg = [sum(ESqInt[i][j][k][l] 
    for j in range(len(wls))
    for k in range(len(angles))
    for l in range(len(pols)))/(len(wls)*len(angles)*len(pols))
    for i in range(len(stack))]
  return (ESqInt, ESqIntAvg)

"""Plot T and R as a function of angle.
Averages over wls are n_angles x n_pols"""
def TRAnglePlot(T,R,wls,angles,save,saveFileName):
  RsWlAvg = [sum(R[j][k][0] 
    for j in range(len(wls)))/len(wls)
    for k in range(len(angles))]
  RpWlAvg = [sum(R[j][k][1] 
    for j in range(len(wls)))/len(wls)
    for k in range(len(angles))]
  TsWlAvg = [sum(T[j][k][0] 
    for j in range(len(wls)))/len(wls)
    for k in range(len(angles))]
  TpWlAvg = [sum(T[j][k][1] 
    for j in range(len(wls)))/len(wls)
    for k in range(len(angles))]
  anglesdeg = (180/pi)*angles
  ax = plt.axes()
  ax.plot(anglesdeg, RsWlAvg, label='R_s')
  ax.plot(anglesdeg, RpWlAvg, label='R_p')
  ax.set_xlabel('Incident angle (deg)')
  ax.legend()
  if save:
    plt.savefig('results/{0}/{0}-RAnglePlot.pdf'.format(saveFileName))
  plt.show()
  ax = plt.axes()
  ax.plot(anglesdeg, TsWlAvg, label='T_s')
  ax.plot(anglesdeg, TpWlAvg, label='T_p')
  ax.set_xlabel('Incident angle (deg)', fontsize=18)
  plt.xticks(fontsize=15)
  plt.yticks(fontsize=15)
  ax.legend()
  if save:
    plt.savefig('results/{0}/{0}-TAnglePlot.pdf'.format(saveFileName))
  plt.show()

"""Plot T and R as a function of wavenumber k (in m).
Averages over angles"""
def TRSpectrumPlot(T,R,wls,angles,save,saveFileName,tPlot=True,
  rPlot=True):
  RsAnglAvg = [sum(R[j][k][0]
    for k in range(len(angles)))/len(angles)
    for j in range(len(wls))]
  RpAnglAvg = [sum(R[j][k][1]
    for k in range(len(angles)))/len(angles)
    for j in range(len(wls))]
  TsAnglAvg = [sum(T[j][k][0]
    for k in range(len(angles)))/len(angles)
    for j in range(len(wls))]
  TpAnglAvg = [sum(T[j][k][1]
    for k in range(len(angles)))/len(angles)
    for j in range(len(wls))]
  RAnglPolAvg = [0.5*(RsAnglAvg[j]+RpAnglAvg[j])
    for j in range(len(wls))]
  TAnglPolAvg = [0.5*(TsAnglAvg[j]+TpAnglAvg[j])
    for j in range(len(wls))]
  k_0 = (1.0/wls) * (1.0e9) #conversion included from nm -> m

  ax = plt.axes()
  if tPlot:
    ax.plot(wls, TsAnglAvg, label='T_s')
    ax.plot(wls, TpAnglAvg, label='T_p')
    ax.plot(wls, TAnglPolAvg, label='T_avg')
  if rPlot:
    ax.plot(wls, RsAnglAvg, label='R_s')
    ax.plot(wls, RpAnglAvg, label='R_p')
    ax.plot(wls, RAnglPolAvg, label='R_avg')
  ax.legend()
  ax.set_xlabel('Wavelength (nm)', fontsize=18)
  plt.xticks(fontsize=15)
  plt.yticks(fontsize=15)
  if save:
    plt.savefig('results/{0}/{0}-TRWlsSpectrum.pdf'.format(saveFileName))
  plt.show()

  ax = plt.axes()
  if tPlot:
    ax.plot(0.01*k_0, TsAnglAvg, label='T_s')
    ax.plot(0.01*k_0, TpAnglAvg, label='T_p')
    ax.plot(0.01*k_0, TAnglPolAvg, label='T_avg')
  if rPlot:
    ax.plot(0.01*k_0, RsAnglAvg, label='R_s')
    ax.plot(0.01*k_0, RpAnglAvg, label='R_p')
    ax.plot(0.01*k_0, RAnglPolAvg, label='R_avg')
  ax.legend(loc='upper left')
  ax.set_xlabel('k (cm^-1)', fontsize=18)
  ax.set_ylim([0.0,1.0])
  plt.xticks(fontsize=15)
  plt.yticks(fontsize=15)
  if save:
    plt.savefig('results/{0}/{0}-TRWaveNumsSpectrum.pdf'.format(saveFileName))
  plt.show()

  return (k_0, RAnglPolAvg, TAnglPolAvg)

"""Plot spectrum of the E^2 integral"""
def ESqIntSpectrumPlot(ESqInt, stack, wls, angles, pols, indices, save, saveFileName):
  ax = plt.axes()

  wnums = (1.0e7/wls) #wavenumbers in cm^-1 to be used in plotting
  print "hello"

  for i in range(len(stack)):
    if stack[i].active:
      ESqIntWls = [sum(ESqInt[i][j][k][l] 
        for k in range(len(angles))
        for l in range(len(pols)))/(len(angles)*len(pols))
        for j in range(len(wls))]
      ax.plot(wnums, ESqIntWls, label=stack[i].name)

  ax.legend()
  ax.set_xlabel('Wavenumber (cm$^-1$)', fontsize=18)
  ax.set_ylabel('$E^2$ integral (arb. units)', fontsize=18)
  plt.xticks(fontsize=15)
  plt.yticks(fontsize=15)
  if save:
    plt.savefig('results/{0}/{0}-ESqIntSpectrum.pdf'.format(saveFileName))
  plt.show()

"""Plot E^2 throughout the device"""
def ESqPlot(E_i, stack, wls, angles, pol, indices, save, saveFileName, pointsPerLayer = 500):
  #Evaluate E^2 across each layer
  #ESqVals is (n_layers x n_wls x n_angles x n_pols x n_pointsPerLayer)
  #xPoints is (n_layers x n_PointsPerLayer)
  #stackXPos is n_layers - start position of each layer
  #phasePoints is n_layers x n_wls x n_angles x n_pointsPerLayer
  stackXPos = [0.0]+[sum([stack[index].depth 
    for index in range(i)])
    for i in range(1, len(stack))]
  xPoints = [[stackXPos[i] + stack[i].depth*l
    for l in linspace(0,1.0,pointsPerLayer)]
    for i in range(len(stack))]
  xEvalPoints = [[stack[i].depth*l
    for l in linspace(0,1.0,pointsPerLayer)]
    for i in range(len(stack))]
  ESq = [[[[[ESqEval(E_i[i][j][k][l],indices[i][j],wls[j],xEvalPoints[i][n],angles[k], pol[l])
    for n in range(pointsPerLayer)] 
    for l in range(len(pol))]
    for k in range(len(angles))]
    for j in range(len(wls))]
    for i in range(len(stack))]
  ESqAvg = [[sum(ESq[i][j][k][l][n]
    for j in range(len(wls))
    for k in range(len(angles))
    for l in range(len(pol)))/(len(pol)*len(angles)*len(wls))
    for n in range(pointsPerLayer)]
    for i in range(len(stack))]
  ax = plt.axes()
  for i in range(len(stack)):
    ax.plot(xPoints[i],ESqAvg[i],label=stack[i].name)
    ax.legend()
  if save:
    plt.savefig('results/{0}/{0}-ESqPlot.pdf'.format(saveFileName))
  ax.set_xlabel('x-position (nm)', fontsize=18)
  ax.set_ylabel('E^2 (arbitrary units)', fontsize=18)
  plt.xticks(fontsize=15)
  plt.yticks(fontsize=15)
  ax.grid(True)
  plt.show()

if __name__ == '__main__':
  main()