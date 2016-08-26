###########################################################
# An updated and revised version of analytic wave propagator.
# Zane Rossi, University of Chicago, 2016.
# Based on code by John Roberts in attempt to minimize and
# steamline code for optimization and debugging
###########################################################

from numpy import *
from collections import *
from itertools import chain
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import os
from scipy.interpolate import interp1d
import math as math

###########################################################
# Primary data structures
###########################################################

# named tuple to hold layer information
# depth (and wavelength) in nm
# index is either index constant, name of lookup file, or information
# for Drude/Lorentz models
# active indicates whether (whether integrated E^2 is calculated
film = namedtuple('layer',['depth','index','name','active'])

# WARNING: activating this will not mesh with current
# optimized normFactor in normPwr and ODPlot
peakWavelength = 4000
blackbody = False

###########################################################
# Main loop
###########################################################

def main():

    ###########################################################
    # Primary options and input parameters
    ###########################################################

    """ENTER SAVE AND OUTPUT OPTIONS HERE"""

    saveFileName = 'test'
    save = False

    # options to plot T/R as function of angle or incident wavelength
    plotTRangle = False
    plotTRspectrum = False

    # option to evaluate E^2 integral in active layer and plot spectrum
    evalESQint = False

    # option to plot normalized power absorbed in active layer across spectrum
    # as well as optical density (axes with method defined limits)
    normPwr = False
    ODPlot  = True

    # option to plot E^2 throughout structure
    plotESQ = False

    """ENTER INCIDENT LIGHT PARAMS AND ENVIRONMENT INDICES HERE"""

    # measurements in nm
    # note that use of linespace is required currently for spectrum plots
    # this can be circumvented.
    wls = linspace(2000, 5000, 500)
    # wls = linspace(4000, 4000, 1)
    # currently unnacepting of angle variation
    # angle represents CCW rotation from direction of propagation
    angles = (pi/180)*array([0.0])
    # pols: 0 is s- [normal to plane of incidence], 1 is p- [parallel]
    pols = [0,1]
    # incident and final indices (currently free space)
    # n_i = 3.43 given incidence through silicon
    n_i = 3.43 # indicative of CaF2 in the average 
    n_f = 1.0

    """ENTER STACK AND MATERIAL PROPERTIES HERE"""

    # list of the parameters for the ITO Drude model
    # these will be used by indexLookup
    # these parameters are: restricted fit, both constraints (see notes 3/05)
    ito_omega_p = 2.73649315189e+14 #s^-1
    ito_epsilon_inf = 3.24
    ito_gamma = 5.96450643459e+12 #s^-1
    itoDrudeParams = [ito_omega_p, ito_epsilon_inf, ito_gamma]

    # list of parameters for erf bandgap model HgTe
    base_index = 2.4  # dimensionless
    nu_0       = 2200 # wavenumber
    delta_nu   = 300  # wavenumber
    erf_amp    = 0.11  # cm^-1
    erfParams  = [base_index, nu_0, delta_nu, erf_amp]

    """ENTER MAIN STACK HERE"""

    # computer optimized stack (courtesy of autoOptimize)
    # stack = [
    #         film(95, 1.399, 'SiO2', True),
    #         film(130, itoDrudeParams, 'ITO', False),
    #         film(640, 2.4 + 0.106j, 'CQD', True),
    #         film(200, 'au', 'Gold', False)
    #         ]
    # Average E^2 integral in active layer (CQD): 
    # [82.01202593364928, 245.38468092862215, 3504.612761494427, 0.9749937541718504]

    # original hand-optimized stack
    # stack = [
    #         film(745, 1.399, 'SiO2', True),
    #         film(10, itoDrudeParams, 'ITO', False),
    #         film(410, 2.4 + 0.106j, 'CQD', True),
    #         film(50, 'au', 'Gold', False)
    #         ]
    # Average E^2 integral in active layer (CQD): 

    # model of the old device using a gold back contact
    # stack = [
    #        film(6.0, 'NiCr', 'NiCr approximation', False),
    #        film(500.0, erfParams, 'CQD', True),
    #        film(150.0, 'au', 'Gold', False)
    #        ]

    # for the determination of optical density
    stack = [
           film(400.0, erfParams, 'CQD', True),
           ]

    # input from experimental data
    # stack = [
    #         # film(200, 1.399, 'SiO2', False),
    #         # film(50, itoDrudeParams, 'ITO', False),
    #         film(10000, 3.43+0.4j, 'CQD', True),
    #         film(150, 'au', 'Gold', False)
    #         ]

    ###########################################################
    # Main processes and evaluation
    ###########################################################

    # save procedure for file storage
    if save:
        os.makedirs('results/{0}'.format(saveFileName))
        saveFile = open('results/{0}/{0}.txt'.format(saveFileName), 'w')
        saveFile.write('stack:\n{0}\n'.format(stack))
        saveFile.write('wavelengths (nm): {0}\n'.format(wls))
        saveFile.write('angles (rad): {0}\n'.format(angles))
        saveFile.write('n_i (initial index): {0}\n'.format(n_i))
        saveFile.write('n_f (final index): {0}\n'.format(n_f))
        saveFile.write('polarizations (0=normal [s], 1=parallel [p]): {0}\n'.format(pols))

    # for each wavelength and each film
    indices = [[indexLookup(layer.index, wl) for wl in wls] for layer in stack]
    # snell's law applied for incident angle in all layers
    # traversal over films, wavelength, and (of course) angle
    t_angles = snell(indices, angles, n_i, n_f)
    (I,P) = genMatrices(stack, wls, angles, n_i, n_f, indices, t_angles)
    
    """The goal of this method is to break it up into two or more blocks
    of return. Namely, that it calculate T and R deriavtes separately from
    E_i, which is arguably important enough in its own right.
    """

    # field calculations
    (E_0, E_f, E_i) = evalField(stack, wls, angles, pols,
    n_i, n_f, indices, t_angles, P, I)

    # transmission and reflection calculations
    (T, R, TAvg, RAvg) = evalTR(stack, E_f, E_0, angles, n_i, n_f, 
        t_angles, addBulkT=True)

    # useful scalar quantities
    print 'TAvg = {0}'.format(TAvg)
    print 'RAvg = {0}'.format(RAvg)
    print 'RAvg + TAvg = {0}'.format(RAvg+TAvg)

    if save:
        saveFile.write('TAvg = {0}\n'.format(TAvg))
        saveFile.write('RAvg = {0}\n'.format(RAvg))
        saveFile.write('RAvg + TAvg = {0}\n'.format(RAvg+TAvg))

    if plotTRangle:
        TRAnglePlot(T,R,wls,angles,save,saveFileName)

    # do not use with standard list wls
    if plotTRspectrum:
        TRSpectrumPlot(T,R,wls,angles,save,saveFileName)

    # do not use with standard list wls
    if evalESQint:
        (ESqInt,ESqIntAvg, PAbsd, PAbsdAvg) = ESqIntEval(stack,wls,angles,pols,indices,E_i)
        PAbsdTot = 0.0
        PAbsdActive = 0.0

        for i in range(len(stack)):
            PAbsdTot += PAbsdAvg[i]

            if stack[i].active:
                print 'Average E^2 integral in active layer ({0}): {1}'.format(
                    stack[i].name,ESqIntAvg[i])
                print 'Average absorbed power in active layer ({0}): {1}'.format(
                    stack[i].name,PAbsdAvg[i])
                PAbsdActive += PAbsdAvg[i]
            if save:
                saveFile.write('Average E^2 integral in active layer ({0}): {1}\n'.format(
                    stack[i].name,ESqIntAvg[i]))
                saveFile.write('Average absorbed power in active layer ({0}): {1}\n'.format(
                    stack[i].name,PAbsdAvg[i]))

        # currently this is rather meaningless without blackbody weighting
        print 'EQE/IQE = {0}'.format((PAbsdActive/PAbsdTot)*(1-TAvg-RAvg))
        ESqIntSpectrumPlot(ESqInt, PAbsd, stack, wls, angles, pols, indices, save, saveFileName)

    if normPwr:
        (ESqInt,ESqIntAvg, PAbsd, PAbsdAvg) = ESqIntEval(stack,wls,angles,pols,indices,E_i)
        normPwrPlot(PAbsd, stack, wls, angles, pols, n_i, n_f, indices, save, saveFileName)

    if ODPlot:
        (ESqInt,ESqIntAvg, PAbsd, PAbsdAvg) = ESqIntEval(stack,wls,angles,pols,indices,E_i)
        ODNormPlot(PAbsd, stack, wls, angles, pols, n_i, n_f, indices, save, saveFileName)

    if plotESQ:
        ESqPlot(E_i, stack, wls, angles, pols, indices, save, saveFileName)

"""Look up the index of refraction for a given wavelength reference file.
May be real, complex, or interpolated from list by wavelength. 
If the index is a list of Drude model parameters, feed them to drude() 
to get the index.
Note on units: the save files store wavelengths in microns so there is a conversion factor"""
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
        # for the modeling of HgTe extinction coefficient via
        # bandgap model, in order:base_index, nu_0, delta_nu, erf_amp
        if len(index) == 4:
            base_index = index[0]
            nu_0       = index[1]
            delta_nu   = index[2]
            erf_amp    = index[3]
            wavenum    = 1.0e7/wl
            if wavenum <= nu_0:     # below bandgap (possible smoothing needed)
                return base_index
            else:                   # above bandgap
                index = base_index + (1j)*(erf_amp*math.erf((1.0e7/wl - nu_0)/delta_nu))
                return index
        else:    
            # in this case the argument is parameters for the Drude model
            # in this case expect the parameters to be in order:
            # omega_p, epsilon_inf, gamma
            omega_p = index[0]
            epsilon_inf = index[1]
            gamma = index[2]
            k = (1.0e9/wl) #wavenumber in m^-1
            (_, n_r, n_i) = drude(k, wl, omega_p, epsilon_inf, gamma)
            return n_r + (1j)*n_i
    else:
        # call based on index name of film
        ntable = loadtxt('{0}-n.csv'.format(index), skiprows=1, unpack=True, delimiter=',')
        ktable = loadtxt('{0}-k.csv'.format(index), skiprows=1, unpack=True, delimiter=',')
        layerIndices = []
        # conversion from microns to nm
        nInterp = interp1d(1.0e3*ntable[0], ntable[1])
        kInterp = interp1d(1.0e3*ktable[0], ktable[1])
        return nInterp(wl)+(1j)*kInterp(wl)

"""Function to give the index of refraction from the Drude model
for a given omega_p (plasma frequency) and gamma (=1/tau the relaxation
time). This function takes redundant information (wavelengths as well
as wavenumbers) since k is used to calculate while wls is used by the
indexLookup function to interpolate. There is a term for the high-frequency
dielectric constant."""
def drude(k, wls, omega_p, epsilon_inf, gamma):
    # speed of light in m/s
    c = 3.0e8
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
                    # print "Initial to first layer: " + str(n_i) + " to " + str(indices[i][j])
                    # print "incident: " + str(180/pi*angles[k]) + " refracted: " + str(180/pi*t_angles[i][j][k])
                elif (i == len(indices)):
                    t_angles[i][j][k] = arcsin(indices[i-1][j]*sin(
                        t_angles[i-1][j][k])/n_f)
                    # print "Interior to final layer: " + str(indices[i-1][j]) + " to " + str(n_f)
                    # print "incident: " + str(180/pi*t_angles[i-1][j][k]) + " refracted: inconsequential" 
                else:
                    t_angles[i][j][k] = arcsin(indices[i-1][j]*sin(
                        t_angles[i-1][j][k])/indices[i][j])
                    # print "Interior to interior: " + str(indices[i-1][j]) + " to " + str(indices[i][j])
                    # print "incident: " + str(180/pi*t_angles[i-1][j][k]) + " refracted: " + str(180/pi*t_angles[i][j][k])
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
    # generate the I matrices (n_layers+1 x n_wls x n_angles x 2)

    I = [[[[
        IMatrix(n_i, indices[i][j], angles[k], t_angles[i][j][k], l) if i == 0
        else IMatrix(indices[i-1][j], n_f, t_angles[i-1][j][k], t_angles[i][j][k], l) if i == len(stack)
        else IMatrix(indices[i-1][j], indices[i][j], t_angles[i-1][j][k], t_angles[i][j][k], l)
        for l in [0,1]] 
        for k in range(len(angles))] 
        for j in range(len(wls))] 
        for i in range(len(stack)+1)]

    #generate the P matrices (n_layers x n_wls x n_angles)
    P = [[[PMatrix(indices[i][j], wls[j], stack[i].depth, t_angles[i][j][k]) 
        for k in range(len(angles))] 
        for j in range(len(wls))] 
        for i in range(len(stack))]
    return (I,P)

"""method to compute matrices for transmission across interfaces"""
def IMatrix(n_1, n_2, theta_1, theta_2, pol):
    # Uses t,r consistent with Petersson et al 1999
    # The indices n_1 and n_1 input will be complex.
    if (pol==0): # perpendicular (s-) polarization
        r = (((n_1*cos(theta_1))-(n_2*cos(theta_2)))/
          ((n_1*cos(theta_1))+(n_2*cos(theta_2))))
        t = ((2*n_1*cos(theta_1))/
          ((n_1*cos(theta_1))+(n_2*cos(theta_2))))
    elif (pol==1): # parallel (p-) polarization
        r = (((n_1*cos(theta_2))-(n_2*cos(theta_1)))/
          ((n_2*cos(theta_1))+(n_1*cos(theta_2))))
        t = ((2*n_1*cos(theta_1))/
          ((n_2*cos(theta_1))+(n_1*cos(theta_2))))
    return array([[1/t, r/t],[r/t, 1/t]])

"""method to compute matrices for transmission through layers
(addition of a phase from propagation)"""
def PMatrix(n, wl, d, theta): 
    delta = (2*pi*d*n*cos(theta))/(wl)
    return array([[exp((-1j)*delta),0],[0,exp(1j*delta)]])

"""Given the modeling parameters and the transmission matrices generate
the E-field distribution throughout the device and return the integral
of E^2 in each layer.  Average over all given parameters (wavelengths,
angles, polarizations [expected to be a list of the form [0] or [0,1] where
0 is normal to the plane of incidence and 1 is parallel to the plane of 
incidence])."""

def evalField(stack, wls, angles, pol, n_i, n_f, indices, t_angles, P, I):
    # generate set of evolution matrices for each layer
    M = [[[[
        matrixCreate(i, len(stack),l, k, j, P, I)
        for l in range(len(pol))] 
        for k in range(len(angles))] 
        for j in range(len(wls))] 
        for i in range(len(stack)+1)]
    
    #Calculate the total E_f and E_r for the stack using M_0
    # currently set up backwards, seemingly according to literature
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

    return (E_0, E_f, E_i)

def matrixCreate(startLayer, endLayer, pol, angle, wLength, P, I):
    M_partial = identity(2)
    for index in range(startLayer, endLayer):
        M_partial = dot(M_partial,dot(P[index][wLength][angle], I[index+1][wLength][angle][pol]))
    return M_partial

def evalTR(stack, E_f, E_0, angles, n_i, n_f, t_angles, addBulkT=False):
    # option to add a bulk transmission coefficient
    # for a thick substrate. In this case the bulk
    # coefficient added will be for the final index (n_f) to air (n=1).
    # this expression is again Hecht 4.68 for normal incidence
    bulkT = (4.0*n_f)/((n_f+1.0)**2)
    if addBulkT:
        T = [[[real(bulkT*(n_f/n_i)*(cos(t_angles[len(stack)][j][k])/cos(angles[k]))*abs(
            E_f[j][k][l][0])**2)
            for l in range(len(E_f[0][0]))]
            for k in range(len(E_f[0]))]
            for j in range(len(E_f))]
    else:
        T = [[[real((n_f/n_i)*(cos(t_angles[len(stack)][j][k])/cos(angles[k]))*abs(
            E_f[j][k][l][0])**2)
            for l in range(len(E_f[0][0]))]
            for k in range(len(E_f[0]))]
            for j in range(len(E_f))]
    R = [[[(abs(E_0[j][k][l][1])/abs(E_0[j][k][l][0]))**2
        for l in range(len(E_f[0][0]))]
        for k in range(len(E_f[0]))]
        for j in range(len(E_f))]

    #averaging for both of these lists
    TAvg = mean(list(chain.from_iterable(chain.from_iterable(T))))
    RAvg = mean(list(chain.from_iterable(chain.from_iterable(R))))
    
    return (T, R, TAvg, RAvg)

"""Weight wavelength of incident light to closer resemble blackbody or
otherwise specially normalized incident radiation. Method will take index 
and wavelength linespace, and return blackbody curve in terms of weight 
on given index. Note: this produces new normalization! peak is in nanometers
currently relevant in ESqIntEval and ESqPlot"""
def weightSpectrum(index, wls, peak, active):
    h = 6.626*10**(-34)         # planck constant
    c = 3.0*10**(8)             # speed of light
    kb = 1.380*10**(-23)        # boltzmann constant
    wien = 2.89777*10**(-3)     # wien's displacement constant
    wavelength = peak*10**(-9)  # conver peak to meters
    T = wien/wavelength         # temperature of blackbody in kelvin

    if peak > wls[-1] or peak < wls[0]:
        raise ValueError("Peak wavelength not permissible in specified wls range")
    elif index < 0 or index > len(wls)-1:
        raise ValueError("Index not permissible in specified wls range")
    elif not active:
        return 1.0
    else:
        wv = wls[index]*10**(-9)    # convert wavelength to meters
        norm = 1/(2*h*c/(wavelength**5)*(1/(exp(h*c/(kb*T*wavelength))-1)))
        return norm*2*h*c/(wv**5)*(1/(exp(h*c/(kb*T*wv))-1))

"""Return the integral of Re(E)^2 in any layer given the initial E-field
vector.  This function accounts for absorption and evaluates the integral
numerically."""
def ESqIntegral(E_0, index, wl, d, theta, pol):
    return integrate.quad(lambda x: ESqEval(E_0, index, wl, x, theta, pol),0,d)[0]

"""Evaluate E^2 at a point (given by x) within a desired layer using the
initial E-field within that layer"""
def ESqEval(E_0, index, wl, x, theta, pol):
    E_f = E_0[0]
    E_r = E_0[1]
    delta = (2*pi*index*x*cos(theta))/wl

    # note reliance on modulus of field, not its real part.
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
and polarizations.  Use the initial E-field vectors E_i from evalField()"""
def ESqIntEval(stack, wls, angles, pols, indices, E_i):
    ESqInt = [[[[ESqIntegral(E_i[i][j][k][l], indices[i][j], 
        wls[j], stack[i].depth, angles[k], pols[l])
        for l in range(len(pols))] 
        for k in range(len(angles))]
        for j in range(len(wls))]
        for i in range(len(stack))]
    #value proportional to power absorbed in layer
    #use proportionality to real index, imag index, and 1/lambda from Petersson eqn. 22
    PAbsd = [[[[real(indices[i][j])*imag(indices[i][j])*ESqInt[i][j][k][l]*(10000.0/wls[j])
        for l in range(len(pols))] 
        for k in range(len(angles))]
        for j in range(len(wls))]
        for i in range(len(stack))]
    # for use with blackbody source
    ESqIntAvg = [sum(weightSpectrum(j, wls, peakWavelength, blackbody)*ESqInt[i][j][k][l] 
    #ESqIntAvg = [sum(ESqInt[i][j][k][l] 
        for j in range(len(wls))
        for k in range(len(angles))
        for l in range(len(pols)))/(len(wls)*len(angles)*len(pols))
        for i in range(len(stack))]
    # for use with blackbody source
    PAbsdAvg = [sum(weightSpectrum(j, wls, peakWavelength, blackbody)*PAbsd[i][j][k][l] 
    #PAbsdAvg = [sum(PAbsd[i][j][k][l] 
        for j in range(len(wls))
        for k in range(len(angles))
        for l in range(len(pols)))/(len(wls)*len(angles)*len(pols))
        for i in range(len(stack))]
    return (ESqInt, ESqIntAvg, PAbsd, PAbsdAvg)

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
    k_0 = 1.0e9/wls #conversion included from nm -> m

    # for use in porting to excel
    # print ",".join(map(lambda elem: str(elem),1.0e7/wls)) 
    # print ",".join(map(lambda elem: str(elem),TAnglPolAvg))

    # ax = plt.axes()
    # if tPlot:
    #     ax.plot(wls, TsAnglAvg, label='T_s')
    #     ax.plot(wls, TpAnglAvg, label='T_p')
    #     ax.plot(wls, TAnglPolAvg, label='T_avg')
    # if rPlot:
    #     ax.plot(wls, RsAnglAvg, label='R_s')
    #     ax.plot(wls, RpAnglAvg, label='R_p')
    #     ax.plot(wls, RAnglPolAvg, label='R_avg')
    # ax.legend()
    # ax.set_xlabel('Wavelength (nm)', fontsize=18)
    # plt.xticks(fontsize=15)
    # plt.yticks(fontsize=15)
    # if save:
    #     plt.savefig('results/{0}/{0}-TRWlsSpectrum.pdf'.format(saveFileName))
    # plt.show()

    ax = plt.axes()
    if tPlot:
        # ax.plot(0.01*k_0, TsAnglAvg, label='T_s')
        # ax.plot(0.01*k_0, TpAnglAvg, label='T_p')
        ax.plot(0.01*k_0, TAnglPolAvg, label='T_avg')

    if rPlot:
        # ax.plot(0.01*k_0, RsAnglAvg, label='R_s')
        # ax.plot(0.01*k_0, RpAnglAvg, label='R_p')
        ax.plot(0.01*k_0, RAnglPolAvg, label='R_avg')
    ax.legend(loc='upper left')
    ax.set_xlabel('k (cm^-1)', fontsize=18)
    # ax.set_ylim([2,5])
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    if save:
        plt.savefig('results/{0}/{0}-TRWaveNumsSpectrum.pdf'.format(saveFileName))
    plt.show()

    return (k_0, RAnglPolAvg, TAnglPolAvg)

"""Plot spectrum of the E^2 integral"""
def ESqIntSpectrumPlot(ESqInt, PAbsd, stack, wls, angles, pols, indices, save, saveFileName):
    ax = plt.axes()

    wnums = (1.0e7/wls) #wavenumbers in cm^-1 to be used in plotting
    for i in range(len(stack)):
        if stack[i].active:
            ESqIntWls = [sum(ESqInt[i][j][k][l] 
                for k in range(len(angles))
                for l in range(len(pols)))/(len(angles)*len(pols))
                for j in range(len(wls))]
            PAbsdWls = [sum(PAbsd[i][j][k][l]
                for k in range(len(angles))
                for l in range(len(pols)))/(len(angles)*len(pols))
                for j in range(len(wls))]
            ax.plot(wnums, ESqIntWls, label='{0} layer $\int E^2$'.format(stack[i].name))
            ax.plot(wnums, PAbsdWls, label='{0} layer power absorbed'.format(stack[i].name))

            plt.ylim([0,4000])

    ax.legend()
    ax.set_xlabel('Wavelength (nm)', fontsize=18)
    ax.set_xlabel('Wavenumber (cm$^-1$)', fontsize=18)
    ax.set_ylabel('(arb. units)', fontsize=18)
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
    # for use of blackbody radiation source
    # ESqAvg = [[sum(weightSpectrum(j, wls, peakWavelength, blackbody)*ESq[i][j][k][l][n]
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

def normPwrPlot(PAbsd, stack, wls, angles, pols, n_i, n_f, indices, save, saveFileName):
    # highly attenuating configuration for normalization factor
    # note that initial film layer should match n_i
    # au may be altered to pure reflector if sufficient accuracy not
    # achieved. All references are to tempStack until noted.
   
    # tempStack = [
    #        film(10000.0, 3.43 + 0.1j, 'abs', True),
    #        film(150.0, 'au', 'Gold', False)
    #        ]

    # indices = [[indexLookup(layer.index, wl) for wl in wls] for layer in tempStack]
    # t_angles = snell(indices, angles, n_i, n_f)
    # (I,P) = genMatrices(tempStack, wls, angles, n_i, n_f, indices, t_angles)

    # (E_0, E_f, E_i) = evalField(tempStack, wls, angles, pols,
    # n_i, n_f, indices, t_angles, P, I)

    # (ESqInt,ESqIntAvg, tempPAbsd, PAbsdAvg) = ESqIntEval(tempStack,wls,angles,pols,indices,E_i)

    # normFactorTable = [sum(tempPAbsd[0][j][k][l] # arbitrary normalization placeholder
    #              for k in range(len(angles))
    #              for l in range(len(pols)))/(len(angles)*len(pols))
    #              for j in range(len(wls))]

    # for generic sped up use (assuming constant absorption hereon)
    normFactor = 2660

    # stack should not be used, without reference to tempStack
    ax = plt.axes()
    wnums = (1.0e7/wls) #wavenumbers in cm^-1 to be used in plotting

    for i in range(len(stack)):
        if stack[i].active:
            # PAbsdWls = [sum(PAbsd[i][j][k][l]/normFactorTable[j] # for explicit calculation use

            PAbsdWls = [sum(PAbsd[i][j][k][l]/normFactor
                for k in range(len(angles))
                for l in range(len(pols)))/(len(angles)*len(pols))
                for j in range(len(wls))]
            ax.plot(wnums, PAbsdWls, label='{0} layer power absorbed'.format(stack[i].name))

            plt.ylim([0,1])

    ax.legend()
    ax.set_xlabel('Wavenumber (cm$^-1$)', fontsize=18)
    ax.set_ylabel('(normalized to source)', fontsize=18)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    if save:
        plt.savefig('results/{0}/{0}-ESqIntSpectrum.pdf'.format(saveFileName))
    plt.show()

def ODNormPlot(PAbsd, stack, wls, angles, pols, n_i, n_f, indices, save, saveFileName):
    # for generic sped up use (assuming constant absorption hereon)
    normFactor = 2660

    # stack should not be used, without reference to tempStack
    ax = plt.axes()
    wnums = (1.0e7/wls) #wavenumbers in cm^-1 to be used in plotting

    for i in range(len(stack)):
        if stack[i].active:
            ODAbsdWls = [-1.0*log10(1.0 - sum(PAbsd[i][j][k][l]/normFactor
                for k in range(len(angles))
                for l in range(len(pols)))/(len(angles)*len(pols)))
                for j in range(len(wls))]
            ax.plot(wnums, ODAbsdWls, label='{0} layer power absorbed'.format(stack[i].name))

            plt.ylim([0,0.2])

    ax.legend()
    ax.set_xlabel('Wavenumber (cm$^-1$)', fontsize=18)
    ax.set_ylabel('OD (unitless)', fontsize=18)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    if save:
        plt.savefig('results/{0}/{0}-ESqIntSpectrum.pdf'.format(saveFileName))
    plt.show()

if __name__ == "__main__":
    main()