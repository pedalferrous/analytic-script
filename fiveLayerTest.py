###########################################################
# Zane Rossi, University of Chicago
# June, 2016, PGS Lab
###########################################################

from numpy import *
from collections import *
from itertools import *
import copy
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from scipy.interpolate import interp1d
import time
# outside file
import thinFilmClean


###########################################################
# Indicate materials in use and inherent limitations.
# i.e. min depth, max depth, depth discretization
###########################################################

# constants required for the operation of thinFilmClean
wls = [4000]
# angle variation deprecated
angles = (pi/180)*array([0.0])
#0 represents s-polarization while 1 is p-polarization
pols = [0,1]
#incident and final indices (currently free space)
n_i = 3.43
n_f = 1.0

# initial stack, with limitation parameters imposed
ito_omega_p = 2.73649315189e+14 #s^-1
ito_epsilon_inf = 3.24
ito_gamma = 5.96450643459e+12 #s^-1
itoDrudeParams = [ito_omega_p, ito_epsilon_inf, ito_gamma]

# prototypical data structure to populate stack
film = namedtuple('layer',['depth','index','name',
    'min_depth', 'max_depth', 'discrete_depth', 'active'])
# modified film to be ported into thinFilmClean.py
portFilm = namedtuple('layer',['depth','index','name','active'])

stack = [
        film(0, 1.399, 'SiO2_0', 10, 150, 10.0, False),
        film(0, itoDrudeParams, 'ITO', 60, 160, 10.0, False),
        film(0, 2.4 + 0.106j, 'CQD', 540, 840, 20.0, False),
        film(0, 1.399, 'SiO2_1', 10, 150, 10.0, False),
        film(0, 'au', 'Gold', 200, 200, 2.0, False)
        ]
# number of layers in stack
layerNum = len(stack)

###########################################################
# various helper methods for debugging
###########################################################

# prints elements of search space serially
# edit to include names of materials
def printSearchSpace(searchSpace, layerNum):
    temp = copy.deepcopy(searchSpace)
    # reduce to flat list of tuples
    for n in range(layerNum-1):
        temp = chain.from_iterable(temp)
    temp = list(temp)
    # print entire list, or truncated version for debugging
    if len(list(temp)) <= 8:       
        for index in range(layerNum):
            print "Layer " + str(index+1) + "\t\t",
        print "Average E^2 in Active Layer"
        for elem in temp:
            depths = str(elem[0]).split(":")
            for layer in depths:
                print str(layer) + "\t\t",
            print str(elem[1])
    else:
        for index in range(layerNum):
            print "Layer " + str(index+1) + "\t\t",
        print "Average E^2 in Active Layer"
        for index in range(4):
            depths = str(temp[index][0]).split(":")
            for layer in depths:
                print str(layer) + "\t\t",
            print str(temp[index][1]) 
        print "\n.  .  .  " + str(len(temp)-8) + " elements " + " .  .  .\n"
        for index in range(4):
            depths = str(temp[index-4][0]).split(":")
            for layer in depths:
                print str(layer) + "\t\t",
            print str(temp[index-4][1])      


# method for finding maximum E^2 avg given search space
def getOptimum(searchSpace, layerNum):
    temp = copy.deepcopy(searchSpace)
    for n in range(layerNum-1):
        temp = chain.from_iterable(temp)
    return max(temp, key=lambda elem: elem[1])

###########################################################
# naive approach, populating a search space and traversing
###########################################################

def main():

    # initialization of proper size multidimensional numpy array
    # inserted values are initialized to type 'O' (object) given eventual replacement by arbitrary tuples
    searchSpace = full([int((stack[i].max_depth-stack[i].min_depth)/stack[i].discrete_depth)+1 
        for i in range(layerNum)],None,dtype='O')


    print "\nCharacter of search space:"
    print "1st dim of size: " + str(len(searchSpace))
    print "2nd dim of size: " + str(len(searchSpace[0]))
    print "3rd dim of size: " + str(len(searchSpace[0][0]))
    print "4th dim of size: " + str(len(searchSpace[0][0][0]))
    print "5th dim of size: " + str(len(searchSpace[0][0][0][0]))
    print ""

    # for extremely provisional use
    # innefiecient means of searchSpace filling (no heuristic)
    # indicating CQD layer in our simulation
    t0 = time.clock()
    activeLayerIndex = 2
    for dim0 in range(len(searchSpace)):
        for dim1 in range(len(searchSpace[0])):
            for dim2 in range(len(searchSpace[0][0])):
                for dim3 in range(len(searchSpace[0][0][0])):
                    for dim4 in range(len(searchSpace[0][0][0][0])):
                        # populating layer depths
                        init0  = stack[0].min_depth + dim0*stack[0].discrete_depth
                        init1  = stack[1].min_depth + dim1*stack[1].discrete_depth
                        init2  = stack[2].min_depth + dim2*stack[2].discrete_depth
                        init3  = stack[3].min_depth + dim3*stack[3].discrete_depth
                        init4  = stack[4].min_depth + dim4*stack[4].discrete_depth
                        initialParameters = [[init0, init1, init2, init3, init4]]
                        for i in range(len(initialParameters[0])):
                            # there is a way to do this more efficiently (do we need non-mutable named tuples)
                            stack[i] = film(initialParameters[0][i], stack[i].index, stack[i].name, stack[i].min_depth, stack[i].max_depth, stack[i].discrete_depth, stack[i].active)

                        portStack = [portFilm(layer.depth, layer.index, layer.name, layer.active) for layer in stack]

                        # required pre-processing for field calculations
                        indices = [[thinFilmClean.indexLookup(layer.index, wl) for wl in wls] for layer in portStack]
                        t_angles = thinFilmClean.snell(indices, angles, n_i, n_f)
                        (I,P) = thinFilmClean.genMatrices(portStack, wls, angles, n_i, n_f, indices, t_angles)
                        # calculate E-field, its average square
                        (E_0, E_f, E_i) = thinFilmClean.evalField(portStack, wls, angles, pols, n_i, n_f, indices, t_angles, P, I)
                        (ESqInt,ESqIntAvg) = thinFilmClean.ESqIntEval(portStack,wls,angles,pols,indices,E_i)
                        # indicates E^2 field in CQD
                        searchSpace[dim0][dim1][dim2][dim3][dim4] = (str(init0) + ":" + str(init1) + ":" + str(init2) + ":" + str(init3) + ":" + str(init4) , ESqIntAvg[activeLayerIndex])
    t1 = time.clock()

    print "\nin " + str(t1-t0) + " seconds"

    printSearchSpace(searchSpace, layerNum)
    print "\nOptimum found with node: "
    print(getOptimum(searchSpace, layerNum))

if __name__ == "__main__":
    main()