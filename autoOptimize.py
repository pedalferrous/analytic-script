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
        film(0, 1.399, 'SiO2', 90, 90, 10.0, False),
        film(0, itoDrudeParams, 'ITO', 130, 130, 5.0, False),
        film(0, 2.4 + 0.106j, 'CQD', 440, 840, 10.0, False),
        film(0, 'au', 'Gold', 200, 200, 2.0, False)
        ]
# number of layers in stack
layerNum = len(stack)

###########################################################
# methods for heuristic searching
###########################################################

# array of same size as searchSpace
# 0: unvisited, 1: visited + inactive
# 2: visited + active.
currentState = zeros([int((stack[i].max_depth
    -stack[i].min_depth)/stack[i].discrete_depth)+1 
        for i in range(layerNum)])
# index of starting point (permits multuple seed locations)
# note that you have to generate 
initialPosition = [zeros(len(stack))]
# list of coordinates for active cells
activeCells = copy.deepcopy(initialPosition)

# takes tuple, returns a list of tuple coordinates of neighbor cells
# ignoring boundaries and negatives
def neighborCells(coords):
    arrayCoords = array(coords)
    dim = len(arrayCoords)
    basis = identity(dim)
    neighborsAdd = [add(arrayCoords,      basis[i,:]) for i in range(dim)]
    neighborsSub = [add(arrayCoords, (-1)*basis[i,:]) for i in range(dim)]
    return tuple(concatenate((neighborsAdd, neighborsSub),axis=0))

# traverse cell along dim by step (integer)
# work horse for flood fill (fed hasTraversed checked cells)
def traverseDim(currentState, activeCell, dim, step):
    # will return new activeCell or None, not modifying original 
    # will return new currentState, not modifying original
    newState  = copy.deepcopy(currentState)
    # copy cell to be modified
    # instantiate changed index value
    newcoord = copy.copy(activeCell)
    newIndex = newcoord[dim] + step 

    # check if valid traversal dim
    if dim < 0 or dim > len(newcoord) - 1:
        raise ValueError("Improper dimension. Check size of search space.")
    # check if new coord is valid in grid
    elif newIndex < 0 or newIndex > currentState.shape[dim]:
        raise ValueError("Cell position invalid. Did you feed into hasTraversed?")
    # both are valid and modification is possible
    else:
        # make current cell visited + inactive
        newState[activeCell] = 1
        # create new active coordinate according to dim and step
        newcoord = newcoord[:dim] + (newIndex,) + newcoord[dim+1:]
        # check if cell is already visited
        if currentState[newcoord] != 0:
            raise ValueError("Cell has already been visited. Did you feed into hasTraversed?")
        # make active said new coord and return results
        newState[newcoord] = 2
        return (newState, newcoord)

# determine whether specified coords (a tuple) have been traversed
# OR if cell is not contained within grid
def hasTraversed(currentState, coords):
    # returns a boolean if invalid for any reason
    for index in range(len(coords)):
        if coords[index] < 0 or coords[index] > currentState.shape[index] - 1:
            return True
        else:
            continue
    if currentState[coords] == 0:
        return False
    else:
        return True

"""
   ________  ______  ____  _______   ________   _       ______  ____  __ __
  / ____/ / / / __ \/ __ \/ ____/ | / /_  __/  | |     / / __ \/ __ \/ //_/
 / /   / / / / /_/ / /_/ / __/ /  |/ / / /     | | /| / / / / / /_/ / ,<
/ /___/ /_/ / _, _/ _, _/ /___/ /|  / / /      | |/ |/ / /_/ / _, _/ /| |
\____/\____/_/ |_/_/ |_/_____/_/ |_/ /_/       |__/|__/\____/_/ |_/_/ |_|
"""

# method will expand around active cells according to
# fitnessMap, returning copies of new state and active
def propagate(currentState, activeCells, fitnessMap):
    newState  = copy.deepcopy(currentState)
    newActive = []
    for cell in activeCells:
        newState[cell] = 1
        neighbors = neighborCells(cell)
        for auxCell in neighbors:
            if not hasTraversed(auxCell):
                # remember to add this cell to newActive
                continue
            else:
                continue
    return (newState, newActive)

# gives barycenter minus geometric center
def centerOffset(fitnessMap):
    return subtract(coordinateCenter(fitnessMap, True), 
        coordinateCenter(fitnessMap, False))

# finds weghted or unweighted coordinate center
def coordinateCenter(fitnessMap, isWeighted):
    if isWeighted:
        # sum coords as vectors multiplied by fitness
        weightedSum = [sum([elem[0][i]*elem[1] for elem in fitnessMap]) for i in range(len(fitnessMap[0][0]))]
        # sum all fitness to give barycentric coordinate
        norm = sum([elem[1] for elem in fitnessMap])
        weightedSum = array([(1.0/norm)*elem for elem in weightedSum])
        return weightedSum
    else:
        # same as above with weight unity on all vectors
        nonWeightedSum = [sum([elem[0][i] for elem in fitnessMap]) for i in range(len(fitnessMap[0][0]))]
        # sum of unity becomes len
        norm = len(fitnessMap)
        nonWeightedSum = array([(1.0/norm)*elem for elem in nonWeightedSum])
        return nonWeightedSum

# returns list of tuples with coordinates and fitness to be used by propagate
def evaluateCells(activeCells, activeLayer):
    fitnessMap = []
    for elem in activeCells:
        # populate parameters in physical measurement from activeCell coordinates
        parameters = [stack[i].min_depth + elem[i]*stack[i].discrete_depth for i in range(len(elem))]
        for i in range(len(parameters)):
            # find a more efficient way to populate
            stack[i] = film(parameters[i], 
                stack[i].index, 
                stack[i].name, 
                stack[i].min_depth, 
                stack[i].max_depth, 
                stack[i].discrete_depth, 
                stack[i].active)
        # convert to portable format given thinFilm.py
        portStack = [portFilm(
            layer.depth, 
            layer.index, 
            layer.name, 
            layer.active) for layer in stack]
        # required pre-processing for field calculations
        indices = [[thinFilmClean.indexLookup(
            layer.index, 
            wl) for wl in wls] for layer in portStack]
        t_angles = thinFilmClean.snell(indices, angles, n_i, n_f)
        (I,P) = thinFilmClean.genMatrices(
            portStack, 
            wls, 
            angles, 
            n_i, 
            n_f, 
            indices, 
            t_angles)
        # calculate E-field, its average square
        (E_0, E_f, E_i) = thinFilmClean.evalField(
            portStack, 
            wls, 
            angles, 
            pols, 
            n_i, 
            n_f, 
            indices, 
            t_angles, 
            P, 
            I)
        (ESqInt,ESqIntAvg) = thinFilmClean.ESqIntEval(
            portStack,
            wls,
            angles,
            pols,
            indices,
            E_i)
        fitnessElem = (elem, ESqIntAvg[activeLayer])
        fitnessMap.append(fitnessElem)
    return fitnessMap

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
    # data identification #:#:#:#... first element of stored tuple
    # in populated search space
    searchSpace = []
    # indicating CQD layer in our simulation
    activeLayerIndex = 2

    # initialization of proper size multidimensional numpy array
    # inserted values are initialized to type 'O' (object) given eventual replacement by arbitrary tuples
    searchSpace = full([int((stack[i].max_depth-stack[i].min_depth)/stack[i].discrete_depth)+1 
        for i in range(layerNum)],None,dtype='O')

    initialParameters = [[800.0,10.0,400.0,50.0]]

    # print "\nCharacter of search space:"
    # print "1st dim of size: " + str(len(searchSpace))
    # print "2nd dim of size: " + str(len(searchSpace[0]))
    # print "3rd dim of size: " + str(len(searchSpace[0][0]))
    # print "4th dim of size: " + str(len(searchSpace[0][0][0]))
    # print ""

    # for extremely provisional use
    # innefiecient means of searchSpace filling (no heuristic)
    # for dim0 in range(len(searchSpace)):
    #     for dim1 in range(len(searchSpace[0])):
    #         for dim2 in range(len(searchSpace[0][0])):
    #             for dim3 in range(len(searchSpace[0][0][0])):
    #                 # populating layer depths
    #                 # delta0 = dim1*(stack[0].max_depth-stack[0].min_depth)/stack[0].discrete_depth
    #                 # delta1 = dim2*(stack[1].max_depth-stack[1].min_depth)/stack[1].discrete_depth
    #                 # delta2 = dim3*(stack[2].max_depth-stack[2].min_depth)/stack[2].discrete_depth
    #                 # delta3 = dim4*(stack[3].max_depth-stack[3].min_depth)/stack[3].discrete_depth
    #                 init0  = stack[0].min_depth + dim0*stack[0].discrete_depth
    #                 init1  = stack[1].min_depth + dim1*stack[1].discrete_depth
    #                 init2  = stack[2].min_depth + dim2*stack[2].discrete_depth
    #                 init3  = stack[3].min_depth + dim3*stack[3].discrete_depth
    #                 initialParameters = [[init0, init1, init2, init3]]
    #                 for i in range(len(initialParameters[0])):
    #                     # there is a way to do this more efficiently (do we need non-mutable named tuples)
    #                     stack[i] = film(initialParameters[0][i], stack[i].index, stack[i].name, stack[i].min_depth, stack[i].max_depth, stack[i].discrete_depth, stack[i].active)

    #                 portStack = [portFilm(layer.depth, layer.index, layer.name, layer.active) for layer in stack]

    #                 # required pre-processing for field calculations
    #                 indices = [[thinFilmClean.indexLookup(layer.index, wl) for wl in wls] for layer in portStack]
    #                 t_angles = thinFilmClean.snell(indices, angles, n_i, n_f)
    #                 (I,P) = thinFilmClean.genMatrices(portStack, wls, angles, n_i, n_f, indices, t_angles)
    #                 # calculate E-field, its average square
    #                 (E_0, E_f, E_i) = thinFilmClean.evalField(portStack, wls, angles, pols, n_i, n_f, indices, t_angles, P, I)
    #                 (ESqInt,ESqIntAvg) = thinFilmClean.ESqIntEval(portStack,wls,angles,pols,indices,E_i)
    #                 # indicates E^2 field in CQD
    #                 searchSpace[dim0][dim1][dim2][dim3] = (str(init0) + ":" + str(init1) + ":" + str(init2) + ":" + str(init3) , ESqIntAvg[activeLayerIndex])

    # printSearchSpace(searchSpace, layerNum)
    # print "\nOptimum found with node: "
    # print(getOptimum(searchSpace, layerNum))

    #######################################################
    # testing heuristic methods
    #######################################################  

    testState = zeros([2,2,2,2])
    testActive = [(0,0,0,0)]
    print "\n__________ test hasTraversed\n"
    print testState
    print "\nstatus of cell (0,0,0,0) -- hasTraversed: " + str(hasTraversed(testState, (0,0,0,0)))
    print "\nstatus of impossible cell (-1,0,0,0) -- hasTraversed: " + str(hasTraversed(testState, (-1,0,0,0)))
    print "\nstatus of impossible cell (0,0,3,0) -- hasTraversed: " + str(hasTraversed(testState, (0,0,3,0)))
    print "\nmodifying (0,0,0,0) from value 0 to value 1\n"
    testState[(1,1,1,1)] = 1
    print testState
    print "\nstatus of cell (1,1,1,1) -- hasTraversed: " + str(hasTraversed(testState, (1,1,1,1)))
    print "\nstatus of cell (1,1,1,0) -- hasTraversed: " + str(hasTraversed(testState, (1,1,1,0)))

    print "\n__________ test traverseDim\n"
    print "\nlist of active cells by coordinate\n"
    print testActive
    for elem in testActive:
        testState[elem] = 2
    print "\ntestState given activated cells\n"
    print testState

    # should leave originals unaltered while returning new active and state
    (newState, activeCell) = traverseDim(testState, testActive[0], 0, 1)
    testActive = [activeCell]
    testState = newState
    print "\nnew state given traversal from (0,0,0,0) to (1,0,0,0)\n"
    print testState
    print "\nnew activeCell list given movement from (0,0,0,0) to (1,0,0,0)\n"
    print testActive

    (newState, activeCell) = traverseDim(testState, testActive[0], 1, 1)
    testActive = [activeCell]
    testState = newState
    print "\nnew state given traversal from (1,0,0,0) to (1,1,0,0)\n"
    print testState
    print "\nnew activeCell list given movement from (1,0,0,0) to (1,1,0,0)\n"
    print testActive

    print "\n__________ test evaluateCells\n"
    # calls with active cells and layer (CQD)
    fitnessMap = evaluateCells(testActive, 2)
    print "\nfitnessMap given active CQD layer and (1,1,0,0) active cell\n"
    print fitnessMap
    fitnessMap = [((1, 1, 0, 0), 200),((0, 1, 1, 0), 100)]
    # baryCenter = coordinateCenter(fitnessMap, True)
    # geomCenter = coordinateCenter(fitnessMap, False)
    # print "\nbarycenter of fitnessMap\n"
    # print baryCenter
    # print "\ngeometric center of fitnessMap\n"
    # print geomCenter
    print "\ndifference between geometric- and bary-center\n"
    print centerOffset(fitnessMap)
    print "\nneighbors of coodinate (0,1,0,0)\n"
    print neighborCells((0,1,0,0))
    
    #######################################################
    # dirty display for 1D slices only
    #######################################################    

    # space = ndarray.flatten(searchSpace)
    # print space

    # porting to mathematica for dirty plots
    # print "{",
    # for elem in space:
    #     print str(elem[1]) + ",",
    # print "}"
    # poster child thus far
    # ('300.0:80.0:600.0:200.0', 3414.849686094519, [447.17418606338254, 286.86860412575743, 3414.849686094519, 0.9938275291899248])
    # ('90.0:130.0:640.0:200.0', 3504.612761494427, [82.01202593364928, 245.38468092862215, 3504.612761494427, 0.9749937541718504])



if __name__ == "__main__":
    main()