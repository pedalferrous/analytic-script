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
# n_i = 3.43
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
        film(0, 1.399, 'SiO2', 95, 95, 10.0, False),
        film(0, itoDrudeParams, 'ITO', 5, 450, 10.0, False),
        film(0, 2.4 + 0.106j, 'CQD', 350, 850, 10.0, False),
        film(0, 'au', 'Gold', 200, 200, 10.0, False)
        ]
# number of layers in stack

# layer arrangement from Bouchon 2012
# stack = [
#         film(0, 'au', 'Gold_1', 50, 50, 5.0, False),
#         film(0, 'ZnS', 'Zinc Sulfide', 800, 1800, 1.0, False),
#         film(0, 'au', 'Gold_2', 200, 200, 10.0, False)
#         ]


layerNum = len(stack)

###########################################################
# methods for heuristic searching
###########################################################

# takes tuple, returns a list of tuples of tuples, 
# first element is tuple of (dim, step)
# second element is a tuple of coordinates for neighbor cell (dim, step) away
# method ignores boundaries and negatives (checked by hasTraversed)
def getNeighborCells(coords):
    arrayCoords = array(coords)
    dim = len(arrayCoords)
    basis = identity(dim,dtype=int)
    # orthogonal neighbors
    neighborsAdd = [((i, 1),tuple(add(arrayCoords,      basis[i,:]))) for i in range(dim)]
    neighborsSub = [((i,-1),tuple(add(arrayCoords, (-1)*basis[i,:]))) for i in range(dim)]
    # final concatenation and return
    return concatenate((neighborsAdd, neighborsSub),axis=0)

# traverse cell along dim by step (integer)
# work horse for flood fill (fed hasTraversed checked cells)
def traverseDim(currentState, activeCell, dim, step):
    # will return new activeCell or None, not modifying original 
    # will return new currentState, not modifying original
    newState  = copy.deepcopy(currentState)
    # copy cell to be modified
    # instantiate changed index value
    newCoord = copy.copy(activeCell)
    newIndex = newCoord[dim] + step 

    # check if valid traversal dim
    if dim < 0 or dim > len(newCoord) - 1:
        raise ValueError("Improper dimension. Check size of search space.")
    # check if new coord is valid in grid
    elif newIndex < 0 or newIndex > currentState.shape[dim]:
        raise ValueError("Cell position invalid. Did you feed into hasTraversed?")
    # both are valid and modification is possible
    else:
        # make current cell visited + inactive
        newState[activeCell] = 1
        # create new active coordinate according to dim and step
        newCoord = newCoord[:dim] + tuple([newIndex]) + newCoord[dim+1:]
        # check if cell is already visited
        if currentState[newCoord] != 0:
            raise ValueError("Cell has already been visited. Did you feed into hasTraversed?")
        # make active said new coord and return results
        newState[newCoord] = 2
        return (newState, newCoord)

# determine whether specified coords (a tuple) have been traversed
# OR if cell is not contained within search space
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

# method will expand around active cells according to
# fitnessMap, returning copies of new state and active
def propagate(currentState, activeCells, fitnessMap, optimum):
    newState   = copy.deepcopy(currentState)
    newActive  = copy.deepcopy(activeCells)
    newOptimum = max(optimum, getFittestPercent(fitnessMap, 0), key=lambda elem: elem[0][1])
    ###################################
    # code for cell choice heuristic
    # possibly deactivating or ignoring them
    ###################################
    # selectedCells = activeCells # default

    # begin search with utilization of all active cells
    # then grade to more selective
    if len(activeCells) < 24:
        percentage = 70
        # print "level 1"
    elif len(activeCells) < 100:
        percentage = 15
        # print "level 2"        
    elif len(activeCells) < 200:
        percentage = 10
        # print "level 3"        
    else:
        percentage = 5
        # print "level 4"    
    # print "active cells: " + str(len(activeCells))    
    fittestPercent = getFittestPercent(fitnessMap, percentage)
    selectedCells  = map(lambda elem: elem[0], fittestPercent) 

    for cell in selectedCells: # loop over most promising
        neighbors      = getNeighborCells(cell)
        newState[cell] = 1 # mark current cell as visited + inactive
        # SLOW: consider replacement by portion of split list from selectedCells
        newActive.remove(cell)
        
        ###############################
        # code for neighbor choice heuristic
        ###############################
        # selectedNeighbors = neighbors # default
        
        # determine most promising neighbors of fittest cells via barycenter
        selectedNeighbors = []
        offset            = getCenterOffset(fittestPercent)
        # threshold for barycenter and minimum active cells for heuristic
        threshold = 0.003
        minActive = 10
        properDim = map(lambda elem: True if abs(elem) >= threshold else False, offset)
        signDim   = map(lambda elem: 1 if elem >= 0 else -1, offset)
        # if no clear advantageous direction, or lack of activity, simple floodfill
        if not any(properDim) or len(activeCells) < minActive:
            properDim = [True]*len(offset)
            signDim   = [0]*len(offset)

        for elem in neighbors:
            # check if dim is proper and movement in right direction
            # elem[0][0] is dim, elem[0][1] is step
            if properDim[elem[0][0]] and signDim[elem[0][0]]*elem[0][1] >= 0:
                selectedNeighbors.append(elem)

        ###############################
        # code for basic propagation (independent of heuristc)
        ###############################
        for auxCell in selectedNeighbors:
            # second element in tuple gives coordinates
            # checks against newState, being continuously altered
            if not hasTraversed(newState, auxCell[1]):
                # traverse dim alters newState correctly and continuously
                (newState, newCoord) = traverseDim(newState, 
                    cell, auxCell[0][0], auxCell[0][1])
                newActive.append(newCoord)
            else:
                continue

    # keep track of state, active, and fittest cell of this iteration
    return (newState, newActive, newOptimum)

# gives barycenter minus geometric center
def getCenterOffset(fitnessMap):
    return subtract(getCoordinateCenter(fitnessMap, True), 
        getCoordinateCenter(fitnessMap, False))

# returns top % performing active cells
def getFittestPercent(fitnessMap, percent):
    if percent < 0 or percent > 100:
        raise ValueError("Invalid percentage choice: must be within [0,100].")
    # revere sort fitnessMap by fitness
    sortedMap = sorted(fitnessMap, key=lambda elem: elem[1], reverse=True)
    lenMap    = int((percent/100.0)*len(sortedMap))
    # always return at least one element
    if lenMap == 0:
        return sortedMap[0:1]
    else:
        return sortedMap[0:lenMap]

# finds weghted or unweighted coordinate center
def getCoordinateCenter(fitnessMap, isWeighted):
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
# this method alone calls thinFilm
def evaluateCells(activeCells, activeLayer):
    fitnessMap = []
    for elem in activeCells:
        # populate parameters in physical measurement from activeCell coordinates
        parameters = [stack[i].min_depth + elem[i]*stack[i].discrete_depth for i in range(len(elem))]
        for i in range(len(parameters)):
            # find a more efficient way to populate
            stack[i] = film(
                parameters[i], 
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

        # for use in minimizing T + R
        # (T, R, TAvg, RAvg) = thinFilmClean.evalTR(portStack, E_f, E_0, angles, n_i, n_f, 
        # t_angles, addBulkT=True)
        # # work to minimize their sum by maximizing its reciprocal
        # # no use of active layer
        # fitnessElem = (elem, 1.0/(TAvg + RAvg))

        fitnessElem = (elem, ESqIntAvg[activeLayer])
        fitnessMap.append(fitnessElem)
    return fitnessMap

###########################################################
# various helper methods for display
###########################################################

# produces prettry planar graph of currentState in dim1-dim2 plane
# currently does not slice automatically (hard-coded)
def prettyPrint(currentState, dim1, dim2):
    sliceParam  = [slice(None) if index == dim1 or index == dim2 else 0 for index in range(len(currentState.shape))]
    stateSlice  = currentState[sliceParam]
    prettySlice = map(lambda elem: map(lambda elem: "#" if elem == 2 else "." if elem == 1 else " ", elem), stateSlice)
    for row in prettySlice:
        for element in row:
            print element,
        print ""

###########################################################
# naive approach, populating a search space and traversing
###########################################################

def main():

    MAX_STEP     = 400 # hard cutoff for propagations
    ACTIVE_LAYER = 2  # which layer in stack active (zero-index)
    DIVISIONS    = 4  # active cells per dimension (cartesian product base)
    # same size as searchSpace
    # 0: unvisited, 1: visited + inactive
    # 2: visited + active.
    currentState = zeros([int((stack[i].max_depth
        -stack[i].min_depth)/stack[i].discrete_depth)+1 
            for i in range(layerNum)])
    # index of starting point (permits multuple seed locations)
    # list of tuples: ALL cell coordinate are tuples
    initialPosition = []
    stateShape      = currentState.shape
    dimIntervals    = [range(0,elem,elem/DIVISIONS-1) if elem > DIVISIONS else range(0,elem) for elem in stateShape]
    print "\ninterval division of search space: \n" + str(dimIntervals) + "\n"
    initialPosition     = array(meshgrid(*dimIntervals)).T.reshape(-1,len(stateShape))
    initialPosition     = map(lambda elem: tuple(elem), initialPosition)

    # list of coordinates for active cells
    activeCells = copy.deepcopy(initialPosition)
    # initialize currentState according to activeCells
    for elem in activeCells:
        currentState[elem] = 2
    # store global optimum thusfar
    optimum = [((0,0,0,0),0.0)]

    # main propagation loop
    print "\nBest Performing Cell per Each Iteration:"
    print "________________________________________\n"
    print "position in space\t\tfitness"
    count = 0
    while len(activeCells) != 0 and count < MAX_STEP:
        count = count + 1
        fitnessMap = evaluateCells(activeCells, ACTIVE_LAYER)
        (newState, newActive, newOptimum) = propagate(currentState, activeCells, fitnessMap, optimum)
        currentState = newState
        activeCells  = newActive
        optimum      = newOptimum
        print str(optimum[0][0]) + "\t\t\t" + "%.3f" % optimum[0][1]
        prettyPrint(newState, 1, 2)

    # prints statistics on search space size and percentage actually searched
    searchCount = 0
    for elem in currentState.flatten():
        if elem == 2 or elem == 1:
            searchCount = searchCount + 1
    totalLen = len(currentState.flatten())
    print "________________________________________"
    print "\ntotal cells searched: " + str(searchCount)
    print "total cells in space: " + str(totalLen)
    print "traversal percentage: "  + "%.3f" % (searchCount*100.0/totalLen)

    #######################################################
    # unit tests
    #######################################################  

    # testState = zeros([2,2,2,2])
    # testActive = [(0,0,0,0)]
    # print "\n__________ test hasTraversed\n"
    # print testState
    # print "\nstatus of cell (0,0,0,0) -- hasTraversed: " + str(hasTraversed(testState, (0,0,0,0)))
    # print "\nstatus of impossible cell (-1,0,0,0) -- hasTraversed: " + str(hasTraversed(testState, (-1,0,0,0)))
    # print "\nstatus of impossible cell (0,0,3,0) -- hasTraversed: " + str(hasTraversed(testState, (0,0,3,0)))
    # print "\nmodifying (1,1,1,1) from value 0 to value 1\n"
    # testState[(1,1,1,1)] = 1
    # print testState
    # print "\nstatus of cell (1,1,1,1) -- hasTraversed: " + str(hasTraversed(testState, (1,1,1,1)))
    # print "\nstatus of cell (1,1,1,0) -- hasTraversed: " + str(hasTraversed(testState, (1,1,1,0)))

    # print "\n__________ test traverseDim\n"
    # print "\nlist of active cells by coordinate\n"
    # print testActive
    # for elem in testActive:
    #     testState[elem] = 2
    # print "\ntestState given activated cells\n"
    # print testState

    # # should leave originals unaltered while returning new active and state
    # (newState, activeCell) = traverseDim(testState, testActive[0], 0, 1)
    # testActive = [activeCell]
    # testState = newState
    # print "\nnew state given traversal from (0,0,0,0) to (1,0,0,0)\n"
    # print testState
    # print "\nnew activeCell list given movement from (0,0,0,0) to (1,0,0,0)\n"
    # print testActive

    # (newState, activeCell) = traverseDim(testState, testActive[0], 1, 1)
    # testActive = [activeCell]
    # testState = newState
    # print "\nnew state given traversal from (1,0,0,0) to (1,1,0,0)\n"
    # print testState
    # print "\nnew activeCell list given movement from (1,0,0,0) to (1,1,0,0)\n"
    # print testActive

    # print "\n__________ test evaluateCells\n"
    # # calls with active cells and layer (CQD)
    # fitnessMap = evaluateCells(testActive, 2)
    # print "\nfitnessMap given active CQD layer and (1,1,0,0) active cell\n"
    # print fitnessMap
    
    # fitnessMap = [((1, 1, 0, 0), 200),((0, 1, 1, 0), 100)]

    # print "\ndifference between geometric- and bary-center\n"
    # print getCenterOffset(fitnessMap)
    # print "\neighbors of coodinate (0,1,0,0)\n"
    # print getNeighborCells((0,1,0,0))
    
if __name__ == "__main__":
    main()