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

# stack = [
#         film(0, 1.399, 'SiO2', 500, 2500, 10.0, False),
#         film(0, itoDrudeParams, 'ITO', 50, 50, 5.0, False),
#         film(0, 2.4 + 0.106j, 'CQD', 300, 800, 5.0, False),
#         film(0, 'au', 'Gold', 100, 100, 10.0, False)
#         ]
# number of layers in stack

# stack = [
#         film(0, itoDrudeParams, 'ITO', 50, 50, 5.0, False),
#         film(0, 2.4 + 0.106j, 'CQD', 300, 900, 5.0, False),
#         film(0, 'au', 'Gold', 100, 100, 10.0, False)
#         ]

stack = [
       film(6.0, 'NiCr', 'NiCr approximation',6.0, 6.0, 5.0, False),
       film(300.0, 2.4 + 0.106j, 'CQD', 200.0, 2000.0, 5.0, True),
       film(100.0, 'au', 'Gold', 100.0, 100.0, 5.0, False)
       ]


layerNum = len(stack)

###########################################################
# methods for heuristic searching
###########################################################

# method to return percentage of cells chosen to be active
# on a given propagation (returns values between 0 and 100 inclusive)
def getPercentage(population, selectivity):
    # return exp(-1.0*population/(selectivity))*100
    # print population/((1.0 - exp(-1.0*population**2/selectivity**2))*0.1*population + 1.0)
    return 100.0/((1.0 - exp(-1.0*population**2/selectivity**2))*0.1*population + 1.0)

# takes tuple, returns a list of tuples of tuples, 
# first element is tuple of (dim, step)
# second element is a tuple of coordinates for neighbor cell (dim, step) away
# method ignores boundaries and negatives (checked by hasTraversed)
def getNeighborCells(coords):
    arrayCoords = array(coords)
    dim = len(arrayCoords)
    basis = identity(dim,dtype=int)
    # orthogonal neighbors
    # tuple of ((dim, step), coord)
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
def propagate(currentState, activeCells, fitnessMap, optimum, ACTIVE_LAYER, selectivity):
    
    # threshold for barycenter and minimum active cells for heuristic
    threshold = 0.003
    minActive = 10
    
    # population selection threshold for cuttoff of greediness
    greedPopulation = 50

    # selectivity for percentage of population propagated per cycle
    # currently an argument of propagate for ease of testing
    # selectivity = 1000

    bestCell   = getFittestPercent(fitnessMap, 0)
    newState   = copy.deepcopy(currentState)
    newActive  = copy.deepcopy(activeCells)
    newOptimum = max(optimum, bestCell, key=lambda elem: elem[0][1])
    
    ###################################
    # code for cell choice heuristic
    # possibly deactivating or ignoring them
    ###################################

    # selectedCells = activeCells # default

    # begin search with utilization of all active cells
    # then grade to more selective
    population     = len(activeCells)
    percentage     = getPercentage(population, selectivity)
    fittestPercent = getFittestPercent(fitnessMap, percentage)

    # print "cells selected  : " + str(len(fittestPercent))
    # print "greed population: " + str(greedPopulation)

    ###################################################
    # Code for direct hill climbing mode
    # percentage required to activate defined above
    ###################################################
    # selectedCells  = map(lambda elem: elem[0], fittestPercent) # default

    # determine if few enough cells for greed
    if len(fittestPercent) < greedPopulation:
        # if only one cell selected, bypass below and simply select it
        if len(fittestPercent) > 1:
            # pass tail of fittest cells to following standard heuristic
            selectedCells  = map(lambda elem: elem[0], fittestPercent[1:]) 
            bestGroup   = [fittestPercent[0]]
            fittestCells = map(lambda elem: elem[0], bestGroup)
           
            directNeighbors = []

            # remember that fittestcell is a list of the one (1) fittest cell
            fittestCell = fittestCells[0]
            # remember to update grid !
            newState[fittestCell] = 1
            newActive.remove(fittestCell)

            neighbor = getNeighborCells(fittestCell)

            # populate with all valid neighbors
            for elem in neighbor:
                coords = elem[1]
                if not hasTraversed(newState, coords):
                    directNeighbors.append(coords)
                else:
                    continue

            # may eventually be altered to remove this call, given
            # redundancy of ACTIVE_LAYER and selection of all neighbors
            # can currently be used for selectivity of greed
            pseudoFitness = evaluateCells(directNeighbors, ACTIVE_LAYER)
            greedyCoords   = getFittestPercent(pseudoFitness, 100)

            for greedyCoord in greedyCoords:
                vdiff = subtract(greedyCoord[0], fittestCell)
                dim  = filter(lambda elem: elem != 0, [i if vdiff[i] != 0 else 0 for i in range(len(vdiff))])
                # dirty fix for zero-th dimension check
                if dim == []:
                    dim = 0
                else:
                    dim = dim[0]
                step = 1 if vdiff[dim] > 0 else -1 if vdiff[dim] < 0 else 0
                bestCoord = ((dim, step),greedyCoord[0])
                if not hasTraversed(newState, bestCoord[1]):
                    # traverse dim alters newState correctly and continuously
                    (newState, newCoord) = traverseDim(newState, 
                        fittestCell, dim, step)
                    newActive.append(newCoord)
        else:
            selectedCells  = map(lambda elem: elem[0], fittestPercent)
    else:
        selectedCells  = map(lambda elem: elem[0], fittestPercent)


    for cell in selectedCells:
        neighbors      = getNeighborCells(cell)
        newState[cell] = 1 # mark current cell as visited + inactive
        # SLOW: consider replacement by portion of split list from selectedCells
        newActive.remove(cell)
        
        ###################################################
        # code for neighbor choice heuristic
        ###################################################
        # selectedNeighbors = neighbors # default
        
        # determine most promising neighbors of fittest cells via barycenter
        selectedNeighbors = []
        offset            = getCenterOffset(fittestPercent)
    
        properDim = map(lambda elem: True if abs(elem) >= threshold else False, offset)

        signDim   = map(lambda elem: 1 if elem >= 0 else -1, offset)
        # if no clear advantageous direction, or lack of activity, simple floodfill
        if not any(properDim) or len(activeCells) < minActive:
            properDim = [True]*len(offset)
            signDim   = [0]*len(offset)

        for elem in neighbors:
            # check if dim is proper and movement in right direction
            dim  = elem[0][0]
            step = elem[0][1]
            if properDim[dim] and signDim[dim]*step >= 0:
                selectedNeighbors.append(elem)

        ###################################################
        # code for basic propagation (independent of heuristc)
        ###################################################

        for auxCell in selectedNeighbors:
            # second element in tuple gives coordinates
            # checks against newState, being continuously altered
            coord = auxCell[1]
            dim   = auxCell[0][0]
            step  = auxCell[0][1]
            if not hasTraversed(newState, coord):
                # traverse dim alters newState correctly and continuously
                (newState, newCoord) = traverseDim(newState, 
                    cell, dim, step)
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

    # used to keep population of chosen cells high
    # if len(sortedMap) < 4:
    #     return sortedMap
    # else:
    #     return sortedMap[0:4]

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
        # populate parameters in physical measurement from activeCell coords
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
        (ESqInt,ESqIntAvg, PAbsd, PAbsdAvg) = thinFilmClean.ESqIntEval(
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
    copyState = copy.deepcopy(currentState)
    sliceParam  = [slice(None) if index == dim1 or index == dim2 else 0 for index in range(len(copyState.shape))]
    stateSlice  = copyState[sliceParam]
    specified = 9999
    stateSlice[25][29] = 9999  
    prettySlice = map(lambda elem: map(lambda elem: "#" if elem == 2 else "." if elem == 1 else "O" if elem == specified else " ", elem), stateSlice)
    for row in prettySlice:
        for element in row:
            print element,
        print ""

###########################################################
# naive approach, populating a search space and traversing
###########################################################

def main():
    # cycle through selectivity (for testing purposes only!)
    for selectivity in range(850,900,50):
        print "selectivity: " + str(selectivity)

        MAX_STEP     = 500 # hard cutoff for propagations
        ACTIVE_LAYER = 1 # which layer in stack active (zero-index)
        DIVISIONS    = 5  # active cells per dimension (cartesian product base)
        THRESHOLD    = 3506
        
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

        # single cell initial position
        initialPosition = [(0,0,0)] # heuristic behavior

        # list of coordinates for active cells
        activeCells = copy.deepcopy(initialPosition)
        # initialize currentState according to activeCells
        for elem in activeCells:
            currentState[elem] = 2
        # store global optimum thusfar
        optimum = [((0,0,0),0.0)]

        # main propagation loop
        print "\nBest Performing Cell per Each Iteration:"
        print "________________________________________\n"   
        print "position in space\t\tfitness"
        
        count = 0

        print "{"

        while len(activeCells) != 0 and count < MAX_STEP and optimum[0][1] < THRESHOLD:
            count = count + 1
            fitnessMap = evaluateCells(activeCells, ACTIVE_LAYER)
            (newState, newActive, newOptimum) = propagate(currentState, activeCells, 
                fitnessMap, optimum, ACTIVE_LAYER, selectivity)
            currentState = newState
            activeCells  = newActive
            optimum      = newOptimum

            # to be used for cell dump
            # if optimum[0][1] > 3505:
            #     print "for selectivity: " + str(selectivity)
            #     print str(optimum[0][0]) + "\t\t\t" + "%.3f" % optimum[0][1] + " at count: " + str(count)
            #     print "________________________________________"
            #     print "{",


            #     for i in range(len(newState)):
            #         for j in range(len(newState[i])):
            #             for k in range(len(newState[i][j])):
            #                 for l in range(len(newState[i][j][k])):
            #                     if newState[i][j][k][l] != 0:
            #                         print "{" + str(i) + "," + str(j) + "," + str(k) + "},",

            #     # for elem in newActive:
            #     #     print "{" + ",".join(map(lambda elem: str(elem),elem[0:3])) + "},",
            #     print "}"
            #     print "________________________________________"
            #     break

            # print str(optimum[0][0]) + "\t\t\t" + "%.3f" % optimum[0][1] + " at count: " + str(count)
            print "{" + str(fitnessMap[0][0][1]) + "," + str(fitnessMap[0][1]) + "},"

        print "}"

        # prints statistics on search space size and percentage actually searched
        searchCount = 0
        for elem in currentState.flatten():
            if elem == 2 or elem == 1:
                searchCount = searchCount + 1
        totalLen = len(currentState.flatten())
        print "\ntotal cells searched: " + str(searchCount)
        print "total cells in space: " + str(totalLen)
        print "traversal percentage: "  + "%.3f" % (searchCount*100.0/totalLen)
        print "________________________________________"


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
    
if __name__ == "__main__":
    main()