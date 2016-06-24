#######################################
# Zane Rossi
# University of Chicago
# PGS Lab: June, 2016
#######################################

from numpy import *
from collections import *
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import os
from scipy.interpolate import interp1d

#######################################
# Instantiation of grid structure and duration
# as wel as various other electromagnetic constants
#######################################

# size of total simulation in nm
FIELD_SIZE = 5000
# arbitrary dimension determined by other constants
DURATION = 2000
# wavelength in nm; goal is acceptance of spectrum
WAVELENGTH = 1000
# fundamntal mesh dz in nm (currently uniform)
DZ = 5
# barrier is a list of layers
BARRIER = []
# E and H fields given FIELD_SIZE and DZ
MESH_SIZE = FIELD_SIZE/DZ
E = [0.0]*MESH_SIZE
H = [0.0]*MESH_SIZE
# electromagnetic constants
MU = 1.0
EP = 1.0

#######################################
# Data Structures
#######################################

# outside specification of layer sturcture
# start position, end position, index of refraction (possible lookup), and
# whether the layer E^2 is to be actively calculated
layer = namedtuple('layer',['start','stop','index','active'])


def main():
	print "There is nothing here yet."


if __name__ == "__main__":
	main()