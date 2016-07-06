from numpy import *
from collections import *
from scipy import loadtxt
import matplotlib.pyplot as plt
import thinFilm

"""This is a script to find results for a range of layer
thicknesses or indices using thinFilm.py.  It varies device parameters,
not the 'external' parameters (i.e., incident light wavelengths, angles,
polarizations)."""

#film is defined the same way as in thinFilm.py
film = namedtuple('layer',['depth','index','name','active'])

def main():

  save = False

  #list of the parameters for the ITO Drude model
  #these will be used by indexLookup
  #these parameters are: restricted fit, both constraints 
  #(see notes 3/05)
  ito_omega_p = 2.73649315189e+14 #s^-1
  ito_epsilon_inf = 3.24
  ito_gamma = 5.96450643459e+12 #s^-1
  itoDrudeParams = [ito_omega_p, ito_epsilon_inf, ito_gamma]

  """
  stack = [film(1000, 1.399, 'SiO$_2$', False),
    film(50, itoDrudeParams, 'ITO', False),
    film(400,2.4+0.106j,'CQD',True),
    film(15,'au','Au',False)
    ]
  """

  stack = [film(6.0, 4+13.5j, 'NiCr approximation', False),
    film(300.0, 2.4+0.046j, 'HgTe CQD film', True),
    film(15.0,'au','Au', False)
    ]

  wls = array([4000]) #nm
  angles = (pi/180)*array([0])
  n_i = 1.4097
  n_f = 1.0
  pols = [0,1]

  paramSweep(wls, angles, n_i, n_f, pols, stack, 'depth', 2, 
    linspace(5.0,30.0,1000), save, 'dots.pdf')

"""this method will take a template stack and sweep a stack
parameter value.  It calls methods in thinFilm.py directly"""
def paramSweep(wls, angles, n_i, n_f, pols, stackTemplate,
  paramType, paramFilmNo, paramVals, save=False, saveFileName='default'):
  #here can select the parameter that will be swept
  #thicknesses of a film (nm)
  #filmDepths = linspace(0.0,500,1000)
  #range of indices to sweep
  #indexVals = linspace(2.0,2.8,1000)
  #range of wavelengths to sweep
  #plotwls = linspace(3000.0,5000.0,1000)
  #this will hold the integral of E^2 for each depth in filmDepths
  ESqResults = []

  #find the E-square integral for each depth using thinFilm
  for paramVal in paramVals:

    stack = stackTemplate

    if paramType == 'depth':
      stack[paramFilmNo] = film(paramVal, stack[paramFilmNo].index,
        stack[paramFilmNo].name, stack[paramFilmNo].active)
    elif paramType == 'index-real':
      stack[paramFilmNo] = film(stack[paramFilmNo].depth, 
        paramVal+(imag(stack[paramFilmNo].index))*(1j),
        stack[paramFilmNo].name, stack[paramFilmNo].active)
    elif paramType == 'index-imag':
      stack[paramFilmNo] = film(stack[paramFilmNo].depth, 
        (paramVal*(1j))+real(stack[paramFilmNo].index),
        stack[paramFilmNo].name, stack[paramFilmNo].active)
    # populates with list of (possibly) complex indices of refraction
    indices = [[thinFilm.indexLookup(layer.index, wl) for wl in wls]
     for layer in stack]

    t_angles = thinFilm.snell(indices, angles, n_i, n_f)
    (I,P) = thinFilm.genMatrices(stack, wls, angles, n_i,
      n_f, indices, t_angles)

    (T,R,TAvg,RAvg,E_i) = thinFilm.evalMatrices(stack, wls, angles, pols,
      n_i, n_f, indices, t_angles, P, I, addBulkT=True)

    (ESqInt,ESqIntAvg) = thinFilm.ESqIntEval(stack,wls,angles,
      pols,indices,E_i)

    for i in range(len(stack)):
      if stack[i].active == True:
        ESqResults.append(ESqIntAvg[i])

  ax = plt.axes()
  ax.plot(paramVals, ESqResults)
  ax.set_xlabel(stack[paramFilmNo].name+' layer '+paramType,fontsize=18)
  ax.set_ylabel('CQD layer $E^2$ integral (arbitrary units)',fontsize=18)
  ax.grid(True)
  plt.xticks(fontsize=15)layer
  plt.yticks(fontsize=15)
  if save:
    plt.savefig(saveFileName)
  plt.show()

if __name__ == '__main__':
  main()
