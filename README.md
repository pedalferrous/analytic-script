# analytic-script
**PROJECT:** Analytic Solver for Thin Film Photovoltaics <br/>
**PROJECT PET NAME:** Every Love Story is a Ghost Story <br/>
**AUTHOR:** Zane Rossi with thinFilmClean.py based in code by John Roberts <br/>
**INSTIUTION:** University of Chicago, PGS Lab, 2016. <br/>

> Giustizia mosse il mio alto fattore;
> fecemi la divina podestate,
> la somma sapïenza e'l primo amore.
> Dinanzi a me non fuor cose create
> se non etterne, e io etterno duro.
> Lasciate ogne speranza, voi ch'intrate'

#######################################
## PRIMARY USES & SUMMARY: 
* Contents of this project may be used in tandem to
  optimize thin film configurations as well as plot
  their internal characteristics, including explicit
  and average e-fields, transmission and reflection
  coefficients, and e-field integrals across structures
  or particular features.
* File `thinFilmClean.py` contains all methods relating
  to analytic thin film solver, using methods defined in
  various photonics papers indicated in code. This file
  may be called explicitly to plot various collections of
  `film`s defined in head of code (`stack`). This 'stack'
  is easily adjustable, with measurements in nm, and
  directions for input of refractive indices, allowing
  constant real, constant complex, csv-format, and Drude
  Model.
* File `autoOptimize.py` calls various methods in the
  aforementioned file, but does not use `stack` defined
  in file (if such `stack` exists). Rather
  `autoOptimize.py` defines its own `stack` with auxiliary
  parameters for minimum, maximum, and smallest  discrete
  step possible to fabricate for each layer. This `stack`
  is modified quietly by the code (specifically in method
  `evaluateCells`, modified into `portStack` to be readable by 
  `thinFilmClean.py` methods). Code
  uses heuristic to populate space defined as n-dimensional and
  discrete, consisting of all possible combinations of layer thicknesses. 
  This space, often very large, is seeded and pseudo-greedily
  searched so that optimal or near optimal configuration
  may be found, although direct calculation is performed for only a
  fraction of the entire space. Details of this optimization, its
  current parameters, and its advances use are defined below.
  `autoOptimize.py` does not use `thinFilmClean.py` graphic methods
  currently, but optimal `stack`s may be either ported 
  by hand, or called with relative ease.
* file `paramVary.py` is deprecated.
* file `thinFilm.py` is deprecated and will eventually be replaced by
  `thinFilmClean.py`.
* auxiliary files contain optical index information,
  or notes on code performance and routes for improvement.

#######################################
## THINFILMCLEAN METHODS:

### main:
  * the following descriptors appear in order seen in method
  * options to save, as well as plot various characteristics
    of device, clearly named and labeled
  * various alterable constants are defined, with their uses
    and dimensions defined in code. these will not
    need to be altered often if one is testing differing
    configurations under same external conditions
      *  `wls` is simply a linspace of wavelengths (in nm)
         to be iterated over. Average plots will average unweighted over all
         specified wavelengths. Blackbody radiation weighting
         may be turned on with specified maximum wavelength in code preamble.
         be aware these results are not directly commensurable with non
         blackbody results, as overall scale is arbitrary.
      *  `angles` is a non weighted list in **RADIANS**
  * fundamental constants are defined
  * main stack is defined, comprising layers, each a `film`
    data type defined early in code (a named tuple). A `stack`
    is simply a list of `film`s, and only one `stack` may currently
    be run at a time (although alterations to code would not
    make the sequential run of `stack`s difficult. Film options are:
      * `layer` : generic name of layer, do not worry about this.
      * `depth` : depth in nm of said layer, integer or float.
      * `index` : index of refraction, either constant, name of external file
                  or Drude parameters defined earlier in code (currently only
                  valid for ITO. 
      * `name`  : unique name for layer, for display and reference purposes.
                  That is, element of named tuple may be referenced by index 
                  or name of field, in this case 'name.'
      * `active`: determines whether this layer's individual e-field wil be 
                  integrated, and thus plotted or otherwise displayed. leave
                  boolean True if this is desired.
  * rest of methods need not be altered, and result is deterministic
    based on parameters supplied above. These methods simply call the
    appropriate helper methods to do dirty work in traversing your
    specified structure.

### evalField:
  * given defined stack and spectral range, this method returns all information
    associated with the e-field resultant in the device.
  * `INPUT`:
    * `stack`: defined above; a list of the `film` data type.
    * `wls`: defined above; a linspace of wavelengths in nm.
    * `angles`: defined above; a list of angles in **RADIANS**.
    * `pol`: polarization; 0 for p-polarization; 1 for s-polarization.
    * `n_i`: optical index of material preceding stack
    * `n_f`: optical index of material following stack
    * `indices`: for each wavelength in `wls`, and for each layer in `stack`,
      a two dimensional list giving (possibly complex, i.e. attenuation and 
      phase coefficient) optical index.
    * `t_angles`: for each layer, wavelength, and angle from their respective
      data structures, a three dimensional list of resultant angles following
      Snell's law via the `snell` method.
    * `P`: propagation matrix; a two dimensional square matrix indicating phase
      and amplitude change of incident waves. Indexed in multi-dimensional
      array by layer, wavelength, angle, and polarization in that order.
    * `I`: interface matrix: a two dimensional square matrix indicating
      reflection and transmission between layers. Same indexing as `P`.
  * `RETURN`:
    * `E_0`: calculated incident and reflected electric field, given as a two
      dimensional vector whose elements indicate magnitude of right and left 
      moving plane-waveforms. 
    * `E_f`: calculated transmitted electric field, using same protocol as
      above.
    * `E_i`: list of e-field present at forward edge of indexed layer, which
      may then be propagated forward by methods like EsqPlot for display.
  * this method as well as `evalTR` follow closely the matrix transfer method
    defined in Pettersson et al. (JAP 1999). 

#######################################
## AUTOOPTIMIZE METHODS:

### main:
  * easily adjustable constants are defined
    * `MAX_STEP`: number of propagations desired (similar to iterations 
      of flood fill, save with accompanying heuristic that allows non linear
      scaling)
    * `ACTIVE_LAYER`: index of layer in stack, defined in autoOptimize.py, that
      one wishes to maximize the e-field integral over. Currently this is 
      implemented for one active layer only, although this can be changed
      with relative ease in accompanying methods (namely, evaluate Cells)
    * `DIVISIONS`: seeds space with (`DIVISIONS` + 1) active cells evenly spaced
      along each dimension of the search space, taking the cartesian product
      of these one dimensional intervals, and thus seeding a lattice. Higher
      value increases search time per propagation, but give shorter path to
      possible optima
  * search space is then instantiated to proper size according to stack, and
    values specified above.
  * code propagates either `MAX_STEP` times, or until there are no active cells
    (if you feel you are running low on active cells too soon, see 
    selectedNeighbors subsection of propagate method). Output is list of best
    cells per propagation, listing coordinates first, and then fitness of
    said cell. Final output give numerical information on traversal of space.

### propagate:
  * this is the workhorse method of the heuristic search algorithm. all
    methods called via this method (save `evaluateCells`) will be described
    in short for now, and given full entries later. They will not in standard
    cases need alteration, serving primarily as helpers.
  * `INPUT`: 
    * `currentState`: numpy array containing all information on unvisited (0)
      cells, visited but inactive (1) cells, and visited active (2) cells.
      Dimension of array is length of stack, and is automatically generated
      in main. Solely an integer array. 
    * `activeCells` : list of tuples denoting coordinates of all cells in
      `currentState` with value 2. Dimension is same as n discussed above.
    * `fitnessMap`  : list of tuples. Each tuple contains two (2) elements,
      namely a tuple and a float. This interior tuple is a coordinate of the
      same form populating `activeCells`. Float is accompanying fitness of
      cell. All heuristic decisions are made off this data structure.
    * `optimum`     : single element of fitnessMap, largest float parameter.
  * `RETURN`:
    * tuple comprising `newState`, `newActive`, and `newOptimum`.
    * these are modified **COPIES** of original input, allowing for sequential
      feeding of method's output into itself. See INPUT for data types.
  * `propagate` makes copies of its inputs, and should **NEVER** modify them. when
    modifying propagate, two main section should be accessed first, forming
    the core of the heuristic:
      * `selectedCells`    : those cells which are chosen to determine the
        fitness of some of their neighbors. these cells are chosen from the
        current active cells. Note, even if cells are not chosen to be within
        `selectedCells`, they should **REMAIN** active for possible future use. 
        currently, the algorithm is set to select varying percentages of the
        active cells based on fitness alone. This selection criteria can be
        modified as seen fit.
      * `selectedNeighbors`: those cells which are traveled to from above
        `selectedCells`. Currently the code chooses which direction to travel
        based on barycentric methods and threshold values. Note, **ALWAYS** check
        `hasTraversed` of cell before propagating to it. This incorporates both
        boundaries and previously traversed cells. The code will throw errors
        if you do not, and you will lament your previous blasé existence.
        Currently, simply modifying how `properDim` is selected is an easy path
        for modifying the heuristic.
          - `properDim` is boolean n-element array detailing which of the n
            dimensions are optimal to propagate towards
          - `signDim` is an integer n-element array indicating whether forward
            backward, or both in specified optimal dimension should be chosen
      * Varying these two modes of selection, a majority of heuristic
        alterations can be comfortably made. If you choose to go beyond these
        methods, be aware of variable mutation, checking cell validity, and
        maintaining active cells not otherwise being accessed. Future room
        for modifications centers in better maintenance of active cells.

### evaluateCells:
  * this is the workhorse method for porting into `thinFilmClean.py` in order
    to fill `fitnessMap`. all extra-file communication should be contained in
    this method (promise me).
  * `INPUT`:
    * `activeCells`: defined above.
    * `activeLayer`: layer in stack by zero index whose value is to be somehow
      manipulated by evaluateCells to be largest if said coordinate is deemed
      optimal.
  * `RETURN`:
    * `fitnessMap`: defined above.
  * `evaluateCells` calls on elements in stack to create a new list of films
    compatible with the `portFilm` data type defined in the preamble. This data
    type corresponds to the `film` data type of `thinFilmClean.py`. This new
    stack is then ported through `thinFilmClean.py` methods to return requisite
    integral for e-field in `activeLayer` (current setup for optimization).
    If one for example wished to minimize transmission, one could feed some
    monotone decreasing function of transmission into the second element of
    the tuples then fed into fitnessMap.

### prettyPrint:
  * for quick debugging of the heuristic, this allows the display of the
    `dim1`-`dim2` plane with `0 -> " ", 1 -> ".", 2 -> "#"`. When space is taken
    to be two dimensional (e.g. all but two layers of stack are given set
    measurement), this allows visualization of whole space and activity of 
    heuristic in traversal. Method stands on its own due to clumsiness of
    slicing ndarray data type.
  * `INPUT`:
    * `currentState`: defined above.
    * `dim1`, `dim2`: integers detailing slicing dimensions. They should be
      **DIFFERENT**, and between 0 and n - 1 inclusive for n-layer stack.
  * `RETURN`: IO, prints to terminal.

