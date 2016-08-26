# analytic-script
**PROJECT:** Analytic Solver for Thin Film Photovoltaics <br/>
**PROJECT PET NAME:** Every Love Story is a Ghost Story <br/>
**AUTHOR:** Zane Rossi with thinFilmClean.py based in code by John Roberts <br/>
**INSTIUTION:** University of Chicago, PGS Lab, 2016. <br/>

> *Giustizia mosse il mio alto fattore;* <br/>
> *fecemi la divina podestate,* <br/>
> *la somma sapïenza e 'l primo amore.* <br/><br/>
> *Dinanzi a me non fuor cose create* <br/>
> *se non etterne, e io etterno duro.* <br/>
> *Lasciate ogne speranza, voi ch'intrate.* <br/>

#######################################
## PRIMARY USES & SUMMARY:
* Contents of this project may be used in tandem to
  optimize thin film configurations as well as plot
  internal characteristics. This includes explicit
  and average e-fields, transmission and reflection
  coefficients, and e-field integrals across structures
  or particular features.
* File `thinFilmClean.py` comprises all methods relating
  to analytic thin film solver, using techniques defined in
  various photonics papers indicated in code. This file
  may be called explicitly to plot various collections of
  `film`s defined in head of code (i.e. in the `stack`). 
  This `stack` is easily adjustable, with measurements in nm, 
  and directions for input of refractive indices, allowing
  constant real, constant complex, csv-format, and Drude
  Model.
* File `autoOptimize.py` calls various methods in the
  aforementioned file, but does not use `stack` defined
  in file (if such `stack` exists). Rather
  `autoOptimize.py` defines its own `stack` with auxiliary
  parameters for minimum, maximum, and smallest  discrete
  step possible to fabricate for each layer. This `stack`
  is modified quietly by the code (specifically in method
  `evaluateCells`, where it is modified to `portStack` and made readable by 
  `thinFilmClean.py` methods). Code
  uses heuristic to populate space defined as n-dimensional and
  discrete, each lattice point a different combination of layer thicknesses. 
  This space, often very large, is seeded and quasi-greedily
  searched so that optimal or near optimal configurations
  may be found, in spite of direct calculation is performed for only a
  fraction of the entire space (in short, a metaheuristic). 
  Details of this optimization, its
  current parameters, and its recent advances use are defined below.
  `autoOptimize.py` does not use `thinFilmClean.py` graphic methods
  currently, but optimal `stack`s may be ported 
  by hand or called with relative ease by a curious user.
* File `paramVary.py` is deprecated.
* File `thinFilm.py` is deprecated and will eventually be superseded by
  `thinFilmClean.py`.
* Auxiliary files contain optical index information,
  or notes on code performance and routes for improvement. Study of their format
  will allow for the inclusion of arbitrary material parameters. 

#######################################
## THINFILMCLEAN METHODS:

### Properties Availible for Display
  * `plotTRangle`   : Displays transmission and reflection coefficients as a function of
  angle (currently not in heavy use).
  * `plotTRspectrum`: Displays transmission and reflection as a function of wavenumber
  (current method restricts this to average over polarizations, with wavelength plot commented out). 
  Settings may be easily manipulated in `TRSpectrumPlot`.
  * `EvalESQint`    : Displays both integral of the square of the e-field in active layers as
  plotted against spectrum (currently wavenumber), and power absorbed.
  * `normPwr`       : Displays normalized power absorbed in active layers, with the normalization
  coefficient currently hard-coded, but with ability to be calculated at cost of time if drastic change
  in representation of e-fields is to occur (not a trivial modification). This is commented out in situ.
  * `ODPlot`        : Plot of optical density of active layers with respect to spectrum (note dependence on
  `normPwr`'s normalization).
  * `plotESQ`       : Plots magnitude of the e-field within the device itself, with respect to distance
  along device. The plot axes are auto-adjusted, with units arbitrary and equal intensity across all 
  wavelengths (unless blackbody weighting is initiated).

### main:
  * The following descriptors appear in line order.
  * Options to save, as well as plot various characteristics
    of device, clearly named and labeled.
  * Various user-mutable constants are defined, with their uses
    and dimensions specified in situ. These will not
    need to be altered often if one is testing differing
    configurations under same external conditions
      *  `wls` is simply a linspace of wavelengths (in nm)
         to be iterated over. Average plots will average unweighted over all
         specified wavelengths. Blackbody radiation weighting
         may be turned on with specified maximum wavelength in code preamble.
         be aware these results are not directly commensurable with non
         blackbody results, as overall scale is arbitrary.
      *  `angles` is a non weighted list in **RADIANS**
  * Fundamental non-user-mutable constants are defined
  * Main stack is defined, comprising layers, each a `film`
    data type defined early in code (a named tuple). A `stack`
    is simply a list of `film`s, and only one `stack` may currently
    be run at a time (although alterations to code would not
    make queuing, or the sequential run of `stack`s, difficult. Film options are:
      * `layer` : generic name of layer; do not worry about this.
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
  * Rest of methods need not be altered, and result is deterministic
    based on parameters supplied above. These methods simply call the
    appropriate helper methods to do dirty work in traversing user's
    specified structure.

### evalField:
  * Given defined stack and spectral range, this method returns all information
    associated with the e-field resultant in device.
  * `INPUT`:
    * `stack`: defined above; a list of the `film` data type.
    * `wls`: defined above; a linspace of wavelengths in nm.
    * `angles`: defined above; a list of angles in **RADIANS**.
    * `pol`: polarization; `0` for p-polarization; `1` for s-polarization.
    * `n_i`: optical index of material preceding stack.
    * `n_f`: optical index of material following stack.
    * `indices`: for each wavelength in `wls`, and for each layer in `stack`,
      a two dimensional list giving (possibly complex, i.e. attenuation and 
      phase coefficient) optical index.
    * `t_angles`: indexed by layer, wavelength, and angle from their respective
      data structures: a three dimensional list of resultant trasnmitted 
      angles following Snell's law via the `snell` method.
    * `P`: propagation matrix: a two dimensional square matrix indicating phase
      and amplitude change of incident waves. Indexed in multi-dimensional
      array by layer, wavelength, angle, and polarization in that order.
    * `I`: interface matrix: a two dimensional square matrix indicating
      reflection and transmission between layers. Same indexing as `P`.
  * `RETURN`:
    * Tuple comprising the values below.
    * `E_0`: calculated incident and reflected electric field, given as a two
      dimensional vector whose elements indicate magnitude of right and left 
      moving plane-waveforms. 
    * `E_f`: calculated transmitted electric field, using same protocol as
      above.
    * `E_i`: list of e-field present at forward edge of indexed layer, which
      may then be propagated forward by methods like `EsqPlot` for display.
  * This method as well as `evalTR` follow closely the matrix transfer method
    defined in Pettersson et al. (JAP 1999). 

### esqIntEval:
  * `INPUT`:
    * `stack`: defined above.
    * `wls`: defined above.
    * `angles`: defined above.
    * `pols`: defined above.
    * `indices`: defined above.
    * `E_i`: defined above in `evalField` return.
  * `RETURN`:
    * Tuple comprising the values below.
    * `ESqInt`: Numerical integral from interpolated function given by 
      parameters contained in e-field definition. Still indexed by layer number,
      wavelength, angle, and polarization in that order. 
    * `ESqIntAvg`: Average of above data structure across all indices save
      layer number. That is, the general performance of each layer.
    * `PAbsd`: Power absorbed as defined by external technique referenced in 
      situ. Indexed through the same means as `ESquInt`.
    * `PAbsdAvg`: Average of above data structure in the same manner as
      `ESqIntAvg`.

#######################################
## AUTOOPTIMIZE METHODS:

### main:
  * Easily adjustable constants are defined
    * `MAX_STEP`: number of propagations desired (similar to iterations 
      of flood fill, save with accompanying heuristic that allows non linear
      scaling)
    * `ACTIVE_LAYER`: index of layer in stack, defined in autoOptimize.py, that
      one wishes to maximize the e-field integral over. Currently this is 
      implemented for one active layer per run, although this can be changed
      with relative ease in accompanying methods (namely, `evaluateCells`)
    * `DIVISIONS`: seeds space with (`DIVISIONS + 1`) active cells evenly spaced
      along each dimension, taking in essence the cartesian product
      of these one dimensional intervals, and thus seeding a lattice. Higher
      value increases search time per propagation, but gives shorter path to
      possible optima. The best heuristic is to have a displacement between 
      far apart seeds be not more than an order of magnitude larger than 
      `MAX_STEP`.
    * `THRESHOLD`: Cutoff value at which iteration will halt in a manner 
      analogous to `MAX_STEP`, or the space being exhausted of unvisited cells.
  * Search space is then instantiated to proper size according to stack, and
    values specified above.
  * Code propagates either `MAX_STEP` times, or until there are no active cells
    (if you feel you are running low on active cells too soon, see 
    `selectedNeighbors` subsection of propagate method). Output is list of the 
    best cell per propagation, listing coordinates first and then fitness of
    said cell. Final output give numerical information on traversal of space.

### propagate:
  * Welcome (no please, stay)! This is the workhorse method of the heuristic 
    search algorithm. All methods called via this method (save `evaluateCells`) 
    will be described in short for now, and given full entries later. They will 
    not in standard cases need alteration, serving primarily as helpers.
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
    * Tuple comprising `newState`, `newActive`, and `newOptimum`.
    * These are modified **COPIES** of original input, allowing for sequential
      feeding of method's output into itself. See INPUT for data types.
  * `propagate` makes copies of its inputs, and should **NEVER** modify them. 
    when modifying propagate, two main section should be accessed first, forming
    the core of the heuristic:
      * heuristic variables:
        * `threshold`: magnitude of a coordinate difference between centroid 
          and barycenter which necessitates traveling along said coordinate
          in signed direction.
        * `minActive`: minimum required cells needed so that `threshold` comes
          into effect. Thus, when there are few active cells, propagation
          direction will become non-selective, resulting in a flood-fill until
          `minActive` is exceeded.
        * `selectivity`: argument for the method `getPercentage`, and whose
           sensible range will obviously vary with the definition of the
           function called by that method.
        * `greedPopulation`: maximum population during which the greedy
          hill-climbing protocol remains in effect. Should be set to something
          on the same order of magnitude as `minActive`, with the preference
          being for slightly larger. If this is set to `0`, greediness halts,
          while if it is set too high, ridges and local maxima will slow
          propagation speed dramatically.
        * `percentage`: variable fed into `getFittestPercent`, with a range
          between `0` and `100` inclusive. Currently set by a call to the
          aforementioned `getPercentage`, appropriately. 
      * `selectedCells`    : those cells which are chosen to determine the
        fitness of some of their neighbors. these cells are chosen from the
        current active cells. Note, even if cells are not chosen to be within
        `selectedCells`, they should **REMAIN** active for possible future use. 
        currently, the algorithm is set to select varying percentages of the
        active cells based on their fitness alone (the percentage has now been
        updated to reference an external method, `getPercentage`, which allows
        for a smooth function (in this case an overshooting logistic) to
        modulate how many cells are chosen). This selection criteria can be
        modified as seen fit, but as a rule of thumb, the number of cells
        selected per propagation should not dip below 10 or so, lest bottleneck
        effects occur. 
      * Between these two selections of cells to propagate from, and neighbors
        of said cells to pick, an auxiliary process has been inserted, with the
        intent of solving 'gapping' of the space when the barycenter heuristic 
        is too diffuse. When activated, the insert will, when the size of 
        `selectedCells` falls below a certain value (namely `greedPopulation`),
        choose the best cell among `selectedCells` and propagate in every
        available direction from that cell. This will provide a hill climber
        effect best used when algorithm has drawn close to a desired extrema. In
        conjunction with each other, this and the constraints of 
        `selectedNeighbors` may be joined to produce good global *and* local 
        performance. Updates are pending as to the best means of optimizing
        this particular portion of the meta-heuristic.
      * `selectedNeighbors`: those cells which are traveled to from above
        `selectedCells`. Currently the code chooses which direction to travel
        based on barycentric methods and threshold values. Note, **ALWAYS** 
        check `hasTraversed` of cell before propagating to it. This incorporates 
        both boundaries and previously traversed cells. The code will throw 
        errors if you do not, and you will lament your previous blasé existence.
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
  * This is the workhorse method for porting into `thinFilmClean.py` in order
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

### hasTravsered:
  * Checks if input cell has been traversed or is invalid choice.
  * `INPUT`:
    * `currentState`: defined above.
    * `coords`: a coordinate tuple of dimension equal to the length of `stack`.
  * `RETURN`:
    * Boolean: `True` if cell has value `1`,`2`, **OR** coordinate is invalid.
      `False` otherwise (namely cell has value `0`).

### travserseDim:
  * Propagates from cell in specified direction a specified distance, returns
    new state and new active cell.
  * `INPUT`:
    * `currentState`: defined above.
    * `activeCell`: coordinate tuple of central cell to be traveled from.
    * `dim`: between `0` and `len(stack) - 1` inclusive. Specifies axis of
      travel. Method checks if this is a valid choice, and will throw error
      if not.
    * `step`: integer which specifies distance of travel along `dim`. Method
      checks if this is in conjunction with `dim` is a valid position, 
      and will throw error if not.
  * `RETURN`:
    * Tuple consisting of those values found below.
    * `newState`: state with `activeCell` position changed to `1` and new site
      of activity changed to `2`.
    * `newCoord`: tuple coordinate of new active cell from traversal.

### getNeighborCells:
  * Returns **ALL** coordinates of neighboring cells adjacent by a face in 
    n-dimensional search space. This includes **INVALID** coordinates not in 
    space.
  * `INPUT`:
    * `coords`: coordinates of cell around which neighbors are pulled.
  * `RETURN`:
    * List of tuple coordinates of possibly invalid neighbors. This is a lazy
    neighbor search because **EVERYTHING** must be run through `hasTraversed`
    in `propagate`.

### getCoordinateCenter:
  * Returns weighted or unweighted centroid of a list of cell positions and 
    fitnesses.
  * `INPUT`:
    * `fitnessMap`: defined above; a list of tuples of tuple coordinates and 
    fitnesses.
    * `isWeighted`: boolean `True` is barycenter is desired, `False` if centroid
    is desired.
  * `RETURN`:
    * Tuple coordinate of weighted or unweighted center determined by 
    `isWeighted`.

### getCenterOffset:
  * Small helper method to return difference between centroid and barycenter.
  That is, the coordinate `getCoordinateCenter(fitnessMap, True)` less 
  `getCoordinateCenter(fitnessMap, True)`. For ease in heuristic. 
  * `INPUT`:
    * `fitnessMap`: defined above; a list of tuples of tuple coordinates and 
    fitnesses.
  * `RETURN`:
    * Tuple representing vector difference between centroid and barycenter: 
      specifically, barycenter less centroid.

### getFittestPercent:
  * Given `fitnessMap` and a percentage, return a revised `fitnessMap`
    containing only the top percentage of its original elements.
  * `INPUT`:
    * `fitnessMap`: defined above.
    * `percent`: a numerical value between `1` and `100` (complain all you want
      about non-unity percentage base).
  * `RETURN`:
    * List of tuples of tuple coordinates and fitness values. This must 
      **ALWAYS** be unpacked before direct insertion into `selectedCells` in 
      `propagate`.

### prettyPrint:
  * For quick debugging of the heuristic, this allows the display of the
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

