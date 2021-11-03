# Change Log

All notable changes to this project should be documented in this file.

The _Unreleased_ section should become the release once the release is ready
and the text can be used as part of the release description.
The pre-release historical sections are using dates and topics as section
headings instead of release numbers.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project loosely adheres to
[Semantic Versioning](https://semver.org/spec/v2.0.0.html),
i.e., major-minor-patch.
However, since there is _currently_ no need to keep the API stable,
keeping SemVer is focused on conceptual stability rather than code strict API
compatibility.
API changes which require small adjustments in the code can be introduced in
minor versions. For example, a new parameter which is optional from a user
point of view can be added in a minor version, although it may require changes
for the caller of a C++ function which has an additional required parameter
after the change. However, conceptual changes such as complete API rewrites
warant a new major version.
Please, open an issue if you need to change this policy.
For greater stability, consider using the interfaces instead of the core C++
library.

Interfaces to the core library which are in separate repositories all share
the same numbering line with this repository, so increasing version number
there increases version here and a new tag (and possibly also release) is
created in this repository although there are no changes in this repository or
the tagged revision was already tagged before.

## [Unreleased]

### Added

- Network kernel. (Vaclav Petras)
- Established dispersers from a cell are now tracked. (Chris Jones)
- Collection of general-use function for container handling for internal use.
- Template parameter for random number generator to Model and Simulation.

### Changed

- Virtual functions are used to switch beween kernels instead of an
  ever-growing if-else statement. (Vaclav Petras)
  * This is possibly faster and it has additional cost only
    when compared to hardcoded simple kernels.
  * Kernel parameters don't need to be provided anymore when the kernel is
    not used.
  * Opens door for adding many kernels in a flexible way.
- Anthropogenic kernel is now different from natural kernel and contains one
  network kernel. (Vaclav Petras)

### Fixed

- Fixed and completed documentation of function parameters. (Vaclav Petras)

## [1.1.0] - 2021-06-18

### Added

- Pest movement based on overpopulation (Vaclav Petras)
  * When cell contains too many pests, pests leave and move to a different cell.
- Total hosts added to run_step (Chris Jones)
  * Allows overpopulations to be based on total_hosts rather than total_populations
  * Allows movement to be based on total_hosts rather than total_populations
- Mortality timing is now more flexible (Chris Jones)
  * requires 2 new parameters mortality_frequency and mortality_frequency_n
- DOI, including a DOI for each release, is now available through Zenodo.

### Changed

- Model class internal attributes and functions are now protected instead of private
  to allow derived classes to access them for greater flexibility (Vaclav Petras).
- Model class kernels are now created and returned by protected functions (Vaclav Petras).
- Config class has now more defaults and subsequent setup now consistently fails when
  required values were not set (Vaclav Petras).
- Movement module has added parameters (`mortality_tracker_vector`, `exposed`, `resistant`, 
  `total_exposed`, and `suitable_cells`) (Chris Jones).
  * Now moves exposed and resistant populations
  * Now moves total_hosts
  * Now moves moratility tracked hosts
  * Adds a new suitable_cell if a movement creates a new location with hosts.
- Treatments now modify total hosts (Chris Jones)
- Citation updated to PoPS Frontiers paper (Chris Jones)
  * Jones, C., Jones, S., Petrasova, A., Petras, V., Gaydos, D., 
  Skrip, M., Takeuchi, Y., Bigsby, K., and Meentemeyer, R., 2021.
  Iteratively forecasting biological invasions with PoPS and a little help from 
  our friends.
  *Frontiers in Ecology and the Environment* 
  [DOI: 10.1002/fee.2357](https://doi.org/10.1002/fee.2357)

## [1.0.2] - 2020-10-09

- [Patch release of rpops](https://github.com/ncsu-landscape-dynamics/rpops/releases/tag/v1.0.2) (no changes in pops-core)

## [1.0.1] - 2020-09-11

- [Patch release of r.pops.spread](https://github.com/ncsu-landscape-dynamics/r.pops.spread/releases/tag/v1.0.1) (no changes in pops-core)

## [1.0.0] - 2020-09-10

- Susceptible-infected (SI) and susceptible-exposed-infected (SEI) host phases,
- Host mortality tracking,
- Host removal and pesticide application treatments,
- Host resistance based on pesticide application treatments,
- Treatments applied only to a ratio of hosts,
- Seasonal spread limited to certain months,
- Limited host movement (experimental),
- Yearly pest removal based on lethal temperature,
- Two different dispersal kernels (natural and anthropogenic),
- Cauchy and exponential radial dispersal kernels,
- Reduced stochasticity options and deterministic versions of kernels,
- Spread rate measurement,
- Quarantine and study area escape measurements,
- C++ API centered around two classes: pops::Model and pops::Config.

## 2020-08-27 - Version 1 preparations

### Changed

- The C++ library project renamed from PoPS to PoPS Core.
- Only the current weather coefficient is now passed to the model class.
- Raster class uses `int` by default for indexing.

### Fixed

- Numerous compiler warnings fixed. Code now compiles with -Wall -Wextra.

## 2020-08-12

### Fixed

- Exposed class now being removed with treatments
  * Added exposed class treatments to model class

## 2020-05-15 - Reducing stochasticity

### Added

- Stochasticity in Simulation can now be optionally reduced. (Vaclav Petras)
  * Simulation constructor takes additional parameters to disable
    stochasticity in generate and disperse functions, particularly in
    generating and establishing dispersers.
  * The disperse function takes additional parameter which is a fixed
    probability value when disperser is established.

## 2020-04-16 - SEI model

### Added

- Exposed host tracking with latency period to run SEI models. (Vaclav Petras)

### Changed

- Step numbering is now expected to start at zero. (Vaclav Petras)

## 2020-01-04 - Movement Module

### Added

- movement modeule to the simulation class to move hosts from one location to another. (Chris Jones)

## 2019-12-18 - More date functionality

### Added

- functions to check if a date is the last day or week of a month

## 2019-10-29 - More raster generalizations

### Added

- Raster can now hold data from some other object as long as the types
  and layout match. (Vaclav Petras)
- Raster Index and Number types exposed using typedefs.

### Changed

- Raster now has optional template parameter for the type of index used
  in the class. (Vaclav Petras)

### Fixed

- Signed and unsigned types are no longer arbitrarily mixed in the
  interface, e.g., `(int, int)` constructor is now `(Index, Index)`
  being consistent with `cols()` and `rows()` functions. (Vaclav Petras)

## 2019-09-05 - Dispersal kernel rewrite

### Added

- More complete list of operators supported by Raster. (Vaclav Petras)
- Rasters with different numerical types can be used in a single
  expression (Vaclav Petras)

### Changed

- Binary operators for rasters are now function templates
  instead of member functions (Vaclav Petras)

## 2019-08-11 - Dispersal kernel rewrite

### Added

- Radial exponential kernel and a uniform kernel (Vaclav Petras)
  * Any combination of kernels is now possible on the library level.
- Functions to convert from strings to kernel-related enums.

### Changed

- Dispersal kernel are organized as classes (functors). (Vaclav Petras)
  * The disperse function is now a function template in the Simulation.
  * Dispersal kernel is a function (functor) called from the Simulation.
  * The Simulation class does not have any kernel-related parameters
    or members.
- Usage of std::cout related to kernels was replaced by exceptions
  in the new code. (Vaclav Petras)
- Kernel-related enums are not strongly typed using enum class
  (Vaclav Petras)

### Fixed

- A call to abs function on the distance used for spread caused
  conversion to int. Now std::abs is used so double is used throughout.
  This changes simulation results. (Vaclav Petras)
- Resolution is now double instead of int. This changes results in cases
  when the resolution is not a integer number (even when it is
  so close to the closest higher integer that it is printed out without
  any decimal places by software considering it, possibly correctly,
  close enough from user perspective). This was also reported in the
  issue #26. (Vaclav Petras)

## 2019-07-11

### Added

- The Raster class has a new function data() which provides access to
  the underlying array. This follows the example of similar addition to
  std::vector in C++11. The Raster class is now handling the memory
  and provides raster algebra semantics, but its users have
  direct access to the underlying data. (Vaclav Petras)

### Changed

- The optional support for GRASS GIS reading and writing GRASS GIS
  rasters was removed from the Raster class. This needs to be now
  handled in the GRASS GIS module. (Vaclav Petras)
- Formerly private variables in Raster class are now only protected
  so that derived classes can have direct access to them, e.g., in the
  constructor. (Vaclav Petras)
- Private/protected variables in Raster class are now using the rows and
  columns terminology and trailing underscores. (Vaclav Petras)

### Fixed

- In some cases (depending on compiler and use of GRASS GIS), sqrt was
  not accessible in the simulation code leading to compilation error.
  Now std::sqrt is used. (Vaclav Petras)

## 2019-07-09

### Added

- Added ability to set number of days for the timestep instead of
  day, week, or month. This is useful if the pest has multiple generations
  in a year. (Chris Jones)
  
### Changed
  
- Made extension of the 52 week to be 8 or 9 (leap years) days instead of 
  7 to ensure that the next year starts on Jan 1st. This makes forecasting
  into the future more streamlined. (Chris Jones)

## 2019-06-14

### Added

- Added types of treatment application. Removing a portion for both
  susceptible and infected, like in SLF case, and a portion from
  susceptible all all infected in given cells, like in SOD case.
  (Vaclav Petras)

## 2018-11-17

- Added mortality function to simulation class

## 2018-10-31 - Minor Refactoring

### Added

- Added Season class to manage seasonality

### Changed

- Removed resolution from the Raster class

## 2018-09-20 - Variable Rename

### Added

- Created a style guide in a CONTRIBUTING file. (Vaclav Petras)

### Changed

- Renamed variables and functions to make the model more clear and easy
  to understand. (Chris Jones)
- Renaming and partial unification of styles across all files. (Vaclav Petras)
 - All is now in namespace called pops.
 - Private variables use trailing underscore if needed.
- Removed unused and legacy functions and types. (Vaclav Petras)

## 2018-08-02 - PoPS Model Separation

### Changed

- The general part of the simulation is now a separate library called
  PoPS. GRASS GIS module code was removed together with some parts of
  this file completely unrelated to the library. (Vaclav Petras)
- All files which were not supposed to be part of the repository
  were removed and purged from the history. However, all the relevant
  history is preserved. (Vaclav Petras)

## 2018-06-21 - Spotted Lanternfly

### Added

- The date class now supports also month increments. (Vaclav Petras)
- Critical temperature as the lowest temperature spores can survive
  in a provided month (Vaclav Petras)
 - Temperature rasters are used at a specified month to check against
   a provided critical temperature and if the condition is met,
   infected trees become susceptible again.

### Changed

- Img was replaced by a generalized Raster template class which can
  handle both integers and floating point numbers. (Vaclav Petras)

## 2018-06-04 - Mortality Addition

### Added

- Mortality (Vaclav Petras)
- Image class constructor taking another image using its dimensions
  and a provided value (Vaclav Petras)
- Multiply for image class is commutative (Vaclav Petras)

## 2018-03-26 - March 2018 Update

### Changed

- Explicitly include necessary standard C++ headers (Vaclav Petras)

## 2017-09-05 - September 2017 Update

### Added

- Long-range dispersal kernel (Anna Petrasova)
 - Events are recorded.
 - The affected points are exported as a vector map.
- Output probability of cell being infected (Vaclav Petras)
- Optionally output one run for series instead of an average (Vaclav Petras)

### Changed

- Spread of SOD based on a single species (Anna Petrasova)
 - Spread for UMCA and oak replaced by single species, assumed tanoak.

## 2017-01-28 - January 2017 Status

### Added

- Simplified weather inputs (Vaclav Petras)
 - Weather can be supplied as as a text file (spatially constant)
 - Weather can be supplied as one variable (non-spatial and non-temporal)

### Changed

- Efficiency improvements (Vaclav Petras)
 - Enable inlining of size getters of Img class which makes all_infected
   function much faster.
 - Move constructor and assignment operator added for cases when RVO
   is not applied.
 - Internal storage changed to one array (usually faster allocation).
- Code cleanup (Vaclav Petras)
 - Use same API style for Von Mises as for std lib distributions.
 - Creating data for Img outside of the object is avoided.
 - Indexing the Img done using operator ().
 - Using operators for all operations which fit semantically.
 - Remove unused variables from the code.
 - The 'using namespace' statement replaced by explicit use 'using' for
   string and other classes or objects.
- Date class API extended to provide readable comparison operators
  replacing usage of method with unclear name (Anna Petrasova)
- The cpl_string dependency was removed. (Vaclav Petras)
- The sporulation object is now seeded only once in the beginning.
  This changes the stochastic output of one run. (Anna Petrasova)

### Fixed

- Initialize memory for the sporulation object for the cases when it is
  zero cases to fix conditional jump which depends on uninitialised
  value. (Vaclav Petras)
- Memory allocation and deallocation is done by the right pair of
  malloc-free or new-delete. (Vaclav Petras)
- Compilation output and other non-repository files removed from the
  repository. (Vaclav Petras)
- Von Mises distribution concentration parameter is float not integer.
  (Vaclav Petras)
- Copy of GDAL code was removed from the repository, using system GDAL
  includes now. (Vaclav Petras)

[unreleased]: https://github.com/ncsu-landscape-dynamics/pops-core/compare/v1.1.0...HEAD
[1.1.0]: https://github.com/ncsu-landscape-dynamics/pops-core/compare/v1.0.2...v1.1.0
[1.0.2]: https://github.com/ncsu-landscape-dynamics/pops-core/compare/v1.0.1...v1.0.2
[1.0.1]: https://github.com/ncsu-landscape-dynamics/pops-core/compare/v1.0.0...v1.0.1
[1.0.0]: https://github.com/ncsu-landscape-dynamics/pops-core/tree/v1.0.0
