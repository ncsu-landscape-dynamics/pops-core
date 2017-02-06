# Change Log

All notable changes to this project should be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/).

## 2017-01-28 - January 2017 status

### Added

- GRASS GIS interface added (Vaclav Petras)
 - Main spatial inputs handled through GRASS GIS libraries so that
   different extents and resolutions can be used. (Not implemented for
   the weather time series.)
 - Text and numerical inputs hardcoded in defines replaced by
   description using GRASS GIS library. Command line parameters are now
   parsed, checked including dependencies, and it is possible to open
   GUI.
 - All hardcoded parameters replaced by command line parameters.
 - Basic descriptions for all parameters.
 - GDAL input and output is in the code.
- Multiple stochastic runs and parallelization (Vaclav Petras)
 - Optionally outputs also stdandard deviation.
 - Aggregates series output from runs.
 - Creates OpenMP threads for chunks of weeks.
 - User can set number of threads from command line (or GUI).
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
 - Use the interal data directly to write using GDAL output.
- Code cleanup (Vaclav Petras)
 - Time measurement code removed (use external tool such as time instead).
 - Commented-out code for seeding by time removed.
 - Use same API style for Von Mises as for std lib distributions.
 - Creating data for Img outside of the object is avoided.
 - Indexing the Img done using operator ().
 - Using operators for all operations which fit semantically.
 - Remove unused variables from the code.
 - The 'using namespace' statement replaced by explicit use 'using' for
   string and other classes or objects.
- Date class API extended to provide readable comparison operators
  replacing usage of method with unclear name (Anna Petrasova)
- Random seed for generators is required or must be explicitly generated.
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
