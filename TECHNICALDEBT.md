# Technical Debt

All issues and imperfections of the code should be documented in this
file.

The goal is to declare the technical debt as described by
[Easterbrook 2014](http://doi.org/10.1038/ngeo2283) or Wikipedia
contributors in
(https://en.wikipedia.org/wiki/Technical_debt)[Technical debt].

The sections are based on sections in CHANGELOG. Entries should be
removed when resolved. Issue from tracker can be optionally linked
in an entry.

## 2018-06-13 - Spotted Lanternfly

### Design

- The decision if to increase by week or month is done with conditional
  (ternary) operator on one long line. Maybe enum and universal
  functions for increase and end of year would make it shorter.
- Season is just a std::pair, but a function to test if month is in
  range (or a range class) would make shorter code without the access to
  first and second members and two queries for month.

### Documentation

- Add to documentation.
- Update README.

## 2018-06-04 - Mortality Addition

### Efficiency

- Mortality-related objects are always created (and take memory).
- Mortality-related infection cohort always updated in spore spread and
  need an additional parameter. However, the update is just a repeating
  previous line (perhaps some container optionally wrapping two images
  would be useful).

### Design/Layout

- Several new lines of mortality model added to already long main
  function.

### Naming/User Interface

- Mortality cohorts may require better name since the individual can be
  potentially also a cohort or stand.

### Documentation

- Add to documentation.
- Update README.

## 2018-03-26 - March 2018 Update

### Style

- Some of the code in NetCDF ifdefs has bad structure (additional empty
  blocks).

## 2017-09-05 - September 2017 Update

### Dead Code

- Dead code for the specific multiple host species.

### Design/Layout

- Several new lines of non-trivial saving of escaping spores added to
  already long main function.

## 2017-01-28 - January 2017 Status

### Functionality

- Image class is only for integers but the weather is floating point,
  so it needs a different treatment in code than other data.
- The date class needs abstraction for the leap years and days in month.
- Formatting of text and parameters for the date class (leading zeros).
- Tests needed for the image class.

### Layout

- Create dedicated file for the direction enum and helper functions.

### Style

- Some of the Date class methods have still misleading names (e.g.,
  increased versus increase).
- Methods and functions should use underscores not camel case (using
  Python and GRASS GIS convention).
