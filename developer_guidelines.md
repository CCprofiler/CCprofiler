# Developer guidelines

## General styleguide
* Use 2 spaces for indenting. **No tabs**.
* No lines longer than 80 characters.
* Function names: Use camelCaps with initial lowercase, then alternate case between words. Use '.' only at beginning of non-exported functions.
* Variable names: lower case with '\_' to separate words
* Class names: camelCaps
* Use <- not = for assignment.
* Use “##” to start full-line comments. Indent at the same level as surrounding code.
* New classes should be defined in 'AllClasses.R'.

## Package specific guidelines
* Use data.table and not data.frame
* Functions should import tables from files or as R data.table
* Use XXX for parallelization
* To import fuctions from other packages use @importFrom
* Write description, example and test for all new functions
