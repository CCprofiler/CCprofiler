# SECprofiler

## Setup
    $ git clone https://github.com/hafenr/SECprofiler
    $ cd SECprofiler
    $ R
    >> require(devtools)
    >> load_all()
    
## Workflow

After updating the roxygen docs issue the command

    >> document()

## Help

More information on how to use/develop this package can be found under:  [](http://r-pkgs.had.co.nz/intro.html)

## ELUTION PROFILING CONTAINER STRUCTURE
    Trace container: List with:
    Class: "Traces"
    $traces (data.table Column1: id Column 2:end numeric quant matrix)
    $trace_type (character, "peptide" or "protein")
    $traces_annotation (data.table, contains id column as key for merging and associated annotations in separate columns)
    $fraction_annotation
    

