## R CMD check results
There were no ERRORs or WARNINGs.

There were 2 NOTEs:

* checking for future file timestamps ... NOTE
  unable to verify current time
  
I checked online for solutions. It seems that there is an issue with the external web resource which extracts current time.

* checking R code for possible problems ... NOTE
  LinProj: multiple local function definitions for 'varlink' with
    different formal arguments
  VarEst: multiple local function definitions for 'varlink' with
    different formal arguments'

In both functions 'LinProj' and 'VarEst' the local function 'varlink' is defined conditionally based on distributions of interest. This 'multiple definition' helps simplify the code.
