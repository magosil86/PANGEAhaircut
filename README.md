# PANGEAhaircut

Code package to post process IVA contigs before they are used to assemble a de novo reference genome for mapping of short reads.

## Authors:

* Oliver Ratmann <oliver.ratmann@imperial.ac.uk>

# Installation

The easiest way to install `PANGEAhaircut` is to use the `devtools` package:

```r
# install.packages("devtools")
library(devtools)
install_github("olli0601/PANGEAhaircut")
```

To work with this `R` package:
fire up `R`, and type 

```r
library(help=PANGEAhaircut)
```

# Content:

* Function `haircutwrap.get.cut.statistics` to compute descriptive statistics that are used to calculate the call probability. This function can be called from within R, from the command line, and on an HPC system. See `?haircutwrap.get.cut.statistics` for help.

* Function `haircutwrap.get.call.for.PNG_ID` to call 10 bp long chunks of cut/raw contigs, based on descriptive statistics of the contigs and the consensus sequence. This function can be called from within R, from the command line, and on an HPC system. See `?haircutwrap.get.call.for.PNG_ID` for help.
