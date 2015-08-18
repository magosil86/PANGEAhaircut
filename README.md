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

* Function `cmd.haircut.pipeline` that combines the above two steps. This produces a UNIX bash script that can be submitted to an HPC system. E.g. typing in R:

```
#DATA		<- SET THIS DIRECTORY
indir			<- paste(DATA, 'contigs_150408_wref', sep='/' )
outdir		<- paste(DATA, 'contigs_150408_wref_cutstat', sep='/' )		
cat(cmd.haircut.pipeline(indir, outdir, batch.n=200, batch.id=2))
```

produces:


```
#######################################################
#
# start: run haircut.pipeline
#
#######################################################

CWD=$(pwd)
echo $CWD
mkdir -p $CWD/cutstat_15-08-18-13-03-29
mkdir -p $CWD/call_15-08-18-13-03-29

 
#######################################################
# start: run haircut.cutstat.contigs.Rscript
####################################################### 
echo 'run Rscript /Users/Oliver/Library/R/3.1/library/PANGEAhaircut/haircut.cutstat.contigs.Rscript'
 Rscript /Users/Oliver/Library/R/3.1/library/PANGEAhaircut/haircut.cutstat.contigs.Rscript -indir=/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2015/2015_PANGEA_haircut/contigs_150408_wref -outdir=$CWD/cutstat_15-08-18-13-03-29 -batch.n=200 -batch.id=2 
echo 'end Rscript /Users/Oliver/Library/R/3.1/library/PANGEAhaircut/haircut.cutstat.contigs.Rscript'
#######################################################
# end: run haircut.cutstat.contigs.Rscript
#######################################################

 
#######################################################
# start: run haircut.call.contigs.Rscript
####################################################### 
echo 'run Rscript /Users/Oliver/Library/R/3.1/library/PANGEAhaircut/haircut.call.contigs.Rscript'
 Rscript /Users/Oliver/Library/R/3.1/library/PANGEAhaircut/haircut.call.contigs.Rscript -indir.st=$CWD/cutstat_15-08-18-13-03-29 -indir.al=/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2015/2015_PANGEA_haircut/contigs_150408_wref -outdir=$CWD/call_15-08-18-13-03-29 
echo 'end Rscript /Users/Oliver/Library/R/3.1/library/PANGEAhaircut/haircut.call.contigs.Rscript'
#######################################################
# end: run haircut.call.contigs.Rscript
#######################################################

mv $CWD/call_15-08-18-13-03-29/* /Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2015/2015_PANGEA_haircut/contigs_150408_wref_cutstat
rm -d $CWD/call_15-08-18-13-03-29 $CWD/cutstat_15-08-18-13-03-29

#######################################################
#
# end: run haircut.pipeline
#
#######################################################
```
