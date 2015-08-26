# PANGEAhaircut

Code package to post process IVA contigs before they are used to assemble a de novo reference genome for mapping of short reads.

## Authors:

* Oliver Ratmann <oliver.ratmann@imperial.ac.uk>

# Installation

This `R` package requires third party code before installation:

* MAFFT http://mafft.cbrc.jp/alignment/software/

After these dependencies are installed, the easiest way to install `PANGEAhaircut` is via the `devtools` package:

```r
# install.packages("devtools")
library(devtools)
install_github("olli0601/PANGEAhaircut")
```

After installation, you may need to make the installed Rscripts executable:

```
chmod 775 /export71/home/or105/R/x86_64-unknown-linux-gnu-library/3.2/PANGEAhaircut/*Rscript
```

To work with this `R` package:
fire up `R`, and type 

```r
library(help=PANGEAhaircut)
```

# Content:
* Function `cmdwrap.align.contigs.with.ref` to align cut and raw IVA contigs against a set of reference HIV sequences. This function can be called from within R, from the command line, and on an HPC system. See `?cmdwrap.align.contigs.with.ref` for help.

* Function `cmd.haircut.check.alignment` to check the alignment of cut and raw IVA contigs against a set of references. See `?cmd.haircut.check.alignment` for help.

* Function `haircutwrap.get.cut.statistics` to compute descriptive statistics that are used to calculate the call probability. This function can be called from within R, from the command line, and on an HPC system. See `?haircutwrap.get.cut.statistics` for help.

* Function `haircutwrap.get.call.for.PNG_ID` to call 10 bp long chunks of cut/raw contigs, based on descriptive statistics of the contigs and the consensus sequence. This function can be called from within R, from the command line, and on an HPC system. See `?haircutwrap.get.call.for.PNG_ID` for help.

* Function `cmd.haircut.pipeline` that combines the above steps. This produces a UNIX bash script that can be submitted to an HPC system. E.g. typing in R:

```
#DATA		<- SET THIS DIRECTORY
indir.cut		<- paste(DATA, 'contigs_150408_unaligned_cut', sep='/' )
indir.raw		<- paste(DATA, 'contigs_150408_unaligned_raw', sep='/' )
outdir		<- paste(DATA, 'contigs_150408_model150816a', sep='/' )		
cat(cmd.haircut.pipeline(indir.cut, indir.raw, outdir, batch.n=200, batch.id=1))
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
mkdir -p $CWD/algnd_15-08-26-11-52-27
mkdir -p $CWD/cutstat_15-08-26-11-52-27
mkdir -p $CWD/call_15-08-26-11-52-27

#######################################################
# start: run cmdwrap.align.contigs.with.ref
#######################################################
mafft --anysymbol --add /Users/Oliver/Dropbox\ \(Infectious\ Disease\)/OR_Work/2015/2015_PANGEA_haircut/contigs_150408_unaligned_cut/12559_1_1_hiv_cut.fasta --auto $CWD/algnd_15-08-26-11-52-27 > /Users/Oliver/Dropbox\ \(Infectious\ Disease\)/OR_Work/2015/2015_PANGEA_haircut/contigs_150408_unaligned_raw/12559_1_1_cut_wRefs.fasta
...
...
...
mafft --anysymbol --add /Users/Oliver/Dropbox\ \(Infectious\ Disease\)/OR_Work/2015/2015_PANGEA_haircut/contigs_150408_unaligned_cut/13557_1_84_hiv_cut.fasta --auto $CWD/algnd_15-08-26-11-52-27 > /Users/Oliver/Dropbox\ \(Infectious\ Disease\)/OR_Work/2015/2015_PANGEA_haircut/contigs_150408_unaligned_raw/13557_1_84_cut_wRefs.fasta
#######################################################
# end: run cmdwrap.align.contigs.with.ref
#######################################################

 
#######################################################
# start: run haircut.check.alignment.Rscript
####################################################### 
echo 'run Rscript /Users/Oliver/Library/R/3.1/library/PANGEAhaircut/haircut.check.alignment.Rscript'
 Rscript /Users/Oliver/Library/R/3.1/library/PANGEAhaircut/haircut.check.alignment.Rscript -indir=$CWD/algnd_15-08-26-11-52-27 -outdir=$CWD/algnd_15-08-26-11-52-27 
echo 'end Rscript /Users/Oliver/Library/R/3.1/library/PANGEAhaircut/haircut.check.alignment.Rscript'
#######################################################
# end: run haircut.check.alignment.Rscript
#######################################################

 
#######################################################
# start: run haircut.cutstat.contigs.Rscript
####################################################### 
echo 'run Rscript /Users/Oliver/Library/R/3.1/library/PANGEAhaircut/haircut.cutstat.contigs.Rscript'
 Rscript /Users/Oliver/Library/R/3.1/library/PANGEAhaircut/haircut.cutstat.contigs.Rscript -indir=$CWD/algnd_15-08-26-11-52-27 -outdir=$CWD/cutstat_15-08-26-11-52-27 
echo 'end Rscript /Users/Oliver/Library/R/3.1/library/PANGEAhaircut/haircut.cutstat.contigs.Rscript'
#######################################################
# end: run haircut.cutstat.contigs.Rscript
#######################################################

 
#######################################################
# start: run haircut.call.contigs.Rscript
####################################################### 
echo 'run Rscript /Users/Oliver/Library/R/3.1/library/PANGEAhaircut/haircut.call.contigs.Rscript'
 Rscript /Users/Oliver/Library/R/3.1/library/PANGEAhaircut/haircut.call.contigs.Rscript -indir.st=$CWD/cutstat_15-08-26-11-52-27 -indir.al=/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2015/2015_PANGEA_haircut/contigs_150408_wref -outdir=$CWD/call_15-08-26-11-52-27 
echo 'end Rscript /Users/Oliver/Library/R/3.1/library/PANGEAhaircut/haircut.call.contigs.Rscript'
#######################################################
# end: run haircut.call.contigs.Rscript
#######################################################

mv $CWD/call_15-08-26-11-52-27/* /Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2015/2015_PANGEA_haircut/contigs_150408_model150816a
rm -d $CWD/algnd_15-08-26-11-52-27 $CWD/call_15-08-26-11-52-27 $CWD/cutstat_15-08-26-11-52-27

#######################################################
#
# end: run haircut.pipeline
#
#######################################################
```
