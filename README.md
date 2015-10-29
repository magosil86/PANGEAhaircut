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

* Function `cmd.haircut.pipeline` that combines the above steps. This produces a UNIX bash script that can be submitted to an HPC system. There are several options to invoke the pipeline, see `?cmd.haircut.pipeline` for help. One popular option is to invoke the pipeline on a single IVA contig from the command line:
```
Rscript haircut.pipeline.Rscript -in.raw=path/raw_contigs.fasta  -in.cut=path/cut_contigs.fasta -out=path/output.fasta -CNF.contig.idx=1
```
where `haircut.pipeline.Rscript` is found in the R package installation directory and `-in.raw` must be specified at a minimum. Another popular option is to type in R:

```
#DATA		<- SET THIS DIRECTORY
indir.cut		<- paste(DATA, 'contigs_150408_unaligned_cut', sep='/' )
indir.raw		<- paste(DATA, 'contigs_150408_unaligned_raw', sep='/' )
outdir		<- paste(DATA, 'contigs_150408_model150816a', sep='/' )		
cat(cmd.haircut.pipeline(indir.cut, indir.raw, outdir, batch.n=200, batch.id=1))
```

This produces:


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

# How it works:
The haircutting procedures is divided into several separate steps:

## Aligning cut/raw/ref contigs in multiple ways
The procedure starts by aligning cut and raw IVA contigs against a set of reference HIV sequences. The raw contigs are produced by IVA. The cut contigs are spliced or reverse versions of the raw contigs. MAFFT --add is used obtain different alignments: (1) raw+ref, (2) cut+ref, (3) raw+cut+ref, (4) raw+cut+ref with length kept at the length of the cut+ref alignment. The raw+cut+ref alignment is produced by adding the raw contigs to alignment (2). Alignment (4) is created with the additional options --keeplength --op 0.1.

## Keeping one of the cut/raw/ref contig alignments for further analysis
The next step is to decide which of the alignments (1-4) are used in subsequent steps. This procedure aims to handle very large insertions into the reference alignment. These arise occasionally when raw IVA contigs correspond to the reverse of an actual part of the HIV genome. The current rule is: if alignment (3) has length <12000 bp, then use (3). Otherwise, use alignment (2). Potential caveat: this may ignore valid raw contigs if they accompany a reverse raw contig, and if the valid raw contig is not amongst the cut contigs.  

## Calculate descriptive statistics for each site of cut/raw contigs and the consensus in the same alignment
The next step calculates statistics that are used to decide if a 10 base pair chunk of a contig is to be kept (i. e. no haircut). The script determines the site frequency composition and coverage amongst the reference sequences. Next, the consensus sequence amongst the references is determined. For each cut/raw contig and the consensus, 
*  the probability that the sequence agrees with the references is calculated for each site (*pAGR*),
*  and the presence of a gap at a site is evaluated (*GAP*).

Next, *pAGR* and *GAP* are smoothed by evaluating the mean over a sliding window of length 200. This gives two statistics *spAGR* and *sGAP* for each site of cut/raw contigs and the consensus in the same alignment file.

## Describe the probability of calling a site of the cut/raw contigs and the consensus in the same alignment
Denote this probability by *mu*. The model is 

E logit(mu) = b0 + spAGR * b1 + sGAP * b2.

The probability *mu* differs for each site and depends on the site specific values *spAGR* and *sGAP*. The parameters b0, b1, b2 are fitted through Beta Binomial regression on ~3,000 contigs that Chris Wymant curated manually in April 2015. The parameters are learned independently for each consecutive 10 base pair chunk. Thus, the calling probability is quite local. Potential caveat: if the cut/ref/raw alignment contains large insertions, then the site in the alignment do not correspond to the sites of the model parameters b0, b1, b2. The script removes all insertions that are only present amongst the cut/raw contigs and are longer 100 sites. All such cases corresponded to reverse raw contigs. The value 100 is chosen because the longest insertion in the reference HIV compendium is about 100 bp. Once these insertions are removed, the coordinates of the descriptive statistics correspond to the coordinates of the model parameters.

The above two statistics separate called and not called (i. e. haircut) sites amongst the April 2015 curated contigs quite well:
* each panel corresponds to a 10 bp chunk, with starting position indicated
* blue: called site amongst curated contigs
* red: not called site amongst curated contigs
* y-axis: is *spAGR* (not what it says in the label)
* there is considerable variation in the calling region as we move along the genome. This is why the model parameters b0, b1, b2 are separately calculated for 10 bp chunks.

![alt tag](https://github.com/olli0601/PANGEAhaircut/blob/master/inst/man_stats1.png)
![alt tag](https://github.com/olli0601/PANGEAhaircut/blob/master/inst/man_stats2.png)


## Calling sites i. e. no haircut
A particular site of a cut/raw contig is called (that is not deleted), if the calling probability at that site if
* it is larger than 0.80, or
* if it is not more than 10 standard deviations below the calling probability of the consensus (*muc*)

mu >= min(0.8, muc - 10 * std dev (muc))

This rule accounts for heterogeneity across the HIV genome. In env and especially the V loops, we expect substantial site variation. At these sites, the calling probability of the consensus sequences is much lower than in more conserved gene regions. Therefore, the above approach calls sites much less stringently at known sites of the genome that are associated with substantial variation. 

Here is an example to illustrate:
* each panel corresponds to a cut/raw contig
* blue line: is the threshold, min(0.8, muc - 10 * std dev (muc))
* black line: is the calling probability of the cut/raw contig
* the dips in the blue line correspond to sites with high sequence variability

![alt tag](https://github.com/olli0601/PANGEAhaircut/blob/master/inst/callprob_ex.png)

The sensitivity (SENS), specificity (SPEC), false discovery rate (FDR) and false omission rates (FOR) associated with the April 2015 curated contigs is quite good:
* for positions 1,400 - 1,599.
* we evaluated plots as the one shown here below, and decided to set the cut off in the above decision rule to 10 standard deviations.

![alt tag](https://github.com/olli0601/PANGEAhaircut/blob/master/inst/man_stats3.png)


## Curating called sites into long chunks
Across each cut/raw contig, neighbouring sites may be called or not called, if *mu* is close to the threshold. We define
* called regions of a cut/raw contig as any set of called sites that is not separated by more than 300 bp. The value of 300 bp corresponds to the short read length.
* hair as called sites that are either at the start or the end of a called region
* internal calls as called sites that are not hair
* gaps as uncalled sites

Then,
* gaps between internal calls are called if they are shorter than 100 bp.
* gaps between internal calls are called if they occur after position 9,700
* hair is removed if it less than 150 bp long. This is half the length of a short read.
* any remaining called sites of less than 50bp are removed.

Finally, 
* the script checks if call cut/raw contigs correspond to each other. If this is the case, only the longer one is called.
