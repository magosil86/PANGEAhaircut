#
#	within R
#
\dontrun{
#DATA	<- SET THIS DIRECTORY
indir	<- paste(DATA, 'contigs_150408_wref', sep='/' )
outdir	<- paste(DATA, 'contigs_150408_wref_cutstat', sep='/' )		
par		<- c('FRQx.quantile'=NA, 'FRQx.thr'=NA, 'CNS_FRQ.window'=200, 'CNS_AGR.window'=200, 'GPS.window'=200)
haircutwrap.get.cut.statistics(indir, par, outdir=outdir)	
}
