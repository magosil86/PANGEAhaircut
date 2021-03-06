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
#
#	run from command line 
#	this produces a command line string that can be run in UNIX alikes
#
\dontrun{
	
#DATA		<- SET THIS DIRECTORY
indir	<- paste(DATA, 'contigs_150408_wref', sep='/' )
outdir	<- paste(DATA, 'contigs_150408_wref_cutstat', sep='/' )		
cmd		<- cmd.haircut.cutstat(indir, outdir)
cat(cmd)
}
#
#	create multiple runs on HPC using the command line version
#
\dontrun{
	
#DATA		<- SET THIS DIRECTORY
indir		<- paste(DATA, 'contigs_150408_wref', sep='/' )
outdir		<- paste(DATA, 'contigs_150408_wref_cutstat', sep='/' )		
batch.n		<- 200
tmp			<- data.table(FILE=list.files(indir, pattern='fasta$', recursive=T))
tmp[, BATCH:= ceiling(seq_len(nrow(tmp))/batch.n)]
tmp			<- tmp[, max(BATCH)]
for(batch.id in seq.int(1,tmp))
{			
	cmd			<- cmd.haircut.cutstat(indir, outdir, batch.n=batch.n, batch.id=batch.id)
	cmd			<- cmd.hpcwrapper(cmd, hpc.nproc= 1, hpc.q='pqeelab', hpc.walltime=4, hpc.mem="5000mb")
	cat(cmd)		
	cmd.hpccaller(paste(DATA,"tmp",sep='/'), paste("hrct",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.'), cmd)	
}	
quit("no")	
}
