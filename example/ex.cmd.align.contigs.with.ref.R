#
#	create multiple runs on HPC using the command line version
#
\dontrun{
	
#DATA		<- SET THIS DIRECTORY
indir		<- paste(DATA, 'contigs_150408_merged_unaligned', sep='/' )
outdir		<- paste(DATA, 'contigs_150408_wref', sep='/' )
batch.n		<- 200
tmp			<- data.table(FILE=list.files(indir, pattern='fasta$', recursive=T))
tmp[, BATCH:= ceiling(seq_len(nrow(tmp))/batch.n)]
tmp			<- tmp[, max(BATCH)]
for(batch.id in seq.int(1,tmp))
{	
	cmd			<- cmdwrap.align.contigs.with.ref(indir, outdir, batch.n=batch.n, batch.id=batch.id)
	cmd			<- cmd.hpcwrapper(cmd, hpc.nproc= 1, hpc.q='pqeelab', hpc.walltime=1, hpc.mem="5000mb")
	cat(cmd)		
	cmd.hpccaller(paste(DATA,"tmp",sep='/'), paste("hrct",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.'), cmd)	
}
}
