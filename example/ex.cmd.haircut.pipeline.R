#
#	create multiple runs on HPC using the command line version
#
\dontrun{
	
#DATA		<- SET THIS DIRECTORY
indir.cut	<- paste(DATA, 'contigs_150408_unaligned_cut', sep='/' )
indir.raw	<- paste(DATA, 'contigs_150408_unaligned_raw', sep='/' )
outdir		<- paste(DATA, 'contigs_150408_model150816a', sep='/' )		
batch.n		<- 200
tmp			<- data.table(FILE=list.files(indir, pattern='fasta$', recursive=T))
tmp[, BATCH:= ceiling(seq_len(nrow(tmp))/batch.n)]
tmp			<- tmp[, max(BATCH)]
for(batch.id in seq.int(1,tmp))
{	
	
	cmd			<- cmd.haircut.pipeline(indir.cut, indir.raw, outdir, batch.n=batch.n, batch.id=batch.id)
	cmd			<- cmd.hpcwrapper(cmd, hpc.nproc= 1, hpc.q='pqeelab', hpc.walltime=4, hpc.mem="5000mb")
	cat(cmd)		
	#cmd.hpccaller(paste(DATA,"tmp",sep='/'), paste("hrct",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.'), cmd)	
}	
}
