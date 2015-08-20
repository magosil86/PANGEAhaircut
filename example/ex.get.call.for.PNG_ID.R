#
#	within R
#
\dontrun{
	
#DATA					<- SET THIS DIRECTORY		
tmp						<- haircut.get.fitted.model.150816a()
ctrmc					<- tmp$coef		
predict.fun				<- tmp$predict
#	get contigs that were used for training
outfile	<- paste(DATA,'contigs_150408_trainingset_subsets.R',sep='/')
ctrain	<- haircut.get.training.contigs(NULL, outfile, NULL)
set(ctrain, NULL, 'CUT', ctrain[, factor(CUT, levels=c('cut','raw'), labels=c('Y','N'))])
setnames(ctrain, 'CUT', 'BLASTnCUT')		
#	get covariates for all contigs
indir.st<- paste(DATA,'contigs_150408_wref_cutstat',sep='/')
indir.al<- paste(DATA,'contigs_150408_wref',sep='/')
outdir	<- paste(DATA,'contigs_150408_model150816a',sep='/')
par		<- c(	'FRQx.quantile'=NA, 'FRQx.thr'=NA, 'CNS_FRQ.window'=200, 'CNS_AGR.window'=200, 'GPS.window'=200, 
		'PRCALL.thrmax'=0.8, 'PRCALL.thrstd'=10, 'PRCALL.cutprdcthair'=150, 'PRCALL.cutprdctcntg'=50, 'PRCALL.cutrawgrace'=100, 'PRCALL.rmintrnlgpsblw'=100 ,'PRCALL.rmintrnlgpsend'=9700)
haircutwrap.get.call.for.PNG_ID(indir.st,indir.al,outdir,ctrmc,predict.fun,par,ctrain=ctrain)	
}
#
#	run from command line 
#	this produces a command line string that can be run in UNIX alikes
#
\dontrun{

#DATA		<- SET THIS DIRECTORY
indir.st	<- paste(DATA,'contigs_150408_wref_cutstat',sep='/')
indir.al	<- paste(DATA,'contigs_150408_wref',sep='/')
outdir		<- paste(DATA,'contigs_150408_model150816a',sep='/')
cmd			<- cmd.haircut.call(indir.st, indir.al, outdir)
cat(cmd)
}
#
#	create multiple runs on HPC using the command line version
#
\dontrun{
	
#DATA		<- SET THIS DIRECTORY
indir.st	<- paste(DATA,'contigs_150408_wref_cutstat',sep='/')
indir.al	<- paste(DATA,'contigs_150408_wref',sep='/')
outdir		<- paste(DATA,'contigs_150408_model150816a',sep='/')
trainfile	<- paste(DATA,'contigs_150408_trainingset_subsets.R',sep='/')
batch.n		<- 200

tmp			<- data.table(INFILE=list.files(indir.st, pattern='\\.R$', recursive=T))
tmp[, BATCH:= ceiling(seq_len(nrow(tmp))/batch.n)]
tmp			<- tmp[, max(BATCH)]
for(batch.id in seq.int(1,tmp))
{			
	cmd			<- cmd.haircut.call(indir.st, indir.al, outdir, trainfile=trainfile, batch.n=batch.n, batch.id=batch.id, prog=PR.HAIRCUT.CALL )
	cmd			<- cmd.hpcwrapper(cmd, hpc.nproc= 1, hpc.q='pqeelab', hpc.walltime=4, hpc.mem="5000mb")
	cat(cmd)		
	cmd.hpccaller(paste(DATA,"tmp",sep='/'), paste("hrct",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.'), cmd)	
}	
quit("no")
}
