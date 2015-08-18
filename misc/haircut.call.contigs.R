#! /Library/Frameworks/R.framework/Versions/3.1/Resources/bin/Rscript
##! /apps/R/3.2.0/lib64/R/bin/Rscript
###############################################################################
cat("#######################################################
# haircut.call.contig
# version 150816	
#######################################################")
args <- commandArgs()
if(!any(args=='--args'))
	args<- vector("numeric",0)
if(any(args=='--args'))
	args<- args[-(1:match("--args", args)) ]
#	default args
verbose			<- 1
mfile			<- NA
trainfile		<- NA
indir.st		<- NA
indir.al		<- NA
outdir			<- NA
batch.n			<- NA
batch.id		<- NA
if(length(args))
{
	#	args input
	tmp<- na.omit(sapply(args,function(arg)
					{	switch(substr(arg,2,6),
								mfile= return(substr(arg,8,nchar(arg))),NA)	}))
	if(length(tmp)>0) mfile<- tmp[1]		
	tmp<- na.omit(sapply(args,function(arg)
					{	switch(substr(arg,2,10),
								trainfile= return(substr(arg,12,nchar(arg))),NA)	}))
	if(length(tmp)>0) trainfile<- tmp[1]		
	tmp<- na.omit(sapply(args,function(arg)
					{	switch(substr(arg,2,9),
								indir.st= return(substr(arg,11,nchar(arg))),NA)	}))
	if(length(tmp)>0) indir.st<- tmp[1]		
	tmp<- na.omit(sapply(args,function(arg)
					{	switch(substr(arg,2,9),
								indir.al= return(substr(arg,11,nchar(arg))),NA)	}))
	if(length(tmp)>0) indir.al<- tmp[1]
	tmp<- na.omit(sapply(args,function(arg)
					{	switch(substr(arg,2,7),
								outdir= return(substr(arg,9,nchar(arg))),NA)	}))
	if(length(tmp)>0) outdir<- tmp[1]			
	tmp<- na.omit(sapply(args,function(arg)
					{	switch(substr(arg,2,8),
								batch.n= return(as.numeric(substr(arg,10,nchar(arg)))),NA)	}))
	if(length(tmp)>0) batch.n<- tmp[1]
	tmp<- na.omit(sapply(args,function(arg)
					{	switch(substr(arg,2,9),
								batch.id= return(as.numeric(substr(arg,11,nchar(arg)))),NA)	}))
	if(length(tmp)>0) batch.id<- tmp[1]
}
if( is.na(mfile) || is.na(indir.st) || is.na(indir.al) || is.na(outdir))
{
	stop('Usage:
haircut.call.contig -indir.st=INDIRST -indir.al=INDIRAL -mfile=MFILE -outdir=OUTDIR [-trainfile=TRAINFILE -batch.n=BATCHN -batch.id=BATCHID]
where
INDIRST\t
INDIRAL\t
MFILE\t
OUTDIR\t
TRAINFILE\t
BATCHN\t
BATCHID
')
}
if(verbose)
{
	cat('\ninput args\n',paste(mfile, trainfile, indir.st, indir.al, outdir, batch.n, batch.id, sep='\n'))
}	
###############################################################################
#	helper functions for debugging
my.dumpframes<- function()
{
	geterrmessage()
	dump.frames()
	cat(paste("\nmy.dumpframes dump 'last.dump' to file",paste(HOME,paste("debug_",paste(strsplit(date(),' ')[[1]],collapse='_'),".rda\n",sep=''),sep='')))
	save(last.dump, file=paste(HOME,paste("debug_",paste(strsplit(date(),' ')[[1]],collapse='_'),".rda",sep=''),sep=''))
	q()
}
###############################################################################
#	set debugging
if(DEBUG)	
	options(error= my.dumpframes)	
###############################################################################
#	run script
require(PANGEAhaircut)
tmp						<- haircut.get.fitted.model.150816a(NULL, mfile)
ctrmc					<- tmp$coef		
predict.fun				<- tmp$predict
#	get contigs that were used for training
ctrain	<- NULL
if(!is.na(trainfile))
{
	ctrain	<- haircut.get.training.contigs(NULL, trainfile, NULL)
	set(ctrain, NULL, 'CUT', ctrain[, factor(CUT, levels=c('cut','raw'), labels=c('Y','N'))])
	setnames(ctrain, 'CUT', 'BLASTnCUT')				
}
#	get covariates for all contigs
par		<- c(	'FRQx.quantile'=NA, 'FRQx.thr'=NA, 'CNS_FRQ.window'=200, 'CNS_AGR.window'=200, 'GPS.window'=200, 
				'PRCALL.thrmax'=0.8, 'PRCALL.thrstd'=10, 'PRCALL.cutprdcthair'=150, 'PRCALL.cutprdctcntg'=50, 'PRCALL.cutrawgrace'=100, 'PRCALL.rmintrnlgpsblw'=100 ,'PRCALL.rmintrnlgpsend'=9700)
haircutwrap.get.call.for.PNG_ID(indir.st,indir.al,outdir,ctrmc,predict.fun,par,ctrain=ctrain, batch.n=batch.n, batch.id=batch.id)
