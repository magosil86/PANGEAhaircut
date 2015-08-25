cmd.flatten.contigs<- function(indir=NA, outfile=NA, ctrain=NULL)
{
}
##--------------------------------------------------------------------------------------------------------
haircut.QC.divergence.output.curated<- function()
{
	#	need to flatten curated contigs only once
	if(0)	
	{
		#	rm LTRs from curated
		indir		<- '/Users/Oliver/Dropbox (Infectious Disease)/PANGEA_data/contigs_150408/CuratedAlignmentsToRefs'
		outdirc		<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2015/2015_PANGEA_haircut/contigs_150408_curated'		
		infiles		<- data.table(INFILE=list.files(indir, pattern='fasta$',recursive=T))
		infiles[, PNG_ID:= gsub('_wref.*','',gsub('\\.fasta','',basename(INFILE)))]
		infiles[, OUTFILE:= gsub('_wref','',gsub('\\.fasta','_nLTR\\.fasta',basename(INFILE)))]	
		invisible(infiles[,{
							cr		<- read.dna(paste(indir,'/',INFILE,sep=''),format='fasta')					 
							cr		<- cr[, seq.int(haircut.find.nonLTRstart(cr), ncol(cr))]
							cr		<- cr[, seq.int(1, haircut.find.lastRefSite(cr))]							
							write.dna(cr, file=paste(outdirc,'/',OUTFILE,sep=''), format='fasta', colsep='', nbcol=-1)
							save(cr, file=paste(outdirc,'/',gsub('fasta','R',OUTFILE),sep=''))
							NULL
						},by='PNG_ID'])
		outdircf	<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2015/2015_PANGEA_haircut/contigs_150408_curatedflat'
		cmd			<- cmdwrap.flatten.contigs(outdirc, outdircf)	
		outfile		<- paste("hrqc",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
		cmd.hpccaller(paste(DATA,"tmp",sep='/'), outfile, cmd)
	}
	#	flatten automated contigs
	indir		<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2015/2015_PANGEA_haircut/contigs_150408_model150816a'
	outdira		<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2015/2015_PANGEA_haircut/contigs_150408_model150816aflat'
	cmd			<- cmdwrap.flatten.contigs(indir, outdira)	
	#	align against automatically created contigs
	outdircf	<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2015/2015_PANGEA_haircut/contigs_150408_curatedflat'
	infiles		<- data.table(CFILE=list.files(outdircf, pattern='\\.fasta$', recursive=T))
	infiles[, PNG_ID:= gsub('_nLTR_flat\\.fasta','',basename(CFILE))]	
	tmp		 	<- data.table(AFILE=list.files(outdira, pattern='\\.fasta$', recursive=T))
	tmp[, PNG_ID:= gsub('_nohair_flat\\.fasta','',basename(AFILE))]
	infiles		<- merge(infiles, tmp, by='PNG_ID')
	infiles[, OUTFILE:= paste(PNG_ID,'_nohair_flat_wcu.fasta',sep='')]	
	tmp			<- infiles[, {
				#AFILE		<- infiles[1,AFILE]
				#CFILE		<- infiles[1,CFILE]
				#OUTFILE		<- infiles[1,OUTFILE]		
				list(CMD=cmd.align.contigs.with.ref(paste(outdircf,'/',CFILE,sep=''), paste(outdira,'/',AFILE,sep=''), paste(outdira,'/',OUTFILE,sep='')))			
			}, by='PNG_ID']
	tmp			<- paste(tmp$CMD, collapse='\n')	
	cmd			<- paste(cmd, tmp, sep='\n')
	#cmd			<- cmd.hpcwrapper(cmd, hpc.nproc= 1, hpc.q='pqeelab', hpc.walltime=1, hpc.mem="5000mb")
	outfile		<- paste("hrqc",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
	cmd.hpccaller(paste(DATA,"tmp",sep='/'), outfile, cmd)	
	#	count differences	
	infiles		<- data.table(FILE=list.files(outdira, pattern='_nohair_flat_wcu\\.fasta$', recursive=T))
	infiles[, PNG_ID:= gsub('_nohair_flat_wcu\\.fasta','',basename(FILE))]		
	crd			<- infiles[, {
				
				#FILE		<- subset(infiles, PNG_ID=='12559_1_5')[,FILE]				
				cr			<- read.dna(file=paste(outdira,FILE,sep='/'),format='fasta')
				if(!is.matrix(cr))
				{
					cat('\ncannot evaluate', FILE)
					tmp		<- rep(NA_real_,3)
				}
				if(is.matrix(cr))
				{
					tmp		<- c(dist.dna(cr, model='N'), dist.dna(cr, model='indel'), ncol(cr))
				}
				list(DS= tmp[1], DG=tmp[2], N=tmp[3])
			}, by='PNG_ID']
	#	merge with conf scores
	infiles		<- data.table(SFILE=list.files(indir, pattern='\\.csv$', recursive=T))
	tmp			<- do.call('rbind',lapply(seq_len(nrow(infiles)), function(i){
				as.data.table(read.csv(file=paste(indir,infiles[i,SFILE],sep='/')))
			}))
	crd			<- merge(crd, tmp, by='PNG_ID', all.x=TRUE) 
	crd			<- subset(crd, !is.na(N))
	crd[, DGp:= -DG/N]
	setkey(crd, DGp)
	#	write to csv
	write.csv(crd, file=paste(indir,'/','model150816a_discrepancies.csv',sep=''), row.names=F)
	
	ggplot(crd, aes(x=-DGp, y=X0)) + geom_point(alpha=0.2)
	ggplot(crd, aes(x=DG, y=X0)) + geom_point(alpha=0.2)
}
##--------------------------------------------------------------------------------------------------------
haircut.QC.flatten.curated<- function()
{
	indir		<- paste(DATA,'contigs_150408_curated',sep='/')
	outdir		<- paste(DATA,'contigs_150408_curatedflat',sep='/')
	infiles 	<- data.table(FILE=list.files(indir, pattern='\\.fasta$', recursive=T, include.dirs=T))
	infiles[, PNG_ID:= gsub('_nLTR','',gsub('\\.fasta','',basename(FILE)))]
	infiles[, CFILE:= paste(PNG_ID,'_curated.fasta',sep='')]
	#	extract just the curated contigs from file
	invisible(infiles[,{
						cat('\nprocess',FILE)
						#FILE	<- infiles[1,FILE]
						#CFILE	<- infiles[1,CFILE]
						#PNG_ID	<- infiles[1,PNG_ID] 
						cr				<- read.dna(paste(indir,'/',FILE,sep=''), format='fasta')
						cr				<- cr[grepl(PNG_ID,rownames(cr)),]
						rownames(cr)	<- paste(rownames(cr),'_cur',sep='')
						cr				<- seq.rmgaps(cr, rm.only.col.gaps=0, rm.char='-', verbose=0)		
						write.dna(cr, file=paste(outdir,CFILE,sep='/'), format='fasta', colsep='', nbcol=-1)
						NULL
					}, by='FILE'])
}
##--------------------------------------------------------------------------------------------------------
haircut.QC.align.curated<- function()
{
	indir		<- paste(DATA,'contigs_150408_curated',sep='/')
	outdir		<- paste(DATA,'contigs_150408_curatedflat',sep='/')
	infiles 	<- data.table(FILE=list.files(indir, pattern='\\.fasta$', recursive=T, include.dirs=T))
	infiles[, PNG_ID:= gsub('_nLTR','',gsub('\\.fasta','',basename(FILE)))]
	infiles[, CFILE:= paste(PNG_ID,'_curated.fasta',sep='')]
	#	align against automatically created contigs
	indir		<- paste(DATA,'/contigs_150408_model150816a',sep='/')
	tmp		 	<- data.table(AFILE=list.files(indir, pattern='\\.fasta$', recursive=T, include.dirs=T))
	tmp[, PNG_ID:= gsub('_wref_nohair\\.fasta','',basename(AFILE))]		
	cat('\nNO AUTOMATED CONTIGS?\n', paste(setdiff(infiles[, PNG_ID], tmp[,PNG_ID]), collapse='\n'))
	infiles		<- merge(infiles, tmp, by='PNG_ID')
	infiles[, OUTFILE:= paste(PNG_ID,'_nohair_wref_wcu.fasta',sep='')]
	
	tmp			<- infiles[, {
				#AFILE		<- infiles[1,AFILE]
				#CFILE		<- infiles[1,CFILE]
				#OUTFILE		<- infiles[1,OUTFILE]		
				list(CMD=cmd.align.contigs.with.ref(paste(outdir,'/',CFILE,sep=''), paste(indir,'/',AFILE,sep=''), paste(outdir,'/',OUTFILE,sep='')))			
			}, by='PNG_ID']
	cmd			<- paste(tmp$CMD, collapse='\n')
	cmd			<- cmd.hpcwrapper(cmd, hpc.nproc= 1, hpc.q='pqeelab', hpc.walltime=1, hpc.mem="5000mb")
	cat(cmd)		
	outfile		<- paste("hrct",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
	cmd.hpccaller(paste(DATA,"tmp",sep='/'), outfile, cmd)
}