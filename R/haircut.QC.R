##--------------------------------------------------------------------------------------------------------
haircut.QC.flatten.automated<- function()
{
	#	flatten automated contigs
	indir		<- paste(DATA,'contigs_150408_model150816a',sep='/')
	outdira		<- paste(DATA,'contigs_150408_model150816aflat',sep='/')
	cmd			<- cmdwrap.flatten.contigs(indir, outdira)	
	outfile		<- paste("hrqc",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
	cmd.hpccaller(paste(DATA,"tmp",sep='/'), outfile, cmd)		
}
##--------------------------------------------------------------------------------------------------------
haircut.QC.align.curated.automated<- function()
{		
	outdira		<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2015/2015_PANGEA_haircut/contigs_150408_model150816aflat'		
	outdircf	<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2015/2015_PANGEA_haircut/contigs_150408_curatedflat'
	infiles		<- data.table(CFILE=list.files(outdircf, pattern='\\.fasta$', recursive=T))
	infiles[, PNG_ID:= gsub('_nLTR_flat\\.fasta','',basename(CFILE))]	
	tmp		 	<- data.table(AFILE=list.files(outdira, pattern='\\.fasta$', recursive=T))
	tmp[, PNG_ID:= gsub('_nohair_flat\\.fasta','',basename(AFILE))]
	infiles		<- merge(infiles, tmp, by='PNG_ID')
	infiles[, OUTFILE1:= paste(PNG_ID,'_tmp.fasta',sep='')]
	infiles[, OUTFILE2:= paste(PNG_ID,'_nohair_flat_wcu.fasta',sep='')]	
	tmp			<- infiles[, {
				#AFILE		<- infiles[1,AFILE]
				#CFILE		<- infiles[1,CFILE]
				#OUTFILE1	<- infiles[1,OUTFILE1]
				#OUTFILE2	<- infiles[1,OUTFILE2]				
				tmp			<- c(	gsub(' ','\\ ',gsub('(','\\(',gsub(')','\\)',paste(outdircf,'/',CFILE,sep=''), fixed=T), fixed=T), fixed=T),
						gsub(' ','\\ ',gsub('(','\\(',gsub(')','\\)',paste(outdira,'/',AFILE,sep=''), fixed=T), fixed=T), fixed=T),
						gsub(' ','\\ ',gsub('(','\\(',gsub(')','\\)',paste(outdira,'/',OUTFILE1,sep=''), fixed=T), fixed=T), fixed=T),
						gsub(' ','\\ ',gsub('(','\\(',gsub(')','\\)',paste(outdira,'/',OUTFILE2,sep=''), fixed=T), fixed=T), fixed=T)	)
				cmd			<- paste('echo ',PNG_ID,'\n',sep='')
				cmd			<- paste(cmd, 'cat ',tmp[1],' ',tmp[2],' > ',tmp[3],'\n',sep='')			
				cmd			<- paste(cmd, PR.MUSCLE,' -quiet -in ',tmp[3],' -out ',tmp[4],'\n','rm ',tmp[3],'\n',sep='')				
				list(CMD=cmd)			
			}, by='PNG_ID']
	cmd			<- paste(tmp$CMD, collapse='\n')	
	outfile		<- paste("hrqc",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
	cmd.hpccaller(paste(DATA,"tmp",sep='/'), outfile, cmd)	
}
##--------------------------------------------------------------------------------------------------------
haircut.QC.divergence.curated.automated<- function()
{
	indir		<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2015/2015_PANGEA_haircut/contigs_150408_model150816a'
	outdira		<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2015/2015_PANGEA_haircut/contigs_150408_model150816aflat'		
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
	write.csv(crd, file=paste(outdira,'/','model150816a_discrepancies.csv',sep=''), row.names=F)
	
	ggplot(crd, aes(x=-DGp, y=X0)) + geom_point(alpha=0.2)
	ggplot(crd, aes(x=DG, y=X0)) + geom_point(alpha=0.2)
}
##--------------------------------------------------------------------------------------------------------
haircut.QC.rmLTR.curated<- function()
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
##--------------------------------------------------------------------------------------------------------
haircut.QC.flatten.curated<- function()
{
	indir		<- paste(DATA,'contigs_150408_curated',sep='/')
	outdir		<- paste(DATA,'contigs_150408_curatedflat',sep='/')
	infiles 	<- data.table(FILE=list.files(indir, pattern='\\.fasta$', recursive=T, include.dirs=T))
	infiles[, PNG_ID:= gsub('_nLTR','',gsub('\\.fasta','',basename(FILE)))]
	infiles[, CFILE:= gsub('\\.fasta','_flat\\.fasta',FILE)]
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
	batch.n		<- 200
	indirc		<- paste(DATA,'contigs_150408_curatedflat',sep='/')
	indira		<- paste(DATA,'contigs_150408_model150816a',sep='/')
	outdir		<- paste(DATA,'contigs_150408_model150816a_cf',sep='/')
	infiles 	<- data.table(FILE=list.files(indirc, pattern='\\.fasta$', recursive=T))
	infiles[, PNG_ID:= gsub('_flat','',gsub('_nLTR','',gsub('\\.fasta','',basename(FILE))))]	
	#	align against automatically created contigs
	
	tmp		 	<- data.table(AFILE=list.files(indira, pattern='\\.fasta$', recursive=T))
	tmp[, PNG_ID:= gsub('_wref_nohair\\.fasta','',basename(AFILE))]		
	cat('\nNO AUTOMATED CONTIGS?\n', paste(setdiff(infiles[, PNG_ID], tmp[,PNG_ID]), collapse='\n'))
	infiles		<- merge(infiles, tmp, by='PNG_ID')
	infiles[, OUTFILE:= paste(PNG_ID,'_nohair_wref_wcu.fasta',sep='')]
	
	tmp			<- infiles[, {
				#AFILE		<- infiles[1,AFILE]
				#CFILE		<- infiles[1,CFILE]
				#OUTFILE		<- infiles[1,OUTFILE]		
				list(CMD=cmd.align.contigs.with.ref(paste(indirc,'/',FILE,sep=''), paste(indira,'/',AFILE,sep=''), paste(outdir,'/',OUTFILE,sep='')))			
			}, by='PNG_ID']
	tmp[, BATCH:= ceiling(seq_len(nrow(tmp))/batch.n)]
	for(batch.id in seq.int(1,tmp[, max(BATCH)]))
	{			
		cmd			<- subset(tmp, BATCH==batch.id)[, paste(CMD, collapse='\n')] 				
		cmd			<- cmd.hpcwrapper(cmd, hpc.nproc= 1, hpc.q='pqeelab', hpc.walltime=4, hpc.mem="5000mb")
		cat(cmd)		
		cmd.hpccaller(paste(DATA,"tmp",sep='/'), paste("hrct",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.'), cmd)	
	}
}
