dev.AC.data<- function()
{
	require(data.table)
	indir		<- "/work/or105/PANGEA_AC/150921/"
	infiles		<- data.table(FASTQ=list.files(indir, pattern='fastq\\.gz$|fastq$', recursive=1))
	infiles[, RUN:= sapply(strsplit(FASTQ,'/',fixed=1),'[[',1)]
	infiles[, ID:= gsub('_S[0-9]+_.*','',basename(FASTQ))]
	infiles[, IDS:= gsub('_S|_','',regmatches(basename(FASTQ),regexpr('_S[0-9]+_',basename(FASTQ))))]
	infiles[, IDL:= gsub('_L|_','',regmatches(basename(FASTQ),regexpr('_L[0-9]+_',basename(FASTQ))))]
	infiles[, IDR:= gsub('_R|_','',regmatches(basename(FASTQ),regexpr('_R[0-9]+_',basename(FASTQ))))]
	
	ctfiles		<- data.table(IVA=list.files(indir, recursive=1))
	ctfiles		<- subset(ctfiles, grepl('contigs_hit_ref',IVA)) 
	ctfiles[, RUN:= sapply(strsplit(IVA,'/',fixed=1),'[[',1)]
	ctfiles[, ID:= gsub('.fasta','',gsub('.assembly_contigs_hit_ref','',sapply(strsplit(IVA,'/',fixed=1),'[[',2)))]	
	infiles		<- merge(infiles, ctfiles, by=c('RUN','ID'), all=1)
	
	ctfiles		<- data.table(KRAKEN=list.files(indir, pattern='kraken.report$', recursive=1))	
	ctfiles[, ID:= gsub('.kraken.report','',sapply(strsplit(KRAKEN,'/',fixed=1),'[[',2))]	
	infiles		<- merge(infiles, ctfiles, by=c('ID'), all=1)
	
	write.csv( infiles, row.names=FALSE, file= paste(indir, '/infiles_151020.csv',sep='') )
	save(infiles, file=paste(indir, '/infiles_151020.R',sep=''))
	if(1)
	{
		#	only consider files with RESXxx, PRESxxx and TASPxxx
		pngfiles	<- subset(infiles, !grepl('^UID|^Undetermined',ID))		
		pngfiles[, FASTQ_NEW:= paste(RUN,'_',basename(FASTQ),sep='')]
		tmp			<- pngfiles[, which(!is.na(IVA))]
		pngfiles[, IVA_NEW:= NA_character_]		
		set(pngfiles, tmp, 'IVA_NEW', pngfiles[tmp, paste(RUN,'_',gsub('\\.assembly_contigs_hit_ref','',basename(IVA)),sep='')])
		#	move files
		tmp			<- subset(pngfiles, !is.na(IVA) & IDR==1)
		setkey(tmp, IVA)
		tmp			<- unique(tmp)
		tmp[, file.rename(paste(indir,'/',IVA,sep=''), paste(indir,'/IVA/',IVA_NEW,sep='')), by='IVA']		
		pngfiles[, file.rename(paste(indir,'/',FASTQ,sep=''), paste(indir,'/FASTQ/',FASTQ_NEW,sep='')), by='FASTQ']
		
		write.csv( pngfiles, row.names=FALSE, file= paste(indir, '/pngfiles_151020.csv',sep='') )
		#	check pngfiles
		stopifnot(nrow(subset(pngfiles, is.na(FASTQ)))==0)	
		#	no IVA contigs --> 124
		write.csv( subset(pngfiles, IDR==1 & is.na(IVA)), row.names=FALSE, file= paste(indir, '/infiles_151020_noIVA.csv',sep='') )
		#	check: same S for each ID? -->  Run2_96_600V3_20140814 RES320
		tmp	<- subset(pngfiles, !is.na(IVA))[, list(IDS_N=length(unique(IDS))), by=c('RUN','ID')]
		subset(tmp, IDS_N>1)
		#	check: more than two R files ?
		tmp	<- subset(pngfiles, !is.na(IVA))[, list(IDR_N=length(IDR)), by=c('RUN','ID')]
		subset(tmp, IDR_N>2)
		#	total fastq 853
		nrow(subset(pngfiles, !is.na(FASTQ) & IDR==1))	
		#	total contigs 729
		write.table( subset(pngfiles, !is.na(IVA) & IDR==1, select=IVA), quote=FALSE, col.names=FALSE, row.names=FALSE, file= paste(indir, '/infiles_151020_yesIVA.csv',sep='') )
		nrow(subset(pngfiles, !is.na(IVA) & IDR==1))
		#	number contigs per RUN
		subset(pngfiles, !is.na(IVA) & IDR==1)[, table(RUN)]
		#	number fastq per RUN
		subset(pngfiles, !is.na(FASTQ) & IDR==1)[, table(RUN)]
		#	no kraken report
		write.csv( subset(pngfiles, IDR==1 & is.na(KRAKEN)), row.names=FALSE, file= paste(indir, '/infiles_151020_noKraken.csv',sep='') )
	}
	if(0)
	{
		#	no FASTQ 
		stopifnot(nrow(subset(infiles, is.na(FASTQ)))==0)	
		#	no IVA contigs
		subset(infiles, IDR==1 & is.na(IVA))
		#	check: same S for each ID?
		tmp	<- subset(infiles, !is.na(IVA))[, list(IDS_N=length(unique(IDS))), by=c('RUN','ID')]
		subset(tmp, IDS_N>1)
		#	check: more than two R files ?
		tmp	<- subset(infiles, !is.na(IVA))[, list(IDR_N=length(IDR)), by=c('RUN','ID')]
		subset(tmp, IDR_N>2)
		#	total fastq
		nrow(subset(infiles, !is.na(FASTQ) & IDR==1))	
		#	total contigs
		nrow(subset(infiles, !is.na(IVA) & IDR==1))
		#	number contigs per RUN
		subset(infiles, !is.na(IVA) & IDR==1)[, table(RUN)]
		#	number fastq per RUN
		subset(infiles, !is.na(FASTQ) & IDR==1)[, table(RUN)]
	}
	
}

dev.align.haircutcontigs<- function()
{
	if(0)
	{
		#	extract the haircut contigs and prepare alignment files
		indir1		<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2015/2015_PANGEA_haircut/contigs_150902_model150816b'
		indir2		<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2015/2015_PANGEA_haircut/contigs_150902_model150816b_manual'
		outdir		<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2015/2015_PANGEA_haircut/contigs_150903_model150816b'
		infiles 	<- rbind( 	data.table(FILE=list.files(indir2, pattern='\\.fasta$', recursive=T, full.name=T)),
								data.table(FILE=list.files(indir1, pattern='\\.fasta$', recursive=T, full.name=T))	)		
		infiles[, PNG_ID:= gsub('_nohair','',gsub('_wref','',gsub('\\.fasta','',basename(FILE))))]
		infiles[, MNL:= factor(grepl('manual',FILE), levels=c(TRUE,FALSE), labels=c('Y','N'))]
		infiles		<- dcast.data.table(infiles, PNG_ID~MNL, value.var='FILE')
		tmp			<- infiles[, which(is.na(Y))]
		set(infiles, tmp, 'Y', infiles[tmp, N])
		infiles		<- melt(infiles, id.vars='PNG_ID', measure.vars='Y', value.name='FILE')
		#	extract just the curated contigs from file
		cs			<- sapply(seq_len(nrow(infiles)), function(i)
				{
					cat('\nprocess',infiles[i,FILE])
					#FILE	<- infiles[i,FILE]
					#OFILE	<- infiles[1,OFILE]
					#PNG_ID	<- infiles[1,PNG_ID] 
					cr				<- read.dna(infiles[i,FILE], format='fasta')
					cr				<- cr[grepl(infiles[i,PNG_ID],rownames(cr)),]	
					as.list(cr)
				})
		#	reformat to list of sequences
		cs			<- unlist(cs, recursive=F)
		names(cs)	<- NULL
		csi			<- data.table(TAXON= sapply(cs, rownames), LEN= sapply(cs, ncol))
		csi[, PNG_ID:= sapply(strsplit(gsub('^\\.','',TAXON),'.',fixed=T),'[[',1)]
		cs			<- do.call(c,lapply(cs, as.list))	
		#	extract all contigs with min length into single alignment and save
		tmp			<- merge(infiles, subset(csi, LEN==min(LEN), PNG_ID)[1,], by='PNG_ID')
		cr			<- read.dna(tmp[1,FILE], format='fasta')
		cr			<- as.list(cr[!grepl(tmp[1,PNG_ID],rownames(cr)),])
		cr			<- c(cr, cs[ subset(csi, LEN==min(LEN))[, TAXON] ])	
		write.dna(cr, file= paste(outdir, '/', 'contigs_minlen.fasta', sep=''), format='fasta', colsep='', nbcol=-1)
		write.dna(cr, file= paste(outdir, '/', 'contigs_stub.fasta', sep=''), format='fasta', colsep='', nbcol=-1)
		save(cr, csi, cs, file= paste(outdir, '/', 'contigs.R', sep=''))
		#	extract all others 
		tmp			<- subset(csi, LEN>min(LEN))
		setkey(tmp, LEN)
		tmp[, BATCH_ID:= ceiling(seq_len(nrow(tmp))/200) ]
		invisible(sapply(tmp[, unique(BATCH_ID)], function(b)
						{
							write.dna(cs[ subset(tmp, BATCH_ID==b)[, TAXON] ], file= paste(outdir, '/', 'contigs_batch',b,'.fasta', sep=''), format='fasta', colsep='', nbcol=-1)
						}))				
	}
	if(1)	#check if we missed any contigs
	{
		indir		<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2015/2015_PANGEA_haircut/contigs_150903_model150816b'
		infile		<- 'contigs_0_20.fasta'
		seq			<- read.dna(paste(indir,infile,sep='/'), format='fasta')
		cpr			<- data.table(TAXON=rownames(seq))
		cmissed		<- merge(data.table(TAXON=setdiff( csi[,TAXON], cpr[,TAXON] )), tmp, by='TAXON')
		cmissed		<- cs[ cmissed[,TAXON] ]
		
		cr			<- read.dna( '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2015/2015_PANGEA_haircut/contigs_150902_model150816b/12559_1_1_wref_nohair.fasta', format='fasta')
		cr			<- cr[1:200,] 
		write.dna(c(cmissed, as.list(cr)), file= paste(indir, '/', 'contigs_missed.fasta', sep=''), format='fasta', colsep='', nbcol=-1)		
	}
	if(1)	#clean up alignment
	{
		indir		<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2015/2015_PANGEA_haircut/contigs_150903_model150816b'
		infile		<- 'contigs_0_20b.fasta'
		seq			<- read.dna(paste(indir,infile,sep='/'), format='fasta')
		tmp			<- which( duplicated(rownames(seq)) )
		rownames(seq)[tmp] 
		seq			<- seq[ setdiff( seq_len(nrow(seq)), tmp ), ]
		tmp			<- which(!grepl('^\\.[0-9]+_', rownames(seq)))		
		seq			<- rbind(seq[tmp,], seq[ setdiff( seq_len(nrow(seq)), tmp ), ])
		write.dna(seq, file= paste(indir, '/', 'contigs_cnsalign_PNGIDn3366_CNTGSn6120.fasta', sep=''), format='fasta', colsep='', nbcol=-1)
		#	even more missed?
		tmp			<- setdiff(csi[,TAXON],rownames(seq)[201:nrow(seq)])
		cr			<- read.dna( '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2015/2015_PANGEA_haircut/contigs_150902_model150816b/12559_1_1_wref_nohair.fasta', format='fasta')
		cr			<- cr[1:200,] 		
		write.dna(c(cs[tmp], as.list(cr)), file= paste(indir, '/', 'contigs_missed2.fasta', sep=''), format='fasta', colsep='', nbcol=-1)
		#	--> turns out these are just crap
	
		indir	 	<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2015/2015_PANGEA_haircut/contigs_150903_model150816b'
		infiles		<- list.files(indir, pattern='*stripped*',full.names=TRUE)		
		invisible(lapply(seq_along(infiles), function(i)
						{
							infile			<- infiles[i]
							seq				<- read.dna(infile, format='fasta')
							tmp				<- which( duplicated(rownames(seq)) )
							rownames(seq)[tmp] 
							seq				<- seq[ setdiff( seq_len(nrow(seq)), tmp ), ]
							rownames(seq)	<- gsub(' (stripped)','',rownames(seq), fixed=TRUE)							
							tmp				<- which(!grepl('^\\.[0-9]+_', rownames(seq)))		
							seq				<- rbind(seq[tmp,], seq[ setdiff( seq_len(nrow(seq)), tmp ), ])
							tmp				<- tail(strsplit(basename(infile),'_')[[1]],1)
							write.dna(seq, file= paste(indir, '/', 'contigs_cnsalign_PNGIDn3366_CNTGSn6120_',tmp, sep=''), format='fasta', colsep='', nbcol=-1)							
						}))
	}
	if(1)
	{
		indir		<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2015/2015_PANGEA_haircut/contigs_150903_model150816b'
		load( paste(indir, '/', 'contigs.R', sep='') )
		#	
		#	add all others to alignment in batches of 200
		batch.id	<- ".14683_1_18.2.2_cut"
		write.dna( cs[ subset(tmp, TAXON==".14683_1_18.2.2_cut")[, TAXON] ], file= paste(outdir, '/', 'contigs_BATCH',batch.id,'.fasta', sep=''), format='fasta', colsep='', nbcol=-1)
		cmd			<- cmd.align.contigs.with.ref(paste(outdir, '/', 'contigs_BATCH',batch.id,'.fasta', sep=''), paste(outdir, '/', 'contigs_minlen.fasta', sep=''), paste(outdir, '/', 'contigs_minlen_BATCH',batch.id,'.fasta', sep=''), options='')		
		
		batch.id	<- "14683_1_18"
		write.dna( cs[ subset(tmp, PNG_ID=="14683_1_18")[, TAXON] ], file= paste(outdir, '/', 'contigs_BATCH',batch.id,'.fasta', sep=''), format='fasta', colsep='', nbcol=-1)
		cmd			<- cmd.align.contigs.with.ref(paste(outdir, '/', 'contigs_BATCH',batch.id,'.fasta', sep=''), paste(outdir, '/', 'contigs_minlen.fasta', sep=''), paste(outdir, '/', 'contigs_minlen_BATCH',batch.id,'.fasta', sep=''), options='')	
		
		batch.id	<- 20
		write.dna( cs[ subset(tmp, BATCH_ID==batch.id)[, TAXON] ], file= paste(outdir, '/', 'contigs_BATCH',batch.id,'.fasta', sep=''), format='fasta', colsep='', nbcol=-1)	
		cmd			<- cmd.align.contigs.with.ref(paste(outdir, '/', 'contigs_BATCH',batch.id,'.fasta', sep=''), paste(outdir, '/', 'contigs_minlen.fasta', sep=''), paste(outdir, '/', 'contigs_minlen_BATCH.fasta', sep=''), options='')
		
	}
		
}


dev.haircut<- function()	
{
	if(0)	#check which files not aligned 
	{
		dscp	<- as.data.table(read.csv("/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2015/2015_PANGEA_haircut/pngfiles_151020.csv", stringsAsFactors=FALSE))
		dscp[, ID:= paste(RUN,'_',ID,sep='')]
		
		dr		<- data.table(BLAST=list.files("~/Dropbox (Infectious Disease)/OR_Work/2015/2015_PANGEA_haircut/contigs_151026_unaligned_raw"))
		dr[, ID:= gsub('_hiv','',sapply(strsplit(BLAST,'.',fixed=1),'[[',1))]
		dscp	<- merge(dscp, dr, by='ID', all.x=1)
		subset(dscp, is.na(BLAST) & IDR==1 & !is.na(IVA))	
		
		subset(dscp, IDR==1)
	}
	if(0)	#check alignment: code
	{
		indir		<- paste(DATA, 'contigs_150408_wref', sep='/' )
		indir		<- paste(DATA, 'contigs_150902_wref', sep='/' )				
		infiles		<- data.table(FILE=list.files(indir, pattern='fasta$', recursive=T))
		infiles[, PNG_ID:= gsub('\\.fasta','',gsub('_frclen|_refc|_refr|_wRefs','',FILE))]
		infiles[, AL_TYPE:= gsub('_*','',gsub('\\.fasta','',gsub('[0-9]*','',FILE)))]				
		tmp			<- infiles[, {
					#FILE		<- infiles[1,FILE]
					tmp			<- paste(indir, FILE, sep='/')
					tmp			<- gsub(' ','\\ ',gsub('(','\\(',gsub(')','\\)',tmp,fixed=T),fixed=T),fixed=T)					
					cmd			<- paste("awk '/^>/ {if (seqlen) print seqlen;print;seqlen=0;next} {seqlen+=length($0)}END{print seqlen}' ", tmp, sep='')
					tmp			<- system(cmd, intern=TRUE)					
					list(LEN=as.numeric(tmp[2]))
				}, by='FILE']
		infiles		<- merge(infiles, tmp, by='FILE')		
		tmp			<- dcast.data.table(infiles, PNG_ID~AL_TYPE, value.var='LEN')
		#	suspect that raw is in reverse if alignment length is larger than 12k		
		stopifnot( nrow(subset(tmp, wRefs>12e3 & refc>12e3))==0 )
		#	mv refc to wRefs in case wRef length > 12e3
		alfixup		<- merge(subset(tmp, wRefs>12e3 & refc<=12e3, PNG_ID), infiles, by='PNG_ID')
		alfixup		<- dcast.data.table(alfixup, PNG_ID~AL_TYPE, value.var='FILE')
		alfixup[, {
						file.rename( paste(indir,wRefs,sep='/'), paste(indir,gsub('wRefs\\.fasta','fndrvs\\.fasta',wRefs),sep='/'))
						file.copy( paste(indir,refc,sep='/'), paste(indir,gsub('refc\\.fasta','wRefs\\.fasta',wRefs),sep='/'), overwrite=TRUE)
						NULL
				}, by='PNG_ID']		
		
		cat(cmd.haircut.check.alignment(indir, indir))
	}
	if(0)	#check alignment: cmd
	{
		indir		<- paste(DATA, 'contigs_151026_wref', sep='/' )		
		argv		<<- cmd.haircut.check.alignment(indir, indir)
		argv		<<- paste('-',unlist(strsplit(argv,' -|\n')),sep='')
	}
	if(0)
	{
		indir		<- paste(DATA, 'contigs_151026_wref', sep='/' )
		outdir		<- paste(DATA, 'contigs_151026_wref_cutstat', sep='/' )				
		argv		<<- cmd.haircut.cutstat(indir, outdir, batch.n=200, batch.id=1)
		argv		<<- paste('-',unlist(strsplit(argv,' -|\n')),sep='')
	}
	if(1)
	{
		indir.st	<- paste(DATA,'contigs_150408_wref_cutstat',sep='/')
		indir.al	<- paste(DATA,'contigs_150408_wref',sep='/')
		outdir		<- paste(DATA,'contigs_150408_model150816a',sep='/')
		indir.st	<- paste(DATA,'contigs_151026_wref_cutstat',sep='/')
		indir.al	<- paste(DATA,'contigs_151026_wref',sep='/')
		outdir		<- paste(DATA,'contigs_151026_model150816b',sep='/')
		
		argv		<<- cmd.haircut.call(indir.st, indir.al, outdir)
		argv		<<- paste('-',unlist(strsplit(argv,' -|\n')),sep='')
		
	}
	if(0)
	{
		#fixup mafft
		indir		<- paste(DATA, 'contigs_150408_merged_unaligned', sep='/' )
		outdir		<- paste(DATA, 'contigs_150408_wref2', sep='/' )
		batch.n		<- 200
		
		tmp			<- data.table(INFILE=list.files(indir, pattern='fasta$', recursive=T))
		tmp[, PNG_ID:= gsub('\\.fasta','',gsub('_cut|_raw','',INFILE))]
		tmp[, OUTFILE:= gsub('\\.fasta','_wRefs\\.fasta', gsub('_hiv|_HIV','',basename(INFILE)))]
		reffile		<- system.file(package="PANGEAhaircut", "HIV1_COM_2012_genome_DNA_WithExtraA1UG.fasta")
		
		i			<- tmp[, which(PNG_ID=='12559_1_13')]
		INFILE		<- tmp[i, INFILE]
		OUTFILE		<- tmp[i, OUTFILE]
		reffile		<- system.file(package="PANGEAhaircut", "HIV1_COM_2012_genome_DNA_WithExtraA1UG.fasta")
		cat(cmd.align.contigs.with.ref(paste(indir,'/',INFILE,sep=''), reffile, paste(outdir,'/',OUTFILE,sep='')))
		#
		indir.cut	<- '/Users/Oliver/Dropbox (Infectious Disease)/PANGEA_data/contigs_150408/HIVcontigs_cut'
		indir.raw	<- '/Users/Oliver/Dropbox (Infectious Disease)/PANGEA_data/contigs_150408/HIVcontigs_Raw'
		outdir		<- paste(DATA, 'contigs_150408_wref2', sep='/' )
		infiles		<- data.table(INFILECUT=list.files(indir.cut, pattern='fasta$', recursive=T))
		infiles[, PNG_ID:= gsub('_hiv','',gsub('\\.fasta','',gsub('_cut|_raw','',INFILECUT)))]
		tmp			<- data.table(INFILECUT=list.files(indir.raw, pattern='fasta$', recursive=T))
		tmp			<- data.table(INFILERAW=list.files(indir.raw, pattern='fasta$', recursive=T))
		tmp[, PNG_ID:= gsub('_hiv','',gsub('\\.fasta','',gsub('_cut|_raw','',INFILERAW)))]
		infiles		<- merge(infiles, tmp, all=TRUE, by='PNG_ID')
		infiles[, OUTFILE1:= paste(PNG_ID,'_c.fasta',sep='')]
		infiles[, OUTFILE2:= paste(PNG_ID,'_refc.fasta',sep='')]
		infiles[, OUTFILE3:= paste(PNG_ID,'_wRefs.fasta',sep='')]
		infiles[, OUTFILE4:= paste(PNG_ID,'_refr.fasta',sep='')]
		infiles[, OUTFILE5:= paste(PNG_ID,'_frclen.fasta',sep='')]
		
		i			<- infiles[, which(PNG_ID=='13554_1_14')]
		i			<- infiles[, which(PNG_ID=='13554_1_17')]
		#i			<- infiles[, which(PNG_ID=='12559_1_13')]
		
		INFILECUT	<- infiles[i, INFILECUT]
		INFILERAW	<- infiles[i, INFILERAW]
		OUTFILE1	<- infiles[i, OUTFILE1]
		OUTFILE2	<- infiles[i, OUTFILE2]
		OUTFILE3	<- infiles[i, OUTFILE3]
		OUTFILE4	<- infiles[i, OUTFILE4]
		OUTFILE5	<- infiles[i, OUTFILE5]
		reffile		<- system.file(package="PANGEAhaircut", "HIV1_COM_2012_genome_DNA_WithExtraA1UG.fasta")
		#reffilen7	<- '/Users/Oliver/Dropbox (Infectious\ Disease)/OR_Work/2015/2015_PANGEA_haircut/contigs_150408_wref2/HIV1_COM_n7.fasta'
		cmd			<- cmd.add.tag.to.fasta.names( paste(indir.cut,'/',INFILECUT,sep=''), paste(outdir,'/',OUTFILE1,sep=''), tag='_cut')
		cmd			<- paste(cmd, cmd.align.contigs.with.ref(paste(outdir,'/',OUTFILE1,sep=''), reffile, paste(outdir,'/',OUTFILE2,sep='')), sep='\n')
		cmd			<- paste(cmd, cmd.align.contigs.with.ref(paste(indir.raw,'/',INFILERAW,sep=''), paste(outdir,'/',OUTFILE2,sep=''), paste(outdir,'/',OUTFILE3,sep='')), sep='\n')
		cmd			<- paste(cmd, cmd.align.contigs.with.ref(paste(indir.raw,'/',INFILERAW,sep=''), paste(outdir,'/',OUTFILE2,sep=''), paste(outdir,'/',OUTFILE5,sep=''), options='--keeplength --op 0.1'), sep='\n')
		cmd			<- paste(cmd, cmd.align.contigs.with.ref(paste(indir.raw,'/',INFILERAW,sep=''), reffile, paste(outdir,'/',OUTFILE4,sep='')), sep='\n')		
		tmp			<- paste(outdir,'/',OUTFILE1,sep='')
		tmp			<- gsub(' ','\\ ',gsub('(','\\(',gsub(')','\\)',tmp,fixed=T),fixed=T),fixed=T)
		cmd			<- paste(cmd, '\n','rm ',tmp,sep='')		
		cat(cmd)		
	}
	if(0)
	{
		#	read just one file
		#	determine statistics for each contig after LTR		
		indir	<- '/Users/Oliver/Dropbox\ (Infectious Disease)/PANGEA_data/InterestingContigAlignments'
		outdir	<- '/Users/Oliver/Dropbox\ (Infectious Disease)/OR_Work/2015/2015_PANGEA_haircut'
		infiles	<- data.table(FILE=list.files(indir, pattern='fasta$', recursive=T))
		infiles[, PNG_ID:= gsub('_wRefs\\.fasta','',gsub('_cut|_raw','',FILE))]
		infiles[, BLASTnCUT:= regmatches(FILE,regexpr('cut|raw',FILE))]
		set(infiles, NULL, 'BLASTnCUT', infiles[, factor(BLASTnCUT, levels=c('cut','raw'), labels=c('Y','N'))])
		par		<- c('FRQx.quantile'=0.05, 'CNS_AGR.window'=200, 'GPS.window'=200)
		
		
		file	<- paste(indir, infiles[1, FILE], sep='/')
		#	read Contigs+Rrefs: cr
		cr		<- read.dna(file, format='fasta')			
		#	determine start of non-LTR position and cut 
		tmp		<- haircut.find.nonLTRstart(cr)
		cr		<- cr[, seq.int(tmp, ncol(cr))]
		#	determine reference sequences. 
		#	non-refs have the first part of the file name in their contig name and are at the top of the alignment
		tmp		<- strsplit(basename(file), '_')[[1]][1]
		tx		<- data.table(TAXON= rownames(cr), CONTIG=as.integer(grepl(tmp, rownames(cr))) )
		stopifnot( all( tx[, which(CONTIG==1)] == seq.int(1, tx[, length(which(CONTIG==1))]) ) )
		cat(paste('\nFound contigs, n=', tx[, length(which(CONTIG==1))]))
		#	determine base frequencies at each site amongst references.
		tmp		<- cr[subset(tx, CONTIG==0)[, TAXON],]
		cnsr	<- haircut.getconsensus(tmp, par, bases=c('a','c','g','t','-') )	#	CoNSensus of References: cnsr
		#	for each contig, determine %agreement with consensus on rolling window
		cnsc	<- rbind(cnsr, cr[subset(tx, CONTIG==1)[, TAXON],])
		#	determine first and last non-gap sites
		tx		<- data.table(	TAXON= rownames(cnsc), 
				FIRST= apply( as.character(cnsc), 1, function(x) which(x!='-')[1] ),
				LAST= ncol(cnsc)-apply( as.character(cnsc), 1, function(x) which(rev(x)!='-')[1] )		)
		#	get cut statistics
		cnsc.df	<- haircut.get.cut.statistics(cnsc, tx, par, outdir=NA, file=NA, mode='rolling')
		cnsc.df	<- haircut.get.cut.statistics(cnsc, tx, par, outdir=outdir, file=file, mode='rolling')		
	}
	if(0)	#extract just the curated contigs from file
	{
		indir		<- '/Users/Oliver/Dropbox (Infectious Disease)/PANGEA_data/contigs_150408/CuratedAlignmentsToRefs'
		outdir		<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2015/2015_PANGEA_haircut/contigs_150408_model150816acur'
		infiles 	<- data.table(FILE=list.files(indir, pattern='\\.fasta$', recursive=T, include.dirs=T))
		infiles[, PNG_ID:= gsub('\\.fasta','',basename(FILE))]
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
		#	align against automatically created contigs
		indir		<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2015/2015_PANGEA_haircut/contigs_150408_model150816a'
		tmp		 	<- data.table(AFILE=list.files(indir, pattern='\\.fasta$', recursive=T, include.dirs=T))
		tmp[, PNG_ID:= gsub('_wref_nohair\\.fasta','',basename(AFILE))]		
		cat('\nNo Automated Contigs?\n', paste(setdiff(infiles[, PNG_ID], tmp[,PNG_ID]), collapse='\n'))
		cat('\nNo Curated Contigs?\n', paste(setdiff(tmp[, PNG_ID], tmp[,PNG_ID]), collapse='\n'))
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
		#PNG_ID	<- infiles[1,PNG_ID] 
		
	}
	if(0)
	{
		#	process several files
		indir	<- '/Users/Oliver/Dropbox\ (Infectious Disease)/PANGEA_data/InterestingContigAlignments'
		outdir	<- '/Users/Oliver/Dropbox\ (Infectious Disease)/OR_Work/2015/2015_PANGEA_haircut/interesting_150408'		
		par		<- c('FRQx.quantile'=NA, 'FRQx.thr'=NA, 'CNS_FRQ.window'=200, 'CNS_AGR.window'=200, 'GPS.window'=200)
		haircutwrap.get.cut.statistics(indir, par, outdir=outdir)
	}
	if(0)
	{
		#	run mafft --add to get contigs+ref
		indir		<- '/Users/Oliver/Dropbox\ (Infectious Disease)/OR_Work/2015/2015_PANGEA_haircut/contigs_150408_merged_unaligned'
		outdir		<- '/Users/Oliver/Dropbox\ (Infectious Disease)/OR_Work/2015/2015_PANGEA_haircut/contigs_150408_wref'
		reffile		<- '/Users/Oliver/Dropbox\ (Infectious Disease)/PANGEA_data/HIV1_COM_2012_genome_DNA_WithExtraA1UG.fasta'
		cmd			<- cmdwrap.align.contigs.with.ref(indir, reffile, outdir, batch.n=200, batch.id=2)
		cmd			<- cmdwrap.align.contigs.with.ref(indir, reffile, outdir)
		cmd			<- cmd.hpcwrapper(cmd, hpc.nproc= 1, hpc.q='pqeelab', hpc.walltime=4, hpc.mem="5000mb")
		cat(cmd)		
		outdir		<- paste(DATA,"tmp",sep='/')
		outfile		<- paste("hrct",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
		cmd.hpccaller(outdir, outfile, cmd)	
		#haircutwrap.align.contigs.with.ref(indir, reffile, outdir)
	}	
	if(0)
	{
		indir	<- '/Users/Oliver/Dropbox\ (Infectious Disease)/OR_Work/2015/2015_PANGEA_haircut/contigs_150408_wref'
		outfile	<- '/Users/Oliver/Dropbox\ (Infectious Disease)/OR_Work/2015/2015_PANGEA_haircut/contigs_150408_cutsubraw.R'
		txe		<- haircutwrap.get.subset.among.raw.and.cut.contigs(indir, outfile)
	}
	if(0)	# get training data
	{	
		#	get contigs that are to be used for training: this determines the 1's and excludes some 0's
		#	read curated contigs 
		indir	<- '/Users/Oliver/Dropbox\ (Infectious Disease)/PANGEA_data/contigs_150408/CuratedAlignmentsToRefs'
		outfile	<- '/Users/Oliver/Dropbox\ (Infectious Disease)/OR_Work/2015/2015_PANGEA_haircut/contigs_150408_curated.R'
		ctrain	<- haircut.get.curated.contigs(indir, outfile)		
		setnames(ctrain, c('FILE','LEN','FIRST','LAST'),c('ANS_FILE','ANS_LEN','ANS_FIRST','ANS_LAST'))
		#	remove cut & raw contigs that are identical as double 1's and as 0's
		indir	<- '/Users/Oliver/Dropbox\ (Infectious Disease)/OR_Work/2015/2015_PANGEA_haircut/contigs_150408_wref'
		outfile	<- '/Users/Oliver/Dropbox\ (Infectious Disease)/OR_Work/2015/2015_PANGEA_haircut/contigs_150408_trainingset_subsets.R'
		ctrain	<- haircut.get.training.contigs(indir, outfile, ctrain)
		set(ctrain, NULL, 'CUT', ctrain[, factor(CUT, levels=c('cut','raw'), labels=c('Y','N'))])
		setnames(ctrain, 'CUT', 'BLASTnCUT')		
		#	now create training data set: 
		#	expand the curated contigs in the training data to 1's per site, and
		#	expand all non-existing curated contigs to 0's, provided they are not excluded because identical to a 1
		indir	<- '/Users/Oliver/Dropbox\ (Infectious Disease)/OR_Work/2015/2015_PANGEA_haircut/contigs_150408_wref_cutstat'
		outdir	<- '/Users/Oliver/Dropbox\ (Infectious Disease)/OR_Work/2015/2015_PANGEA_haircut/contigs_150408_train'
		outfile	<- 'contigs_150408_train'
		par		<- c(	'FRQx.quantile'=NA, 'FRQx.thr'=NA, 'CNS_FRQ.window'=200, 'CNS_AGR.window'=200, 'GPS.window'=200, 
						'PRCALL.thrmax'=0.8, 'PRCALL.thrstd'=10, 'PRCALL.cutprdcthair'=150, 'PRCALL.cutprdctcntg'=50, 'PRCALL.cutrawgrace'=100, 'PRCALL.rmintrnlgpsblw'=100 ,'PRCALL.rmintrnlgpsend'=9700)		
		haircut.get.training.data(indir, ctrain, par, outdir, outfile)
	}
	if(0)	#	fit model
	{
		indir	<- '/Users/Oliver/Dropbox\ (Infectious Disease)/OR_Work/2015/2015_PANGEA_haircut/contigs_150408_train'
		outfile	<- '/Users/Oliver/Dropbox\ (Infectious Disease)/OR_Work/2015/2015_PANGEA_haircut/model_150816a.R'
		tmp		<- haircut.get.fitted.model.150816a(indir, outfile)
	}
	if(0)	#	call contigs on training data and plot
	{		
		#mfile	<- '/Users/Oliver/Dropbox\ (Infectious Disease)/OR_Work/2015/2015_PANGEA_haircut/model_150816a.R'
		#	get model coefficients across the chunks
		tmp						<- haircut.get.fitted.model.150816a()
		ctrmc					<- tmp$coef		
		predict.fun				<- tmp$predict
		#	get contigs that were used for training		
		outfile	<- '/Users/Oliver/Dropbox\ (Infectious Disease)/OR_Work/2015/2015_PANGEA_haircut/contigs_150408_trainingset_subsets.R'
		ctrain	<- haircut.get.training.contigs(NULL, outfile, NULL)
		set(ctrain, NULL, 'CUT', ctrain[, factor(CUT, levels=c('cut','raw'), labels=c('Y','N'))])
		setnames(ctrain, 'CUT', 'BLASTnCUT')		
		#	get covariates for all contigs
		indir.st<- '/Users/Oliver/Dropbox\ (Infectious Disease)/OR_Work/2015/2015_PANGEA_haircut/contigs_150408_wref_cutstat'
		indir.al<- '/Users/Oliver/Dropbox\ (Infectious Disease)/OR_Work/2015/2015_PANGEA_haircut/contigs_150408_wref'
		outdir	<- '/Users/Oliver/Dropbox\ (Infectious Disease)/OR_Work/2015/2015_PANGEA_haircut/contigs_150408_model150816a'
		par		<- c(	'FRQx.quantile'=NA, 'FRQx.thr'=NA, 'CNS_FRQ.window'=200, 'CNS_AGR.window'=200, 'GPS.window'=200, 
						'PRCALL.thrmax'=0.8, 'PRCALL.thrstd'=10, 'PRCALL.cutprdcthair'=100, 'PRCALL.cutrawgrace'=100, 'PRCALL.rmintrnlgpsblw'=100 ,'PRCALL.rmintrnlgpsend'=9700)
		haircutwrap.get.call.for.PNG_ID(indir.st,indir.al,outdir,ctrmc,ctrev,predict.fun,par,ctrain=ctrain)
	}
	if(0)	# evaluate training data
	{	
		#	
		indir	<- '/Users/Oliver/Dropbox\ (Infectious Disease)/OR_Work/2015/2015_PANGEA_haircut/contigs_150408_train'
		outfile	<- 'model_150811a.R'
		#
		#	find 20 worst false pos
		#
		ctrfp	<- do.call('rbind', lapply(seq(1,10000,200), function(site)
				{
					ctr		<- haircut.calculate.training.data(indir, site)	
					tmp		<- seq(ctr[, min(SITE)-1],ctr[, max(SITE)],10)
					ctr[, CHUNK:=cut(SITE, breaks=tmp, labels=tmp[-length(tmp)])]	
					do.call('rbind',lapply( ctr[, unique(CHUNK)], function(chunk){
								ctrch	<- subset(ctr, CHUNK==chunk)
								ctrchfp	<- subset( ctrch, ANS_CALL==0 & AGRpc>subset( ctrch, ANS_CALL==1 )[, min(AGRpc)])
								ctrchfp	<- ctrchfp[, list(ANS_CALL=mean(ANS_CALL), AGRpc=mean(AGRpc), GPS=mean(GPS), CNS_FRQr=mean(CNS_FRQr)), by=c('PNG_ID','TAXON')]
								setkey(ctrchfp, AGRpc)
								ctrchfp[ seq.int(max(1,nrow(ctrchfp)-20), nrow(ctrchfp)), ]	
							}))										
				}))		
		save(ctrfp, file= paste(indir,'/fp.R',sep=''))
		setkey(ctrfp, PNG_ID, TAXON)
		ctrfp[, length(unique(TAXON))]
		#
		#	show all three types of information across sites
		#
		invisible(do.call('rbind', lapply(seq(1,10000,200), function(site)
					{
						cat('\nProcess',site)
						ctr		<- haircut.calculate.training.data(indir, site)	
						tmp		<- seq(ctr[, min(SITE)-1],ctr[, max(SITE)],10)
						ctr[, CHUNK:=cut(SITE, breaks=tmp, labels=tmp[-length(tmp)])]
						ggplot(ctr, aes(x=GPS, y=FRQ, colour=factor(ANS_CALL, levels=c(0,1), labels=c('N','Y')))) + 
								geom_point(alpha=0.6, size=1) + 
								facet_wrap(~CHUNK, ncol=5) +
								theme_bw() + labs(x='%gappiness', y='%agreement with consensus\nline: %agreement of consensus with references', colour='In curated contigs') + theme(legend.position='bottom')
						ggsave(file=paste(outdir,'/',outfile,'_ANSCALLBYFRQGPS_SITE',site,'.pdf',sep=''), w=9, h=9)						
					})))
		#
		#	calculate model coefficients across the chunks
		#
		tmp						<- haircut.get.fitted.model.150811a(indir, outfile)
		ctrmc					<- tmp$coef
		ctrev					<- tmp$ev
		model.150811a.predict	<- tmp$predict
	
		
		
		
		#	select threshold
		ctrt		<- subset(ctrev, THR==0.8)
		ggplot(melt(ctrt, id.vars=c('CHUNK'), measure.vars=c('THR','SENS','SPEC','FDR','FOR')), aes(x=CHUNK, y=value, colour=variable)) + geom_line() +
					scale_x_continuous(breaks=seq(0,20e3,1e3), expand=c(0,0)) +
					scale_y_continuous(breaks=seq(0,1,0.2), minor_breaks=seq(0,1,0.05)) +
					theme_bw() + labs(x='base relative to HXB2') + theme(panel.grid.major=element_line(colour="grey50", size=0.2), panel.grid.minor=element_line(colour="grey70", size=0.1))
		ggsave(file=paste(outdir,'/',outfile,'_model.150811a_THR_const80.pdf',sep=''), w=10, h=4)
		#
		setkey(ctrev, CHUNK, FDR)
		ctrt		<- ctrev[, {
					z		<- which(FDR<= max(0.01, min(FDR, na.rm=T)))[1]
					if(!is.finite(z) | THR[z]<0.8)
						z	<- which(THR==0.8)
					list(THR=THR[z], SENS=SENS[z], SPEC=SPEC[z], FDR=FDR[z], FOR=FOR[z])
				}, by='CHUNK']
		ggplot(melt(ctrt, id.vars=c('CHUNK'), measure.vars=c('THR','SENS','SPEC','FDR','FOR')), aes(x=CHUNK, y=value, colour=variable)) + geom_line() +
				scale_x_continuous(breaks=seq(0,20e3,1e3), expand=c(0,0)) +
				scale_y_continuous(breaks=seq(0,1,0.2), minor_breaks=seq(0,1,0.05)) +
				theme_bw() + labs(x='base relative to HXB2') + theme(panel.grid.major=element_line(colour="grey50", size=0.2), panel.grid.minor=element_line(colour="grey70", size=0.1))
		ggsave(file=paste(outdir,'/',outfile,'_model.150811a_THR_BYFDR01.pdf',sep=''), w=10, h=4)
		#	just const threshold??



		tmp		<- ctrev[, list(FDR=min(FDR)), by='CHUNK']
		ggplot(tmp, aes(x=FDR)) + geom_histogram()
		subset(tmp, FDR>0.05)
		
		ctr		<- haircut.calculate.training.data(indir, 1)	
		tmp		<- seq(ctr[, min(SITE)-1],ctr[, max(SITE)],10)
		ctr[, CHUNK:=cut(SITE, breaks=tmp, labels=tmp[-length(tmp)])]
		
		ctrp	<- merge(ctr, ctrmc, by='CHUNK')
		ctrp[, PR_CALL:=model.150811a.predict(AGRpc, GPS, BETA0, BETA1, BETA2)]
		setkey(ctrp, CHUNK)
		ctrev	<- as.data.table(expand.grid(THR= seq(0.05,0.95,0.01), CHUNK= as.character(ctrp[, unique(CHUNK)]), stringsAsFactors=FALSE))
		ctrev	<- ctrev[, {
					tmp	<- ctrp[CHUNK,][,table(ANS_CALL, factor(PR_CALL>=THR, levels=c('TRUE','FALSE'), labels=c('TRUE','FALSE')))]
					list(TP= tmp['1','TRUE'], FP= tmp['0','TRUE'], FN= tmp['1','FALSE'], TN= tmp['0','FALSE'] )					
				}, by=c('CHUNK','THR')]
		ctrev	<- merge(ctrev, ctrev[, list(SENS= TP/(TP+FN), SPEC=TN/(TN+FP), FDR=FP/(FP+TP), FOR=FN/(FN+TN)), by=c('CHUNK','THR')], by=c('CHUNK','THR'))
		#	plot Sensitivity & Specificity
		ggplot(melt(ctrev, id.vars=c('CHUNK','THR'), measure.vars=c('SENS','SPEC','FDR','FOR')), aes(x=THR, y=100*value, colour=CHUNK)) + 
				geom_line() + labs(x='threshold on predicted call probability',y='%', colour='site\n(base relative to HXB2)') +
				theme(legend.position='bottom') + facet_wrap(~variable, ncol=2, scales='free') +
				guides(col = guide_legend(ncol=5, byrow=TRUE))
		ggsave(file=paste(outdir,'/',outfile,'_model.150811a_SensSpec_SITE',site,'.pdf',sep=''), w=9, h=9)
		#
		ctrev
		
		ctrp[CHUNK, table(ANS_CALL, PR_CALL>=THR)]
		

		#	calculate model coefficients across the chunks		
		ctrmc	<- do.call('rbind',lapply(ctr[, unique(CHUNK)], function(chunk)
				{
					ctrch	<- subset(ctr, CHUNK==chunk)
					ctrchm	<- gamlss(ANS_CALL~AGRpc+GPS, data=ctrch, family=BI())				#as good as 'AGRpc+GPS+CNS_FRQr' in terms of FN, FP
					ctrchmc	<- data.table(CHUNK=chunk, BETA0=coef(ctrchm)[1], BETA1=coef(ctrchm)[2], BETA2=coef(ctrchm)[3])					
				}))
		#	evaluate FN FP
		#	return the whole lot
		#	determine threshold based on FN FP, and make call
	
		
		
		#	add rolling mean by AGRpc
		site	<- 1610; chunk<- '1610'
		site	<- 9350; chunk<- '9350'
		site	<- 3650; chunk<- '3650'
		ctr		<- haircut.load.training.data(indir, site)		
		tmp		<- seq(ctr[, min(SITE)-1],ctr[, max(SITE)],10)
		ctr[, CHUNK:=cut(SITE, breaks=tmp, labels=tmp[-length(tmp)])]
		ctrch	<- subset(ctr, CHUNK==chunk)
		setkey(ctrch, FRQ)
		ctrch[, ANS_CALLrm:=ctrch[, rollapply(ANS_CALL, width=100, FUN=mean, align="center", partial=TRUE)]]	
		
		#setkey(ctrch, GPSc, AGRpc)
		#ctrch	<- merge(ctrch, ctrch[, list(TAXON=TAXON, SITE=SITE, BLASTnCUT=BLASTnCUT, ANS_CALLrm2=rollapply(ANS_CALL, width=100, FUN=mean, align="center", partial=TRUE)), by='GPSc'] , by=c('TAXON','SITE','BLASTnCUT','GPSc'))		
		#tmp		<- seq(0,1.01,0.01)
		#ctr[, GPSC:= cut(GPS, breaks=tmp, labels=tmp[-length(tmp)])]
		#ctr[, AGRpcC:= cut(AGRpc, breaks=tmp, labels=tmp[-length(tmp)])]		
		
		ctrchm	<- gamlss(ANS_CALL~FRQ, data=ctrch, family=BI())
		ctrch[, PR_CALL:= predict(ctrchm, type='response', what='mu')]
		ggplot(melt(ctrch, id.vars=c('FRQ','GPS'), measure.vars=c('ANS_CALLrm','PR_CALL')), aes(x=FRQ, y=value, colour=variable)) + geom_line()
		ctrchm2	<- gamlss(ANS_CALL~FRQ, sigma.formula=~FRQ, data=ctrch, family=BB(), i.control=glim.control(cyc=100,cc=1e-4))
		ctrch[, PR_CALL2:= predict(ctrchm2, type='response', what='mu')]
		ggplot(melt(ctrch, id.vars=c('FRQ','GPS'), measure.vars=c('ANS_CALLrm','PR_CALL','PR_CALL2')), aes(x=FRQ, y=value, colour=variable)) + geom_line()
		
		ctrchs	<- ctrch[, {
					z	<- sample(length(AGRpc), min(length(AGRpc),50))
					list(AGRpc=AGRpc[z], GPS=GPS[z], ANS_CALL=ANS_CALL[z], CHUNK=CHUNK[z])	
				}, by=c('AGRpcC','GPSC')]
		
		ctrchms	<- gamlss(ANS_CALL~AGRpc+GPS, data=ctrchs, family=BI())
		ctrchm2	<- gamlss(ANS_CALL~AGRpc+GPS, data=ctrchs, family=BB())
		ctrchm3	<- gamlss(ANS_CALL~AGRpc+GPS, sigma.formula=~AGRpc+GPS, data=ctrchs, family=BB())
		#	predict sigma: how does it look?
		ctrchsp	<- as.data.table(expand.grid(AGRpc=seq(0.01,0.99,0.01), GPS=seq(0.1,0.9,0.1)))
		mu		<- predict(ctrchm2, newdata=ctrchsp, what='mu', type='response', se.fit=TRUE)
		
		predict(ctrchm2, what='mu', type='response', se.fit=TRUE)
		max(predict(ctrchm2, what='mu', type='response', se.fit=TRUE)$se.fit)
		
		ctrchsp[, mu:=mu]
		ctrchsp[, TYPE:='s']
		ctrchp	<- as.data.table(expand.grid(AGRpc=seq(0.01,0.99,0.01), GPS=seq(0.1,0.9,0.1)))
		mu		<- predict(ctrchm, newdata=ctrchp, what='mu', type='response')
		ctrchp[, mu:=mu]
		ctrchp[, TYPE:='a']
		ctrchp	<- rbind(ctrchp, ctrchsp)
		ggplot(ctrchp, aes(x=AGRpc, y=mu, group=TYPE, colour=TYPE)) + geom_line() + facet_grid(GPS~.)
		sigma	<- predict(ctrchm3, newdata=ctrchp, what='sigma', type='response')				
		ctrchp[, mu:=mu]
		ctrchp[, sigma:=sigma]
		ctrchp[, varn1:=mu*(1-mu)]
		
		
		
		#ctrchm	<- gamlss(ANS_CALL~AGRpc+bs(GPSc,df=3), data=ctrch, family=BI())
		
		
		
		
		ctrchm	<- gamlss(ANS_CALL~AGRpc, data=ctrch, family=BI())				#as good as 'AGRpc+GPS+CNS_FRQr' in terms of FN, FP
		ctrchmc	<- data.table(CHUNK=chunk,BETA0=coef(ctrchm)[1], BETA1=coef(ctrchm)[2])
		#ctrch[, PR_CALL:= predict(ctrchm, type='response', what='mu')]
		ctrch	<- merge(ctrch, ctrchmc, by='CHUNK')
		ctrch[, PR_CALL2:= exp(BETA0+BETA1*AGRpc)/(exp(BETA0+BETA1*AGRpc)+1)]
		
		
		ggplot(melt(ctrch, id.vars='AGRpc', measure.vars=c('ANS_CALLrm','PR_CALL')), aes(x=AGRpc, y=value, colour=variable)) + geom_line()
		ggsave(file=paste(outdir,'/',outfile,'_SITE',chunk,'_ANS_CALLrmBYAGRpc.pdf',sep=''), w=6, h=4)
		
		
		plot(ctrchm)	#looks pretty good
		Rsq(ctrchm)		#60%
		AIC(ctrchm)
		
		
		ggplot(ctrch, aes(x=PR_CALL, y=ANS_CALL)) + geom_point()
		#	get sens and spec
		ctrev	<- data.table(THR= seq(0.05,0.95,0.01))
		ctrev	<- ctrev[, {
					tmp	<- ctrch[, table(ANS_CALL, PR_CALL>=THR)]
					list(TP= tmp['1','TRUE'], FP= tmp['0','TRUE'], FN= tmp['1','FALSE'], TN= tmp['0','FALSE'] )					
				}, by='THR']
		ctrev	<- merge(ctrev, ctrev[, list(SENS= TP/(TP+FN), SPEC=TN/(TN+FP), FDR=FP/(FP+TP), FOR=FN/(FN+TN)), by='THR'], by='THR')
		ggplot(melt(ctrev, id.vars='THR', measure.vars=c('SENS','SPEC','FDR','FOR')), aes(x=THR, y=value, colour=variable)) + geom_line() + theme(legend.position='bottom')
		ggsave(file=paste(outdir,'/',outfile,'_SITE',chunk,'_SensSpec.pdf',sep=''), w=4, h=4)
		ggplot(ctrev, aes(x=1-SPEC, y=SENS)) + geom_line() + scale_x_continuous(limit=c(0,1),expand=c(0,0)) + scale_y_continuous(limit=c(0,1),expand=c(0,0))
		ggsave(file=paste(outdir,'/',outfile,'_SITE',chunk,'_ROC.pdf',sep=''), w=4, h=4)
		#ctrchm2	<- gamlss(ANS_CALL~AGRpc+GPS+CNS_FRQr, data=ctrch, family=BB())	#BB does not always converge
	}
}
##--------------------------------------------------------------------------------------------------------
pipeline.various<- function()
{
	if(0)
	{		
		indir.cut	<- paste(DATA, 'contigs_150902_unaligned_cut', sep='/' )
		indir.raw	<- paste(DATA, 'contigs_150902_unaligned_raw', sep='/' )
		outdir		<- paste(DATA,'contigs_150902_model150816b',sep='/')
		batch.n		<- 200
		tmp			<- data.table(FILE=list.files(indir.raw, pattern='fasta$', recursive=T))
		tmp[, BATCH:= ceiling(seq_len(nrow(tmp))/batch.n)]
		tmp			<- tmp[, max(BATCH)]
		for(batch.id in seq.int(1,tmp))
		{	
			cmd			<- cmd.haircut.pipeline(indir.cut, indir.raw, outdir, batch.n=batch.n, batch.id=batch.id)
			cmd			<- cmd.hpcwrapper(cmd, hpc.nproc= 1, hpc.q='pqeelab', hpc.walltime=1, hpc.mem="5000mb")
			#cat(cmd)		
			outfile		<- paste("hrct",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
			cmd.hpccaller(paste(DATA,"tmp",sep='/'), outfile, cmd)	
		}
	}
	if(0)
	{		
		indir.cut	<- paste(DATA, 'contigs_150408_unaligned_cut', sep='/' )
		indir.raw	<- paste(DATA, 'contigs_150408_unaligned_raw', sep='/' )
		outdir		<- paste(DATA,'contigs_150408_model150816b',sep='/')
		batch.n		<- 200
		tmp			<- data.table(FILE=list.files(indir.raw, pattern='fasta$', recursive=T))
		tmp[, BATCH:= ceiling(seq_len(nrow(tmp))/batch.n)]
		tmp			<- tmp[, max(BATCH)]
		for(batch.id in seq.int(1,tmp))
		{	
			cmd			<- cmd.haircut.pipeline(indir.cut, indir.raw, outdir, batch.n=batch.n, batch.id=batch.id)
			cmd			<- cmd.hpcwrapper(cmd, hpc.nproc= 1, hpc.q='pqeelab', hpc.walltime=1, hpc.mem="5000mb")
			#cat(cmd)		
			outfile		<- paste("hrct",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
			cmd.hpccaller(paste(DATA,"tmp",sep='/'), outfile, cmd)	
		}
	}
	if(0)
	{		
		indir.cut	<- paste(DATA, 'contigs_150408_unaligned_cut', sep='/' )
		indir.raw	<- paste(DATA, 'contigs_150408_unaligned_raw', sep='/' )
		outdir		<- paste(DATA, 'contigs_150408_wref', sep='/' )
		indir.cut	<- paste(DATA, 'contigs_150902_unaligned_cut', sep='/' )
		indir.raw	<- paste(DATA, 'contigs_150902_unaligned_raw', sep='/' )
		outdir		<- paste(DATA, 'contigs_150902_wref', sep='/' )		
		indir.cut	<- paste(DATA, 'contigs_151026_unaligned_cut', sep='/' )
		indir.raw	<- paste(DATA, 'contigs_151026_unaligned_raw', sep='/' )
		outdir		<- paste(DATA, 'contigs_151026_wref', sep='/' )		
		
		batch.n		<- 100
		tmp			<- data.table(FILE=list.files(indir.raw, pattern='fasta$', recursive=T))
		tmp[, BATCH:= ceiling(seq_len(nrow(tmp))/batch.n)]
		tmp			<- tmp[, max(BATCH)]
		for(batch.id in seq.int(1,tmp))
		{	
			cmd			<- cmdwrap.align.contigs.with.ref(indir.cut, indir.raw, outdir, batch.n=batch.n, batch.id=batch.id)
			cmd			<- cmd.hpcwrapper(cmd, hpc.nproc= 1, hpc.q='pqeph', hpc.walltime=1, hpc.mem="1000mb")
			cat(cmd)		
			outfile		<- paste("hrct",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
			cmd.hpccaller(DATA, outfile, cmd)	
		}
	}
	if(1)
	{
		indir		<- paste(DATA, 'contigs_150408_wref', sep='/' )
		outdir		<- paste(DATA, 'contigs_150408_wref_cutstat', sep='/' )
		indir		<- paste(DATA, 'contigs_150902_wref', sep='/' )
		outdir		<- paste(DATA, 'contigs_150902_wref_cutstat', sep='/' )		
		indir		<- paste(DATA, 'contigs_151026_wref', sep='/' )
		outdir		<- paste(DATA, 'contigs_151026_wref_cutstat', sep='/' )		
		
		batch.n		<- 200
		tmp			<- data.table(FILE=list.files(indir, pattern='wRefs\\.fasta$', recursive=T))
		tmp[, BATCH:= ceiling(seq_len(nrow(tmp))/batch.n)]
		tmp			<- tmp[, max(BATCH)]
		for(batch.id in seq.int(1,tmp))
		{			
			cmd			<- cmd.haircut.cutstat(indir, outdir, batch.n=batch.n, batch.id=batch.id)
			cmd			<- cmd.hpcwrapper(cmd, hpc.nproc= 1, hpc.q='pqeelab', hpc.walltime=4, hpc.mem="5000mb")
			cat(cmd)		
			cmd.hpccaller(DATA, paste("hrct",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.'), cmd)	
		}	
	}
	if(0)
	{
		indir.st	<- paste(DATA,'contigs_150408_wref_cutstat',sep='/')
		indir.al	<- paste(DATA,'contigs_150408_wref',sep='/')
		outdir		<- paste(DATA,'contigs_150408_model150816a',sep='/')
		trainfile	<- paste(DATA,'contigs_150408_trainingset_subsets.R',sep='/')
		indir.st	<- paste(DATA,'contigs_150902_wref_cutstat',sep='/')
		indir.al	<- paste(DATA,'contigs_150902_wref',sep='/')
		outdir		<- paste(DATA,'contigs_150902_model150816b',sep='/')
		trainfile	<- NA
		batch.n		<- 200
		
		tmp			<- data.table(INFILE=list.files(indir.al, pattern='wRefs\\.fasta$', recursive=T))
		tmp[, BATCH:= ceiling(seq_len(nrow(tmp))/batch.n)]
		tmp			<- tmp[, max(BATCH)]
		for(batch.id in seq.int(1,tmp))
		{			
			cmd			<- cmd.haircut.call(indir.st, indir.al, outdir, trainfile=trainfile, batch.n=batch.n, batch.id=batch.id)
			cmd			<- cmd.hpcwrapper(cmd, hpc.nproc= 1, hpc.q='pqeelab', hpc.walltime=4, hpc.mem="5000mb")
			cat(cmd)		
			cmd.hpccaller(paste(DATA,"tmp",sep='/'), paste("hrct",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.'), cmd)	
		}	
	}
	if(0)
	{
		#haircut.QC.flatten.curated()
		haircut.QC.flatten.automated()
		haircut.QC.align.curated()
	}
}
##--------------------------------------------------------------------------------------------------------
##	HAIRCUT program, version 15086 to: 
##	- align contigs to references
##	- calculate and save haircut statistics
##--------------------------------------------------------------------------------------------------------
prog.haircut.150806<- function()
{	
	if(0)
	{
		indir		<- paste(DATA, 'contigs_150408_wref', sep='/' )
		outdir		<- paste(DATA, 'contigs_150408_wref_cutstat', sep='/' )		
		par			<- c('FRQx.quantile'=NA, 'FRQx.thr'=NA, 'CNS_FRQ.window'=200, 'CNS_AGR.window'=200, 'GPS.window'=200)
		batch.n		<- NA
		batch.id	<- NA
		haircutwrap.get.cut.statistics(indir, par, outdir=outdir, batch.n=batch.n, batch.id=batch.id)
	}
	if(0)
	{
		indir		<- paste(DATA, 'contigs_150408_wref', sep='/' )
		outdir		<- paste(DATA, 'contigs_150408_wref_cutstat', sep='/' )		
		batch.n		<- 200
		tmp			<- data.table(FILE=list.files(indir, pattern='fasta$', recursive=T))
		tmp[, BATCH:= ceiling(seq_len(nrow(tmp))/batch.n)]
		tmp			<- tmp[, max(BATCH)]
		for(batch.id in seq.int(1,tmp))
		{	
			
			cmd			<- cmd.haircut.pipeline(indir, outdir, batch.n=batch.n, batch.id=batch.id)
			cmd			<- cmd.hpcwrapper(cmd, hpc.nproc= 1, hpc.q='pqeelab', hpc.walltime=4, hpc.mem="5000mb")
			cat(cmd)		
			cmd.hpccaller(paste(DATA,"tmp",sep='/'), paste("hrct",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.'), cmd)	
		}	
	}
	if(0)
	{
		indir		<- paste(DATA, 'contigs_151026_wref', sep='/' )
		outdir		<- paste(DATA, 'contigs_151026_wref_cutstat', sep='/' )				
		cat(cmd.haircut.cutstat(indir, outdir, batch.n=200, batch.id=1))		
	}
	if(0)
	{
		#	get model coefficients across the chunks
		mfile					<- paste(DATA,'model_150816a.R',sep='/')		
		tmp						<- haircut.get.fitted.model.150816a(NA, mfile)
		ctrmc					<- tmp$coef		
		predict.fun				<- tmp$predict
		#	get contigs that were used for training
		outfile	<- paste(DATA,'contigs_150408_trainingset_subsets.R',sep='/')
		ctrain	<- haircut.get.training.contigs(NA, outfile, NA)
		set(ctrain, NULL, 'CUT', ctrain[, factor(CUT, levels=c('cut','raw'), labels=c('Y','N'))])
		setnames(ctrain, 'CUT', 'BLASTnCUT')		
		#	get covariates for all contigs
		indir.st<- paste(DATA,'contigs_150408_wref_cutstat',sep='/')
		indir.al<- paste(DATA,'contigs_150408_wref',sep='/')
		outdir	<- paste(DATA,'contigs_150408_model150816a',sep='/')
		par		<- c(	'FRQx.quantile'=NA, 'FRQx.thr'=NA, 'CNS_FRQ.window'=200, 'CNS_AGR.window'=200, 'GPS.window'=200, 
				'PRCALL.thrmax'=0.8, 'PRCALL.thrstd'=10, 'PRCALL.cutprdcthair'=150, 'PRCALL.cutprdctcntg'=50, 'PRCALL.cutrawgrace'=100, 'PRCALL.rmintrnlgpsblw'=100 ,'PRCALL.rmintrnlgpsend'=9700,
				'PRCALL.mxgpinref'=100)
		haircutwrap.get.call.for.PNG_ID(indir.st,indir.al,outdir,ctrmc,predict.fun,par,ctrain=ctrain)
	}
}
##--------------------------------------------------------------------------------------------------------
##	HAIRCUT program to call parts of contigs
##--------------------------------------------------------------------------------------------------------
haircutprog.get.call.for.PNG_ID<- function()
{
	verbose			<- 1
	mfile			<- paste(DATA,'model_150816a.R',sep='/')
	trainfile		<- paste(DATA,'contigs_150408_trainingset_subsets.R',sep='/')
	indir.st		<- paste(DATA,'contigs_150408_wref_cutstat',sep='/')
	indir.al		<- paste(DATA,'contigs_150408_wref',sep='/')
	outdir			<- paste(DATA,'contigs_150408_model150816a',sep='/')
	batch.n			<- NA
	batch.id		<- NA
	
	if(exists("argv"))
	{
		#	args input
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,6),
									mfile= return(substr(arg,8,nchar(arg))),NA)	}))
		if(length(tmp)>0) mfile<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									trainfile= return(substr(arg,12,nchar(arg))),NA)	}))
		if(length(tmp)>0) trainfile<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,9),
									indir.st= return(substr(arg,11,nchar(arg))),NA)	}))
		if(length(tmp)>0) indir.st<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,9),
									indir.al= return(substr(arg,11,nchar(arg))),NA)	}))
		if(length(tmp)>0) indir.al<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									outdir= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) outdir<- tmp[1]			
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,8),
									batch.n= return(as.numeric(substr(arg,10,nchar(arg)))),NA)	}))
		if(length(tmp)>0) batch.n<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,9),
									batch.id= return(as.numeric(substr(arg,11,nchar(arg)))),NA)	}))
		if(length(tmp)>0) batch.id<- tmp[1]
	}
	if(verbose)
	{
		cat('\ninput args\n',paste(mfile, trainfile, indir.st, indir.al, outdir, batch.n, batch.id, sep='\n'))
	}	
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
}