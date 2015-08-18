##--------------------------------------------------------------------------------------------------------
##	get contigs to be used for training. this excludes duplicate raw/cut contigs for which a curated answer is available 
##--------------------------------------------------------------------------------------------------------
haircut.get.training.contigs<- function(indir, outfile, ctrain)
{
	options(show.error.messages = FALSE)		
	readAttempt		<-try(suppressWarnings(load(outfile)))
	options(show.error.messages = TRUE)	
	if( inherits(readAttempt, "try-error")	)
	{
		infiles	<- data.table(FILE=list.files(indir, pattern='fasta$', recursive=T))
		infiles[, PNG_ID:= gsub('_wRefs\\.fasta','',gsub('_cut|_raw','',basename(FILE)))]
		infiles[, BLASTnCUT:= regmatches(FILE,regexpr('cut|raw',FILE))]		 
		#	identify identical raw and cut contigs, and identify cut contigs that are part of the raw contig
		txe		<- infiles[,{
					#png_id	<- '12559_1_11'
					#files	<- subset(infiles, PNG_ID==png_id)[, FILE]
					#blastncut<- subset(infiles, PNG_ID==png_id)[, BLASTnCUT]
					cat('\nProcess', PNG_ID)
					tx	<- haircut.get.subset.among.raw.and.cut.contigs(indir, FILE, PNG_ID, BLASTnCUT)
					tx
				}, by='PNG_ID']
		#	consider identical contigs. merge cut contigs into equal contigs, so we can look at all equal pairs
		txec	<- merge(subset(txe, EQ), ctrain, all.x=TRUE, by=c('PNG_ID','TAXON'))		
		tmp		<- subset(txe, EQ & CUT=='cut', select=c(PNG_ID, OCNTG, CCNTG, TAXON, EQ_FIRST, EQ_LAST))
		setnames(tmp, c('TAXON','CCNTG','EQ_FIRST','EQ_LAST'),c('TAXON_CUT','CCNTG_CUT','EQ_FIRST_CUT','EQ_LAST_CUT'))
		txec	<- merge(txec, tmp, all.x=TRUE, allow.cartesian=TRUE,by=c('PNG_ID','OCNTG'))
		tmp		<- subset(ctrain, select=c(PNG_ID, TAXON, ANS_LEN))
		setnames(tmp, c('TAXON','ANS_LEN'), c('TAXON_CUT','ANS_LEN_CUT'))
		txec	<- merge(txec, tmp, all.x=1, by=c('PNG_ID','TAXON_CUT'))
		#	work out which contigs to use for training
		txec[, USE_IN_TRAIN:='Y']
		#	if raw and cut contigs in curated and identical: don t use raw (double counting).
		tmp		<- txec[, which(CUT=='raw' & is.na(CCNTG_CUT) & !is.na(ANS_LEN) & !is.na(ANS_LEN_CUT))]
		cat(paste('\nFound raw and cut contigs in curated that are identical: don t use raw for training to avoid double counting, n=', length(tmp)))
		set(txec, tmp, 'ANS_FILE', NA_character_)
		set(txec, tmp, c('ANS_FIRST','ANS_LAST'), NA_integer_)	#need to set ANS_LEN to NA later
		set(txec, tmp, 'USE_IN_TRAIN', 'N')
		#			
		#	if raw and concatenated cut contigs identical and only raw in curated: don t use cut in training as 0
		tmp		<- subset(txec, USE_IN_TRAIN=='Y' & CUT=='raw' & is.na(CCNTG_CUT) & !is.na(ANS_LEN) & is.na(ANS_LEN_CUT))[, TAXON_CUT] 
		cat(paste('\nFound raw and cut contigs identical and only raw in curated: don t use cut in training as 0, n=', length(tmp)))
		set(txec, txec[, which(TAXON%in%tmp & CUT=='cut')], 'USE_IN_TRAIN', 'N')
		#	if raw and concatenated cut contigs identical and only cut in curated: don t use raw in training as 0
		tmp		<- txec[, which(USE_IN_TRAIN=='Y' & CUT=='raw' & is.na(CCNTG_CUT) & is.na(ANS_LEN) & !is.na(ANS_LEN_CUT))] 
		cat(paste('\nFound raw and cut contigs identical and only cut in curated: don t use raw in training as 0, n=', length(tmp)))
		set(txec, tmp, 'USE_IN_TRAIN', 'N')
		#
		#	if cut contig subset of raw contig and both in curated: should not happen
		tmp		<- nrow(subset(txec, USE_IN_TRAIN=='Y' & CUT=='raw' & !is.na(CCNTG_CUT) & !is.na(ANS_LEN) & !is.na(ANS_LEN_CUT)))
		stopifnot(tmp==0)
		#	if cut contig subset of raw contig and only cut in curated: use cut only partial
		tmp		<- txec[, which(USE_IN_TRAIN=='Y' & CUT=='cut' & !is.na(CCNTG_CUT) & !is.na(EQ_FIRST) & !is.na(ANS_LEN) & !is.na(ANS_LEN_CUT))]
		cat(paste('\nFound cut contig subset of raw contig and only cut in curated: use cut only partial, n=', length(tmp)))
		set(txec, tmp, 'USE_IN_TRAIN', 'P')
		#	if cut contig subset of raw contig and only cut in curated: dont use raw as 0's
		tmp		<- txec[, which(USE_IN_TRAIN=='Y' & CUT=='raw' & !is.na(CCNTG_CUT) & !is.na(EQ_FIRST_CUT) & is.na(ANS_LEN) & !is.na(ANS_LEN_CUT))]
		cat(paste('\nFound cut contig subset of raw contig and only cut in curated: dont use raw as 0s, n=', length(tmp)))
		set(txec, tmp, 'USE_IN_TRAIN', 'N')
		#	if cut contig subset of raw contig and only cut in curated: use cut only partial
		tmp		<- subset(txec, CUT=='raw' & !is.na(CCNTG_CUT) & !is.na(EQ_FIRST_CUT) & is.na(ANS_LEN) & !is.na(ANS_LEN_CUT))[, TAXON_CUT]
		tmp		<- txec[, which(TAXON%in%tmp & CUT=='cut')]
		cat(paste('\nFound cut contig subset of raw contig and only cut in curated: dont use raw as 0s, n=', length(tmp)))		
		set(txec, tmp, 'USE_IN_TRAIN', 'P')
		#	can now set ANS_LEN to NA from above case 
		set(txec, txec[, which(is.na(ANS_FIRST) & !is.na(ANS_LEN))], 'ANS_LEN', NA_integer_)
		#	add contigs that are in curated and unique amongst cut and raw
		tmp		<- merge(subset(txe, !EQ), ctrain, by=c('PNG_ID','TAXON'))
		tmp[, USE_IN_TRAIN:='Y']
		stopifnot( length(intersect( tmp[, TAXON], txec[,TAXON] ))==0	)
		txec	<- rbind(txec, tmp, use.names=TRUE, fill=TRUE)
		#	if cut contig subset of raw contig and only raw in curated: dont use cut as 0's
		tmp		<- subset(txec, CUT=='raw' & !is.na(CCNTG_CUT) & !is.na(EQ_FIRST_CUT) & !is.na(ANS_LEN) & is.na(ANS_LEN_CUT), c(PNG_ID, OCNTG))
		tmp		<- subset( merge(tmp, txe, by=c('PNG_ID','OCNTG')), CUT=='cut' )[, unique(TAXON)]
		cat(paste('\nFound cut contig subset of raw contig and only raw in curated: dont use cut as 0, n=', length(tmp)))
		tmp		<- txec[, which(TAXON%in%tmp & CUT=='cut')]
		set(txec, tmp, 'USE_IN_TRAIN', 'N')
		#	delete tmp cols
		set(txec, NULL, c('TAXON_CUT','CCNTG_CUT','ANS_LEN_CUT','EQ_FIRST_CUT','EQ_LAST_CUT'), NULL)
		setkey(txec, PNG_ID, TAXON, CUT)
		txec	<- unique(txec)					
		save(txec, file=outfile)
		#print( txec[, table(USE_IN_TRAIN)] )
	}
	txec
}
##--------------------------------------------------------------------------------------------------------
##	get training data set: contig statistics aligned to sites where curated contigs were not cut 
##--------------------------------------------------------------------------------------------------------
haircut.get.training.data<- function(indir, ctrain, par, outdir, outfile)
{
	require(plyr)
	infiles	<- data.table(INFILE=list.files(indir, pattern='\\.R$', recursive=T))
	infiles[, BATCH:= ceiling(seq_len(nrow(infiles))/100)]
	for(b in infiles[, unique(BATCH)])
	{
		#	add curated answers to data by batch
		tmp		<- subset(infiles, BATCH==b)[, {
					cat(paste('\nProcess', INFILE ))
					#	INFILE	<- infiles[12,INFILE]
					load( paste(indir,INFILE,sep='/'))						
					#	select contig+refs for which curated answer available				
					#cm		<- merge( unique(subset(cnsc.df, select=c(PNG_ID, TAXON, BLASTnCUT))), subset( ctrain, select=c(PNG_ID, TAXON, BLASTnCUT, USE_IN_TRAIN, EQ_FIRST, EQ_LAST, ANS_FILE, ANS_LEN, ANS_FIRST, ANS_LAST) ), all.x=1, by=c('PNG_ID','TAXON','BLASTnCUT') )
					cm		<- merge( unique(subset(cnsc.df, TAXON!='consensus',select=c(PNG_ID, TAXON, BLASTnCUT))), subset( ctrain, select=c(PNG_ID, TAXON, BLASTnCUT, USE_IN_TRAIN, EQ_FIRST, EQ_LAST, ANS_FILE, ANS_LEN, ANS_FIRST, ANS_LAST) ), all.x=1, by=c('PNG_ID','TAXON','BLASTnCUT') )
					#	code below identifies ANS_CALL=1
					#	we also need ANS_CALL=0; these are all those that don t have a curated file + all those that have USE_IN_TRAIN=='Y'					
					cm		<- subset(cm, is.na(USE_IN_TRAIN) |  USE_IN_TRAIN=='Y' |  USE_IN_TRAIN=='P')
					if( cm[, all(is.na(ANS_LEN))] )
					{
						cat(paste('\nNo data contigs matched in training data set', INFILE ))
						ca	<- merge(subset(cnsc.df, select=c(TAXON, SITE, FRQ, FRQ_STD, AGRpc, GPS, PNG_ID, BLASTnCUT)), subset(cm, select=TAXON), by='TAXON')
						ca[, ANS_CALL:=0L]							
					}
					if( cm[, !all(is.na(ANS_LEN))] && ncol(cnsc)!=subset(cm, !is.na(ANS_LEN))[1, ANS_LEN] )
					{			
						cat(paste('\nMatching contigs exist in training data set, but are not necessarily aligned', INFILE )) 
						cr		<- read.dna(subset(cm, !is.na(ANS_LEN))[1,  ANS_FILE], format='fasta')
						#	determine start of non-LTR position and cut 
						cr		<- cr[, seq.int(haircut.find.nonLTRstart(cr), ncol(cr))]
						#	determine consensus to calculate offset
						tmp		<- strsplit(basename(subset(cm, !is.na(ANS_LEN))[1,  ANS_FILE]), '_')[[1]][1]
						tx		<- data.table(TAXON= rownames(cr), CONTIG=as.integer(grepl(tmp, rownames(cr))) )			
						tmp		<- cr[subset(tx, CONTIG==0)[, TAXON],]
						rp		<- haircut.get.frequencies(tmp, bases=c('a','c','g','t','-') )
						crc		<- haircut.get.consensus.from.frequencies(rp, par)$DNAbin
						cdc		<- cnsc['consensus', ]
						#	calculate site offset in the curated sequence (crc)
						offset	<- haircut.calculate.offset(crc,cdc)
						#	update ANS_FIRST and ANS_LAST according to offset of curated contig in data contigs
						cm		<- cm[, list(ANS_LEN=ANS_LEN+offset[ANS_LEN], ANS_FIRST=ANS_FIRST+offset[ANS_FIRST], ANS_LAST=ANS_LAST+offset[ANS_LAST]), by=c('PNG_ID','TAXON','BLASTnCUT','USE_IN_TRAIN','EQ_FIRST','EQ_LAST')]
						if(any(offset!=0) && ncol(cnsc)!=subset(cm, !is.na(ANS_LEN))[1, ANS_LEN])
							warning(paste('\nPerhaps check: Found unequal lengths for', INFILE))
					}
					if( cm[, !all(is.na(ANS_LEN))] )
					{
						#	expand to SITEs that were called manually
						ca		<- subset(cm, !is.na(ANS_LEN))[, list(SITE=seq.int(ANS_FIRST,ANS_LAST), ANS_CALL=1L), by='TAXON']
						#	keep only those contigs that are not in USE_IN_TRAIN=='N'
						tmp		<- merge(subset(cnsc.df, select=c(TAXON, SITE, FRQ, FRQ_STD, AGRpc, GPS, PNG_ID, BLASTnCUT)), subset(cm, select=c(TAXON, USE_IN_TRAIN, EQ_FIRST, EQ_LAST)), by='TAXON')
						ca		<- merge(tmp, ca, all.x=1, by=c('TAXON','SITE'))
						ca		<- subset(ca, USE_IN_TRAIN=='Y' | is.na(USE_IN_TRAIN) | (USE_IN_TRAIN=='P' & SITE>=EQ_FIRST & SITE<=EQ_LAST))
						set(ca, ca[,which(is.na(ANS_CALL))],'ANS_CALL',0L)
						set(ca, NULL, c('USE_IN_TRAIN','EQ_FIRST','EQ_LAST'), NULL)
					}
					
					tmp		<- subset(cnsc.df, TAXON=='consensus')
					setnames(tmp, c('FRQ','FRQ_STD','AGRpc','GPS'), c('CNS_FRQ','CNS_FRQ_STD','CNS_AGRpc','CNS_GPS'))
					ca		<- merge(ca, subset(tmp, select=c(SITE, CNS_FRQ, CNS_FRQ_STD, CNS_AGRpc, CNS_GPS)), by='SITE')
					ca
				}, by='INFILE']
		#	save batches in chunks of sites
		cat(paste('\nSave batch to file', b ))
		tmp[, CHUNK:= cut(SITE, breaks=c(-1,seq.int(200, 10000, 200),Inf), labels=c(seq.int(200, 10200, 200) ))]
		for(c in tmp[, unique(CHUNK)])
		{
			ctr	<- subset(tmp, CHUNK==c)				
			save(ctr, file=paste(outdir, '/', outfile,'_sites',c,'_batch',b,'.R',sep=''))
		}		
		cat(paste('\nSaved batch to file', b ))
	}	
	#	load batches and save chunks
	ofiles	<- data.table(FILE=list.files(outdir, pattern='\\.R$', recursive=T))
	ofiles	<- subset(ofiles, !grepl('SITES', FILE))
	set(ofiles, NULL, 'SITES', ofiles[, substring(regmatches(FILE, regexpr('sites[0-9]+',FILE)), 6)])
	ofiles[, {
				ctr	<- do.call('rbind',lapply(FILE, function(x)
								{
									load(paste(outdir,x,sep='/'))
									ctr
								}))
				save(ctr, file=paste(outdir, '/', outfile,'_SITES',SITES,'.R',sep=''))
			}, by='SITES']
	#	rm intermediates files
	invisible( file.remove(ofiles[, paste(outdir,FILE,sep='/')]) )
	#Perhaps check: Found unequal lengths for 15065_1_10_cut_wRefs_HAIRCUTSTAT_thr5_aw200_fw100_gw200.R	1367
	#Perhaps check: Found unequal lengths for 15070_1_9_raw_wRefs_HAIRCUTSTAT_thr5_aw200_fw100_gw200.R	1522
	#Perhaps check: Found unequal lengths for 15099_1_87_raw_wRefs_HAIRCUTSTAT_thr5_aw200_fw100_gw200.R	1924
	#Perhaps check: Found unequal lengths for 15172_1_43_cut_wRefs_HAIRCUTSTAT_thr5_aw200_fw100_gw200.R	2325
	NULL
}
##--------------------------------------------------------------------------------------------------------
##	determine if cut contigs are identical to a stretch or possibly the whole raw contig with the same contig ID 
##--------------------------------------------------------------------------------------------------------
haircutwrap.get.subset.among.raw.and.cut.contigs<- function(indir, outfile)
{
	options(show.error.messages = FALSE)		
	readAttempt		<-try(suppressWarnings(load(outfile)))
	options(show.error.messages = TRUE)	
	if( inherits(readAttempt, "try-error")	)
	{
		infiles	<- data.table(FILE=list.files(indir, pattern='fasta$', recursive=T))
		infiles[, PNG_ID:= gsub('_wRefs\\.fasta','',gsub('_cut|_raw','',basename(FILE)))]
		infiles[, BLASTnCUT:= regmatches(FILE,regexpr('cut|raw',FILE))]		 
		#	identify identical raw and cut contigs, and identify cut contigs that are part of the raw contig
		txe		<- infiles[,{
					#png_id	<- '12559_1_11'
					#files	<- subset(infiles, PNG_ID==png_id)[, FILE]
					#blastncut<- subset(infiles, PNG_ID==png_id)[, BLASTnCUT]
					cat('\nProcess', PNG_ID)
					tx	<- haircut.get.subset.among.raw.and.cut.contigs(indir, FILE, PNG_ID, BLASTnCUT)
					tx
				}, by='PNG_ID']
		save(txe, file=outfile)
	}
	txe
}
##--------------------------------------------------------------------------------------------------------
##	return EQ= 1 and EQ_FIRST EQ_LAST positions, if cut contig are subsets of the raw contig in that: 
##	- 	there is only once cut contig which disagrees on up to x positions, where x is the difference in alignment lengths in the cut and raw files
##	- 	the cut contigs matches the raw contig without gaps, where the raw contig is set to start from the start position of the cut contig
##--------------------------------------------------------------------------------------------------------
haircut.get.subset.among.raw.and.cut.contigs<- function(indir, files, png_id, blastncut)
{
	crs		<- lapply(files, function(x)
			{
				cr		<- read.dna(file=paste(indir,'/',x,sep=''), format='fasta')
				if(is.matrix(cr))
				{
					cr		<- cr[, seq.int(haircut.find.nonLTRstart(cr), ncol(cr))]
					tmp		<- strsplit(basename(x), '_')[[1]][1]
					tx		<- data.table(TAXON= rownames(cr), CONTIG=as.integer(grepl(tmp, rownames(cr))) )
					cr		<- cr[ subset(tx, CONTIG==1)[, TAXON], ]	
				}
				if(!is.matrix(cr))
					cr		<- as.DNAbin(matrix(vector('character',0), nrow=2, byrow=T, dimnames=list(c('DUMMY','DUMMY'),c())))
				cr					
			})
	names(crs)	<- blastncut
	#	get contig table
	tx		<- do.call('rbind',lapply(seq_along(crs), function(i)	data.table(	TAXON=rownames(crs[[i]]), 
								CUT= blastncut[i], 
								FIRST= apply( as.character(crs[[i]]), 1, function(x) which(x!='-')[1] ),
								LAST= ncol(crs[[i]])-apply( as.character(crs[[i]]), 1, function(x) which(rev(x)!='-')[1] ) + 1L,
								CRS_ID=i)	))
	tx		<- subset(tx, !is.na(FIRST) & !is.na(LAST))	#some contigs may just be in LTR
	tx[, CNTG:=tx[, gsub(paste(png_id,'.',sep=''),'',substring(TAXON, regexpr(png_id, TAXON)))]]
	tx[, OCNTG:= tx[, sapply(strsplit(CNTG,'.',fixed=T),'[[',1)]]
	tx[, CCNTG:= NA_character_]		
	tx[, EQ:=FALSE]
	tx[, EQ_FIRST:=NA_integer_]
	tx[, EQ_LAST:=NA_integer_]
	tmp		<- tx[, which(grepl('.',CNTG,fixed=T))]
	if(length(tmp))
		set(tx, tmp, 'CCNTG', tx[tmp, sapply(strsplit(CNTG,'.',fixed=T),'[[',2)])	
	tmp		<- subset(tx, CUT=='cut' & !is.na(CCNTG))[, list(CCNTGn=length(CCNTG)), by='OCNTG']
	tmp		<- subset(tmp, CCNTGn==1)[, OCNTG]	#check for cut contigs that should be present in multiple cuts by naming scheme, but after LTR removal there is only one cut
	if(length(tmp))
	{
		cat('\nFound lone cuts for which multiple cuts are expected by naming scheme, n=', length(tmp))
		tmp	<- tx[, which(CUT=='cut' & OCNTG%in%tmp)]
		set(tx, tmp, 'CCNTG', NA_character_)
		set(tx, tmp, 'CNTG', tx[tmp, OCNTG])
	}#	check if there are contigs with corresponding name in 'cut' and that are identical / subset
	#	differences in gaps are allowed: these will come from different alignment
	txe		<- dcast.data.table(tx, CNTG~CUT, value.var='TAXON')
	txe		<- subset( txe, !is.na(cut) & !is.na(raw) )
	if(nrow(txe) && c('cut','raw')%in%colnames(txe))
	{						
		txe		<- txe[, {
					x							<- sub('^-*','',sub('-*$','',paste(as.vector(as.character(crs[['cut']][cut,])),collapse='')))
					y							<- sub('^-*','',sub('-*$','',paste(as.vector(as.character(crs[['raw']][raw,])),collapse='')))
					#print(CNTG)
					#print(x)
					#print(y)
					z							<- nchar(x)
					if(nchar(x)==nchar(y))
					{						
						z2	<- x==y
						z3	<- 1L
						attr(z3, 'match.length')<- nchar(y)
					}						
					if(nchar(x)!=nchar(y))
					{
						z	<- gsub('-','',x)
						z2	<- z==gsub('-','',y)
						z3	<- haircut.large.regexpr(z,gsub('-','',y))
						z	<- nchar(x)
					}						
					list(EQtmp=z2, EQ_FIRSTtmp=as.integer(z3), EQ_LASTtmp=as.integer(z3+attr(z3, 'match.length')-1L+nchar(x)-z))					
				}, by='CNTG']
		tx		<- merge(tx, txe, all.x=TRUE, by='CNTG')
		tmp		<- tx[, which(!is.na(EQ_FIRSTtmp) & EQ_FIRSTtmp>0)]
		set(tx, tmp, 'EQ_FIRST', tx[tmp,FIRST]+tx[tmp,EQ_FIRSTtmp]-1L)
		set(tx, tmp, 'EQ_LAST', tx[tmp,FIRST]+tx[tmp,EQ_LASTtmp]-tx[tmp,EQ_FIRSTtmp])
		set(tx, tmp, 'EQ', tx[tmp,EQtmp])
		set(tx, NULL, c('EQ_FIRSTtmp','EQ_LASTtmp','EQtmp'),NULL)				
	}
	#	see if cut contigs are subset of the raw contig
	txe		<- subset(tx, !is.na(CCNTG))
	if(nrow(txe))
	{
		tmp		<- subset(tx, CUT=='raw')
		setnames(tmp, colnames(tmp)[ !grepl('OCNTG',colnames(tmp))], paste( colnames(tmp)[ !grepl('OCNTG',colnames(tmp))], '_raw',sep='' ))
		txe		<- merge(txe, tmp, by='OCNTG')
		txe		<- txe[, {
					z		<- sub('-*$','',paste(as.vector(as.character(crs[['cut']][TAXON,])),collapse=''))
					x		<- sub('^-*','',z)
					tmp		<- nchar(z)-nchar(x)	#number '-' clipped up to start of cut contig
					tmp		<- max(1L, 1L+tmp-abs(diff(c(ncol(crs[['cut']]), ncol(crs[['raw']])))))	#read pos in raw contig
					y		<- sub('-*$','',paste(as.vector(as.character(crs[['raw']][TAXON_raw, seq.int(tmp, ncol(crs[['raw']])) ])),collapse=''))
					#	y starts just a few sites before x should fit, need to determine exact position
					#print(CNTG)
					#print(x)
					#print(y)
					#print(tmp)				
					z		<- gsub('-','',x)	#rm internal '-' because cut and raw are not in same alignment
					z2		<- gsub('-','',y)	#rm internal '-' because cut and raw are not in same alignment
					z3		<- regexpr( substr(z,1,min(9,nchar(z))),  z2 )
					#print(z3)					
					if(z3>0)
					{
						z2	<- substring(z2, z3)
						z2	<- haircut.large.regexpr(z, z2)
					}
					if(	z3<0 || z3> (abs(diff(c(ncol(crs[['cut']]), ncol(crs[['raw']]))))+1)  )	#positioned raw to fit cut within alignment accuracy - reject match if that s not the case
					{
						z3	<- -1L					
					}					
					#print(z2)	
					list(EQ_FIRSTtmp=as.integer(z3), EQ_LASTtmp=as.integer(z3+attr(z2, 'match.length')-1L+nchar(x)-nchar(z)) )					
				}, by=c('CNTG','OCNTG')]
		txe		<- subset(txe, EQ_FIRSTtmp>0)
		tx		<- merge(tx,txe,all.x=TRUE,by=c('CNTG','OCNTG'))
		tmp		<- tx[, which(!is.na(EQ_FIRSTtmp) & EQ_FIRSTtmp>0)]
		set(tx, tmp, 'EQ_FIRST', tx[tmp,FIRST]+tx[tmp,EQ_FIRSTtmp]-1L)
		set(tx, tmp, 'EQ_LAST', tx[tmp,FIRST]+tx[tmp,EQ_LASTtmp]-tx[tmp,EQ_FIRSTtmp])
		set(tx, NULL, c('EQ_FIRSTtmp','EQ_LASTtmp'),NULL)
		set(tx, tx[, which( CUT=='raw' & OCNTG%in%txe[, unique(OCNTG)] )], 'EQ', TRUE)
		set(tx, tx[, which( !is.na(EQ_FIRST))], 'EQ', TRUE)
	} 					
	tx		<- merge(tx, data.table(INFILE=files, CUT=blastncut), by='CUT')
	subset(tx, select=c(INFILE, TAXON, CUT, CNTG, OCNTG, CCNTG, FIRST, LAST, EQ, EQ_FIRST, EQ_LAST))
}
##--------------------------------------------------------------------------------------------------------
##	get curated contigs: names of curated contigs and first/last sites where they were not cut 
##--------------------------------------------------------------------------------------------------------
haircut.get.curated.contigs<- function(indir, outfile)
{
	options(show.error.messages = FALSE)		
	readAttempt		<-try(suppressWarnings(load(outfile)))
	options(show.error.messages = TRUE)
	tmp				<- !inherits(readAttempt, "try-error")
	if(!tmp)
	{
		infiles	<- data.table(FILE=list.files(indir, recursive=T, pattern='fasta$'))
		infiles[, PNG_ID:= gsub('\\.fasta','',basename(FILE))]
		stopifnot( nrow(infiles)==infiles[, length(unique(PNG_ID))] )					
		ctrain	<- infiles[, {
					file	<- paste(indir, FILE, sep='/')
					cat(paste('\nProcess', file))
					#	read Contigs+Rrefs: cr
					cr		<- read.dna(file, format='fasta')
					#	determine start of non-LTR position and cut 
					tmp		<- haircut.find.nonLTRstart(cr)
					cat(paste('\nFound end of LTR at=', tmp-1))
					cr		<- cr[, seq.int(tmp, ncol(cr))]
					#	determine reference sequences. 
					#	non-refs have the first part of the file name in their contig name and are at the top of the alignment
					tmp		<- strsplit(basename(file), '_')[[1]][1]
					tx		<- data.table(TAXON= rownames(cr), CONTIG=as.integer(grepl(tmp, rownames(cr))) )
					stopifnot( all( tx[, which(CONTIG==1)] == seq.int(1, tx[, length(which(CONTIG==1))]) ) )
					cat(paste('\nFound contigs, n=', tx[, length(which(CONTIG==1))]))
					cns		<- cr[subset(tx, CONTIG==1)[, TAXON],]
					#	determine first and last non-gap sites
					tx		<- data.table(	TAXON= rownames(cns), PNG_ID=PNG_ID, LEN= ncol(cr),
							FIRST= apply( as.character(cns), 1, function(x) which(x!='-')[1] ),
							LAST= ncol(cns)-apply( as.character(cns), 1, function(x) which(rev(x)!='-')[1] )+1L		)
					subset(tx, !is.na(FIRST) & !is.na(LAST))	#	some contigs only map into LTR
				}, by='FILE']
		set(ctrain, NULL, 'FILE', ctrain[, paste(indir, FILE, sep='/')])
		save(ctrain, file=outfile)
	}
	ctrain
}
##--------------------------------------------------------------------------------------------------------
##	Binomial model using FRQ and GAP
##--------------------------------------------------------------------------------------------------------
#' @import data.table zoo plyr ape reshape2 ggplot2
#' @export
haircut.get.fitted.model.150816a<- function(indir=NA, outfile=NA)
{
	if(is.na(outfile) | is.na(indir))		#	load from R package
	{
		options(show.error.messages = FALSE)		
		readAttempt		<-try(suppressWarnings(load(system.file(package='PANGEAhaircut', "model_150816a.R"))))
		options(show.error.messages = TRUE)	
		if( inherits(readAttempt, "try-error")	)
			stop('Cannot find model file on system location')
	}
	if(!is.na(outfile) & !is.na(indir))		#	load from outfile, and if not there create
	{
		options(show.error.messages = FALSE)		
		readAttempt		<-try(suppressWarnings(load(outfile)))
		options(show.error.messages = TRUE)	
		if( inherits(readAttempt, "try-error")	)
		{
			ctrmc	<- do.call('rbind', lapply(seq(1,10001,200), function(site)
							{
								cat('\nProcess',site,'\n')							
								ctr		<- haircut.load.training.data(indir, site)	
								tmp		<- seq(ctr[, min(SITE)-1],ctr[, max(SITE)+10],10)
								tmp		<- tmp[tmp<=10000 | tmp==floor(ctr[, max(SITE)+10]/10)*10]
								ctr[, CHUNK:=cut(SITE, breaks=tmp, labels=tmp[-length(tmp)])]							
								ctrmc	<- do.call('rbind',lapply(ctr[, unique(CHUNK)], function(chunk)
												{
													ctrch	<- subset(ctr, CHUNK==chunk)
													ctrchm	<- gamlss(ANS_CALL~FRQ+GPS, data=ctrch, family=BI())				#as good as 'AGRpc+GPS+CNS_FRQr' in terms of FN, FP
													ctrchmc	<- data.table(CHUNK=chunk, BETA0=coef(ctrchm)[1], BETA1=coef(ctrchm)[2], BETA2=coef(ctrchm)[3])					
												}))
							}))
			#	deal with end of genome where little data is available: we estimated coefs for all sites > 1e4, now split up using CHUNK notation
			ctr		<- haircut.load.training.data(indir, 10001)
			tmp		<- seq(ctr[, min(SITE)-1],ctr[, max(SITE)+10],10)
			tmp		<- as.data.table(expand.grid(CHUNK=as.character(tmp[tmp>10000]), BETA0=subset(ctrmc, CHUNK==10000)[, BETA0], BETA1=subset(ctrmc, CHUNK==10000)[, BETA1], BETA2=subset(ctrmc, CHUNK==10000)[, BETA2], stringsAsFactors=F))
			ctrmc	<- rbind(ctrmc, tmp)
			#	model predict function, so we save mem by not having to call 'predict'
			model.150816a.predict<- function(frq, gps, b0, b1, b2)
			{	
				stopifnot(all(!is.na(frq)), all(!is.na(gps)))
				b0[which(is.na(b0))]	<- 0
				b1[which(is.na(b1))]	<- 0
				b2[which(is.na(b2))]	<- 0
				exp(b0+b1*frq+b2*gps)/(exp(b0+b1*frq+b2*gps)+1)	
			}
			#	calculate Sensitivity & Specificity on training data
			ctrev	<- do.call('rbind', lapply(seq(1,10200,200), function(site)
							{
								cat('\nProcess',site,'\n')	
								ctr		<- haircut.load.training.data(indir, site)	
								tmp		<- seq(ctr[, min(SITE)-1],ctr[, max(SITE)+10],10)							
								ctr[, CHUNK:=cut(SITE, breaks=tmp, labels=tmp[-length(tmp)])]
								ctrp	<- merge(ctr, ctrmc, by='CHUNK')
								ctrp[, PR_CALL:=model.150816a.predict(FRQ, GPS, BETA0, BETA1, BETA2)]														
								setkey(ctrp, CHUNK)
								if(0)
								{
									ctrev	<- as.data.table(expand.grid(THR= seq(0.05,0.95,0.01), CHUNK= as.character(ctrp[, unique(CHUNK)]), stringsAsFactors=FALSE))
									ctrev	<- ctrev[, {																																
												tmp	<- ctrp[CHUNK,][,table(ANS_CALL, factor(PR_CALL>=THR, levels=c('TRUE','FALSE'), labels=c('TRUE','FALSE')))]
												list(TP= tmp['1','TRUE'], FP= tmp['0','TRUE'], FN= tmp['1','FALSE'], TN= tmp['0','FALSE'] )					
											}, by=c('CHUNK','THR')]
									ctrev	<- merge(ctrev, ctrev[, list(SENS= TP/(TP+FN), SPEC=TN/(TN+FP), FDR=FP/(FP+TP), FOR=FN/(FN+TN)), by=c('CHUNK','THR')], by=c('CHUNK','THR'))
									
								}
								if(1)
								{
									ctrev	<- as.data.table(expand.grid(THR= c(seq(2,10,1), seq(15,30,5)), CHUNK= as.character(ctrp[, unique(CHUNK)]), stringsAsFactors=FALSE))
									ctrev	<- ctrev[, {
												ctrp[, CNS_PR_CALL:=CNS_FRQ-THR*CNS_FRQ_STD]
												set(ctrp, NULL, 'CNS_PR_CALL', ctrp[,model.150816a.predict(CNS_PR_CALL, CNS_GPS, BETA0, BETA1, BETA2)])										
												tmp	<- ctrp[CHUNK,][,table(ANS_CALL, factor(PR_CALL>=CNS_PR_CALL, levels=c('TRUE','FALSE'), labels=c('TRUE','FALSE')))]
												list(TP= tmp['1','TRUE'], FP= tmp['0','TRUE'], FN= tmp['1','FALSE'], TN= tmp['0','FALSE'] )					
											}, by=c('CHUNK','THR')]
									ctrev	<- merge(ctrev, ctrev[, list(SENS= TP/(TP+FN), SPEC=TN/(TN+FP), FDR=FP/(FP+TP), FOR=FN/(FN+TN)), by=c('CHUNK','THR')], by=c('CHUNK','THR'))								
								}
								#	plot Sensitivity & Specificity
								ggplot(melt(ctrev, id.vars=c('CHUNK','THR'), measure.vars=c('SENS','SPEC','FDR','FOR')), aes(x=THR, y=100*value, colour=CHUNK)) +
										scale_x_reverse() +
										geom_line() + labs(x='standard deviations of consensus call prob\nPR_CALL_CNT >= PR_CALL_CNS-x*std(PR_CALL_CNS)',y='%', colour='site\n(base relative to HXB2)') +
										theme_bw() + theme(legend.position='bottom') + facet_wrap(~variable, ncol=2, scales='free') +									
										guides(col = guide_legend(ncol=5, byrow=TRUE)) 
								ggsave(file=paste(outdir,'/model.150816a_SensSpec_SITE',site,'.pdf',sep=''), w=9, h=9)
								#
								ctrev							
							}))		
			set(ctrev, NULL, 'CHUNK', ctrev[, as.numeric(CHUNK)])			
			save(ctrmc, ctrev, model.150816a.predict, file=outfile)
		}
	}
	list(coef=ctrmc, predict=model.150816a.predict)
}
##--------------------------------------------------------------------------------------------------------
##	deal with large expressions by performing them subsequently
##--------------------------------------------------------------------------------------------------------
haircut.large.regexpr<- function(pattern, x, n=500, ...)
{
	tmp		<- seq_len(ceiling( nchar(pattern)/500 ))*500
	tmp		<- matrix( c(1, tmp[-length(tmp)]+1, tmp[-length(tmp)], nchar(pattern), rep(NA, 2*length(tmp))), byrow=T, nrow=4  )
	for(j in seq_len(ncol(tmp)))
	{
		z			<- regexpr(substr(pattern, tmp[1,j], tmp[2,j]), x, ...)
		tmp[4,j]	<- attr(z,'match.length')
		tmp[3,j]	<- as.numeric(z)
	}
	tmp		<- tmp[, as.logical(cummin( tmp[4,]>0 )), drop=FALSE]
	z		<- all( tmp[3, ]+tmp[4, ]-1==tmp[2, ] )		#check if complete match
	if(ncol(tmp) & z)
	{
		ans							<- tmp[3,1]
		attr(ans,'match.length')	<- sum(tmp[4,])
	}
	if(!ncol(tmp) | !z)
	{
		ans							<- -1
		attr(ans,'match.length')	<- -1
	}
	ans	
}