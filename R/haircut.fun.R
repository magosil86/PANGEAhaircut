


##--------------------------------------------------------------------------------------------------------
##	wrapper to call 'haircutwrap.get.call.for.PNG_ID'
##--------------------------------------------------------------------------------------------------------
#' @title Call 10bp chunks of cut/raw contigs with the same PANGEA_ID
#' @import data.table zoo plyr ape reshape2 ggplot2 
#' @export
#' @example example/ex.get.call.for.PNG_ID.R
haircutwrap.get.call.for.PNG_ID<- function(indir.st,indir.al,outdir,ctrmc,predict.fun,par,ctrain=NULL,batch.n=NA,batch.id=NA)
{
	infiles	<- data.table(INFILE=list.files(indir.st, pattern='\\.R$', recursive=T))
	infiles[, PNG_ID:= gsub('_wRefs.*','',INFILE)]
	#infiles[, BLASTnCUT:= regmatches(INFILE,regexpr('cut|raw',INFILE))]
	#set(infiles, NULL, 'BLASTnCUT', infiles[, factor(BLASTnCUT, levels=c('cut','raw'), labels=c('Y','N'))])
	alfiles <- data.table(ALFILE=list.files(indir.al, pattern='\\.fasta$', recursive=T))
	alfiles[, PNG_ID:= gsub('_wRefs.*','',ALFILE)]
	#alfiles[, BLASTnCUT:= regmatches(basename(ALFILE),regexpr('cut|raw',basename(ALFILE)))]
	#set(alfiles, NULL, 'BLASTnCUT', alfiles[, factor(BLASTnCUT, levels=c('cut','raw'), labels=c('Y','N'))])
	infiles	<- merge(infiles, alfiles, by='PNG_ID')
	if(!is.na(batch.n) & !is.na(batch.id))
	{
		infiles[, BATCH:= ceiling(seq_len(nrow(infiles))/batch.n)]
		infiles		<- subset(infiles, BATCH==batch.id)
	}
	#	predict by PANGEA_ID
	cnsc.info	<-  infiles[,
			{
				cat(paste('\nProcess', PNG_ID))
				with.warning	<- 0
				#	catch any warnings that indicate that automatic calling may have failed
				#	if this happens, set confidence scores to 0
				tryCatch(
						{
							if(0)	#devel
							{
								PNG_ID<- png_id	<- '15172_1_32'
								#PNG_ID<- png_id	<- '12559_1_11'
								#PNG_ID<- png_id	<- '14728_1_84'
								#PNG_ID<- png_id	<- '14938_1_10'
								#PNG_ID<- png_id	<- '14728_1_82'
								#PNG_ID<- png_id	<- '12559_1_24'
								#PNG_ID<- png_id	<- '12559_1_81'
								#PNG_ID<- png_id	<- '12559_1_87'
								#PNG_ID<- png_id	<- '12559_1_13'
								#PNG_ID<- png_id	<- '13549_1_74'
								#PNG_ID<- png_id	<- '13554_1_12'
								#PNG_ID<- png_id	<- '13554_1_14'
								#PNG_ID<- png_id	<- '13554_1_27'
								#PNG_ID<- png_id	<- '13554_1_33'
								#PNG_ID<- png_id	<- '14760_1_1'
								#PNG_ID<- png_id	<- '15034_1_75'
								#PNG_ID<- png_id	<- '14944_1_17'
								#PNG_ID<- png_id	<- '15065_1_24'
								files	<- subset(infiles, PNG_ID==png_id)[, INFILE]
								alfiles	<- subset(infiles, PNG_ID==png_id)[, ALFILE]								
								tmp		<- haircut.get.call.for.PNG_ID(indir.st, indir.al, png_id, files, alfiles, bc, par, ctrmc, predict.fun)
							}
							if(1)
								tmp		<- haircut.get.call.for.PNG_ID(indir.st, indir.al, PNG_ID, INFILE, ALFILE, par, ctrmc, predict.fun)							
						}, 
						warning=function(w)
						{ 
							tmp				<<- tmp 
							with.warning	<<- 1							
						})
				crs		<- tmp$crs
				cnsc.df	<- tmp$cnsc.df	
				#	handle output
				if(any(grepl(PNG_ID,rownames(crs[['N']]))))
				{
					tmp		<- paste(outdir,'/',PNG_ID,'_wref_nohaironraw.fasta',sep='')
					cat('\nWrite to file', tmp)
					write.dna(crs[['N']], file=tmp, format='fasta', colsep='', nbcol=-1)							
				}
				if(any(grepl(PNG_ID,rownames(crs[['Y']]))))
				{
					tmp		<- paste(outdir,'/',PNG_ID,'_wref_nohaironcut.fasta',sep='')
					cat('\nWrite to file', tmp)
					write.dna(crs[['Y']], file=tmp, format='fasta', colsep='', nbcol=-1)							
				}
				#	see if there is curated contig available
				if(!is.null(ctrain))
				{
					cnsc.df	<- merge(cnsc.df, subset(ctrain, select=c(PNG_ID, TAXON, BLASTnCUT, ANS_FIRST, ANS_LAST)), all.x=T, by=c('PNG_ID','TAXON','BLASTnCUT'))		
					cnsc.df[, CUR_CALL:=NA_integer_]
					set(cnsc.df, cnsc.df[, which(!is.na(ANS_FIRST) & SITE>=ANS_FIRST & SITE<=ANS_LAST)], 'CUR_CALL', 1L)
					set(cnsc.df, cnsc.df[, which(!is.na(ANS_FIRST) & (SITE<ANS_FIRST | SITE>ANS_LAST))], 'CUR_CALL', 0L)
					set(cnsc.df, cnsc.df[, which(is.na(CUR_CALL))], 'CUR_CALL', 0L)								
				}	
				if(is.null(ctrain))
					cnsc.df[, CUR_CALL:=NA_integer_]
				#	save as R
				tmp	<- paste(outdir, '/', PNG_ID, '_wref_nohaironcutraw.R',sep='')
				cat('\nSave to file', tmp)
				save(cnsc.df, crs, file=tmp)
				#	plot
				cnsc.df[, TAXONnCUT:= paste(TAXON,BLASTnCUT,sep='_BLASTnCUT:')]
				ggplot(cnsc.df, aes(x=SITE, fill=BLASTnCUT, group=TAXONnCUT)) +
						geom_ribbon(aes(ymax=CALL, ymin=0), alpha=0.5) +
						geom_line(aes(y=PR_CALL), colour='black') +
						geom_line(aes(y=CNS_PR_CALL), colour='blue') +
						geom_line(aes(y=CUR_CALL), colour='red') +
						scale_x_continuous(breaks=seq(0,15e3, ifelse(cnsc.df[,max(SITE)]>5e2, 5e2, floor(cnsc.df[,max(SITE)/3])))) + 
						facet_wrap(~TAXONnCUT, ncol=1) + theme_bw() + theme(legend.position='bottom') +
						labs(fill='Contig BLASTnCUT', x='position on consensus w/o LTR', y='fill: predicted call\nblack line: predictive probability\nblue line: threshold\nred line: curated call')
				tmp		<- paste(outdir, '/', PNG_ID, '_wref_nohaironcutraw.pdf',sep='')
				cat('\nPlot to file', tmp)
				ggsave(w=10, h=3*cnsc.df[, length(unique(TAXON))], file=tmp)
				#	report confidence score
				tmp		<- subset(cnsc.df, CALL==1)[, list(QUANTILE=c(0,0.01,0.05,0.1,0.2,0.5), PR_CALL=quantile(PR_CALL, p=c(0,0.01,0.05,0.1,0.2,0.5))), by=c('TAXON','BLASTnCUT')]
				if(with.warning)
					set(tmp, NULL, 'PR_CALL', 0)
				tmp
			}, by='PNG_ID']
	#	write quantiles of PR_CALL to file
	if(is.na(batch.n) || is.na(batch.id))
		file		<- paste(outdir, '/model150816a_QUANTILESofPRCALLbyCONTIG.csv',sep='')
	if(!is.na(batch.n) & !is.na(batch.id))
		file		<- paste(outdir, '/model150816a_QUANTILESofPRCALLbyCONTIG_batchn',batch.n,'_batchid',batch.id,'.csv',sep='')
	cnsc.info	<- dcast.data.table(cnsc.info, PNG_ID+TAXON+BLASTnCUT~QUANTILE, value.var='PR_CALL')
	write.csv(cnsc.info, row.names=FALSE, file=file)	
}	


##--------------------------------------------------------------------------------------------------------
##	predict Calls for contigs with same PANGEA ID, based on fitted model
##	update: 
##		1)	do not return duplicate contigs (ie cut and raw, if both are to be kept)
##		2)	do not return raw contigs if cut exists and if raw extends into LTR
##--------------------------------------------------------------------------------------------------------
haircut.get.call.for.PNG_ID<- function(indir.st, indir.al, png_id, files, alfiles, par, ctrmc, predict.fun)	
{
	#	load covariates
	stopifnot(length(files)==1, length(alfiles)==1)
	load( paste(indir.st, '/', files, sep='') )
	tmp		<- subset(cnsc.df, TAXON=='consensus', c(SITE, FRQ, FRQ_STD, GPS))
	setnames(tmp, c('FRQ','FRQ_STD','GPS'), c('CNS_FRQ','CNS_FRQ_STD','CNS_GPS'))
	cnsc.df	<- merge(subset(cnsc.df, TAXON!='consensus'), tmp, by='SITE')					
	#	load alignment	and cut LTR and anything that extends past references
	cr		<- read.dna(file=paste(indir.al,alfiles,sep='/'), format='fasta')
	cr		<- cr[, seq.int(haircut.find.nonLTRstart(cr), ncol(cr))]
	cr		<- cr[, seq.int(1, haircut.find.lastRefSite(cr))]
	#	get contig table
	tmp		<- rownames(cr)[grepl(png_id,rownames(cr))]						
	tx		<- data.table(		TAXON=tmp, 
								BLASTnCUT= factor(grepl('cut',tmp), levels=c(TRUE,FALSE),labels=c('Y','N')), 
								FIRST= apply( as.character(cr[tmp,,drop=FALSE]), 1, function(x) which(x!='-')[1] ),
								LAST= ncol(cr)-apply( as.character(cr[tmp,,drop=FALSE]), 1, function(x) which(rev(x)!='-')[1] ) + 1L
								)	
	tx		<- subset(tx, !is.na(FIRST) & !is.na(LAST))	#some contigs may just be in LTR
	tx[, CNTG:=tx[, gsub('_cut','',gsub(paste(png_id,'.',sep=''),'',substring(TAXON, regexpr(png_id, TAXON))))]]
	tx[, OCNTG:= tx[, sapply(strsplit(CNTG,'.',fixed=T),'[[',1)]]
	tx[, CCNTG:= NA_character_]		
	tmp		<- tx[, which(grepl('.',CNTG,fixed=T))]
	if(length(tmp))
		set(tx, tmp, 'CCNTG', tx[tmp, sapply(strsplit(CNTG,'.',fixed=T),'[[',2)])	
	tmp		<- subset(tx, BLASTnCUT=='Y' & !is.na(CCNTG))[, list(CCNTGn=length(CCNTG)), by='OCNTG']
	tmp		<- subset(tmp, CCNTGn==1)[, OCNTG]	#check for cut contigs that should be present in multiple cuts by naming scheme, but after LTR removal there is only one cut
	if(length(tmp))
	{
		cat('\nFound lone cuts for which multiple cuts are expected by naming scheme, n=', length(tmp))
		tmp	<- tx[, which(BLASTnCUT=='Y' & OCNTG%in%tmp)]
		set(tx, tmp, 'CCNTG', NA_character_)
		set(tx, tmp, 'CNTG', tx[tmp, OCNTG])
	}
	#	calculate PR_CALL of contig and of consensus
	tmp		<- seq(cnsc.df[, floor(min(SITE)/10)*10-10],cnsc.df[, max(SITE)+10],10)	
	cnsc.df[, CHUNK:=cut(SITE, breaks=tmp, labels=tmp[-length(tmp)])]
	cnsc.df	<- merge(cnsc.df, ctrmc, by='CHUNK', all.x=TRUE)
	if(cnsc.df[, any(is.na(BETA0))])
	{
		warning('Found NA BETA0 for PNG_ID ',PNG_ID,': likely because the cut/raw contig alignment is much longer than expected. Suggests that the site-specific call probability model may not match to the sites in the alignment. Check manually. ')
		tmp	<- cnsc.df[, which(SITE==subset(cnsc.df, !is.na(BETA0))[, max(SITE)])[1]]
		tmp2<- cnsc.df[, which(is.na(BETA0))]
		set(cnsc.df, tmp2, 'BETA0', cnsc.df[tmp, BETA0])
		set(cnsc.df, tmp2, 'BETA1', cnsc.df[tmp, BETA1])
		set(cnsc.df, tmp2, 'BETA2', cnsc.df[tmp, BETA2])
	}
	cnsc.df[, PR_CALL:= predict.fun(FRQ, GPS, BETA0, BETA1, BETA2)]
	cnsc.df[, CNS_PR_CALL:= CNS_FRQ-par['PRCALL.thrstd']*CNS_FRQ_STD]
	set(cnsc.df, cnsc.df[, which(CNS_PR_CALL<0)], 'CNS_PR_CALL', 0)
	set(cnsc.df, NULL, 'CNS_PR_CALL', cnsc.df[,predict.fun(CNS_PR_CALL, CNS_GPS, BETA0, BETA1, BETA2)])
	set(cnsc.df, cnsc.df[, which(CNS_PR_CALL>par['PRCALL.thrmax'])], 'CNS_PR_CALL', par['PRCALL.thrmax']) 
	#	call if call prob of contig is not too low in relation to the call prob of the consensus 
	cnsc.df[, CALL:= as.integer(PR_CALL>=CNS_PR_CALL)]	
	if(0)
	{
		ggplot(subset(cnsc.df, BLASTnCUT=='Y'), aes(x=SITE)) + facet_wrap(~TAXON, ncol=1) +
				geom_line(aes(y=PR_CALL), colour='black') +
				geom_line(aes(y=CNS_FRQ), colour='red') +
				geom_line(aes(y=CNS_PR_CALL), colour='blue') +
				geom_line(aes(y=CNS_GPS), colour='orange') +
				geom_line(aes(y=FRQ), colour='green') +
				geom_line(aes(y=GPS), colour='DarkGreen')
	} 	
	#	determine predicted sites
	cnsc.1s	<- cnsc.df[, {
				z	<- gsub('0*$','',paste(CALL,collapse=''))
				#print(TAXON)
				#print(z)
				z	<- gregexpr('1+',z)[[1]]	
				list(CALL_ID= seq_along(z), CALL_POS=as.integer(z+min(SITE)-1L), CALL_LEN=attr(z, 'match.length'))
			}, by=c('PNG_ID','TAXON','BLASTnCUT')]
	cnsc.1s[, CALL_LAST:= CALL_POS+CALL_LEN-1L]
	cnsc.1s	<- subset(cnsc.1s, CALL_LEN>0)
	#	fill internal predicted gaps		
	if(!nrow(cnsc.1s))
		warning('Could not determine any calls for either the raw or cut contigs with PNG_ID ',PNG_ID,'. Check manually.')
	if((!is.na(par['PRCALL.rmintrnlgpsblw'] | !is.na(par['PRCALL.rmintrnlgpsend']))) && nrow(cnsc.1s))
	{
		cnsc.g	<- cnsc.1s[, {
					if(length(CALL_ID)==1)
						ans	<- NA_integer_
					if(length(CALL_ID)>1)
						ans	<- CALL_POS[seq.int(2,length(CALL_POS))]-CALL_LAST[seq.int(1,length(CALL_LAST)-1)]-1L
					list(CALL_LAST=CALL_LAST[seq.int(1,length(CALL_LAST)-1)], GAP_LEN= ans)
				}, by=c('TAXON','BLASTnCUT')]
		cnsc.g	<- merge(cnsc.1s, cnsc.g,  by=c('TAXON','BLASTnCUT','CALL_LAST'), all.x=1)
		if(!is.na(par['PRCALL.rmintrnlgpsblw']))
			tmp	<- cnsc.g[, which(GAP_LEN<par['PRCALL.rmintrnlgpsblw'])]
		if(!is.na(par['PRCALL.rmintrnlgpsblw']))
			tmp	<- union(tmp, cnsc.g[, which(CALL_LAST>par['PRCALL.rmintrnlgpsend'] & !is.na(GAP_LEN))])		
		for(i in tmp)	#	add ith called region to next call region
		{
			tmp2	<- cnsc.df[, which(TAXON==cnsc.g$TAXON[i] & BLASTnCUT==cnsc.g$BLASTnCUT[i] & SITE>cnsc.g$CALL_LAST[i] & SITE<=cnsc.g$CALL_LAST[i]+cnsc.g$GAP_LEN[i])]
			stopifnot( cnsc.df[tmp2,all(CALL==0)])
			cat('\nFound internal predicted gap and set to CALL=1, taxon=',cnsc.g[i,TAXON],', cut=',cnsc.g[i,as.character(BLASTnCUT)],', pos=',cnsc.g[i,CALL_LAST+1L],', len=', length(tmp2))
			set(cnsc.df, tmp2, 'CALL', 1L)
			set(cnsc.g, i+1L, 'CALL_POS', cnsc.g[i,CALL_POS])
			set(cnsc.g, i+1L, 'CALL_LEN', cnsc.g[i+1L,CALL_LEN]+cnsc.g[i,CALL_LEN]+cnsc.g[i,GAP_LEN])
			set(cnsc.g, i, 'CALL_ID', NA_integer_)
		}
		cnsc.1s		<- subset(cnsc.g, !is.na(CALL_ID))
	}
	#	check if called contig has gaps of CALL=='0': if yes, return last non-gap before first CALL=='0
	if(!is.na(par['PRCALL.cutprdcthair']) && nrow(cnsc.1s))
	{
		setkey(cnsc.1s, TAXON, BLASTnCUT, CALL_POS)
		tmp		<- cnsc.1s[, which(CALL_LEN<par['PRCALL.cutprdcthair'])]	
		for(i in tmp)	#keep raw
		{
			if( (i-1)>0	  &&  cnsc.1s[i-1,TAXON]==cnsc.1s[i,TAXON] &&  cnsc.1s[i-1,BLASTnCUT]==cnsc.1s[i,BLASTnCUT]  &&  cnsc.1s[i-1,GAP_LEN]>2*par['PRCALL.cutprdcthair'])
			{
				cat('\nFound predicted extra hair of length <',par['PRCALL.cutprdcthair'],'delete, n=',cnsc.1s$CALL_LEN[i])
				set(cnsc.df, cnsc.df[, which(TAXON==cnsc.1s$TAXON[i] & BLASTnCUT==cnsc.1s$BLASTnCUT[i] & SITE>=cnsc.1s$CALL_POS[i] & SITE<=cnsc.1s$CALL_LAST[i])], 'CALL', 0L)
				set(cnsc.1s, i, 'CALL_ID', NA_integer_)
				set(cnsc.1s, i, 'GAP_LEN', cnsc.1s[i,GAP_LEN]+cnsc.1s[i-1,GAP_LEN]+cnsc.1s[i,CALL_LEN])
			}				
		}
		cnsc.1s	<- subset(cnsc.1s, !is.na(CALL_ID))								
	}
	#	check if all called chunks in cut and raw contigs correspond to each other
	if(!is.na(par['PRCALL.cutrawgrace']) && nrow(cnsc.1s))
	{
		cnsc.1s	<- merge(cnsc.1s, tx, by=c('TAXON','BLASTnCUT'))
		tmp		<- subset(cnsc.1s, BLASTnCUT=='Y', select=c(OCNTG, TAXON, CALL_ID, CALL_POS, CALL_LEN ))
		setnames(tmp, c('TAXON','CALL_ID','CALL_POS','CALL_LEN'),c('TAXON_CUT','CALL_ID_CUT','CALL_POS_CUT','CALL_LEN_CUT'))
		tmp		<- merge(subset(cnsc.1s, BLASTnCUT=='N'), tmp, by='OCNTG', allow.cartesian=TRUE)
		tmp		<- subset(tmp, abs(CALL_POS-CALL_POS_CUT)<par['PRCALL.cutrawgrace'] )	
		#	of the corresponding calls, keep the longer one 
		#	dont extend shorter for now as alignments may not match
		tmp2	<- subset(tmp, CALL_LEN>CALL_LEN_CUT)
		for(i in seq_len(nrow(tmp2)))	#keep raw
		{			
			cat('\nkeep only raw:', tmp2[i,TAXON])
			z	<- cnsc.df[, which( TAXON==tmp2$TAXON_CUT[i] & BLASTnCUT=='Y' & SITE>=tmp2$CALL_POS_CUT[i] & SITE<(tmp2$CALL_POS_CUT[i]+tmp2$CALL_LEN_CUT[i]))]
			stopifnot( cnsc.df[z,all(CALL==1)])
			set(cnsc.df, z, 'CALL', 0L)
			set(cnsc.1s, cnsc.1s[, which(TAXON==tmp2$TAXON_CUT[i] & BLASTnCUT=='Y' & CALL_ID==tmp2$CALL_ID_CUT[i])],'CALL_ID',NA_integer_)
		}
		tmp2	<- subset(tmp, CALL_LEN<=CALL_LEN_CUT)
		for(i in seq_len(nrow(tmp2)))	#keep cut
		{	
			cat('\nkeep only cut:', tmp2[i,TAXON_CUT])
			z	<- cnsc.df[, which( TAXON==tmp2$TAXON[i] & BLASTnCUT=='N' & SITE>=tmp2$CALL_POS[i] & SITE<=tmp2$CALL_LAST[i])]
			stopifnot( cnsc.df[z,all(CALL==1)])
			set(cnsc.df, z, 'CALL', 0L)
			set(cnsc.1s, cnsc.1s[, which(TAXON==tmp2$TAXON[i] & BLASTnCUT=='N' & CALL_ID==tmp2$CALL_ID[i])],'CALL_ID',NA_integer_)
		}
		cnsc.1s	<- subset(cnsc.1s, !is.na(CALL_ID))		
	}
	#	check that remaining contigs have sufficient length
	if(!is.na(par['PRCALL.cutprdctcntg']) && nrow(cnsc.1s))
	{		
		tmp		<- cnsc.1s[, which(CALL_LEN<par['PRCALL.cutprdctcntg'])]	
		set(cnsc.1s, tmp, 'CALL_ID', NA_integer_)
		for(i in tmp)	#keep raw
		{
			cat('\nFound predicted contig of length <',par['PRCALL.cutprdctcntg'],'delete, n=',cnsc.1s$CALL_LEN[i])
			set(cnsc.df, cnsc.df[, which( TAXON==cnsc.1s$TAXON[i] & BLASTnCUT==cnsc.1s$BLASTnCUT[i] & SITE>=cnsc.1s$CALL_POS[i] & SITE<=cnsc.1s$CALL_LAST[i])], 'CALL', 0L)
		}										
		cnsc.1s	<- subset(cnsc.1s, !is.na(CALL_ID))								
	}
	#	check there is no dust
	tmp			<- subset(cnsc.df, CALL==1)[, list(CALL_N= length(CALL)), by=c('TAXON','BLASTnCUT')][, CALL_N]
	if(length(tmp))
		stopifnot(min(tmp)>40)
	#	produce fasta output:
	#	select cut and raw contigs with a call, then set all characters with CALL==0 to -
	crs			<- lapply(crs, as.character)
	tmp			<- subset(cnsc.df, BLASTnCUT=='N' & CALL==1 )[, unique(TAXON)]
	tmp2		<- rownames(crs[['N']])[ !grepl(png_id,rownames(crs[['N']])) | rownames(crs[['N']])%in%tmp ] 
	crs[['N']]	<- crs[['N']][tmp2,]
	for(tx in tmp)
		crs[['N']][tx, subset(cnsc.df, BLASTnCUT=='N' & TAXON==tx & CALL==0 & SITE<=ncol(crs[['N']]))[, SITE]]	<- '-'
	
	tmp			<- subset(cnsc.df, BLASTnCUT=='Y' & CALL==1 )[, unique(TAXON)]
	tmp2		<- rownames(crs[['Y']])[ !grepl(png_id,rownames(crs[['Y']])) | rownames(crs[['Y']])%in%tmp ] 
	crs[['Y']]	<- crs[['Y']][tmp2,]
	for(tx in tmp)
		crs[['Y']][tx, subset(cnsc.df, BLASTnCUT=='Y' & TAXON==tx & CALL==0 & SITE<=ncol(crs[['Y']]))[, SITE]]	<- '-'
	crs			<- lapply(crs, as.DNAbin)			
	list(crs=crs, cnsc.df=cnsc.df)
}
##--------------------------------------------------------------------------------------------------------
##	process all files in indir with 'haircut.get.cut.statistics'
##--------------------------------------------------------------------------------------------------------
#' @title Compute descriptive statistics that are used to calculate the call probability
#' @import data.table zoo plyr ape reshape2 ggplot2
#' @export
#' @example example/ex.get.cut.statistics.R
haircutwrap.get.cut.statistics<- function(indir, par, outdir=indir, batch.n=NA, batch.id=NA)	
{
	require(zoo)
	#	read just one file
	#	determine statistics for each contig after LTR
	infiles		<- data.table(FILE=list.files(indir, pattern='fasta$', recursive=T))
	infiles[, PNG_ID:= gsub('_wRefs\\.fasta','',gsub('_cut|_raw','',FILE))]
	if(!is.na(batch.n) & !is.na(batch.id))
	{
		infiles[, BATCH:= ceiling(seq_len(nrow(infiles))/batch.n)]
		infiles		<- subset(infiles, BATCH==batch.id)
	}
	
	#	check which contigs not yet processed
	infiles		<- merge(infiles, infiles[, {
						file	<- paste(indir, FILE, sep='/')
						tmp		<- paste(outdir, '/', gsub('\\.fasta',paste('_HAIRCUTSTAT_thr',100*par['FRQx.quantile'],'_aw',par['CNS_AGR.window'],'_fw',par['CNS_FRQ.window'],'_gw',par['GPS.window'],'.R',sep=''),basename(file)), sep='')
						options(show.error.messages = FALSE)		
						readAttempt		<-try(suppressWarnings(load(tmp)))
						options(show.error.messages = TRUE)	
						list(	DONE=!inherits(readAttempt, "try-error")	)			
					}, by='FILE'], by='FILE')
	cat(paste('\nFound processed files, n=', infiles[, length(which(DONE))]))
	infiles		<- subset(infiles, !DONE)
	#
	#	infiles[, which(grepl('15172_1_32',FILE))]	fls<- 1431
	#	process files
	for(fls in infiles[, seq_along(FILE)])
	{
		file	<- paste(indir, infiles[fls, FILE], sep='/')
		cat(paste('\nProcess', file))
		#	read Contigs+Rrefs: cr
		cr		<- read.dna(file, format='fasta')
		if(!is.matrix(cr) || nrow(cr)==0)
			warning('Found unexpected file format for file ', file)
		if(is.matrix(cr) && nrow(cr)>0)
		{
			#	determine start of non-LTR position and cut 
			cr		<- cr[, seq.int(haircut.find.nonLTRstart(cr), ncol(cr))]
			#	cut at last site of references
			cr		<- cr[, seq.int(1, haircut.find.lastRefSite(cr))]		
			#	determine reference sequences. 
			#	non-refs have the first part of the file name in their contig name and are at the top of the alignment
			tmp		<- strsplit(basename(file), '_')[[1]][1]
			tx		<- data.table(TAXON= rownames(cr), CONTIG=as.integer(grepl(tmp, rownames(cr))) )
			tx		<- tx[order(-CONTIG, na.last=TRUE)]
			cr		<- cr[tx[,TAXON],]			
			cat(paste('\nFound contigs, n=', tx[, length(which(CONTIG==1))]))
			#	determine base frequencies at each site amongst references.
			tmp		<- cr[subset(tx, CONTIG==0)[, TAXON],]
			rp		<- haircut.get.frequencies(tmp, bases=c('a','c','g','t','-') )
			tmp		<- haircut.get.consensus.from.frequencies(rp, par)
			cnsr	<- tmp$DNAbin
			cnsr.df	<- tmp$DATATABLE
			#	for each contig, determine %agreement with consensus on rolling window
			cnsc	<- rbind(cnsr, cr[subset(tx, CONTIG==1)[, TAXON],])
			#	determine first and last non-gap sites
			tmp		<- as.character(cnsc)
			tx		<- data.table(	TAXON= rownames(cnsc), 
					FIRST= apply( tmp, 1, function(x) which(x!='-')[1] ),
					LAST= ncol(cnsc)-apply( tmp, 1, function(x) which(rev(x)!='-')[1] )+1L		)
			tx		<- subset(tx, !is.na(FIRST) & !is.na(LAST))	#	some contigs only map into LTR
			#	determine all cut statistics
			cnsc.df	<- haircut.get.cut.statistics(cnsc, rp, tx, par)		
			#ggplot(cnsc.df, aes(x=SITE)) +facet_wrap(~TAXON, ncol=1) + geom_line(aes(y=FRQ), colour='black') + geom_line(aes(y=AGRpc), colour='blue') + geom_line(aes(y=GPS), colour='red') + geom_line(aes(y=FRQ-2*FRQ_STD), colour='DarkGreen')
			cnsc.df[, PNG_ID:= infiles[fls, PNG_ID]]
			cnsc.df[, BLASTnCUT:= cnsc.df[, factor(grepl('cut',TAXON),levels=c(TRUE,FALSE),labels=c('Y','N'))]]
			cat(paste('\nSave contigs, n=', cnsc.df[, length(unique(TAXON))]))
			#	save
			file	<- paste(outdir, '/', gsub('\\.fasta',paste('_HAIRCUTSTAT_thr',100*par['FRQx.quantile'],'_aw',par['CNS_AGR.window'],'_fw',par['CNS_FRQ.window'],'_gw',par['GPS.window'],'.R',sep=''),basename(file)), sep='')
			cat(paste('\nSave to', file))
			save(cnsc, cnsc.df, file=file)	
		}				
	}
}
##--------------------------------------------------------------------------------------------------------
##	get cut statistics:
##		- FRQ, FRQ_STD, AGRpc, GPS
##--------------------------------------------------------------------------------------------------------
haircut.get.cut.statistics<- function(cnsc, rp,  tx, par)
{
	require(zoo)
	cnsc.df	<- tx[,{					
				tmp		<- cnsc[ c('consensus',TAXON), seq.int(FIRST,LAST)]
				agrpc	<- 1-rollapply( seq_len(LAST-FIRST+1), width=par['CNS_AGR.window'], FUN= function(z) dist.dna(tmp[,z], model='indel' )/length(z), align='center', partial=T )					
				tmp		<- data.table(	BASE=as.vector(as.character(cnsc[TAXON, seq.int(FIRST,LAST)])), SITE=seq.int(FIRST,LAST))
				tmp		<- merge(rp,tmp,by=c('SITE','BASE'))
				set(tmp, NULL, 'FRQ_STD', tmp[, sqrt(FRQ*(1-FRQ)/COV)])
				freqr	<- rollapply( seq_len(LAST-FIRST+1), width=par['CNS_AGR.window'], FUN= function(z) mean(tmp$FRQ[z]), align='center', partial=T )
				freqsr	<- rollapply( seq_len(LAST-FIRST+1), width=par['CNS_AGR.window'], FUN= function(z) mean(tmp$FRQ_STD[z]), align='center', partial=T )
				tmp		<- as.character( cnsc[ TAXON, seq.int(FIRST,LAST)] )=='-'
				gps		<- rollapply( seq_len(LAST-FIRST+1), width=par['GPS.window'], FUN= function(z) mean(tmp[z]), align='center', partial=T )
				list( SITE=seq.int(FIRST,LAST),  FRQ=freqr, FRQ_STD=freqsr, AGRpc=agrpc, GPS=gps   )
			},by='TAXON']	
	cnsc.df
}

##--------------------------------------------------------------------------------------------------------
##	load training data set for sites 
##--------------------------------------------------------------------------------------------------------
haircut.load.training.data<- function(indir, site)
{
	tmp		<- cut(site, breaks=c(-1,seq.int(200, 10200, 200),Inf), labels=c(paste('<',seq.int(200, 10200, 200),sep=''),'>10200'))
	tmp		<- gsub('<|>','',tmp)
	tmp		<- list.files(indir, pattern=paste('SITES',tmp,'\\.R',sep=''), full.names=T)
	stopifnot(length(tmp)==1)
	load(tmp)
	ctr
}
##--------------------------------------------------------------------------------------------------------
##	calculate offset of the first sequence to the second sequence, assuming that the only difference are gap characters
##--------------------------------------------------------------------------------------------------------
haircut.calculate.offset<- function(x,y)
{
	x		<- as.character(x)	
	y		<- as.character(y)		
	stopifnot( gsub('-','',paste(as.vector(x), collapse=''))==gsub('-','',paste(as.vector(y), collapse='')) )
	z		<- rbind.fill.matrix(x,y)
	offset	<- rep(0, ncol(z))			
	k		<- seq_len(ncol(z))+offset
	k		<- k[ k>0 & k<=ncol(z)]
	k		<- which( z[1, seq_along(k)]!=z[2, k] )[1]
	while(!is.na(k))
	{
		if( z[1,k]!='-' )
			offset[ seq.int(k,length(offset)) ] <- offset[ seq.int(k,length(offset)) ]+1
		if( z[1,k]=='-' )
			offset[ seq.int(k,length(offset)) ] <- offset[ seq.int(k,length(offset)) ]-1
		k		<- seq_len(ncol(z))+offset
		k		<- k[ k>0 & k<=ncol(z)]
		k		<- which( z[1, seq_along(k)]!=z[2, k] )[1]		
	}
	if(length(offset)>length(x))
		offset	<- offset[ seq_len(length(x)) ]
	offset
}

##--------------------------------------------------------------------------------------------------------
##	find the first site outside the LTR in the contig + reference alignment 
##--------------------------------------------------------------------------------------------------------
haircut.find.nonLTRstart<- function(cr)
{
	ans				<- seq.find.pos.of.pattern(cr, pattern='a-*t-*g-*g-*g-*t-*g-*c-*g-*a-*g-*a-*g-*c-*g-*t-*c-*a', row.name='B.FR.83.HXB2_LA')
	stopifnot(length(ans)==1, ans>0)
	ans
}
##--------------------------------------------------------------------------------------------------------
##	find the last site inside the reference alignment 
##--------------------------------------------------------------------------------------------------------
haircut.find.lastRefSite<- function(cr)
{
	
	ans				<- seq.find.pos.of.pattern(cr, pattern='t-*t-*t-*t-*a-*g-*t-*c-*a-*g-*t-*g-*t-*g-*g-*a-*a-*a-*a-*t-*c-*t-*c-*t-*a-*g-*c-*a', row.name='B.FR.83.HXB2_LA')
	stopifnot(length(ans)==1, ans>0)
	as.integer(ans+attr(ans,'match.length')-1L)
}
##--------------------------------------------------------------------------------------------------------
##	determine the base frequencies in the reference alignment
##--------------------------------------------------------------------------------------------------------
haircut.get.frequencies	<- function(seq, bases=c('a','c','g','t','-') )
{	
	#	calculate base frequency Profile among References: rp
	rp		<- sapply(seq_len(ncol(seq)), function(i) base.freq(seq[,i], freq=T, all=T))
	rp		<- as.data.table( t(rp[bases,]) )
	rp[, SITE:= seq_len(nrow(rp))]
	setnames(rp, c('a','c','g','t','-'), c('BASEa','BASEc','BASEg','BASEt','GAP'))
	#	calculate first and last positions of each sequence
	tmp		<- as.character(seq)
	tx		<- data.table(	TAXON= rownames(tmp), 
			FIRST= apply( tmp, 1, function(x) which(x!='-')[1] ),
			LAST= ncol(tmp)-apply( tmp, 1, function(x) which(rev(x)!='-')[1] ) + 1L		)
	#	calculate coverage
	tmp		<- data.table(SITE=c(tx[, seq_len(max(FIRST))], tx[, seq.int(min(LAST), ncol(seq))]))
	tmp		<- tmp[, list( COV=tx[, as.numeric(length(which(SITE>=FIRST & SITE<=LAST)))]), by='SITE']
	rp		<- merge(rp, tmp, by='SITE', all.x=TRUE)
	set(rp, rp[, which(is.na(COV))], 'COV', nrow(seq))
	#	correct - count
	set(rp, NULL, 'GAP', rp[, GAP+COV-nrow(seq)]) 
	#	get frequencies
	rp		<- melt(rp, id.vars=c('SITE','COV'), variable.name='BASE', value.name='FRQ')
	set(rp, NULL, 'BASE', rp[, gsub('GAP','-',gsub('BASE','',BASE))])
	set(rp, NULL, 'BASE', rp[, factor(BASE, levels=bases, labels=bases)])
	set(rp, NULL, 'FRQ', rp[, FRQ/COV])
	set(rp, rp[, which(COV==0)], 'FRQ', 0.)
	rp
}	
##--------------------------------------------------------------------------------------------------------
##	calculate consensus from frequency profile
##--------------------------------------------------------------------------------------------------------
haircut.get.consensus.from.frequencies	<- function(rp, par)
{		
	#	get consensus
	crp		<- rp[, {
						z	<- which.max(FRQ)
						list(CNS_BASE= BASE[z], CNS_FRQ=FRQ[z], CNS_COV=COV[z])
					}, by='SITE']
	#	get consensus base above lower quantile. Some consensus bases will be ?
	if(!is.na(par['FRQx.thr']))
		set(crp, crp[, which(CNS_FRQ<par['FRQx.thr'])], 'CNS_BASE', '?')
	#	get DNAbin
	tmp		<- crp[, as.DNAbin(t(as.matrix(CNS_BASE)))]
	rownames(tmp)	<- 'consensus'
	list(DNAbin=tmp, DATATABLE=crp)	
}
##--------------------------------------------------------------------------------------------------------
##	find the first site with pattern in alignment 
##	TODO: put into hivclust
##--------------------------------------------------------------------------------------------------------
seq.find.pos.of.pattern<- function(seq, pattern, row.name)
{
	stopifnot(is.matrix(seq), !is.na(pattern), !is.na(row.name))
	tmp		<- which(grepl(row.name, rownames(seq), fixed=1))
	stopifnot(length(tmp)==1)	
	regexpr(pattern, paste(as.character( seq[tmp, ] ),collapse=''))	
}