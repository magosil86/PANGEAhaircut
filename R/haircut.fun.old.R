##--------------------------------------------------------------------------------------------------------
##	predict Calls for contigs with same PANGEA ID, based on fitted model 
##--------------------------------------------------------------------------------------------------------
haircut.get.call.for.PNG_ID.150811<- function(indir.st, indir.al, png_id, files, alfiles, bc, par, ctrmc, predict.fun)	
{
	#	load covariates
	cnsc.df	<- do.call('rbind',lapply(files, function(x)
					{
						load( paste(indir.st, '/', x, sep='') )
						cnsc.df
					})) 
	#	load alignment
	crs		<- lapply(alfiles, function(x){
				cr	<- read.dna(file=paste(indir.al,x,sep='/'), format='fasta')
				cr[, seq.int(haircut.find.nonLTRstart(cr), ncol(cr))]										
			})
	names(crs)	<- bc
	#	predict
	tmp		<- seq(cnsc.df[, floor(min(SITE)/10)*10],cnsc.df[, max(SITE)+10],10)
	cnsc.df[, CHUNK:=cut(SITE, breaks=tmp, labels=tmp[-length(tmp)])]	
	cnsc.df	<- merge(cnsc.df, ctrmc, by='CHUNK')
	cnsc.df[, PR_CALL:= predict.fun(AGRpc, GPS, BETA0, BETA1, BETA2)]							
	cnsc.df[, CALL:= as.integer(PR_CALL>=par['PRCALL.thr'])]	
	#	check if called contig has gaps of CALL=='0': if yes, return last non-gap before first CALL=='0
	if(!is.na(par['PRCALL.cutprdcthair']))
	{
		tmp		<- cnsc.df[, {
					z	<- gsub('0*$','',paste(CALL,collapse=''))
					#print(TAXON)
					#print(z)
					z	<- gregexpr('1+',z)[[1]]	
					list(CALL_POS=as.integer(z+min(SITE)-1L), CALL_LEN=attr(z, 'match.length'))
				}, by=c('PNG_ID','TAXON','BLASTnCUT')]
		tmp		<- subset( tmp, CALL_LEN<par['PRCALL.cutprdcthair'] )
		if(nrow(tmp))
		{
			cat('\nFound predicted extra hair of length <',par['PRCALL.cutprdcthair'],'delete, n=',tmp[,sum(CALL_LEN)])
			set(tmp, NULL, 'CALL_LEN', tmp[, CALL_POS+CALL_LEN-1])
			for(i in seq_len(nrow(tmp)))
			{				
				set(cnsc.df, cnsc.df[, which(TAXON==tmp$TAXON[i] & BLASTnCUT==tmp$BLASTnCUT[i] & SITE>=tmp$CALL_POS[i] & SITE<=tmp$CALL_LEN[i])], 'CALL', 0L)
			}				
		}										
	}
	#	produce fasta output:
	#	select cut and raw contigs with a call
	crs			<- lapply(crs, as.character)
	tmp			<- subset(cnsc.df, BLASTnCUT=='N' & CALL==1 )[, unique(TAXON)]
	tmp			<- rownames(crs[['N']])[ !grepl(png_id,rownames(crs[['N']])) | rownames(crs[['N']])%in%tmp ] 
	crs[['N']]	<- crs[['N']][tmp,]
	tmp			<- subset(cnsc.df, BLASTnCUT=='Y' & CALL==1 )[, unique(TAXON)]
	tmp			<- rownames(crs[['Y']])[ !grepl(png_id,rownames(crs[['Y']])) | rownames(crs[['Y']])%in%tmp ] 
	crs[['Y']]	<- crs[['Y']][tmp,]
	#	set all characters with CALL==0 to -
	tmp			<- subset(cnsc.df, BLASTnCUT=='N' & CALL==1 )[, unique(TAXON)]
	for(tx in tmp)
	{
		crs[['N']][tx, subset(cnsc.df, BLASTnCUT=='N' & TAXON==tx & CALL==0)[, SITE]]	<- '-'
	}
	tmp			<- subset(cnsc.df, BLASTnCUT=='Y' & CALL==1 )[, unique(TAXON)]
	for(tx in tmp)
	{
		crs[['Y']][tx, subset(cnsc.df, BLASTnCUT=='Y' & TAXON==tx & CALL==0)[, SITE]]	<- '-'
	}
	crs			<- lapply(crs, as.DNAbin)		
	
	list(crs=crs, cnsc.df=cnsc.df)
}
haircut.get.call.for.PNG_ID.150816<- function(indir.st, indir.al, png_id, files, alfiles, bc, par, ctrmc, predict.fun)	
{
	#	load covariates
	cnsc.df	<- do.call('rbind',lapply(files, function(x)
					{
						load( paste(indir.st, '/', x, sep='') )
						tmp	<- subset(cnsc.df, TAXON=='consensus', c(SITE, FRQ, FRQ_STD, GPS))
						setnames(tmp, c('FRQ','FRQ_STD','GPS'), c('CNS_FRQ','CNS_FRQ_STD','CNS_GPS'))
						merge(subset(cnsc.df, TAXON!='consensus'), tmp, by='SITE')
					})) 
	#	load alignment	and cut LTR and anything that extends past references
	crs		<- lapply(alfiles, function(x)
			{
				cr	<- read.dna(file=paste(indir.al,x,sep='/'), format='fasta')
				cr	<- cr[, seq.int(haircut.find.nonLTRstart(cr), ncol(cr))]
				cr[, seq.int(1, haircut.find.lastRefSite(cr))]																
			})
	names(crs)	<- bc
	#
	if( diff(sapply(crs, ncol))>par['PRCALL.cutrawgrace']/2 )
		warning('Found large difference in alignment length for cut/raw contigs with PNG_ID ',PNG_ID,': it may be that identical/subset contigs are not correctly identified. Check manually.')
	#	get contig table
	tx		<- do.call('rbind',lapply(seq_along(crs), function(i)
					{
						tmp		<- rownames(crs[[i]])[grepl(png_id,rownames(crs[[i]]))]						
						data.table(	TAXON=tmp, 
								BLASTnCUT= bc[i], 
								FIRST= apply( as.character(crs[[i]][tmp,,drop=FALSE]), 1, function(x) which(x!='-')[1] ),
								LAST= ncol(crs[[i]])-apply( as.character(crs[[i]][tmp,,drop=FALSE]), 1, function(x) which(rev(x)!='-')[1] ) + 1L,
								CRS_ID=i)
					}))	
	tx		<- subset(tx, !is.na(FIRST) & !is.na(LAST))	#some contigs may just be in LTR
	tx[, CNTG:=tx[, gsub(paste(png_id,'.',sep=''),'',substring(TAXON, regexpr(png_id, TAXON)))]]
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
	if(cnsc.df[, !any(is.na(BETA0))])
	{
		warning('Found NA BETA0 for PNG_ID',PNG_ID,': likely because one cut/raw contig alignment is much longer than expected. Suggests that the site-specific call probability model may not match to the sites in the alignment. Check manually. ')
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
##	predict Calls for contigs with same PANGEA ID, based on fitted model
##	update: 
##		1)	do not return duplicate contigs (ie cut and raw, if both are to be kept)
##		2)	do not return raw contigs if cut exists and if raw extends into LTR
##--------------------------------------------------------------------------------------------------------
haircut.get.call.for.PNG_ID.150814<- function(indir.st, indir.al, png_id, files, alfiles, bc, par, ctrmc, predict.fun)	
{
	#	load covariates
	cnsc.df	<- do.call('rbind',lapply(files, function(x)
					{
						load( paste(indir.st, '/', x, sep='') )
						cnsc.df
					})) 
	#	load alignment	and cut LTR and anything that extends past references
	crs		<- lapply(alfiles, function(x)
			{
				cr	<- read.dna(file=paste(indir.al,x,sep='/'), format='fasta')
				cr	<- cr[, seq.int(haircut.find.nonLTRstart(cr), ncol(cr))]
				cr[, seq.int(1, haircut.find.lastRefSite(cr))]																
			})
	names(crs)	<- bc
	#
	#	get contig table
	tx		<- do.call('rbind',lapply(seq_along(crs), function(i)
					{
						tmp		<- rownames(crs[[i]])[grepl(png_id,rownames(crs[[i]]))]						
						data.table(	TAXON=tmp, 
								BLASTnCUT= factor(bc[i],levels=c('cut','raw'),labels=c('Y','N')), 
								FIRST= apply( as.character(crs[[i]][tmp,,drop=FALSE]), 1, function(x) which(x!='-')[1] ),
								LAST= ncol(crs[[i]])-apply( as.character(crs[[i]][tmp,,drop=FALSE]), 1, function(x) which(rev(x)!='-')[1] ) + 1L,
								CRS_ID=i)
					}))	
	tx		<- subset(tx, !is.na(FIRST) & !is.na(LAST))	#some contigs may just be in LTR
	tx[, CNTG:=tx[, gsub(paste(png_id,'.',sep=''),'',substring(TAXON, regexpr(png_id, TAXON)))]]
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
	#	calculate PR_CALL of contig
	tmp		<- seq(cnsc.df[, floor(min(SITE)/10)*10],cnsc.df[, max(SITE)+10],10)
	cnsc.df[, CHUNK:=cut(SITE, breaks=tmp, labels=tmp[-length(tmp)])]	
	cnsc.df	<- merge(cnsc.df, ctrmc, by='CHUNK')
	cnsc.df[, PR_CALL:= predict.fun(AGRpc, GPS, BETA0, BETA1, BETA2)]
	#	calculate PR_CALL of consensus
	cnsc.df[, CNS_PR_CALL:= predict.fun(CNS_FRQr, CNS_GPSr, BETA0, BETA1, BETA2)]
	
	ggplot(subset(cnsc.df, BLASTnCUT=='Y'), aes(x=SITE)) + facet_wrap(~TAXON, ncol=1) +
			geom_line(aes(y=PR_CALL), colour='black') +
			geom_line(aes(y=CNS_PR_CALL), colour='blue') +
			geom_line(aes(y=CNS_FRQr), colour='red') +
			geom_line(aes(y=CNS_GPSr), colour='orange') +
			geom_line(aes(y=AGRpc), colour='green') +
			geom_line(aes(y=GPS), colour='DarkGreen') 
	
	cnsc.df[, CALL:= as.integer(PR_CALL>=par['PRCALL.thr'])]
	
	
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
	if(!is.na(par['PRCALL.rmintrnlgpsblw'] | !is.na(par['PRCALL.rmintrnlgpsend'])))
	{
		cnsc.g	<- cnsc.1s[, {
					if(length(CALL_ID)==1)
						ans	<- NA_integer_
					if(length(CALL_ID)>1)
						ans	<- CALL_POS[seq.int(2,length(CALL_POS))]-CALL_LAST[seq.int(1,length(CALL_LAST)-1)]-1L
					list(CALL_LAST=CALL_LAST[seq.int(1,length(CALL_LAST)-1)], GAP_LEN= ans)
				}, by=c('TAXON','BLASTnCUT')]
		cnsc.g	<- merge(cnsc.1s, cnsc.g,  by=c('TAXON','BLASTnCUT','CALL_LAST'), all.x=1)
		tmp		<- cnsc.g[, which(GAP_LEN<par['PRCALL.rmintrnlgpsblw'] | (CALL_LAST>par['PRCALL.rmintrnlgpsend'] & !is.na(GAP_LEN)))]
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
	#	check if all called chunks in cut and raw contigs correspond to each other
	if(!is.na(par['PRCALL.cutrawgrace']))
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
		}
		if(length(tmp2))
			set(cnsc.1s, cnsc.1s[, which(TAXON%in%tmp2[, TAXON_CUT] & BLASTnCUT=='Y' & CALL_ID%in%tmp2[, CALL_ID_CUT])],'CALL_ID',NA_integer_)	
		tmp2	<- subset(tmp, CALL_LEN<=CALL_LEN_CUT)
		for(i in seq_len(nrow(tmp2)))	#keep cut
		{	
			cat('\nkeep only cut:', tmp2[i,TAXON_CUT])
			z	<- cnsc.df[, which( TAXON==tmp2$TAXON[i] & BLASTnCUT=='N' & SITE>=tmp2$CALL_POS[i] & SITE<=tmp2$CALL_LAST[i])]
			stopifnot( cnsc.df[z,all(CALL==1)])
			set(cnsc.df, z, 'CALL', 0L)
		}
		if(length(tmp2))
			set(cnsc.1s, cnsc.1s[, which(TAXON%in%tmp2[, TAXON] & BLASTnCUT=='N' & CALL_ID%in%tmp2[, CALL_ID])],'CALL_ID',NA_integer_)
		cnsc.1s	<- subset(cnsc.1s, !is.na(CALL_ID))		
	}	
	#	check if called contig has gaps of CALL=='0': if yes, return last non-gap before first CALL=='0
	if(!is.na(par['PRCALL.cutprdcthair']))
	{
		tmp		<- subset( cnsc.1s, CALL_LEN<par['PRCALL.cutprdcthair'] )
		if(nrow(tmp))
		{
			cat('\nFound predicted extra hair of length <',par['PRCALL.cutprdcthair'],'delete, n=',tmp[,sum(CALL_LEN)])
			set(tmp, NULL, 'CALL_LEN', tmp[, CALL_POS+CALL_LEN-1])
			for(i in seq_len(nrow(tmp)))
			{				
				set(cnsc.df, cnsc.df[, which(TAXON==tmp$TAXON[i] & BLASTnCUT==tmp$BLASTnCUT[i] & SITE>=tmp$CALL_POS[i] & SITE<=tmp$CALL_LEN[i])], 'CALL', 0L)
			}				
		}										
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
##	calculate cut statistics for each contig
##--------------------------------------------------------------------------------------------------------
haircut.get.cut.statistics.v150811<- function(cnsc, tx, par, outdir=NA, file=NA, mode='rolling')
{
	require(zoo)
	stopifnot(mode%in%c('rolling','overall'))
	#	count overall disagreement with consensus		
	if(mode=='overall')
	{
		cnsc.df	<- merge(tx, subset(tx, TAXON!='consensus')[, 
						{
							tmp		<- cnsc[ c('consensus',TAXON), seq.int(FIRST,LAST)]
							
							list( DSGR= dist.dna(tmp, model='indel' ), GPS=mean( as.character( cnsc[ TAXON, seq.int(FIRST,LAST)] )=='-' ), LEN=LAST-FIRST+1 )
						},by='TAXON'], by='TAXON')	
	}
	#	count rolling agreement with consensus
	if(mode=='rolling')
	{
		cnsc.df	<- merge(tx, subset(tx, TAXON!='consensus')[, 
						{
							tmp		<- cnsc[ c('consensus',TAXON), seq.int(FIRST,LAST)]
							agrpc	<- 1-rollapply( seq_len(LAST-FIRST+1), width=par['CNS_AGR.window'], FUN= function(z) dist.dna(tmp[,z], model='indel' )/length(z), align='center', partial=T )
							tmp		<- as.character( cnsc[ TAXON, seq.int(FIRST,LAST)] )=='-'
							gps		<- rollapply( seq_len(LAST-FIRST+1), width=par['GPS.window'], FUN= function(z) mean(tmp[z]), align='center', partial=T )
							list( SITE=seq.int(FIRST,LAST),  AGRpc=agrpc, GPS=gps   )
						},by='TAXON'], by='TAXON')		
		if(!is.na(file) & !is.na(outdir))
		{
			ggplot(cnsc.df, aes(x=SITE, ymax=AGRpc, ymin=0, y=GPS, fill=TAXON, group=TAXON)) + 
					geom_ribbon(alpha=0.5) + geom_line(colour='grey30') +
					scale_x_continuous(breaks=seq(0,15e3, 5e2)) + 
					facet_grid(TAXON~.) + theme_bw() + theme(strip.text=element_blank(), legend.position='bottom') +
					labs(fill='Contig', x='position on consensus w/o LTR', y='line: % gaps\nfill: % agreement with consensus')
			ggsave(w=10, h=2*cnsc.df[, length(unique(TAXON))], file= paste(outdir, '/', gsub('\\.fasta',paste('_CNSAGR_aw',par['CNS_AGR.window'],'_gw',par['GPS.window'],'.pdf',sep=''),basename(file)), sep=''))	
		}
	}
	cnsc.df
}
##--------------------------------------------------------------------------------------------------------
##	fit simple Binomial regression model to training data 
##--------------------------------------------------------------------------------------------------------
haircut.get.fitted.model.150811a<- function(indir, outfile)
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
												ctrchm	<- gamlss(ANS_CALL~AGRpc+GPS, data=ctrch, family=BI())				#as good as 'AGRpc+GPS+CNS_FRQr' in terms of FN, FP
												ctrchmc	<- data.table(CHUNK=chunk, BETA0=coef(ctrchm)[1], BETA1=coef(ctrchm)[2], BETA2=coef(ctrchm)[3])					
											}))
						}))
		#	deal with end of genome where little data is available: we estimated coefs for all sites > 1e4, now split up using CHUNK notation
		ctr		<- haircut.load.training.data(indir, 10001)
		tmp		<- seq(ctr[, min(SITE)-1],ctr[, max(SITE)+10],10)
		tmp[tmp>10000]
		tmp		<- as.data.table(expand.grid(CHUNK=as.character(tmp[tmp>10000]), BETA0=subset(ctrmc, CHUNK==10000)[, BETA0], BETA1=subset(ctrmc, CHUNK==10000)[, BETA1], BETA2=subset(ctrmc, CHUNK==10000)[, BETA2], stringsAsFactors=F))
		ctrmc	<- rbind(ctrmc, tmp)
		#	model predict function, so we save mem by not having to call 'predict'
		model.150811a.predict<- function(agrpc, gps, b0, b1, b2)
		{	
			stopifnot(all(!is.na(agrpc)), all(!is.na(gps)))
			b0[which(is.na(b0))]	<- 0
			b1[which(is.na(b1))]	<- 0
			b2[which(is.na(b2))]	<- 0
			exp(b0+b1*agrpc+b2*gps)/(exp(b0+b1*agrpc+b2*gps)+1)	
		}
		#	calculate Sensitivity & Specificity on training data
		ctrev	<- do.call('rbind', lapply(seq(1,11000,200), function(site)
						{
							cat('\nProcess',site,'\n')							
							ctr		<- haircut.load.training.data(indir, site)	
							tmp		<- seq(ctr[, min(SITE)-1],ctr[, max(SITE)+10],10)
							ctr[, CHUNK:=cut(SITE, breaks=tmp, labels=tmp[-length(tmp)])]
							print(ctr)
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
						}))		
		set(ctrev, NULL, 'CHUNK', ctrev[, as.numeric(CHUNK)])
		
		save(ctrmc, ctrev, model.150811a.predict, file=outfile)
	}
	list(coef=ctrmc, ev=ctrev, predict=model.150811a.predict)
}
##--------------------------------------------------------------------------------------------------------
##	fit Beta Binomial regression model to training data 
##--------------------------------------------------------------------------------------------------------
haircut.get.fitted.model.150816b<- function(indir, outfile)
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
												ctrch		<- subset(ctr, CHUNK==chunk)
												tmp			<- tryCatch( gamlss(ANS_CALL~FRQ+GPS, sigma.formula=~FRQ+GPS, data=ctrch, family=BB(), control=gamlss.control(trace=FALSE), i.control=glim.control(cyc=100,cc=1e-4)), warning=function(w) w, error=function(e) e)
												if(!is(tmp,'warning') & !is(tmp,'error'))
												{
													ctrchm	<- tmp
													tmp		<- tryCatch(vcov(ctrchm),warning=function(w) w, error=function(e) e)
													tmp2	<- !is(tmp,'warning') & !is(tmp,'error')
												}
												if(is(tmp,'warning') | is(tmp,'error'))
												{
													ctrchm	<- gamlss(ANS_CALL~FRQ+GPS, sigma.formula=~FRQ+GPS, data=ctrch, family=BI(), control=gamlss.control(trace=FALSE))
													tmp2	<- TRUE
												}													
												ctrchmc	<- data.table(CHUNK=chunk, BETA0=coef(ctrchm)[1], BETA1=coef(ctrchm)[2], BETA2=coef(ctrchm)[3], FAMILY=family(ctrchm)[1], CONVERGED=tmp2	)					
											}))
						}))
		#	deal with end of genome where little data is available: we estimated coefs for all sites > 1e4, now split up using CHUNK notation
		ctr		<- haircut.load.training.data(indir, 10001)
		tmp		<- seq(ctr[, min(SITE)-1],ctr[, max(SITE)+10],10)
		tmp		<- as.data.table(expand.grid(CHUNK=as.character(tmp[tmp>10000]), BETA0=subset(ctrmc, CHUNK==10000)[, BETA0], BETA1=subset(ctrmc, CHUNK==10000)[, BETA1], BETA2=subset(ctrmc, CHUNK==10000)[, BETA2], stringsAsFactors=F))
		ctrmc	<- rbind(ctrmc, tmp)
		#	model predict function, so we save mem by not having to call 'predict'
		model.150816b.predict<- function(frq, gps, b0, b1, b2)
		{	
			stopifnot(all(!is.na(frq)), all(!is.na(gps)))
			b0[which(is.na(b0))]	<- 0
			b1[which(is.na(b1))]	<- 0
			b2[which(is.na(b2))]	<- 0
			exp(b0+b1*frq+b2*gps)/(exp(b0+b1*frq+b2*gps)+1)	
		}
		#	calculate Sensitivity & Specificity on training data
		ctrev	<- do.call('rbind', lapply(seq(1,11000,200), function(site)
						{
							cat('\nProcess',site,'\n')							
							ctr		<- haircut.load.training.data(indir, site)	
							tmp		<- seq(ctr[, min(SITE)-1],ctr[, max(SITE)+10],10)
							ctr[, CHUNK:=cut(SITE, breaks=tmp, labels=tmp[-length(tmp)])]
							print(ctr)
							ctrp	<- merge(ctr, ctrmc, by='CHUNK')
							ctrp[, PR_CALL:=model.150816b.predict(FRQ, GPS, BETA0, BETA1, BETA2)]
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
						}))		
		set(ctrev, NULL, 'CHUNK', ctrev[, as.numeric(CHUNK)])
		
		save(ctrmc, ctrev, model.150816b.predict, file=outfile)
	}
	list(coef=ctrmc, ev=ctrev, predict=model.150816b.predict)
}
##--------------------------------------------------------------------------------------------------------
##	return EQ= 1 if raw and contigs are identical in that: 
##	- 	there is only once cut contig which disagrees on up x positions, where x is the difference in alignment lengths in the cut and raw files
##	- 	the concatenated cut contigs disagree on up x positions, where x is the difference in alignment lengths in the cut and raw files
##		gaps between cut contigs are not counted
##--------------------------------------------------------------------------------------------------------
haircut.get.identical.among.raw.and.cut.contigs.v150811 <- function(indir, files, png_id, blastncut)
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
	}
	#	check if there are contigs with corresponding name in 'cut' and that are identical
	#	differences in gaps are allowed: these will come from different alignment
	txe		<- dcast.data.table(tx, CNTG~CUT, value.var='TAXON')
	txe		<- subset( txe, !is.na(cut) & !is.na(raw) )
	if(nrow(txe) && c('cut','raw')%in%colnames(txe))
	{						
		txe		<- txe[, {
					x<- sub('^-*','',sub('-*$','',paste(as.vector(as.character(crs[['cut']][cut,])),collapse='')))
					y<- sub('^-*','',sub('-*$','',paste(as.vector(as.character(crs[['raw']][raw,])),collapse='')))
					if(nchar(x)==nchar(y))
						z	<- x==y
					if(nchar(x)!=nchar(y))
						z	<- gsub('-','',x)==gsub('-','',y)
					list(EQ=z)					
				}, by='CNTG']
		set(tx, which( tx[, CNTG]%in%subset(txe, EQ)[, CNTG] ), 'EQ', TRUE)						
	}
	#	concatenate cut contigs 
	txe		<- subset(tx, !is.na(CCNTG))
	if(nrow(txe))
	{
		txe		<- txe[, {
					z	<- sub('^-*','',paste(as.vector(as.character(crs[['cut']][TAXON[1],])),collapse=''))
					tmp	<- ncol(crs[['cut']])-nchar(z)	#number of initial gaps removed
					#print(tmp)
					z	<- sub('-*$','',z)
					tmp	<- tmp+nchar(z)					#number of sites to remove from next contig
					#print(tmp)
					#print(seq_len(length(TAXON))[-1])
					for(k in seq_len(length(TAXON))[-1])
					{
						if( tmp+1 < ncol(crs[['cut']])	)	#in some cases eg 15065_1_10, reversed contigs and non-reversed contigs are prodived. these cannot be concatenated to reconstitute the original raw contig
						{
							z2	<- paste(as.vector(as.character(crs[['cut']][TAXON[k], seq.int(tmp+1, ncol(crs[['cut']]))])), collapse='')
							z2	<- sub('-*$','',z2)
							tmp	<- tmp+nchar(z2)	
							#print(tmp)
							z	<- paste(z, z2, sep='')	
						}										
					}
					#print(nchar(z))
					list( CSEQ=z )	
				}, by='OCNTG']
		#	check if concatenated contig equals raw contig
		txe		<- subset(merge(tx, txe, by='OCNTG'), CUT=='raw')
		txe		<- txe[, {
					x<- sub('^-*','',sub('-*$','',paste(as.vector(as.character(crs[['raw']][TAXON,])),collapse='')))
					z	<- FALSE
					if(nchar(x)==nchar(CSEQ))
						z	<- x==CSEQ
					#print(abs(diff(c(nchar(x),nchar(CSEQ))))<=abs(diff(sapply(crs,ncol))))
					#print( gsub('-','',x) )
					#print( gsub('-','',CSEQ) )
					if(	abs(diff(c(nchar(x),nchar(CSEQ))))<=abs(diff(sapply(crs,ncol)))	)	#contigs may have cut and been put elsewhere, so only allow for a difference in # gaps that is not larger than the difference in alignment length
						z	<- gsub('-','',x)==gsub('-','',CSEQ)
					list(EQ=z)					
				}, by='OCNTG']
		set(tx, which( tx[, OCNTG]%in%subset(txe, EQ)[, OCNTG] ), 'EQ', TRUE)						
	}					
	tx		<- merge(tx, data.table(INFILE=files, CUT=blastncut), by='CUT')
	tx
}
haircut.get.training.contigs.byidentical.v150811<- function(indir, outfile, ctrain)
{
	options(show.error.messages = FALSE)		
	readAttempt		<-try(suppressWarnings(load(outfile)))
	options(show.error.messages = TRUE)	
	if( inherits(readAttempt, "try-error")	)
	{
		infiles	<- data.table(FILE=list.files(indir, pattern='fasta$', recursive=T))
		infiles[, PNG_ID:= gsub('_wRefs\\.fasta','',gsub('_cut|_raw','',basename(FILE)))]
		infiles[, BLASTnCUT:= regmatches(FILE,regexpr('cut|raw',FILE))]		 
		#	identify identical raw and cut contigs
		txe		<- infiles[,{
					#PNG_ID	<- '12559_1_11'
					#files	<- subset(infiles, PNG_ID=='12559_1_10')[, FILE]
					#blastncut<- subset(infiles, PNG_ID=='12559_1_10')[, BLASTnCUT]
					cat('\nProcess', PNG_ID)
					tx	<- haircut.get.identical.among.raw.and.cut.contigs.v150811(indir, FILE, PNG_ID, BLASTnCUT)
					tx
				}, by='PNG_ID']
		#	consider identical contigs. merge cut contigs into equal contigs, so we can look at all equal pairs
		txec	<- merge(subset(txe, EQ), ctrain, all.x=TRUE, by=c('PNG_ID','TAXON'))		
		tmp		<- subset(txe, EQ & CUT=='cut', select=c(PNG_ID, OCNTG, CCNTG, TAXON))
		setnames(tmp, c('TAXON','CCNTG'),c('TAXON_CUT','CCNTG_CUT'))
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
		#	if raw and concatenated cut contigs in curated and identical: should not happen
		tmp		<- txec[, which(CUT=='raw' & !is.na(CCNTG_CUT) & !is.na(ANS_LEN) & !is.na(ANS_LEN_CUT))]
		stopifnot(length(tmp)==0)
		#	if raw and concatenated cut contigs identical and raw in curated: don t use cut in training as 0
		tmp		<- subset(txec, CUT=='raw' & !is.na(CCNTG_CUT) & !is.na(ANS_LEN) & is.na(ANS_LEN_CUT))[, TAXON_CUT] 
		cat(paste('\nFound raw and concatenated cut contigs identical and raw in curated: don t use cut in training as 0, n=', length(tmp)))
		set(txec, txec[, which(TAXON%in%tmp & CUT=='cut')], 'USE_IN_TRAIN', 'N')
		#	if raw and concatenated cut contigs identical and cut in curated: don t use raw in training as 0
		tmp		<- txec[, which(CUT=='raw' & !is.na(CCNTG_CUT) & is.na(ANS_LEN) & !is.na(ANS_LEN_CUT))] 
		cat(paste('\nFound raw and concatenated cut contigs identical and cut in curated: don t use raw in training as 0, n=', length(tmp)))
		set(txec, tmp, 'USE_IN_TRAIN', 'N')		
		#	if raw and cut contigs identical and cut in curated: don t use raw in training as 0
		tmp		<- txec[, which(CUT=='raw' & is.na(CCNTG_CUT) & is.na(ANS_LEN) & !is.na(ANS_LEN_CUT))] 
		cat(paste('\nFound raw and cut contigs identical and cut in curated: don t use raw in training as 0, n=', length(tmp)))
		set(txec, tmp, 'USE_IN_TRAIN', 'N')
		#	if raw and concatenated cut contigs identical and raw in curated: don t use cut in training as 0
		tmp		<- subset(txec, CUT=='raw' & is.na(CCNTG_CUT) & !is.na(ANS_LEN) & is.na(ANS_LEN_CUT))[, TAXON_CUT] 
		cat(paste('\nFound raw and cut contigs identical and raw in curated: don t use cut in training as 0, n=', length(tmp)))
		set(txec, txec[, which(TAXON%in%tmp & CUT=='cut')], 'USE_IN_TRAIN', 'N')
		#	can now set ANS_LEN to NA from above case 
		set(txec, txec[, which(is.na(ANS_FIRST) & !is.na(ANS_LEN))], 'ANS_LEN', NA_integer_)
		#	delete tmp rows for cut-cut
		#subset(txec, CUT=='raw')
		#tmp		<- subset(txec, CUT=='cut' )[, {
		#			stopifnot( all(USE_IN_TRAIN==USE_IN_TRAIN[1])	)
		#			list(CUT=CUT[1], FIRST=FIRST[1], LAST=LAST[1], CRS_ID=CRS_ID[1], CNTG=CNTG[1], OCNTG=OCNTG[1], CCNTG=CCNTG[1], EQ=EQ[1], INFILE=INFILE[1], ANS_FILE=ANS_FILE[1], ANS_LEN=ANS_LEN[1], ANS_FIRST=ANS_FIRST[1], ANS_LAST=ANS_LAST[1], TAXON_CUT=TAXON_CUT[1], CCNTG_CUT=CCNTG_CUT[1], ANS_LEN_CUT=ANS_LEN_CUT[1], USE_IN_TRAIN=USE_IN_TRAIN[1])
		#		}, by=c('PNG_ID','TAXON')]
		#subset(txec, CUT=='raw')		
		#	delete tmp cols
		set(txec, NULL, c('TAXON_CUT','CCNTG_CUT','ANS_LEN_CUT'), NULL)
		setkey(txec, PNG_ID, TAXON, CUT)
		txec	<- unique(txec)
		#
		#	add contigs that are in curated and unique amongst cut and raw
		#
		tmp		<- merge(subset(txe, !EQ), ctrain, by=c('PNG_ID','TAXON'))
		tmp[, USE_IN_TRAIN:='Y']
		stopifnot( length(intersect( tmp[, TAXON], txec[,TAXON] ))==0	)
		txec	<- rbind(txec, tmp, use.names=TRUE)
		save(txec, file=outfile)
		#print( txec[, table(USE_IN_TRAIN)] )
	}	
	txec
}
##--------------------------------------------------------------------------------------------------------
##	determine the consensus sequence in the reference alignment
##	variable sites may have an NA consensus base call, depending on the consensus quantile cutoff 'FRQx.quantile'
##--------------------------------------------------------------------------------------------------------
haircut.getconsensus.v150811	<- function(seq, par, bases=c('a','c','g','t','-') )
{
	stopifnot('FRQx.quantile'%in%names(par))
	#	calculate base frequency Profile among References: rp
	rp		<- sapply(seq_len(ncol(seq)), function(i) base.freq(seq[,i], freq=F, all=T))
	rp		<- as.data.table( t(rp[bases,]) )
	rp[, SITE:= seq_len(nrow(rp))]
	rp		<- melt(rp, id.vars='SITE', variable.name='BASE', value.name='FRQ')
	set(rp, NULL, 'BASE', rp[, factor(BASE, levels=bases, labels=bases)])
	#	renormalize bases and ignore all other calls amongst references	
	rp		<- rp[ , list(BASE=BASE, FRQ=FRQ/sum(FRQ), FRQxu= max(FRQ), FRQx= max(FRQ/sum(FRQ))), by='SITE']
	#	get consensus threshold as lower quantile of a one-inflated Beta distribution, par['FRQx.quantile']	
	setkey(rp, SITE)
	tmp		<- unique(rp)
	#	ggplot( tmp, aes(x=FRQx)) + geom_histogram()	#	looks like BEINF
	#set(rp, NULL, 'DUMMY', rp[, as.integer(ceiling(SITE/500))])
	#ggplot(subset(rp, DUMMY<=2 ), aes(x=SITE, y=FRQ, colour=BASE, group=BASE)) + geom_line() + facet_wrap(~DUMMY, ncol=1, scales='free')
	#ggplot(subset(rp, BASE=='a'), aes(x=SITE, y=FRQx, colour=BASE, group=BASE)) + geom_line() + facet_wrap(~DUMMY, ncol=1, scales='free') + theme_bw() + theme(strip.text = element_blank())
	if(is.na(par['FRQx.thr']))
	{
		qu.m	<- gamlss(FRQx~1, data=tmp, family=BEINF)
		qu.par	<- c('mu'=as.double(predict(qu.m, type='response', what='mu')[1]), 'sigma'=as.double(predict(qu.m, type='response', what='sigma')[1]), 'nu'=as.double(predict(qu.m, type='response', what='nu')[1]), 'tau'=as.double(predict(qu.m, type='response', what='tau')[1]) ) 
		qu.par	<- qBEINF( par['FRQx.quantile'], mu=qu.par['mu'], sigma=qu.par['sigma'], nu=qu.par['nu'], tau=qu.par['tau'])		
	}
	if(!is.na(par['FRQx.thr']))
		qu.par	<- par['FRQx.thr']
	cat(paste('\nMinimum base frequency to call a consensus base=', round(qu.par, d=3)))
	set(rp, rp[, which(FRQx<qu.par)], 'BASE', NA_character_)
	set(rp, NULL, 'BASE', rp[, factor(as.character(BASE), levels=bases, labels=bases)])
	#	get consensus base above lower quantile. Some consensus bases will be NA
	crp		<- rp[, list(CNS_BASE= BASE[ which.max(FRQ) ], CNS_FRQ=FRQxu[1], CNS_AGRpc=FRQx[1]), by='SITE']
	set(crp, crp[, which(is.na(CNS_BASE))], 'CNS_BASE','?')	
	tmp		<- crp[, as.DNAbin(t(as.matrix(CNS_BASE)))]
	rownames(tmp)	<- 'consensus'
	list(DNAbin=tmp, DATATABLE=crp)
}
##--------------------------------------------------------------------------------------------------------
##	wrapper to call 'haircutwrap.get.call.for.PNG_ID.150811'
##--------------------------------------------------------------------------------------------------------
haircutwrap.get.call.for.PNG_ID.150811<- function(indir.st,indir.al,outdir,ctrmc,predict.fun,par,ctrain=NULL)
{
	infiles	<- data.table(INFILE=list.files(indir.st, pattern='\\.R$', recursive=T))
	infiles[, PNG_ID:= gsub('_wRefs.*','',gsub('_cut|_raw','',INFILE))]
	infiles[, BLASTnCUT:= regmatches(INFILE,regexpr('cut|raw',INFILE))]
	set(infiles, NULL, 'BLASTnCUT', infiles[, factor(BLASTnCUT, levels=c('cut','raw'), labels=c('Y','N'))])
	alfiles <- data.table(ALFILE=list.files(indir.al, pattern='\\.fasta$', recursive=T))
	alfiles[, PNG_ID:= gsub('_wRefs.*','',gsub('_cut|_raw','',basename(ALFILE)))]
	alfiles[, BLASTnCUT:= regmatches(basename(ALFILE),regexpr('cut|raw',basename(ALFILE)))]
	set(alfiles, NULL, 'BLASTnCUT', alfiles[, factor(BLASTnCUT, levels=c('cut','raw'), labels=c('Y','N'))])
	infiles	<- merge(infiles, alfiles, by=c('PNG_ID','BLASTnCUT'))
	
	#	predict by PANGEA_ID
	invisible( infiles[,
					{
						cat(paste('\nProcess', PNG_ID))
						#PNG_ID<- png_id	<- '13548_1_21'
						#files	<- subset(infiles, PNG_ID==png_id)[, INFILE]
						#alfiles	<- subset(infiles, PNG_ID==png_id)[, ALFILE]
						#bc		<- subset(infiles, PNG_ID==png_id)[, BLASTnCUT]			
						#tmp		<- haircut.get.call.for.PNG_ID.v150811(indir.str, indir.al, png_id, files, alfiles, bc, par, ctrmc, predict.fun)
						tmp		<- haircut.get.call.for.PNG_ID.v150811(indir.str, indir.al, PNG_ID, INFILE, ALFILE, BLASTnCUT, par, ctrmc, predict.fun)
						crs		<- tmp$crs
						cnsc.df	<- tmp$cnsc.df	
						#	handle output
						tmp		<- paste(outdir,'/',PNG_ID,'_wref_nohaironraw.fasta',sep='')
						cat('\nWrite to file', tmp)
						write.dna(crs[['N']], file=tmp, format='fasta', colsep='', nbcol=-1)
						tmp		<- paste(outdir,'/',PNG_ID,'_wref_nohaironcut.fasta',sep='')
						cat('\nWrite to file', tmp)
						write.dna(crs[['Y']], file=tmp, format='fasta', colsep='', nbcol=-1)
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
								geom_line(aes(y=CUR_CALL), colour='red') +
								scale_x_continuous(breaks=seq(0,15e3, ifelse(cnsc.df[,max(SITE)]>5e2, 5e2, floor(cnsc.df[,max(SITE)/3])))) + 
								facet_wrap(~TAXONnCUT, ncol=1) + theme_bw() + theme(legend.position='bottom') +
								labs(fill='Contig BLASTnCUT', x='position on consensus w/o LTR', y='fill: predicted call\nblack line: predictive probability\nred line: curated call')
						tmp	<- paste(outdir, '/', PNG_ID, '_wref_nohaironcutraw.pdf',sep='')
						cat('\nPlot to file', tmp)
						ggsave(w=10, h=3*cnsc.df[, length(unique(TAXON))], file=tmp)
						NULL
					}, by='PNG_ID'] )
}	
##--------------------------------------------------------------------------------------------------------
##	wrapper to call 'haircutwrap.get.call.for.PNG_ID'
##--------------------------------------------------------------------------------------------------------
haircutwrap.get.call.for.PNG_ID.150814<- function(indir.st,indir.al,outdir,ctrmc,predict.fun,par,ctrain=NULL)
{
	infiles	<- data.table(INFILE=list.files(indir.st, pattern='\\.R$', recursive=T))
	infiles[, PNG_ID:= gsub('_wRefs.*','',gsub('_cut|_raw','',INFILE))]
	infiles[, BLASTnCUT:= regmatches(INFILE,regexpr('cut|raw',INFILE))]
	set(infiles, NULL, 'BLASTnCUT', infiles[, factor(BLASTnCUT, levels=c('cut','raw'), labels=c('Y','N'))])
	alfiles <- data.table(ALFILE=list.files(indir.al, pattern='\\.fasta$', recursive=T))
	alfiles[, PNG_ID:= gsub('_wRefs.*','',gsub('_cut|_raw','',basename(ALFILE)))]
	alfiles[, BLASTnCUT:= regmatches(basename(ALFILE),regexpr('cut|raw',basename(ALFILE)))]
	set(alfiles, NULL, 'BLASTnCUT', alfiles[, factor(BLASTnCUT, levels=c('cut','raw'), labels=c('Y','N'))])
	infiles	<- merge(infiles, alfiles, by=c('PNG_ID','BLASTnCUT'))
	
	#	predict by PANGEA_ID
	cnsc.info	<-  infiles[,
			{
				cat(paste('\nProcess', PNG_ID))
				if(0)	#devel
				{
					PNG_ID<- png_id	<- '15172_1_32'
					PNG_ID<- png_id	<- '12559_1_11'
					PNG_ID<- png_id	<- '14728_1_84'
					PNG_ID<- png_id	<- '14938_1_10'
					PNG_ID<- png_id	<- '14728_1_82'
					PNG_ID<- png_id	<- '12559_1_24'
					PNG_ID<- png_id	<- '12559_1_81'
					#PNG_ID<- png_id	<- '12559_1_87'
					files	<- subset(infiles, PNG_ID==png_id)[, INFILE]
					alfiles	<- subset(infiles, PNG_ID==png_id)[, ALFILE]
					bc		<- subset(infiles, PNG_ID==png_id)[, BLASTnCUT]
					tmp		<- haircut.get.call.for.PNG_ID.150814(indir.str, indir.al, png_id, files, alfiles, bc, par, ctrmc, predict.fun)
				}
				#
				tmp		<- haircut.get.call.for.PNG_ID.150814(indir.str, indir.al, PNG_ID, INFILE, ALFILE, BLASTnCUT, par, ctrmc, predict.fun)
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
						geom_line(aes(y=CNS_FRQr), colour='blue') +
						geom_line(aes(y=CUR_CALL), colour='red') +
						scale_x_continuous(breaks=seq(0,15e3, ifelse(cnsc.df[,max(SITE)]>5e2, 5e2, floor(cnsc.df[,max(SITE)/3])))) + 
						facet_wrap(~TAXONnCUT, ncol=1) + theme_bw() + theme(legend.position='bottom') +
						labs(fill='Contig BLASTnCUT', x='position on consensus w/o LTR', y='fill: predicted call\nblack line: predictive probability\nred line: curated call')
				tmp	<- paste(outdir, '/', PNG_ID, '_wref_nohaironcutraw.pdf',sep='')
				cat('\nPlot to file', tmp)
				ggsave(w=10, h=3*cnsc.df[, length(unique(TAXON))], file=tmp)
				#	report confidence score
				subset(cnsc.df, CALL==1)[, list(QUANTILE=c(0,0.01,0.05,0.1,0.2,0.5), PR_CALL=quantile(PR_CALL, p=c(0,0.01,0.05,0.1,0.2,0.5))), by=c('PNG_ID','TAXON','BLASTnCUT')]
			}, by='PNG_ID']
}	
##--------------------------------------------------------------------------------------------------------
##	process all files in indir with 'haircut.get.cut.statistics'
##--------------------------------------------------------------------------------------------------------
haircutwrap.get.cut.statistics.v150811<- function(indir, par, outdir=indir)	
{
	require(zoo)
	#	read just one file
	#	determine statistics for each contig after LTR
	infiles		<- data.table(FILE=list.files(indir, pattern='fasta$', recursive=T))
	infiles[, PNG_ID:= gsub('_wRefs\\.fasta','',gsub('_cut|_raw','',FILE))]
	infiles[, BLASTnCUT:= regmatches(FILE,regexpr('cut|raw',FILE))]
	set(infiles, NULL, 'BLASTnCUT', infiles[, factor(BLASTnCUT, levels=c('cut','raw'), labels=c('Y','N'))])
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
	#	infiles[, which(grepl('12559_1_81',FILE))]	fls<- 51
	#	process files
	for(fls in infiles[, seq_along(FILE)])
	{
		file	<- paste(indir, infiles[fls, FILE], sep='/')
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
		#	determine base frequencies at each site amongst references.
		tmp		<- cr[subset(tx, CONTIG==0)[, TAXON],]
		tmp		<- haircut.getconsensus.v150811(tmp, par, bases=c('a','c','g','t','-') )	#	CoNSensus of References: cnsr
		cnsr	<- tmp$DNAbin
		cnsr.df	<- tmp$DATATABLE
		#	for each contig, determine %agreement with consensus on rolling window
		cnsc	<- rbind(cnsr, cr[subset(tx, CONTIG==1)[, TAXON],])
		#	determine first and last non-gap sites
		tx		<- data.table(	TAXON= rownames(cnsc), 
				FIRST= apply( as.character(cnsc), 1, function(x) which(x!='-')[1] ),
				LAST= ncol(cnsc)-apply( as.character(cnsc), 1, function(x) which(rev(x)!='-')[1] )+1L		)
		tx		<- subset(tx, !is.na(FIRST) & !is.na(LAST))	#	some contigs only map into LTR
		#	get cut statistics
		cnsc.df	<- haircut.get.cut.statistics.v150811(cnsc, tx, par, outdir=NA, file=NA, mode='rolling')
		#	get rolling CNS_FRQ
		cnsr.df[, CNS_FRQr:= rollapply( seq_len(nrow(cnsr.df)), width=par['CNS_FRQ.window'], FUN= function(z) mean(cnsr.df$CNS_FRQ[z]), align='center', partial=T )]		
		#	get rolling CNS_GPS
		tmp		<- subset(tx, TAXON=='consensus')[, {
					tmp		<- as.character( cnsr.df$CNS_BASE[seq.int(FIRST,LAST)] )=='-'
					list(SITE=seq.int(FIRST,LAST), CNS_GPSr=rollapply( seq_len(LAST-FIRST+1), width=par['GPS.window'], FUN= function(z) mean(tmp[z]), align='center', partial=T ))
				}, by='TAXON']
		cnsr.df	<- merge(cnsr.df, subset(tmp, select=c(SITE,CNS_GPSr)), all.x=1, by='SITE')
		cnsc.df	<- merge(cnsc.df, subset(cnsr.df, select=c(SITE, CNS_FRQr,CNS_GPSr)), all.x=TRUE, by='SITE')
		cnsc.df[, PNG_ID:= infiles[fls, PNG_ID]]
		cnsc.df[, BLASTnCUT:= infiles[fls, BLASTnCUT]]
		cat(paste('\nSave contigs, n=', cnsc.df[, length(unique(TAXON))]))
		#	save
		file	<- paste(outdir, '/', gsub('\\.fasta',paste('_HAIRCUTSTAT_thr',100*par['FRQx.quantile'],'_aw',par['CNS_AGR.window'],'_fw',par['CNS_FRQ.window'],'_gw',par['GPS.window'],'.R',sep=''),basename(file)), sep='')
		cat(paste('\nSave to', file))
		save(cnsc, cnsc.df, file=file)		
	}	
	#	plot by PANGEA_ID
	invisible( infiles[,
					{
						cat(paste('\nPlot', PNG_ID))
						tmp		<- sapply(FILE, function(x)
								{
									paste(outdir, '/', gsub('\\.fasta',paste('_HAIRCUTSTAT_thr',100*par['FRQx.quantile'],'_aw',par['CNS_AGR.window'],'_fw',par['CNS_FRQ.window'],'_gw',par['GPS.window'],'.R',sep=''),basename(x)), sep='')
								}) 
						tmp		<- do.call('rbind',lapply(tmp, function(x)
										{
											load(x)
											cnsc.df
										}))
						tmp[, TAXONnCUT:= paste(TAXON,BLASTnCUT,sep='_BLASTnCUT:')]
						ggplot(tmp, aes(x=SITE, ymax=AGRpc, ymin=0, fill=BLASTnCUT, group=TAXONnCUT)) + 
								geom_ribbon(alpha=0.5) + 
								geom_line(aes(y=GPS), colour='grey30') +
								geom_line(aes(y=CNS_FRQr), colour='#FF7F00') +
								scale_x_continuous(breaks=seq(0,15e3, ifelse(tmp[,max(SITE)]>5e2, 5e2, floor(tmp[,max(SITE)/3])))) + 
								facet_wrap(~TAXONnCUT, ncol=1) + theme_bw() + theme(legend.position='bottom') +
								labs(fill='Contig BLASTnCUT', x='position on consensus w/o LTR', y='black line: % gaps\norange line: % consensus call amongst references\nfill: % agreement with consensus')
						ggsave(w=10, h=2*cnsc.df[, length(unique(TAXON))], file= paste(outdir, '/', PNG_ID, '_CNSAGR_aw', par['CNS_AGR.window'],'_fw',par['CNS_FRQ.window'], '_gw',par['GPS.window'],'.pdf',sep=''))
						NULL
					}, by='PNG_ID'] )
}