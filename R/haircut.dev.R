


dev.haircut<- function()	
{
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
	if(0)
	{
		#	process several files
		indir	<- '/Users/Oliver/Dropbox\ (Infectious Disease)/PANGEA_data/InterestingContigAlignments'
		outdir	<- '/Users/Oliver/Dropbox\ (Infectious Disease)/OR_Work/2015/2015_PANGEA_haircut/interesting_150408'		
		par		<- c('FRQx.quantile'=0.05, 'FRQx.thr'=0.566, 'CNS_FRQ.window'=100, 'CNS_AGR.window'=200, 'GPS.window'=200)
		haircutwrap.get.cut.statistics(indir, par, outdir=outdir)
	}
	if(0)
	{
		#	run mafft --add to get contigs+ref
		indir	<- '/Users/Oliver/Dropbox\ (Infectious Disease)/PANGEA_data/contigs_150408'
		outdir	<- '/Users/Oliver/Dropbox\ (Infectious Disease)/OR_Work/2015/2015_PANGEA_haircut/contigs_150408_wref'
		reffile	<- '/Users/Oliver/Dropbox\ (Infectious Disease)/PANGEA_data/HIV1_COM_2012_genome_DNA_WithExtraA1UG.fasta'
		haircutwrap.align.contigs.with.ref(indir, reffile, outdir)
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
	if(1)	#	fit model
	{
		indir	<- '/Users/Oliver/Dropbox\ (Infectious Disease)/OR_Work/2015/2015_PANGEA_haircut/contigs_150408_train'
		outfile	<- '/Users/Oliver/Dropbox\ (Infectious Disease)/OR_Work/2015/2015_PANGEA_haircut/model_150816a.R'
		tmp		<- haircut.get.fitted.model.150816a(indir, outfile)
	}
	if(0)	#	call contigs on training data and plot
	{		
		mfile	<- '/Users/Oliver/Dropbox\ (Infectious Disease)/OR_Work/2015/2015_PANGEA_haircut/model_150816a.R'
		#	get model coefficients across the chunks
		tmp						<- haircut.get.fitted.model.150816a(NULL, mfile)
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
		haircutwrap.get.call.for.PNG_ID(indir.st,indir.al,outdir,ctrmc,ctrev,predict.fun,txe,par,ctrain=ctrain)
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
##	HAIRCUT program, version 15086 to: 
##	- align contigs to references
##	- calculate and save haircut statistics
##--------------------------------------------------------------------------------------------------------
prog.haircut.150806<- function()
{
	if(0)
	{		
		indir	<- paste(DATA, 'contigs_150408', sep='/' )
		outdir	<- paste(DATA, 'contigs_150408_wref', sep='/' )
		outdir	<- paste('/Users/Oliver/Dropbox\ (Infectious Disease)/OR_Work/2015/2015_PANGEA_haircut', 'contigs_150408_wref', sep='/' )
		reffile	<- paste(DATA, 'HIV1_COM_2012_genome_DNA_WithExtraA1UG.fasta', sep='/' )	
		haircutwrap.align.contigs.with.ref(indir, reffile, outdir)	
	}
	if(0)
	{
		indir	<- paste(DATA, 'contigs_150408_wref', sep='/' )
		outdir	<- paste(DATA, 'contigs_150408_wref_cutstat', sep='/' )		
		par		<- c('FRQx.quantile'=0.05, 'FRQx.thr'=0.566, 'CNS_FRQ.window'=100, 'CNS_AGR.window'=200, 'GPS.window'=200)
		haircutwrap.get.cut.statistics(indir, par, outdir=outdir)
	}
	if(0)
	{
		indir	<- paste(DATA, 'contigs_150408_wref', sep='/' )
		outdir	<- paste(DATA, 'contigs_150408_wref_cutstat', sep='/' )		
		par		<- c('FRQx.quantile'=NA, 'FRQx.thr'=NA, 'CNS_FRQ.window'=200, 'CNS_AGR.window'=200, 'GPS.window'=200)
		haircutwrap.get.cut.statistics.150815(indir, par, outdir=outdir)
	}
	if(0)
	{
		indir	<- paste(DATA, 'contigs_150408_wref', sep='/' )
		outdir	<- paste(DATA, 'contigs_150408_wref_cutstat', sep='/' )		
		par		<- c('FRQx.quantile'=NA, 'FRQx.thr'=NA, 'CNS_FRQ.window'=200, 'CNS_AGR.window'=200, 'GPS.window'=200)
		haircutwrap.get.cut.statistics.150815(indir, par, outdir=outdir)
	}
	if(0)
	{
		indir.st	<- paste(DATA,'contigs_150408_wref_cutstat',sep='/')
		indir.al	<- paste(DATA,'contigs_150408_wref',sep='/')
		outdir		<- paste(DATA,'contigs_150408_model150816a',sep='/')
		cmd			<- cmd.haircut.call(indir.st, indir.al, outdir)		
	}
	if(1)
	{
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
			cmd			<- cmd.haircut.call(indir.st, indir.al, outdir, trainfile=trainfile, batch.n=batch.n, batch.id=batch.id)
			cmd			<- cmd.hpcwrapper(cmd, hpc.nproc= 1, hpc.q='pqeelab', hpc.walltime=4, hpc.mem="5000mb")
			cat(cmd)		
			outdir		<- paste(DATA,"tmp",sep='/')
			outfile		<- paste("hrct",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
			cmd.hpccaller(outdir, outfile, cmd)	
			stop()
		}	
		quit("no")
	}
	if(0)
	{
		#	get model coefficients across the chunks
		mfile					<- paste(DATA,'model_150816a.R',sep='/')		
		tmp						<- haircut.get.fitted.model.150816a(NULL, mfile)
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