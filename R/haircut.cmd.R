PR.PACKAGE					<- "PANGEAhaircut"
PR.STARTME					<- system.file(package=PR.PACKAGE, "misc", "PANGEAhaircut.startme.R")
#PR.STARTME					<- '/Users/Oliver/git/HPTN071sim/source/rPANGEAHIVsim/misc/rPANGEAHIV.startme.R'
#PR.STARTME					<- '/work/or105/libs/HPTN071sim/source/rPANGEAHIVsim/misc/rPANGEAHIV.startme.R'
PR.VARIOUS					<- paste(PR.STARTME," -exe=VARIOUS",sep='')

#' @export
PR.HAIRCUT.CALL				<- paste('Rscript',system.file(package=PR.PACKAGE, "haircut.call.contigs.Rscript"))
#' @export
PR.HAIRCUT.CUTSTAT			<- paste('Rscript',system.file(package=PR.PACKAGE, "haircut.cutstat.contigs.Rscript"))
#' @export
HPC.MPIRUN					<- {tmp<- c("mpirun","mpiexec"); names(tmp)<- c("debug","cx1.hpc.ic.ac.uk"); tmp}
#' @export
HPC.CX1.IMPERIAL			<- "cx1.hpc.ic.ac.uk"		#this is set to system('domainname',intern=T) for the hpc cluster of choice
#' @export
HPC.MEM						<- "1750mb"
#' @export
HPC.CX1.IMPERIAL.LOAD		<- "module load intel-suite mpi mafft/7 R/3.2.0"


######################################################################################
#' @export
cmd.hpcsys<- function()
{
	tmp<- system('domainname',intern=T)
	if(!nchar(tmp))	tmp<- "debug"
	tmp
}
##--------------------------------------------------------------------------------------------------------
##	call to various 
##	olli originally written 06-08-2015
##--------------------------------------------------------------------------------------------------------
cmd.various<- function(prog= PR.VARIOUS)
{
	cmd		<- "#######################################################
# start: run VARIOUS
#######################################################"
	cmd		<- paste(cmd, '\n', prog, '\n', sep='')
	cmd		<- paste(cmd,"#######################################################
# end: run VARIOUS
#######################################################\n",sep='')
	cmd
}
##--------------------------------------------------------------------------------------------------------
##	process all files in indir with 'haircut.align.contigs.with.ref'
##--------------------------------------------------------------------------------------------------------
#' @title Align cut and raw contigs to set of references
#' @import data.table zoo plyr ape reshape2 ggplot2
#' @export
#' @example example/ex.cmd.align.contigs.with.ref.R
cmdwrap.align.contigs.with.ref<- function(indir, outdir, reffile=NA, batch.n=NA, batch.id=NA)
{
	if(is.na(reffile))
		reffile	<- system.file(package="PANGEAhaircut", "HIV1_COM_2012_genome_DNA_WithExtraA1UG.fasta")
	infiles	<- data.table(INFILE=list.files(indir, pattern='fasta$',recursive=T))
	infiles	<- subset(infiles, !grepl('Curated', INFILE))
	infiles[, OUTFILE:= gsub('\\.fasta','_wRefs\\.fasta', gsub('_hiv|_HIV','',basename(INFILE)))]
	if(!is.na(batch.n) & !is.na(batch.id))
	{
		infiles[, BATCH:= ceiling(seq_len(nrow(infiles))/batch.n)]
		infiles		<- subset(infiles, BATCH==batch.id)
	}	
	#INFILE<- '12559_1_1.fasta'
	#OUTFILE<- '12559_1_1_wRefs.fasta'
	#cat(cmd.align.contigs.with.ref(paste(indir,'/',INFILE,sep=''), reffile, paste(outdir,'/',OUTFILE,sep='')))
	tmp		<- infiles[, {
				list(CMD=cmd.align.contigs.with.ref(paste(indir,'/',INFILE,sep=''), reffile, paste(outdir,'/',OUTFILE,sep='')))
			}, by='INFILE']
	cmd		<- "\n#######################################################
# start: run cmdwrap.align.contigs.with.ref
#######################################################"
	cmd		<- paste(cmd,paste(tmp$CMD, collapse='\n'),sep='\n')
	cmd		<- paste(cmd,"\n#######################################################
# end: run cmdwrap.align.contigs.with.ref
#######################################################\n",sep='')
	cmd
}
##--------------------------------------------------------------------------------------------------------
##	call to MAFFT to align contigs with reference compendium
##--------------------------------------------------------------------------------------------------------
cmd.align.contigs.with.ref<- function(infile, reffile, outfile)
{
	#mafft --reorder --anysymbol --add new_sequences --auto input
	tmp		<- c( 	gsub(' ','\\ ',gsub('(','\\(',gsub(')','\\)',infile,fixed=T),fixed=T),fixed=T),
			gsub(' ','\\ ',gsub('(','\\(',gsub(')','\\)',reffile,fixed=T),fixed=T),fixed=T),
			gsub(' ','\\ ',gsub('(','\\(',gsub(')','\\)',outfile,fixed=T),fixed=T),fixed=T)
	)
	cmd		<- paste('mafft --anysymbol --add ',tmp[1],' --auto ',tmp[2],' > ',tmp[3], sep='')
	cmd
}
##--------------------------------------------------------------------------------------------------------
##	command line generator for 'haircut.call.contigs.Rscript'
##--------------------------------------------------------------------------------------------------------
#' @export
cmd.haircut.call<- function(indir.st, indir.al, outdir, mfile=NA, trainfile=NA, batch.n=NA, batch.id=NA, prog=PR.HAIRCUT.CALL )	
{
	cmd<- "\n#######################################################
# start: run haircut.call.contigs.Rscript
#######################################################"
	cmd		<- paste(cmd, paste("\necho \'run ",prog,"\'\n",sep=''))
	cmd		<- paste(cmd, paste(prog,' -indir.st=',indir.st,' -indir.al=',indir.al,' -outdir=',outdir, sep=''))
	if(!is.na(mfile))
		cmd	<- paste(cmd, ' -mfile=',mfile, sep='')	
	if(!is.na(trainfile))
		cmd	<- paste(cmd, ' -trainfile=',trainfile, sep='')
	if(!is.na(batch.n) & !is.na(batch.id))
		cmd	<- paste(cmd, ' -batch.n=',batch.n, ' -batch.id=',batch.id, sep='')
	cmd		<- paste('\n',cmd,paste("\necho \'end ",prog,"\'\n",sep=''))
	cmd		<- paste(cmd,"#######################################################
# end: run haircut.call.contigs.Rscript
#######################################################\n",sep='')
	cmd
}
##--------------------------------------------------------------------------------------------------------
##	command line generator to run 'haircut.cutstat.contigs.Rscript' and 'haircut.call.contigs.Rscript' one after each other
##--------------------------------------------------------------------------------------------------------
#' @title Command line generator to run the Haircut pipeline
#' @description The haircut pipeline involves two steps. 
#' First, descriptive statistics are calculated for each contig and the consensus sequence of that associated references.
#' Second, these descriptive statistics are used to calculate a call probability for 10 base pair long contig chunks.
#' The call probability is modeled as a function of the descriptive statistics. 
#' The underlying statistical model is pre-computed and supplied with the R package.
#' Based on the call probabilites, 10bp chunks of each contig are called (yes=1, no=0). 
#' If cut/raw contigs correspond to each other, only one of both is returned.
#' 	 
#' @example example/ex.cmd.haircut.pipeline.R
#' @export
cmd.haircut.pipeline<- function(indir, outdir, batch.n=NA, batch.id=NA)
{
	#create specific outdir
	cmd				<- "#######################################################
#
# start: run haircut.pipeline
#
#######################################################\n"	
	#create temporary directories
	aldir		<- paste('algnd_',format(Sys.time(),"%y-%m-%d-%H-%M-%S"),sep='')
	cutdir		<- paste('cutstat_',format(Sys.time(),"%y-%m-%d-%H-%M-%S"),sep='')
	outdir.lcl	<- paste('call_',format(Sys.time(),"%y-%m-%d-%H-%M-%S"),sep='')
	cmd			<- paste(cmd,"CWD=$(pwd)\n",sep='\n')
	cmd			<- paste(cmd,"echo $CWD\n",sep='')
	aldir		<- paste("$CWD/",aldir,sep='')
	cutdir		<- paste("$CWD/",cutdir,sep='')
	outdir.lcl	<- paste("$CWD/",outdir.lcl,sep='')
	cmd			<- paste(cmd,"mkdir -p ",aldir,'\n',sep='')
	cmd			<- paste(cmd,"mkdir -p ",cutdir,'\n',sep='')
	cmd			<- paste(cmd,"mkdir -p ",outdir.lcl,'\n',sep='')
	#run alignment of references on batch into aldir
	cmd			<- paste(cmd, cmdwrap.align.contigs.with.ref(indir, aldir, batch.n=batch.n, batch.id=batch.id), sep='')
	#run cutstat on all seqs in aldir
	cmd			<- paste(cmd, cmd.haircut.cutstat(aldir, cutdir), sep='')
	#run call on all seqs in cutdir
	cmd			<- paste(cmd, cmd.haircut.call(cutdir, indir, outdir.lcl), sep='')
	#copy to destination
	cmd			<- paste(cmd, "\nmv ",outdir.lcl,"/* ",outdir,"\n",sep='')
	cmd			<- paste(cmd, "rm -d ",aldir," ",outdir.lcl," ",cutdir,"\n",sep='')
	cmd			<- paste(cmd, "\n#######################################################
#
# end: run haircut.pipeline
#
#######################################################\n",sep='')
	cmd
}
##--------------------------------------------------------------------------------------------------------
##	command line generator for 'haircut.cutstat.contigs.Rscript'
##--------------------------------------------------------------------------------------------------------
#' @export
cmd.haircut.cutstat<- function(indir, outdir, batch.n=NA, batch.id=NA, prog=PR.HAIRCUT.CUTSTAT )	
{
	cmd<- "\n#######################################################
# start: run haircut.cutstat.contigs.Rscript
#######################################################"
	cmd		<- paste(cmd, paste("\necho \'run ",prog,"\'\n",sep=''))
	cmd		<- paste(cmd, paste(prog,' -indir=',indir,' -outdir=',outdir, sep=''))
	if(!is.na(batch.n) & !is.na(batch.id))
		cmd	<- paste(cmd, ' -batch.n=',batch.n, ' -batch.id=',batch.id, sep='')
	cmd		<- paste('\n',cmd,paste("\necho \'end ",prog,"\'\n",sep=''))
	cmd		<- paste(cmd,"#######################################################
# end: run haircut.cutstat.contigs.Rscript
#######################################################\n",sep='')
	cmd
}
##--------------------------------------------------------------------------------------------------------
##	batch file wrapper
##	olli originally written 26-08-2014
##--------------------------------------------------------------------------------------------------------
#' @export
cmd.hpcwrapper<- function(cmd, hpcsys= cmd.hpcsys(), hpc.walltime=24, hpc.mem="1750mb", hpc.nproc=1, hpc.q=NA)
{
	wrap<- "#!/bin/sh"
	#hpcsys<- HPC.CX1.IMPERIAL
	if(hpcsys%in%c(HPC.CX1.IMPERIAL,'(none)'))
	{				
		tmp	<- paste("#PBS -l walltime=",hpc.walltime,":59:59,pcput=",hpc.walltime,":45:00",sep='')
		wrap<- paste(wrap, tmp, sep='\n')		
		tmp	<- paste("#PBS -l select=1:ncpus=",hpc.nproc,":mem=",hpc.mem,sep='')
		wrap<- paste(wrap, tmp, sep='\n')
		wrap<- paste(wrap, "#PBS -j oe", sep='\n')
		if(!is.na(hpc.q))
			wrap<- paste(wrap, paste("#PBS -q",hpc.q), sep='\n\n')
		wrap<- paste(wrap, HPC.CX1.IMPERIAL.LOAD, sep='\n')
	}
	else if(hpcsys=='debug')
		cat(paste("\ndetected no HPC system and no hpcwrapper generated, domain name is",hpcsys))
	else
		stop(paste("unknown hpc system with domain name",hpcsys))
	
	cmd<- lapply(seq_along(cmd),function(i){	paste(wrap,cmd[[i]],sep='\n')	})
	if(length(cmd)==1)
		cmd<- unlist(cmd)
	cmd	
}
##--------------------------------------------------------------------------------------------------------
##	batch file caller
##	olli originally written 26-08-2014
##--------------------------------------------------------------------------------------------------------
#' @export
cmd.hpccaller<- function(outdir, outfile, cmd)
{
	if( nchar( Sys.which("qsub") ) )
	{
		file	<- paste(outdir,'/',gsub(':','',outfile),'.qsub',sep='')
		cat(paste("\nwrite HPC script to",file,"\n"))
		cat(cmd,file=file)
		cmd		<- paste("qsub",file)
		cat( cmd )
		cat( system(cmd, intern=TRUE) )
		Sys.sleep(1)
	}
	else
	{
		file	<- paste(outdir,'/',gsub(':','',outfile),'.sh',sep='')
		cat(paste("\nwrite Shell script to\n",file,"\nStart this shell file manually\n"))
		cat(cmd,file=file)
		Sys.chmod(file, mode = "777")	
		Sys.sleep(1)
	}
	file
}


