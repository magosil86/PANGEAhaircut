PR.PACKAGE					<- "PANGEAhaircut"
PR.STARTME					<- system.file(package=PR.PACKAGE, "misc", "PANGEAhaircut.startme.R")
#PR.STARTME					<- '/Users/Oliver/git/HPTN071sim/source/rPANGEAHIVsim/misc/rPANGEAHIV.startme.R'
#PR.STARTME					<- '/work/or105/libs/HPTN071sim/source/rPANGEAHIVsim/misc/rPANGEAHIV.startme.R'
PR.VARIOUS					<- paste(PR.STARTME," -exe=VARIOUS",sep='')
PR.HAIRCUT.CALL				<- paste(PR.STARTME," -exe=HAIRCUT.CALL",sep='')

HPC.MPIRUN					<- {tmp<- c("mpirun","mpiexec"); names(tmp)<- c("debug","cx1.hpc.ic.ac.uk"); tmp}
HPC.CX1.IMPERIAL			<- "cx1.hpc.ic.ac.uk"		#this is set to system('domainname',intern=T) for the hpc cluster of choice
HPC.MEM						<- "1750mb"
HPC.CX1.IMPERIAL.LOAD		<- "module load intel-suite mpi mafft/7 R/3.2.0"


######################################################################################
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
##	command line generator for 'prog.HPTN071.input.parser.v4'
##	olli originally written 16-08-2015
##--------------------------------------------------------------------------------------------------------
cmd.haircut.call<- function(indir.st, indir.al, outdir, mfile, trainfile=NA, batch.n=NA, batch.id=NA, prog=PR.HAIRCUT.CALL )	
{
	cmd<- "#######################################################
# start: run haircutprog.get.call.for.PNG_ID
#######################################################"
	cmd		<- paste(cmd, paste("\necho \'run ",prog,"\'\n",sep=''))
	cmd		<- paste(cmd, paste(prog,' -mfile=', mfile,' -indir.st=',indir.st,' -indir.al=',indir.al,' -outdir=',outdir, sep=''))
	if(!is.na(trainfile))
		cmd	<- paste(cmd, ' -trainfile=',trainfile, sep='')
	if(!is.na(batch.n) & !is.na(batch.id))
		cmd	<- paste(cmd, ' -batch.n=',batch.n, ' -batch.id=',batch.id, sep='')
	cmd		<- paste('\n',cmd,paste("echo \'end ",prog,"\'\n",sep=''))
	cmd		<- paste(cmd,"#######################################################
# end: run haircutprog.get.call.for.PNG_ID
#######################################################\n",sep='')
	cmd
}
##--------------------------------------------------------------------------------------------------------
##	batch file wrapper
##	olli originally written 26-08-2014
##--------------------------------------------------------------------------------------------------------
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


