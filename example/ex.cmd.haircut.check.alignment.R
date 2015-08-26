#
#	run from command line 
#	this produces a command line string that can be run in UNIX alikes
#
\dontrun{
	
#DATA		<- SET THIS DIRECTORY
ind			<- paste(DATA, 'contigs_150408_wref', sep='/' )
cat(cmd.haircut.check.alignment(indir, indir))
}
