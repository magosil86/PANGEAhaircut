#
#	run from command line 
#	this produces a command line string that can be run in UNIX alikes
#	only one raw IVA infile specified
#
\dontrun{
	 
in.raw	<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2015/2015_PANGEA_haircut/contigs_151026_unaligned_raw/78268.assembly_contigs_hit_ref_hiv.fasta'		
cat(cmd.haircut.pipeline(in.raw))
}
#
#	cut versions of IVA infile are also available and out file is also specified
#
\dontrun{

in.raw	<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2015/2015_PANGEA_haircut/contigs_151026_unaligned_raw/78270.assembly_contigs_hit_ref_hiv.fasta'
in.cut	<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2015/2015_PANGEA_haircut/contigs_151026_unaligned_cut/78270.assembly_contigs_hit_ref_hiv_cut.fasta'
out		<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2015/2015_PANGEA_haircut/78270_nohair.fasta'		
cat(cmd.haircut.pipeline(in.raw, in.cut, out))
}
