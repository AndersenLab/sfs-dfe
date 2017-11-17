


config_dat=Channel.fromPath('directory_config.dat')

spectra_QX1211 = Channel.fromPath("QX1211/*.sfs").combine(['QX1211'])
//spectra_XZ1516 = Channel.fromPath("QX1211/*.sfs").map({ ["QX1211", it]})

spectra_QX1211.into { spectra_in }

// model pop size
spectra_in.combine([0, 1]) // conpop
	   	  .combine([0]) // selmode
	   	  .combine([file('directory_config.dat')])
	   	  .into { spectra }

process run_multi_dfe {
	
	echo true

	tag { "${outgroup} - ${sfs} - ${conpop} - ${selmode}" }

	input:
		set file(sfs), val(outgroup), val(conpop), val(selmode), file('directory_config.dat') from spectra

	output:
		//set file("${sfs}.expSFS.out"), val(outgroup), val(conpop), val(selmode) into spectra
		set val(sfs), file("${sfs}.MAXL.out"), val(outgroup), val(conpop), val(selmode) into parse_maxl
		set val(sfs), file("${sfs}.obs-exp.neu-sel.SFS.out"), val(outgroup), val(conpop), val(selmode) into parse_sfs_file

	script:
		if (selmode == 0 || selmode == 1) {
			nspikes = "-nspikes 5"
			ranrep = "-ranrep 10"
		} else {
			nspikes = ""
			ranrep = ""
		}


	"""
		MultiDFE -conpop ${conpop} \\
				 -sfsfold 0 \\
				 -selmode ${selmode} \\
				 ${nspikes} \\
				 ${ranrep} \\
				 -file ${sfs}
	"""
}

process parse_maxl_proc {
	
	executor 'local'

	input:
		set val(sfs), file("${sfs}.MAXL.out"), val(outgroup), val(conpop), val(selmode) from parse_maxl

	output:
		file("maxml_out.tsv") into maxml_out

	tag { "${outgroup} - ${sfs} - ${conpop} - ${selmode}" }

	"""
	#!/usr/bin/env python 
	import sys
	from collections import OrderedDict

	# Open MAXL file
	with open("${sfs}.MAXL.out", 'r') as f:
		maxl = f.read().strip()
	maxl = OrderedDict([x.split(":") for x in maxl.split("\t")])
	maxl.update({"outgroup": "${outgroup}",
				 "conpop": "${conpop}",
				 "selmode": "${selmode}"})
	with open('maxml_out.tsv', 'w') as f:
		f.write('\\t'.join(maxl.keys()) + "\\n")
		f.write('\\t'.join(maxl.values()))
	"""

}


process merge_maxl_out {

	publishDir 'results/', mode: 'copy'

	input:
		file("maxml_out*.tsv") from maxml_out.toSortedList()

	output:
		file("MAXL.out.tsv")

	"""
	#!/usr/bin/env Rscript --vanilla
	library(tidyverse)

	readr::write_tsv(
		dplyr::bind_rows(
			lapply(list.files(pattern='*.tsv'),
				   function(x) {
				   		readr::read_tsv(x)
				   }
			)),
		path='MAXL.out.tsv'
	)

	"""
}


process parse_sfs_file {

	executor 'local'

	input:
		set val(sfs), file("${sfs}.obs-exp.neu-sel.SFS.out"), val(outgroup), val(conpop), val(selmode) from parse_sfs_file

	output:
		file("sfs_out.txt") into sfs_out

	"""
	cat "${sfs}.obs-exp.neu-sel.SFS.out" | \
	awk 'NR == 1 { print \$0 "sfs\\toutgroup\\tconpop\\tselmode" } NR > 1 && \$0 != "" { print \$0 "\\t${sfs}\\t${outgroup}\\t${conpop}\\t${selmode}" }' > sfs_out.txt
	"""
}

process merge_sfs_file {

	publishDir 'results/', mode: 'copy'

	executor 'local'

	input:
		file("sfs_out*.txt") from sfs_out.toSortedList()

	output:
		file("obs-exp.neu-sel.SFS.out.tsv")

	"""
	cat *.txt | awk 'NR == 1 { print } \$0 !~ "selmode" { print }' > obs-exp.neu-sel.SFS.out.tsv
	"""

}


