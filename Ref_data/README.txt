If you want to make a new reference.....
Copy the fasta into this folder and generate a mock gtf file corresponding to the gene lengths.
Make a new directory for the reference:
	mkdir Ref3

Make an interactive session on Maestro.
	salloc -c 6

Load STAR 
	module load STAR/2.7.8a

	STAR --runMode genomeGenerate --runThreadN 6 --genomeDir <path/to/new/directory/Ref3> --genomeFastaFiles <./Ref3.fasta> --genomeSAindexNbases 5 --sjdbGTFfile <./Ref3_custom.gtf>
