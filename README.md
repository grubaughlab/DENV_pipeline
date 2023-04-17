## Virus agnostic (but mostly dengue) analysis pipeline

Created by Verity Hill and Chrispin Chaguza, Grubaugh Lab

Still a bit beta - please feel free to have a go and report any issues you get! I need to run more tests before I tag a release.

This pipeline takes raw read data in the form of fastq files, maps them against provided bed files and then provides a series of outputs for further analysis including consensus sequences. IMPORTANT: the bed files must correspond to the wet lab protocol that you used and the reference sequence used to generate them otherwise the sequences generated will be incorrect. 

It calls input files as a virus type if it has more than 50% coverage of the reference genome provided.

If running on a server, it is highly recommended to run using screen/tmux or similar.

### Installation instructions

- Clone this repo
- conda env create -f environment.yml
- pip install .

If step 2 fails on a server because of a "bus error", then first run the command "salloc" to request more memory. If this also fails, I've found that mamba works well so if that's installed on your server give that a go


### Main inputs

As a minimum you need fastQ files to analyse. There are two ways you can provide these:

1. If you are running on the Yale cluster using their symlinks, simply provide the symlink emailed to you by YCRC (the second half of the link) using ``--symlink`` and the pipeline will deal with it for you.
2. Otherwise, must be separately in their own folder named by sample. In each file, must have the forward and reverse fastq files, defined by them containing "R1" and "R2" somewhere in the name and with the sample name at the start of the file name. See example input file for more information.


- For the second option, you can use ``--indir`` to indicate where the folders of samples are kept. Default is the same as the output directory (this is for when you already have a input/output folder where you're working)
- NB sample names are found by looping through directories in the input directory. The script ignores temporary and download folders, but will get confused if there are non-sample directories in there other than that.

To get consensus files for dengue, you don't need anything else. If you want to try other viruses you need some bed files to compare the sequencing data to:

- Here, we provide each of the four dengue virus serotypes as package data, and these are used as default. 
- Otherwise, please use the same format and provide the stem of the files in a "refs.txt" text file in the same folder. Use the ``--primer-directory`` option to provide the path to the directory.

### Main outputs

Specify the main output directory using ``--outdir``. Default is the generation of a new folder with today's date.

Within this, there will be:

1. results:
	- DENV.serotype.calls.final.tsv: Contains virus calls per sample which have more than 50% coverage
	- DENV.top.serotype.cals.all.samples.tsv: Contains all top calls per sample, regardless of coverage
	- summary.all.samples.tsv: contains information about all possible options provided per sample
	- alignment - contains alignments by virus type NB not to be used for phylogenetics because it is only rough for estimating coverage. If a trimmed bed file was provided then these are trimmed down, otherwise they are only untrimmed.
	- consensus - consensus sequences of the called virus for each sample
	- depth - text files of each position of the genome and their sequencing depth by sample
	- variants - contains a summary file of the number of variants by sample, and then a file for each sample containing additional information about variants

2.  Within results, there are some QC plots:

	- variant_plot: plots sample id against the number of variants (i.e. mutations compared to the reference genome) in between 20% and 80% of the reads. Useful for detecting co-infections
	- ct_plot: plots genome coverage against Ct value coloured by call.
	
	 To make the Ct plot, you must provide:
	 	- File containing the Ct values using ``--ct-file``
	 	- the name of the column containing Ct values with ``--ct-column``
	 	- the name of the column containing IDs which match the ID names on the fastq files/directories using ``--id-column``

3. Log_files:

	Log files for each step of the pipeline, named by snakemake rule and sample - mostly useful for debugging.


4. temporary_files (if option ``--temp`` is used):

	All files produced in the run of the pipeline that don't make it to results.
	This includes intermediate files produced by software, and bam/fasta/variant files for virus comparisons which did not meet the threshold for the sample to be called as that virus
	
	Mostly for debugging purposes
	
5. downloads (if option ``--download`` is used):

	The same as results, but without the bam files. 
	
	For download and storage in places with less data storage capacity (eg dropbox)
	




### All options

``--outdir`` location where files will be stored

``--indir`` directory containing samples. Each sample must be a folder with the forward and reverse runs in. Default is same as output directory.


``--primer-directory`` or ``-pd`` location where bed files etc for references are. Default is the dengue directory provided as package information

``--depth`` minimum depth to call consensus. Default is 20


``--temp`` keep temporary files

``--tempdir`` location of temporary files. Default is folder in output directory called "temporary_files"

``--download`` produce downloads directory (see above)

``--slurm`` parallelise on a cluster which uses slurm. The outer pipeline will run on the login node, but all actual analysis will be submitted to compute nodes.


``--overwrite`` delete old results files and make new ones in the run

``--verbose`` print more stuff to screen

``--help`` print help message to screen and quit


#### Yale specific options

``--symlink`` for use on Yale's HPC. Provide the second half of the link emailed by genomics core. 

    
    
