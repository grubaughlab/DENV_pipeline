## Virus agnostic (but mostly dengue) analysis pipeline

Created by Verity Hill and Chrispin Chaguza, Grubaugh Lab

Still a bit beta - please feel free to have a go and report any issues you get! I need to run more tests before I tag a release.

This pipeline takes raw read data in the form of fastq files, maps them against provided bed files and then provides a series of outputs for further analysis including consensus sequences. IMPORTANT: the bed files must correspond to the wet lab protocol that you used and the reference sequence used to generate them otherwise the sequences generated will be incorrect. 

It calls input files as a virus type if it has more than 50% coverage of the reference genome provided.

If running on a server, it is highly recommended to run using screen/tmux or similar.

# Installation instructions

1. Clone this repo by going to your command line interface and typing ```git clone https://github.com/grubaughlab/DENV_pipeline.git```
2. Navigate to your local copy of the repo by typing ```cd DENV_pipeline```
3. Create the conda environment by typing ```conda env create -f environment.yml``` You will need conda to be installed first, installation instructions found here: https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html
4. Activate the environment by typing ```conda activate analysis_env```
5. Install packages from the repo by typing ```pip install .```


NB If step 3 fails on a server because of a "bus error", then first run the command "salloc" to request more memory. If this also fails, I've found that mamba works well so if that's installed on your server give that a go

# Running the pipeline

```denv_pipeline -h``` will give you all of the usage instructions as follows:

```    
            DENVseq pipeline
                Version 0.1
        Chrispin Chaguza & Verity Hill
                Grubaugh Lab
    
    
usage: denv_pipeline [--config CONFIG] [--symlink SYMLINK] [--indir INDIR] [--outdir OUTDIR] [--primer-directory PRIMER_DIRECTORY] [--depth DEPTH] [--temp] [--tempdir TEMPDIR] [--download] [--slurm]
                     [--verbose] [--help] [--overwrite] [--ct-file CT_FILE] [--ct-column CT_COLUMN] [--id-column ID_COLUMN]

optional arguments:
  --config CONFIG       config file containing all relevant arguments
  --symlink SYMLINK     argument for generating symlinks
  --indir INDIR         directory containing samples. Each sample must be a folder with the forward and reverse runs in. Default is same as output directory
  --outdir OUTDIR       location where files will be stored.
  --primer-directory PRIMER_DIRECTORY, -pd PRIMER_DIRECTORY
                        location where bed files etc for references are
  --depth DEPTH         depth to map sequences to. Default=20
  --temp                keep intermediate files
  --tempdir TEMPDIR     where the temporary files go
  --download            make a folder without bam files for download
  --slurm               flag for if running on HPC with slurm
  --verbose, -v
  --help, -h
  --overwrite           overwrite current results
  --ct-file CT_FILE     to produce a plot of Ct against coverage, provide a csv file containing Ct information by sample
  --ct-column CT_COLUMN
                        Name of Ct column in Ct file for plot
  --id-column ID_COLUMN
                        Name of ID column in Ct file to make Ct plot
```
### Example commands

If you are using our bed files and reference sequences, the following command is the simplest way to run a sample:

```denv_pipeline --indir <input-directory>```

If you want to use your own bed files and reference sequences:

```denv_pipeline --indir <input-directory> --reference-directory <reference-directory>```

### Main inputs

As a minimum you need fastQ files to analyse. There are two ways you can provide these:

Option A (most people!): The fastq files must be separately in their own folder named by sample. In each file, must have the forward and reverse fastq files, defined by them containing "R1" and "R2" somewhere in the name and with the sample name at the start of the file name. See example input file for more information.

Option B (for Yale users) If you are running on the Yale cluster using their symlinks, simply provide the symlink emailed to you by YCRC (the second half of the link) using ``--symlink`` and the pipeline will deal with it for you.

- For the second option, you can use ``--indir`` to indicate where the folders of samples are kept. Default is the same as the output directory (this is for when you already have a input/output folder where you're working)
- NB sample names are found by looping through directories in the input directory. The script ignores temporary and download folders, but will get confused if there are non-sample directories in there other than that.

To get consensus files for dengue if you are using our sequencing protocol (https://www.protocols.io/view/dengueseq-a-pan-serotype-whole-genome-amplicon-seq-kqdg39xxeg25/v2), you don't need anything else. 

If you want to try other viruses, or use your own reference and bed files:

- Look at our directory "DENV_primers_and_refs" for formatting file names etc
- Provide the stem of each file in the "refs.txt" tesxt file in the same folder
- Use the ``--reference-directory`` option to provide the path to the directory. 


You can also provide all of your arguments using a config file. This is specified using the ```--config``` option. An template can be found in the main repository. Any command line arguments you specify will overwrite config file arguments.


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

    
    
