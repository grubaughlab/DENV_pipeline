## Repo for denv pipeline development

Scripts originally came from Chrispin Chaguza

To do:
- Conda environment
- setup.py
- Actual snakemake
- license

Pipeline steps:
- Detect the run number
- make dirs
- Run the scripts
- sort the files into the alignment/consensus etc
- make QC plots (x2)
- run nextstrain?

End of pip

To install:

- Clone this repo
- conda env create -f environment.yml
- python setup.py install 

If step 2 fails on a server because of a "bus error", then first run the command "salloc" to request more memory
