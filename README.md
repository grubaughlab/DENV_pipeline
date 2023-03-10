## Repo for denv pipeline development

Scripts originally came from Chrispin Chaguza


Pipeline steps:
- Detect the run number
- make dirs
- Run the scripts
- sort the files into the alignment/consensus etc
- make QC plots (x2)

End of pipeline

To install:

- Clone this repo
- conda env create -f environment.yml
- pip install .

If step 2 fails on a server because of a "bus error", then first run the command "salloc" to request more memory
