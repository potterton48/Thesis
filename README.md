# Thesis
These code files were used to produce results in the thesis "Computational methods that predict residence times of GPCR ligands" by Andrew Potterton.

Most of the code files will not work if they are run on other machines, they were specifically developed for the file setup/databases on my machine. The main purpose of sharing the code is for being open how the data was produced. However, one may adapt any parts of the code for use in one's own project. The python libraries used in the scripts are listed at the top of each file.

### Summary of each file:
- Kinetic_data.py (Python3) is a sectioned piece of code (you should run each section using an IDE that supports sectioning - like Spyder). This script follows the steps needed to setup molecular dynamics for GPCR-ligand complexes. Some of the steps automatically run for all receptors, and ligands. Other steps (because of their computational cost) run for all ligands of a single receptor.
- RAMD_Collect_Data.py (Python3) automatically collects tau-RAMD log data, and determines whether extra replicas/ensembles need to be submitted. The tau-RAMD protocol is detailed here (https://doi.org/10.1021/acs.jctc.8b00230). 


All code has a MIT licence.
