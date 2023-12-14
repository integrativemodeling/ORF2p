# Files in this directory:

`README.md` - this file
`thumb_ctd_zn.pdb` - starting structure of the Thumb:CTD ORF2p fragment derived from AlphaFold prediction
`run_stub/` - template folder for a molecular dynamics simulation
`run_stub/prepare_md.bash` - modeling script
`run_*/` - output files
`run_*/md.tpr` - binary run files for production MD runs
`run_*/md_nojump_fit.xtc` - postprocessed trajectories with i) only protein ii) PBC treatment iii) superposition to starting frame 
`start.pdb` - pdb file with protein only that can be used for analysis and visualization of XTC trajectories
`index.ndx` - index file with various groups defined for simulation and analysis

# Running and MD simulation:
1. copy stub directory:
`cp -r run_stub run_X`
2. change directory
`cd run_X`
3. modify path to the GROMACS executable in `prepare_md.bash`
4. start the simulation
`./prepare_md.bash ../thumb_ctd_zn.pdb`
