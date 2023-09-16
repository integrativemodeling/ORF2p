These scripts demonstrate the use of [IMP](https://integrativemodeling.org/) in the modeling of the human LINE1 ORF2p protein using diverse types of data as described in XXXXX.

## List of files and directories:

- `dependencies\environment.yml`  : Files required to deploy a conda environment to reproduce the modeling and analysis

- `data`		            contains all relevant data

  - `ORF2p.fasta`  : Reference protein sequence 
  
  - `af_apo_final_lowpass_4.5.pdb`  : AlphaFold2 model of ORF2p aligned to the APO cryo-em density map
  
- `data_gmm\*`		  GMM representation of the selected fragments of the APO cryo-em density map
- `data_npc\negative_stain`		  2D Negative stain EM images from RNA:DNA and RNA datasets
- `data_npc\xls`		  Chemical cross-links and predicted restraints used for modeling

- `modeling`
  - `run_stub\sample_orf2p_mmcif.py`  Template script to rerun the modeling
  - `run_*`  Modeling output used in the publication

- `analysis`  Output from the [SAMPCON](https://github.com/salilab/imp-sampcon) analysis pipeline

- `validation`		                      Validation of modeling results with EM2D data that were not used for sampling
  - `ORF2p_em2d_validation.ipynb`  Jupyter notebook with the validation code and results
  - `validation\output`  Best-fitting models and various figures used in the publication

## Information

_Author(s)_: Arthur Zalevsky

_License_: [CC-BY-SA-4.0](https://creativecommons.org/licenses/by-sa/4.0/legalcode).
This work is freely available under the terms of the Creative Commons
Attribution-ShareAlike 4.0 International License.
