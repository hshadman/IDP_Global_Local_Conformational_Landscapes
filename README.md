
# Examine Global and Local Conformational Landscapes of Intrinsically Disordered Proteins or Polymers Using Polymer Physics Quantities

## This repository provides code to generate plots of instantaneous shape ratio (_R<sub>s</sub>_) and relative shape anisotropy (_RSA_). These shape parameters have previously been defined in this publication: https://www.cell.com/biophysj/fulltext/S0006-3495(24)00272-8  

For a given intrinsically disordered protein (IDP) or region (IDR) of a protein, the global (i.e. whole chain) conformational ensembles can be examined by generating a scatter plot of _R<sub>s</sub>_ against _RSA_ of the full chain. For the same chain, local conformational ensembles can be examined by using a moving/sliding window across the chain and monitoring _R<sub>s</sub>_ and _RSA_ for each subchain.

The below two jupyter notebook files provide code to generate plots of conformational ensembles.

1. _global_conf_landscape_map.ipynb_: Provides code with examples to generate scatter plots of global ensembles.
2. _subchain_code_with_examples.ipynb_: Provides code with examples to generate _R<sub>s</sub>_ and _RSA_ plots at the local (subchain) level

Two human IDR simulation files are used as examples, taken from https://www.nature.com/articles/s41586-023-07004-5.

_ete1.npy_ and _ete2.npy_ provide end-to-end distance values for the two IDR simulation trajectories.

_rg1.npy_ and _rg2.npy_ provide radius of gyration for the two IDR simulation trajectories.

_traj1.xtc_ and _traj2.xtc_ -> trajectory files for the two IDR simulations.

_top1.pdb_ and _top2.pdb_ -> topology files for the two IDR simulations.

_IDR_fasta_sequences.csv_ provides fasta sequences for all the human IDRs (published by https://www.nature.com/articles/s41586-023-07004-5).

_IDR_nu_values.csv_ provides ν values for all the human IDRs (published by https://www.nature.com/articles/s41586-023-07004-5).

We use the conformational ensembles of the Gaussian Walk (GW) polymer model as a reference for the ensembles of other proteins and polymers. The GW file used in these examples (reference_GW_chainlen_100_for_RSA.csv) can be downloaded from Zenodo. 

The code here can be easily adapted for any protein/polymer chain. 