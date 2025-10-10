![Tests](https://github.com/hshadman/IDP_Global_Local_Conformational_Landscapes/actions/workflows/tests.yml/badge.svg)

# **_PyHeteroMap_:** Resolving Local and Global Conformational Heterogeneity of the Human Intrinsically Disordered Proteome Using (RSA, _R<sub>s</sub>_) Maps

---

## Overview

This repository provides code to generate plots of instantaneous shape ratio (_R<sub>s</sub>_) against relative shape anisotropy (RSA).  
These shape parameters have previously been defined in this publication:  
[Biophysical Journal (2024)](https://www.cell.com/biophysj/fulltext/S0006-3495(24)00272-8)

---

## Description

For a given intrinsically disordered protein (IDP) or region (IDR) of a protein:

- **Global ensembles:** can be examined by generating a scatter plot of _R<sub>s</sub>_ against RSA of the full chain.  
- **Local ensembles:** can be examined using a moving/sliding window across the chain to monitor _R<sub>s</sub>_ and RSA for each subchain.

---

## Features

For a given peptide/protein trajectory, **_PyHeteroMap_** can generate:

1. **(RSA, _R<sub>s</sub>_) scatter plots**  
   - Compare a trajectory against a Gaussian Walk (GW) reference.  
   - Compute metrics such as the f<sub>C_shape</sub> score that quantify its conformational diversity.

2. **Local polymer property plots**  
   - Display how polymer properties such as ⟨RSA⟩, ⟨Rₛ⟩, and others vary at the subchain level.

---

## Applications

**_PyHeteroMap_** mainly targets the ~28,000 human IDR simulations published by  
[Tesei et al. (2024), *Nature*](https://www.nature.com/articles/s41586-023-07004-5)

- Two such human IDR simulations are included as examples in the `examples/` folder.  
- The source code can easily be adapted to read other types of protein trajectories.  
- Trajectory analysis is performed using **MDTraj**.  
- **_PyHeteroMap_** can also simulate new **Gaussian Walk (GW)** chains of any chain length and number of snapshots.

---

## Included Data

- **Tesei_2024_IDR-ome_fasta_sequences.csv**  
  Provides fasta sequences for all human IDRs published by Tesei et al. (2024).

- **reference_GW_chainlength_100.csv**  
  Gaussian Walk (GW) reference ensemble used in the examples.  
  Located at:  src/pyheteromap/reference_GW_chainlength_100.csv


---

## References

- [Biophysical Journal (2024)](https://www.cell.com/biophysj/fulltext/S0006-3495(24)00272-8)  
- [Tesei et al. (2024), *Nature*](https://www.nature.com/articles/s41586-023-07004-5)

---

