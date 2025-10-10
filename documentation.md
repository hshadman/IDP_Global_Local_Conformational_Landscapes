# PyHeteroMap Documentation

- [Overview](#overview)
- [Repository Structure](#repository-structure)
- [Class Description](#pyheteromap---a-python-class)
- [Outputs](#outputs)
- [Developer Notes](#developer-notes)

## Overview
**PyHeteroMap** is a python package that analyzes single protein/peptide MD simulation trajectories. [Tesei et al. (2024), *Nature*](https://www.nature.com/articles/s41586-023-07004-5) published simulations of 28058 Intrinsically Disordered Regions (IDRs) from the human proteome. **PyHeteroMap** is a tool that helps analyze their global and local conformational landscapes, directly from trajectory. However, the source code can be minimally adjusted to analyze other types of trajectories. It mainly has three functionalities. It can:

1. Generate a map of the conformational landscape of a given IDR/peptide chain in the form of an (RSA, _R<sub>s</sub>_) scatter plot against that of a GW reference, directly from trajectory. The GW is a polymer chain model that provides a reference landscape for other proteins/polymers. Data for GW is already provided.

2. Generate subchain plots (see examples) illustrating how RSA, _R<sub>s</sub>_, ν and other polymer properties vary along all moving windows (subchains) of the IDR/peptide chain, given its trajectory and fasta sequence. 

3. Simulate a Gaussian Walk (GW) polymer model chain of any chain length and any number of snapshots. Note: extremely high chain lengths or snapshot counts might require cautious interpretation. 

The package uses MDTraj to analyze trajectories. The examples provided are two human IDR trajectories (Tesei et al. 2024).

The plots generated are not saved automatically, those need to be saved by the user. The source code behind PyHeteroMap can be customized as needed. 

---

## Repository Structure

```
github_subchain_code/
│
├── src/
│   └── pyheteromap/
│       ├── __init__.py
│       ├── pyheteromap.py          # The Main PyHeteroMap class 
│       ├── reference_GW_chainlength_100.csv  # Default Gaussian reference data
│       └── Tesei_2024_IDR-ome_fasta_sequences.csv  # Human IDR-ome fasta sequences (Tesei et al.)
│
├── examples/
│   ├── example_usage.ipynb         # Demonstration Jupyter notebook
│   └── test_data/                  # Example trajectory and topology files (two separate IDRs from human IDR-ome data, Tesei et al. 2024).
│       ├── traj1.xtc
│       ├── traj2.xtc
│       ├── top1.pdb
│       ├── top2.pdb    
│       └── test_idr.csv
├── tests/
│   ├── test_import.py
│   ├── test_default_gw_loads.py
│   ├── test_override_csv_path.py
│   ├── test_requires_trajectory.py
│   └── conftest.py
│
├── LICENSE
├── README.md
├── documentation.md                # (This file)
└── pyproject.toml
```

---

## PyHeteroMap - a python class

### 1. Class Initialization
```python
PyHeteroMap(
    seq_name,
    gw_reference_csv="reference_GW_chainlength_100.csv",
    traj_file_dir=None,
    prmtop_file_dir=None,
    skip_frames=10,
    afrc_returns_angstrom=True
)
```

**Purpose:**  
Initializes the PyHeteroMap object and optionally loads molecular dynamics (MD) trajectories.

**Inputs:**
- `seq_name`: (str) Identifier for the protein/peptide sequence analyzed.
- `gw_reference_csv`: (str) Path to GW reference chain data (default: reference_GW_chainlength_100.csv).
- `traj_file_dir`: (str, optional) Path to trajectory file (e.g., `.xtc`).
- `prmtop_file_dir`: (str, optional) Path to topology file (e.g., `.pdb`).
- `skip_frames`: (int) Number of frames to skip from the start of trajectory.
- `afrc_returns_angstrom`: (bool) If True, AnalyticalFRC (AFRC) results are returned in Ångström units.

**Behavior:**  
If trajectory and topology paths are provided, initializes MDTraj objects.

---

### 2. Gaussian Walk (GW) Reference Handling

#### `_load_gw_reference()`  
Loads the GW reference chain from the specified CSV. Raises an error if the file is missing or empty.

#### `reinitialize_gw_reference()`  
Reloads the reference CSV into memory, replacing the current GW dataframe.

#### `regenerate_GW_chain(chain_length, nosnaps, interval=1, mu=0, sigma=1)`  
Simulates a new GW chain, with number of monomers =  chain_length, number of snapshots = nosnaps, and stores it as `self.gw_df`. Interval specifies snapshots to skip. For a GW chain, the displacement of one monomer to the next (in the x-, y-, and z-directions) is randomly picked from a gaussian distribution with mean = mu and standard deviation = sigma. By default mu is 0 and sigma is 1. 

#### `export_gw_csv(out_path)`  
Exports the loaded or regenerated GW data to a CSV file to out_path. out_path is a string. 


The GW data is stored in the attribute gw_df. Each row is a snapshot of the GW simulation. For the default GW reference file, the columns are provided:

- `ratio` – instantaneous shape ratio of the chain at that snapshot.
- `chain_length` – number of monomers in the GW chain.
- `RSA` – RSA of the GW chain at that snapshot.

However, when a GW chain is regenerated, additional columns are provided, including radius of gyration (Rg2 is radius of gyration squared), end-to-end distance (Rend2 is end-to-end distance squared) and (x, y, z) coordinates of every monomer of the chain at every snapshot. 

---

### 3. Trajectory Management

#### `set_trajectory(traj_file_dir, prmtop_file_dir)`  
Loads a molecular dynamics (MD) trajectory using MDTraj, computes end-to-end distances and verifies topology consistency. If skip_frames was not None, skips the first few frames as specified by skip_frames (at the init stage).

---

### 4. Subchain-Level Analysis

#### `initialize_30mer_subchain(fasta_source, k_frac=3)`  
Extracts subchains of the peptide/protein (from its trajectory) and computes per-subchain properties (RSA, _R<sub>s</sub>_, ν, etc.). Supports input via FASTA string, CSV file, or FASTA file. A sliding/moving window is used to select subchains. For peptides/proteins <=60 residues, the moving window size is 1/k_frac of the number of residues. For >60 residues, the moving window size is fixed at 30 residues. A fasta_source (containing fasta residue sequence of the peptide) must be provided. It can be a plain sequence string, .fasta or .fa file, or extracted from the Tesei_2024_IDR-ome_fasta_sequences.csv file (provided with the package). If fasta_source=Tesei_2024_IDR-ome_fasta_sequences.csv, this is for the human IDR-ome data (Tesei et al. 2024), the sequence is obtained automatically from the provided seq_name. 

#### Subchain Plot Functions  
NOTE: All subchain plots (below) require prior execution of `initialize_30mer_subchain()`. If using fasta_source = Tesei_2024_IDR-ome_fasta_sequences.csv, the path to the Tesei_2024_IDR-ome_fasta_sequences.csv file should be noted.

All subchain plots will show the corresponding global (whole chain) polymer property for the peptide/protein in the form of a dotted horizontal line or text at the top of the plot. The x-axis shows the mid-residue of the subchain. Gray region shows standard deviation (or error in the case of ν). Y-axis shows local values of the polymer property. The subchain plots are colored by amino acid residue type. 

These are the subchain plots possible:

- `plot_subchain_RSA()` – plots ⟨RSA⟩ where RSA is the relative shape anisotropy.  
- `plot_subchain_Rs()` – plots ⟨Rₛ⟩ (instantaneous shape ratio).  
- `plot_subchain_Rg()` – plots $⟨R_g⟩$ (radius of gyration).  
- `plot_subchain_Rg_over_Rgtheta()` – plots $⟨R_g / ⟨R_g^\theta⟩⟩$ where $⟨R_g^\theta⟩$ is obtained from the AFRC package (Alston et al. 2023). $⟨R_g^\theta⟩$ is the radius of gyration of a peptide chain if it behaved as an ideal chain, computed directly from sequence.  
- `plot_subchain_nu()` – plots ν (Flory scaling exponent), using formula used by Tesei et al. 2024.

The output is stored in the attribute _subchain_df. The columns are:

- `fasta_sequence` – fasta sequence of the subchain.
- `start_res_zeroindex` – start index of the subchain (zero index format).
- `end_res_zeroindex` – end index of the subchain (zero index format).
- `full_protein_fasta` – fasta sequence of the full IDR/chain.
- `full_protein_rgyr` – mean radius of gyration of the full IDR/chain.
- `full_protein_ratio` – mean instantaneous shape ratio of the full IDR/chain.
- `full_protein_nu_recompute` – ν (Flory scaling exponent) of the full IDR/chain.
- `full_protein_nu_recompute_err` – ν (Flory scaling exponent) error of the full IDR/chain.
- `full_protein_RSA` – mean RSA of the full IDR/chain.
- `full_protein_rg_rg_theta_mean` – $⟨R_g / ⟨R_g^\theta⟩⟩$ of the full IDR/chain (note use of mean symbol).
- `seq_name` – seq_name provided.
- `Rg/nm` – instantaneous Rg of the subchain.
- `nu_recompute` – ν (Flory scaling exponent) of the subchain.
- `nu_recompute_err` – ν (Flory scaling exponent) error of the subchain.
- `ete` – instantaneous end-to-end distance of the subchain.
- `inst_ratio` – instantaneous shape ratio of the subchain.
- `RSA` – instantaneous RSA of the subchain.
- `AFRC_mean_rg_theta` – $⟨R_g^\theta⟩$ of the subchain (note use of mean symbol).
- `rg_rg_theta_mean` – $R_g / ⟨R_g^\theta⟩$ of the subchain (note use of mean symbol).
- `mid_residue` – the mid residue of the subchain.






---






### 5. (RSA, _Rₛ_) Scatter Plots of Peptide and GW Reference

#### `mod_RSA_Rs_compute_3dplot_from_seq_name(provided_color='magenta')`  
Computes and visualizes (RSA, _Rₛ_) scatter plots of the protein and the GW reference. Calculates and displays ν and computes f<sub>C_shape</sub> scores. The f<sub>C_shape</sub> score computes how many of GW points are close to at least one protein/peptide point on the (RSA, _Rₛ_) scatter plot. It is a quantitative measure of the conformational diversity of that protein/peptide. 

#### `mod_RSA_Rs_protein_3dplot_against_GW(protein_var, protein_label, second_obj, provided_color)`  
Plots and displays the (RSA, _Rₛ_) scatter plot and associated quantities. 

#### `mod_RSA_Rs_compute_3dplot_from_csv(self, csv_path,label, rsa_col='RSA', ratio_col='ratio',provided_color='magenta')`  
Plots and displays the (RSA, _Rₛ_) scatter plot and associated quantities directly from a csv file, no trajectory required. rsa_col is the column header for the RSA data. ratio_col is the column header for the instantaneous shape ratio (_Rₛ_) data. csv_path is the path to the csv file.


#### `RSA_based_fC(...)`  
Computes the f<sub>C_shape</sub> score. 

---

## Outputs

After running `PyHeteroMap`, key results are stored in instance variables:

| Attribute | Description |
|------------|--------------|
| `gw_df` | GW reference DataFrame (either loaded or regenerated). |
| `protein_df` | Frame-wise protein trajectory data with polymer properties. |
| `_traj_full` | Full MDTraj trajectory. |
| `_traj` | Processed trajectory (after skipping frames). |
| `_subchain_df` | Subchain-level data including $R_g$, ν, RSA, etc. |
| `fC_value` | f<sub>C_shape</sub> score: a quantitative measurement of the conformational diversity of the protein/peptide. |


All figures are returned via Matplotlib and can be further customized by adjusting the source code.

---


---

## Developer Notes

- **Unit consistency:** MDTraj returns distances in **nm**, AFRC reference in **Å** (scaling handled internally).  
- **Trajectory requirement:** All methods depending on trajectory data must follow a successful call to `set_trajectory()`.  
- **Performance:** For lengthy trajectories, consider increasing skip_frames.  
- **File organization:** GW (reference) CSV is packaged and accessed relative to the installed module directory.  
- **Tests:** Basic smoke and data integrity tests are located under `tests/`. Run with `pytest -q`.

---
