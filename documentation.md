# PyHeteroMap Documentation

## Overview
**PyHeteroMap** is a Python package designed for structural and polymer-physics-based analysis of Intrinsically Disordered Region (IDR) trajectories. It provides functionality for comparing IDR simulations with Gaussian chain (GW) reference ensembles, subchain-level analyses, and extraction of key polymer shape metrics such as the radius of gyration (Rg), end-to-end distance (Ree), asphericity, acylindricity, and ν (Flory exponent).

The package integrates data from coarse-grained (CG) and all-atom simulations using MDTraj and provides utilities to compute metrics, compare to reference Gaussian chains, and visualize per-residue subchain properties in a reproducible framework.

---

## Repository Structure

```
github_subchain_code/
│
├── src/
│   └── pyheteromap/
│       ├── __init__.py
│       ├── pyheteromap.py          # Main PyHeteroMap class implementation
│       └── reference_GW_chainlength_100.csv  # Default Gaussian reference data
│
├── examples/
│   └── example_usage.ipynb         # Demonstration Jupyter notebook
│
├── tests/
│   ├── test_import.py
│   ├── test_default_gw_loads.py
│   ├── test_override_csv_path.py
│   └── test_requires_trajectory.py
│
├── LICENSE
├── README.md
├── documentation.md                # (This file)
└── pyproject.toml
```

---

## PyHeteroMap Class — Detailed Description

### Class Initialization
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
Initializes the PyHeteroMap analysis object and optionally loads molecular dynamics trajectories.

**Inputs:**
- `seq_name`: (str) Identifier for the protein or sequence analyzed.
- `gw_reference_csv`: (str) Path to Gaussian reference chain data (default distributed CSV).
- `traj_file_dir`: (str, optional) Path to trajectory file (e.g., `.xtc`).
- `prmtop_file_dir`: (str, optional) Path to topology file (e.g., `.pdb`).
- `skip_frames`: (int) Number of frames to skip from the start of trajectory.
- `afrc_returns_angstrom`: (bool) If True, AnalyticalFRC results are returned in Ångström units.

**Behavior:**  
If trajectory paths are provided, initializes MDTraj objects, sets up internal trajectory slices, and prepares reference data from the Gaussian chain CSV.

---

### Gaussian Chain (GW) Reference Handling

#### `_load_gw_reference()`  
Loads the Gaussian reference chain from the specified CSV. Raises an error if the file is missing or empty.

#### `reinitialize_gw_reference()`  
Reloads the reference CSV into memory, replacing the current Gaussian chain DataFrame.

#### `regenerate_GW_chain(chain_length, nosnaps, interval=1, mu=0, sigma=1)`  
Generates a new synthetic Gaussian chain ensemble and stores it as `self.gw_df`.

#### `export_gw_csv(out_path)`  
Exports the loaded or regenerated Gaussian chain data to a CSV file.

---

### Trajectory Management

#### `set_trajectory(traj_file_dir, prmtop_file_dir)`  
Loads a molecular dynamics trajectory using MDTraj, computes end-to-end distances, sets up subchain boundaries, and verifies topology consistency.

---

### Subchain-Level Analysis

#### `initialize_30mer_subchain(fasta_source, k_frac=3)`  
Extracts subchains (default 30 residues) from an IDR trajectory and computes per-subchain properties (Rg, ν, RSA, etc.). Supports input via FASTA string, CSV file, or FASTA file.

#### Subchain Plot Functions  
All subchain plots require prior execution of `initialize_30mer_subchain()`.

- `plot_subchain_RSA()` – plots ⟨RSA⟩ (relative shape anisotropy) per residue.  
- `plot_subchain_Rs()` – plots ⟨Rₛ⟩ (shape ratio) per residue.  
- `plot_subchain_Rg()` – plots ⟨R_g⟩ per residue.  
- `plot_subchain_Rg_over_Rgtheta()` – plots ⟨R_g / R_gθ⟩ ratio per residue.  
- `plot_subchain_nu()` – plots ⟨ν⟩ (Flory scaling exponent) per residue.

---

### Protein vs. Gaussian Reference Comparison

#### `mod_RSA_Rs_compute_3dplot_from_seq_name(provided_color='magenta')`  
Computes and visualizes RSA–Rₛ relationships for the protein vs. Gaussian chain ensemble.

#### `mod_RSA_Rs_protein_3dplot_against_GW(protein_var, protein_label, second_obj, provided_color)`  
Core comparison and plotting method combining experimental data (trajectory-derived) and Gaussian reference distributions.

#### `RSA_based_fC(...)`  
Computes the overlap fraction (fC_shape) between the Gaussian ensemble and trajectory-derived conformational space.

---

### Calculation and Utility Methods

#### `calculate_nu_KLL_from_seq_name(start_residue, end_residue)`  
Calculates ν (Flory scaling exponent) following Tesei et al. 2024, fitting ⟨Rij⟩ = R₀·i^ν.

#### `Rij(traj)` (static)  
Pairwise distance calculation and scaling exponent estimation.

#### `assign_properties(residue)` (static)  
Assigns colors and categories (polar, apolar, charged) for visualization.

---

## Outputs

After running `PyHeteroMap`, key results are stored in instance variables:

| Attribute | Description |
|------------|--------------|
| `gw_df` | Gaussian reference DataFrame (either loaded or regenerated). |
| `protein_df` | Frame-wise protein trajectory data with shape metrics. |
| `_traj_full` | Full MDTraj trajectory. |
| `_traj` | Processed trajectory (after skipping frames). |
| `_subchain_df` | Subchain-level data including Rg, ν, RSA, etc. |
| `fC_value` | Fractional overlap between GW and protein ensembles. |
| `grid_protein`, `grid_GW` | 2D histogram data used for RSA–Rₛ comparison plots. |

Generated plots include:  
- RSA–Rₛ 2D scatter and histograms  
- Per-residue ⟨Rg⟩, ⟨RSA⟩, ⟨Rₛ⟩, ⟨ν⟩  
- Subchain metrics colored by amino acid type  

All figures are returned via Matplotlib and can be further customized.

---

## Example Usage

```python
from pyheteromap import PyHeteroMap

# Initialize
phm = PyHeteroMap("IDR_Example")

# Load trajectory
phm.set_trajectory("traj1.xtc", "top1.pdb")

# Initialize subchain analysis
phm.initialize_30mer_subchain("IDR_fasta_sequences.csv")

# Generate and visualize results
phm.plot_subchain_RSA(6, 4)
phm.plot_subchain_Rg(6, 4)
phm.mod_RSA_Rs_compute_3dplot_from_seq_name(provided_color="magenta")
```

---

## Developer Notes

- **Unit consistency:** MDTraj returns distances in **nm**, AFRC reference in **Å** (scaling handled internally).  
- **Trajectory requirement:** All methods depending on trajectory data must follow a successful call to `set_trajectory()`.  
- **Performance:** For long IDRs (>300 residues), consider increasing `skip_frames` or sampling fewer subchains.  
- **File organization:** Gaussian reference CSV is packaged and accessed relative to the installed module directory.  
- **Tests:** Basic smoke and data integrity tests are located under `tests/`. Run with `pytest -q`.

---
