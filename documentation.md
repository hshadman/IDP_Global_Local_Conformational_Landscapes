# PyHeteroMap Documentation

## Overview
`PyHeteroMap` is a Python class for **subchain and Gaussian walk (GW) reference analysis** of intrinsically disordered protein (IDR) trajectories. It enables users to compare coarse-grained or all-atom molecular dynamics trajectories to a Gaussian chain reference, compute shape descriptors, and generate subchain-level biophysical profiles.

The implementation is designed for reproducible and modular analyses of IDR trajectories, with support for trajectory-based calculations (`mdtraj`) and analytical Flory random coil (AFRC) comparisons.

---

## Repository Structure
```
github_subchain_code/
├── src/
│   └── pyheteromap/
│       ├── __init__.py
│       └── pyheteromap.py
├── examples/
│   └── example_usage.ipynb
├── tests/
│   ├── test_import.py
│   ├── test_requires_trajectory.py
│   ├── test_override_csv_path.py
│   └── test_default_gw_loads.py
├── reference_GW_chainlength_100.csv
├── README.md
├── documentation.md   ← (this file)
├── LICENSE
└── pyproject.toml
```

---

## Class: `PyHeteroMap`

### Purpose
Perform polymer-physics-based analyses of IDR trajectories:
- Compare radius of gyration (Rg), end-to-end distance (Ree), and shape anisotropy against Gaussian walk (GW) references.
- Compute Flory exponent (ν) following Tesei et al. (2024).
- Calculate local (subchain) compaction and shape metrics in 30-mer windows.
- Visualize shape ratio (Rs), relative shape anisotropy (RSA), and Rg/Rgθ trends.

---

### Initialization
```python
PyHeteroMap(
    seq_name: str,
    gw_reference_csv: str = "reference_GW_chainlength_100.csv",
    traj_file_dir: Optional[str] = None,
    prmtop_file_dir: Optional[str] = None,
    skip_frames: int = 10,
    afrc_returns_angstrom: bool = True
)
```
**Arguments:**
- `seq_name`: name of the sequence or system.
- `gw_reference_csv`: path to Gaussian chain reference CSV (default included in repo).
- `traj_file_dir`, `prmtop_file_dir`: optional trajectory and topology files to load.
- `skip_frames`: number of initial frames to skip from trajectory.
- `afrc_returns_angstrom`: whether AFRC returns are in Å.

**Behavior:**
- Initializes empty dataframes for trajectory- and GW-derived properties.
- If trajectory and topology paths are given, automatically loads the trajectory.

---

### Key Methods

#### 1. Trajectory and Reference Handling
| Method | Description |
|---------|--------------|
| `set_trajectory(traj_file, prmtop_file)` | Loads trajectory and topology (CG or all-atom). |
| `_load_gw_reference()` | Loads the default Gaussian chain CSV reference. |
| `reinitialize_gw_reference()` | Reloads or resets GW reference file. |
| `regenerate_GW_chain(chain_length, nosnaps, interval, mu, sigma)` | Generates a new Gaussian chain simulation in memory. |
| `export_gw_csv(out_path)` | Saves generated Gaussian chain data to CSV. |

#### 2. Core Analysis
| Method | Description |
|---------|--------------|
| `mod_RSA_Rs_compute_3dplot_from_seq_name()` | Computes and visualizes RSA vs Rs for a trajectory vs GW. |
| `calculate_nu_KLL_from_seq_name(start, end)` | Calculates Flory exponent ν using distance scaling. |
| `initialize_30mer_subchain(fasta_source, k_frac=3)` | Generates sliding-window subchain data and computes properties. |

#### 3. Visualization
| Method | Description |
|---------|--------------|
| `plot_subchain_RSA()` | Plot local ⟨RSA⟩ per residue. |
| `plot_subchain_Rs()` | Plot local ⟨Rₛ⟩ per residue. |
| `plot_subchain_Rg()` | Plot local ⟨Rg⟩ per residue. |
| `plot_subchain_Rg_over_Rgtheta()` | Plot ⟨Rg⟩/⟨Rgθ⟩ (α) per residue. |
| `plot_subchain_nu()` | Plot local Flory exponent ν per residue. |

---

### Important Internal Functions
| Function | Purpose |
|-----------|----------|
| `RSA_based_fC()` | Compute fC_shape overlap fraction between trajectory and GW ensembles. |
| `_gw_calculate_tensors()` | Compute gyration tensor elements for GW chain. |
| `_gw_diagonalize_tensor()` | Diagonalize tensor and compute RSA/asphericity/acylindricity. |
| `Rij()` | Estimate ν using pairwise distance scaling (Tesei et al. 2024). |
| `assign_properties()` | Assign residue colors (polar, apolar, positive, negative). |

---

### Outputs
- **DataFrames:** `gw_df`, `protein_df`, `_subchain_df` (depending on method used).
- **Plots:** Generated via Matplotlib; standard 2D and marginal density plots.
- **Numerical metrics:** ν, fC_shape, ⟨RSA⟩, ⟨Rg⟩, ⟨Rs⟩, ⟨α⟩.

---

### Example Usage
```python
from pyheteromap import PyHeteroMap

# Initialize and load trajectory
h = PyHeteroMap("example_IDR", traj_file_dir="traj.xtc", prmtop_file_dir="top.pdb")

# Compare against Gaussian chain
h.mod_RSA_Rs_compute_3dplot_from_seq_name()

# Compute and visualize subchain metrics
sub_df = h.initialize_30mer_subchain("IDR_fasta_sequences.csv")
h.plot_subchain_RSA(4, 3)
```

---

## Developer Notes
- **Python version:** 3.9.x recommended
- **Dependencies:** pinned in `pyproject.toml`
- **Build system:** Hatchling
- **Testing:** pytest-based unit tests located in `/tests`
- **License:** MIT

---

## Citation
If you use `PyHeteroMap` in a publication, please cite:
> Hossain S., *PyHeteroMap: Subchain and GW-reference analysis for IDR trajectories*, GitHub (2025).

