[![DOI](https://zenodo.org/badge/744492651.svg)](https://zenodo.org/doi/10.5281/zenodo.10681431)

# _PyConforMap_: Draw pretty maps of your polymer or disordered protein conformational ensembles!

## This repository provides an easy-to-implement python module called _PyConforMap_ that generates scatter plots of instantaneous shape ratio (_R<sub>s</sub>_) against relative radius of gyration (_R<sub>g</sub>/R<sub>g</sub><sup>mean</sup>_).  

**PLEASE READ ALL DOCUMENTATION**

There are two main main metrics: the relative radius of gyration (_R<sub>g</sub>/R<sub>g</sub><sup>mean</sup>_) and the instantaneous shape ratio (_R<sub>s</sub>_). _R<sub>s</sub>_ is computed as _R<sub>s</sub>_ = _R<sub>ee</sub><sup>2</sup>/R<sub>g</sub><sup>2</sup>_ where _R<sub>ee</sub>_ and _R<sub>g</sub>_ are (instantaneous) end-to-end distance and (instantaneous) radius of gyration respectively.

The _R<sub>g</sub>/R<sub>g</sub><sup>mean</sup>_ is a measure of (relative) size for a protein or polymer chain, and _R<sub>s</sub>_ is a measure of its shape. _R<sub>s</sub>_ is expected to be low (~2 or lower) for compact structures and high for highly extended structures (~12 or higher). A single _R<sub>g</sub>/R<sub>g</sub><sup>mean</sup>_ value and corresponding _R<sub>s</sub>_ value for a polymer together is how we define its instantaneous conformation. When all the _R<sub>g</sub>/R<sub>g</sub><sup>mean</sup>_ and  _R<sub>s</sub>_ values of a polymer are plotted together, they constitute what we call a 2D map of the conformational landscape of that polymer.

## The _PyConforMap_ Module Generates Two-Dimensional Scatter Plots
This module generates 2D scatter plots of _R<sub>s</sub>_ against _R<sub>g</sub>/R<sub>g</sub><sup>mean</sup>_ for a protein/polymer simulation (data and protein label/identity provided by user) and a Gausssian Walk (GW) polymer model simulation (data for 720000 snapshots of a GW model of length 100 included with repository). Each point on the scatter plot (belonging to either GW or a given protein/polymer) represents a conformation snapshot, and has coordinates (_R<sub>g</sub>/R<sub>g</sub><sup>mean</sup>_, _R<sub>s</sub>_). The GW model is intended to be a reference model, whose conformational landscape map (i.e. as represented by all the (_R<sub>g</sub>/R<sub>g</sub><sup>mean</sup>_, _R<sub>s</sub>_) points) serves as a 'universal' or reference map for those of other proteins/polymers. Using the 2D scatter plot, an _f<sub>C</sub>_, representing the fraction of the GW points 'close' (i.e. within a pre-defined radius) to at least one protein/polymer point, is automatically calculated. _f<sub>C</sub>_ is a quantity that represents the conformational diversity of the protein/polymer provided, and can be used to rank the conformational diversities of different proteins/polymers. The included GW file is 'GW_chainlen100.csv.' The python module can be additionally used to conduct a new GW simulation with different chain length and number of snapshots, should the user wish to do so. On the scatter plot, it is important that the protein/polymer points do not significantly exceed the boundaries defined by the reference (GW) points. Most of the protein/polymer points should be 'close' (i.e. within a pre-defined radius) to at least one GW point.

## The Module Code Requires One Input File
The needed input is a csv file (for a given protein/polymer simulation) with 2 columns. The first column contains _R<sub>g</sub><sup>2</sup>_ values and the second column contains _R<sub>ee</sub><sup>2</sup>_ values. In this (user provided) file, each row represents a protein/polymer conformation snapshot from the simulation. An example input is the 'example_protein.csv' csv file (included with repository). 

## Files Included with Repository
The 'code_input_output.md' file provides technical details (input arguments, expected outputs) of the module. The 'pyconformap.py' file contains the source code for the module.  The 'illustrated_example.ipynb' jupyter notebook file shows examples to illustrate implementation of the code. The 'GW_chainlen100.csv' is the reference GW simulation and 'example_protein.csv' is the simulation of an example protein.

## Packages Required
The module requires the _pandas_, _numpy_, _matplotlib_, _scipy_, _itertools_, _more_itertools_ and _random_ python packages. They are automatically loaded when the 'pyconformap.py' file is executed, as shown in the illustrated examples.

## Publication
_PyConforMap_ is companion to this [publication](https://www.cell.com/biophysj/abstract/S0006-3495(24)00272-8).

## How to Cite
If you use this module, please cite us using the provided DOI. 

## Contact Information
If you have comments/suggestions or a bug report, please feel free to email me at hossain.shadman17@gmail.com, or contact me through my social media links provided in the home page. 
