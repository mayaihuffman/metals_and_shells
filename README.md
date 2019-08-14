# README
Maya I. Huffman\
mih3@williams.edu\
August 14, 2019
## add_hydrogens.py
### Overview
add_hydrogens ended up being used only for putting hydrogens on non-metal-binding homologs of copper- and zinc-containing proteins. It puts hypothetical hydrogen atoms on the protein using PyRosetta and then cuts down the PDB and adds the metal and coordinating residues from its metal binding counterpart to prepare for global analysis.
### Inputs
**Nonmetal pdb file.** The names of the nonmetal PDB files were formatted as metal PDB code (lowercase), chain letter (uppercase), underscore, followed by the PDB code (lowercase) and chain letter (uppercase) of the non-metal-containing homolog matched and aligned by the Dali server. Ex. *1ck7A_1qjsA.pdb* It is assumed here that the naming scheme for the original PDBs (with the metals) are the four character PDB codes (no chain letter, lowercase) followed by *_0001.pdb*

**Metal.** This is necessary in order to run *shell_res_num.py** in order to place the proper metals and coordinating residues in the non-metal-binding protein. The metal must be in a list (ex. ["CU"]). It should be the two letter abbreviation with both letters uppercase.
### Outputs
Re-writes PDB with the same file name as the input PDB file, now with hydrogens, only the chain of interest, and a metal and coordinating residues from the metalloprotein homolog.
## atom_walk.py
### Overview
Analyzes "branches" of the first and second shell of a metal site using k-means clustering to divide the shells into branches. First the shells are built just as in global_analysis and shell_res_num, but then the shells are split into as many branches as there are coordinating atoms for that metal.
### Inputs
See run_atom_walk.sh for inputs.
### Outputs
Writes and appends lines to a csv file with information about metal site branches including pKa, partial charge averages, hydropathy averages and distributions, and several distances, angles, and dihedrals for geometric analysis between different levels of the branch.
## constants_dictionaries.py
### Overview
Formerly referred to as a partial charge dictionary, but now has dictionaries of all of the constants referenced in the analysis programs. This includes:
* **Partial charge dictionary** for each atom in each type of amino acid
* **Dielectric constants dictionary** for each amino acid
* **pKa dictionary** for amino acids with listed pKa values for their sidechains
* **Hydropathy dictionary** for each amino acid

## coord_envi.py
### Overview
Analyzes "branches" of the first and second shells of a metal site by keeping the shell members attached to their coordinating atom root rather than combining all shell residues into a full metal shell and subsequently dividing the shell into sections (the way atom_walk does).
### Inputs
See run_coord_envi.sh for inputs.
### Outputs
Write and appends lines to csv file with information about metal site branches including pka, partial charge averages, hydropathy averages and distributions, and several distances, angles, and dihedrals for geometric analysis between different levels of the branch.
## global_analysis.py
### Overview
Defines and analyzes coordination, first, and second shell of a metal site. Every shell is reduced to a single average or sum and direction or distance is recorded.
### Inputs
See run_global.sh for inputs
### Outputs
Write and appends lines to csv file with information about metal site including pKa, partial charge averages, and hydropathy values and distributions.
## graphs.py
### Overview
These graph functions tend to be made very specifically for one purpose (or even a one-time use). Therefore, even though there is a list of csvs for many of these graph functions, some of the graphs may not work if there are fewer than or more than the number of csvs in he default csv list. In addition, different graphs require different types of csvs because some csvs have more information than others. I have left in default csv lists that will not work as you likely will not have these exact csvs. But the default settings present will give an idea of what is necessary in order for these plots to work.
### inter_metal_comparison_histograms_global
#### Overview
Creates histograms based on data from global analysis.
#### Inputs
**CSV list.** These csvs should come from global analysis.

**Metal names.** These should correspond to order of csvs in csv list.

**Color** If highlighting one color at a time, I have been using red to highlight one data set and 0.6 (a medium-light grey color) for all other data sets.

**Bins.** Right now applies the same number of bins to every histogram, but in the future it might be more helpful to take out this input and make predetermined bins to better suit each histogram.
#### Outputs
2d histograms showing distributions of pKa, charges, hydropathy, secondary structure.
### intra_metal_comparison_histograms_global
#### Overview
Plots first shell versus second shell within a data set (metal type).
#### Inputs
**CSV list.** Should be global analysis csvs.

**Metal names.** Should correspond with the order of csvs in csv list.

**Colors.** Always put in two colors, one for first shell and other for second shell.

**Bins.** Right now same bin number applied to all histograms.
#### Outputs
Creates a series of histograms comparing first and second shells within a single metal.
### inter_metal_comparison_histograms_kmeans
#### Overview
Essentially does the same thing as global version but uses each branch as a data point instead of each shell.
#### Inputs
Same inputs as global version of this program, but uses "atom walk" csvs.
#### Outputs
Same outputs as inter metal comparison global except for some extra ones regarding distances and angles.
### intra_metal_comparison_histograms_kmeans
#### Overview
Essentially does the same thing as the global version of intra metal comparison histograms, uses branches as data points instead of shells.
#### Inputs
Same inputs as intra metal comparison but requires "atom walk" csvs.
#### Outputs
Same outputs as intra metal comparison global but there are some extra graphs with distances and angles.
### plot_3d_scatter
#### Overview
Plots various 3D scatterplots from the atom walk csvs.
#### Inputs
**CSV list.** Requires atom walk csvs
#### Outputs
3D scatters for distances, angles, and dihedrals.
### find_avg_and_stdev
#### Overview
Makes new csvs with standard deviations and averages
#### Inputs
**CSV list.** Requires atom walk csvs.

**New CSV names.** Names for the new csv files that will be written.
#### Outputs
CSVs with averages and standard deviations.
### scatter_branch_info_2d
#### Overview
2D scatter plots of various data points such as energy and hydropathy.
#### Inputs
**CSV list.** Intended for atom walk csvs.

**Metal names.** Correspond to order of csvs in csv list.

**Colors.** Color of dots on scatter plot.
#### Outputs
Various 2D scatter plots.
### scatter_branch_info_3d
#### Overview
Same as 2D scatter plots but the new axis is pKa.
#### Inputs
See inputs for scatter branch info 2d.
#### Outputs
Various 3D scatterplots with pKa as z axis.
### scatter_2d_global_vs_redox
#### Overview
Intended for the plastocyanin PDBs with associated redox potentials - right now redox potential is only matched up to specific PDB codes already listed.
#### Inputs
One CSV with plastocyanin global info.
#### Outputs
Various scatter plots of data versus the redox potential.
### scatter_2d_coorenvi_vs_pka
#### Overview
Scatter plots of different branch information versus pKa.
#### Inputs
**CSV list.** Should be coor envi csvs but atom walk might also work.

**Colors.** Works best if you make all of them 0.6 grey except for the top one (the last csv in the list), which will be red.

**Last metal.** The metal of the last csv in the csv list (will be the top plot on the scatter).
#### Outputs
Various scatter plots with one data set plotted in red and everything else in grey.
### scatter_2d_coorenvi_branches
#### Overview
Essentially the same ting as the atom walk 2d scatter.
#### Inputs
See inputs for scatter branch info 2D.
#### Outputs
2D scatterplots but not always against pKa.
### heatmap_vs_pka
#### Overview
Essentially achieves the same end as scatter 2d coor envi vs pKa except it takes into account how many times a value falls into a certain range or bin.
#### Inputs
**CSV list.** Should use coor envi csvs (again, atom walk, might work).

**Subdivisions.** Basically the same thing as "bins" for a histogram except that I created them by hand and this isn't a histogram.
#### Outputs
Heatmaps in the style and representing the smae information as 2D coor envi vs pKa.
## nrpdb.py
### Overview
Start with a list of non-redundant PDB codes, select only those not containing any metals, and then run a form of global analysis on these proteins by seeding random points and building/analyzing shells from there.
### shorten_nrpdb
#### Overview
Takes the list of non-redundant proteins and maeks a new list of only the pdb codes (plus their chain letters) that do not have any metals.
#### Inputs
**Full nrpdb list.** Full nrpdb list can be taken from the [Vector Alignment Search Tool (VAST)](https://structure.ncbi.nlm.nih.gov/Structure/VAST/nrpdb.html) at non-redundancy level 10e-7, giving upwards of 14,000 sequence non-redundant PDB codes. This will have to be reformatted first.

**Short nrpdb list name** Name what will become the new file of only non-metal binding PDBs.
#### Outputs
The new list of PDB codes that does not contain any proteins with metals (should be slightly more than 10,000 PDB codes).
### start_csv_no_ss
#### Overview
Makes the header for the new csv. This layout is the same as global except that it does not contain any secondary structure analysis because that tended to cause too much failure and took a long time.
#### Inputs
Name of new csv.
#### Outputs
New csv with header labels.
### prep_and_run_nrpdbs
#### Overview
Reads, cleans, poses, shortens, and analyzes PDBs based on the PDB codes (no pre-existing files needed)
#### Inputs
**Start and end indices.** In place so that the program can be run in chunks rather than all at once. The default indices are in place to run the whole thing at once.

**Shortened nrpdb list.** A list of PDB codes and chain letters allows Rosetta to pose the proteins to be analyzed.

**Output file.** Should already be created with start_csv_no_ss.
#### Outputs
Appends lines to the csv output file.
## run_atom_walk.sh
### Overview
Shell script for atom walk analysis.
### Inputs
**Metal PDB list.** This hsould be formatted in a way that the path to the metal pdb folder is in place and also the pdb itself is formatted as the four letter pdb code and followed by "0001.pdb". A member of this list might look like: "./fe3_pdbs/1dv6_0001.pdb" if there is a fe3_pdbs folder in the working directory.

**Metal list.** this list should have the brackets and quotes in it and contain the two-letter abbreviation for a metal with both letters uppercase. Ex. ["CU"]

**CSV file name.** The csv does not have to exist before running the program, but the desired name of the csv that will contain all of the data is necessary.

Example: bash run.sh metal_list ["CU"] output.csv
### Outputs
See atom_walk.py for output.
## run_coord_envi.sh
### Overview
Shell script for coord envi analysis.
### Inputs
See run_atom_walk.sh for inputs.

Example: bash run_coord_envi.sh metal_list ["CU"] output.csv
### Outputs
See coord_envi.py for outputs.
## run_global.sh
### Overview
Shell script for global analysis.
### Inputs
See run_atom_walk.sh for inputs.

Example: bash run_global.sh metal_list ["CU"] output.csv
### Outputs
See global_analysis.py for outputs.
## shell_res_num.py
### Overview
This is an isolated version of the code that appears in global_analysis.py and also appears in modified staets in atom walk and coord envi.
### Inputs
**PDB file name or PDB code.** If PDB is to be taken from file (pdb_from_text = True), then the pdb_file_name has to be a string of the entire file name including ".pdb" at the end. If PDB is to be taken from the RCSB url (pdb_from_text = False), then you just need to put in the four-character pdb code (uppercase) as the input. The rest of the url is part of the function that reads in the PDB from the URL.

**PDB from text.** If using a downloaded PDB file, mark this as true. If reading from a URL, mark as false.

**Metal.** Need list of metals ex. ["CU"] so that only certain metals are remembered and returned.

Examples: find_shells_info("4dyz_0001.pdb", metals = ["CU"], pdb_from_text = True)\
find_shells_info("2VB2", pdb_from_text = False)
### Outputs
Returns all of the different lists and dictionaries that are used in other codes such as global. The most useful of these are probably the residue numbers lists because it allows for easy visualization of shells in PyMol.
## split_pdb.py
### Overview
Not of any particular use in its own right; it is just used in order to understand how a what each part a PDB line slice means and how these line slices are made into dictionaries.
## mih_transform.py
### Overview
Transform was written by William Hansen and was sent to me on June 18, 2019. I had to make some minor changes in order for it not to generate error messages when using in Python 3 (the program was written in Python 2). Other minor changes were made by either William Hansen or myself. The Transform class is not used in any of my own code. The actual version of Transform has since been updated, making this version both out of date and edited to serve a different purpose.
