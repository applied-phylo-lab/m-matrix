This folder contains code and data files for scenarios where a focal trait is neutral while other traits are under stabilizing seletion.

Structure of data files (files whose names start with "sim_out_") is explained below.

Column 2: Mean of a trait.
Column 3: Mean within-population variance of a trait.

Numbers in column 1 correspond to populations in output files of 2-traits simulations and correspond to traits in other files.
However, in all files, rows are ordered in such a way that rows from the same population are clustered together.

For 2-traits simulations, the ordering is:
"population 1, trait 1",
"population 1, trait 2",
"population 2, trait 1",
"population 2, trait 2",
...,
"population 50, trait 1",
"population 50, trait 2".

For 5-traits simulations, the ordering is:
"population 1, trait 1",
"population 1, trait 2",
"population 1, trait 3",
"population 1, trait 4",
"population 1, trait 5",
...,
"population 50, trait 1",
"population 50, trait 2",
"population 50, trait 3",
"population 50, trait 4",
"population 50, trait 5",

10-traits simulations' output files follow the same rule.
