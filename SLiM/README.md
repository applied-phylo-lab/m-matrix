This folder contains codes of SLiM simulations and output files.

# stabilizing/
Simulations where all traits under consideration are under stabilizing selection

# neutral/
Simulations where a focal trait is under no selection, whereas other traits are under stabilizing selection

# directional_WF/
Simulations where a focal trait is under directional selection, whereas other traits are under stabilizing selection

# directional_nWF/
Simulations where a focal trait is under directional selection, whereas other traits are under stabilizing selection. Non-Wright-Fisher populations are simulated.

Each of the folders contains the following files, as explained below. See README files within each directory for descriptions of files specific to each of them.

# 2_traits.txt
SLiM code for simulating 2 traits. 

# 5_traits.txt
SLiM code for simulating 5 traits.

# 10_traits.txt
SLiM code for simulating 10 traits.

# sim_out_end_<number_of_traits>t_<number_of_pleiotropic_loci>.txt or sim_qn_out_end_<number_of_traits>t_<number_of_pleiotropic_loci>.txt or  sim_dir_out_end_<number_of_traits>t_<number_of_pleiotropic_loci>.txt or sim_dir_nwf_out_end_<number_of_traits>t_<number_of_pleiotropic_loci>.txt
Data files containing output at the end of the simulations.

# sim_out_all_<number_of_traits>t_<number_of_pleiotropic_loci>.txt
Data files containing simulation output over time.

# sim_out_all_<number_of_traits>t_<number_of_pleiotropic_loci>_processed.txt
Data files containing simulation output over time. The same as the above, but with top lines removed such that they are ready for making figures.


