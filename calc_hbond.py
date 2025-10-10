"""
Hydrogen Bond Analysis and Visualization Script

This script performs hydrogen bond analysis on molecular dynamics trajectories
for both wild-type (WT) and mutant (Mut) systems using MDAnalysis. It:
1. Calculates hydrogen bond occupancy maps for both systems
2. Performs a two-proportion Z-test to identify statistically significant 
   differences in hydrogen bonding patterns between WT and Mut
3. Generates a VMD TCL script to visualize significant differences:
   - Blue cylinders: hydrogen bonds more frequent in WT
   - Red cylinders: hydrogen bonds more frequent in Mut
4. Creates a histogram of Z-scores for quality control

The visualization uses CA atom coordinates from a reference structure to 
draw cylinders representing residue-residue hydrogen bond differences.

Author: Jonas Paulus usign Cline and QWEN3 235B
Date: 2025-10-09
"""
import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA
import numpy as np
import matplotlib.pyplot as plt

# Configuration parameters
# Trajectory files
top_wt = "WT.tpr"
traj_wt = "WT_0_1_fit1.xtc"
top_mut = "md_Mut_0_1.tpr"
traj_mut = "fitted_mutant.xtc"

# visualisation structure for coordinate lookup
plot_structure = "vis_structure.gro"

# Analysis parameters
n_frames_wt = 2001  # Number of frames in WT trajectory
n_frames_mut = 2001  # Number of frames in Mut trajectory
z_score_threshold = 3  # Z-score threshold for significant differences
se_factor = 0.3  # SE multiplier for threshold adjustment

# Output files
tcl_output = "draw_hbonds_z_test.tcl"
histogram_output = "z_test_histogram.png"

#%% ================== process WT data =======================================
# Initialize universe with WT files
u = mda.Universe(top_wt, traj_wt)

hbonds = HBA(universe=u)#, between=["protein", "protein"]
hbonds.hydrogens_sel = hbonds.guess_hydrogens("protein")
hbonds.acceptors_sel = hbonds.guess_acceptors("protein")
hbonds.run(backend = "multiprocessing", n_workers = 11)
results = hbonds.results["hbonds"]


protein = u.select_atoms("protein")
hb_map = np.zeros([protein.atoms.n_residues, protein.atoms.n_residues])
for i in range(results.shape[0]):
    resid_don = u.atoms.resids[int(results[i,1])]
    resid_acc = u.atoms.resids[int(results[i,3])]
    #rint(resid_acc, ", ", resid_don)
    if resid_don < 520 and resid_acc < 520:
        hb_map[resid_don-1, resid_acc-1] += 1

WT_map = hb_map/ u.trajectory.n_frames
#%% ========================== Process Mutation Data ================================================
u = MDAnalysis.Universe(top_mut, traj_mut)

hbonds = HBA(universe=u)#, between=["protein", "protein"]
hbonds.hydrogens_sel = hbonds.guess_hydrogens("protein")
hbonds.acceptors_sel = hbonds.guess_acceptors("protein")
hbonds.run(backend = "multiprocessing", n_workers = 11)
results = hbonds.results["hbonds"]

protein = u.select_atoms("protein")
hb_map = np.zeros([protein.atoms.n_residues, protein.atoms.n_residues])
for i in range(results.shape[0]):
    resid_don = u.atoms.resids[int(results[i,1])]
    resid_acc = u.atoms.resids[int(results[i,3])]
    #rint(resid_acc, ", ", resid_don)
    if resid_don < 520 and resid_acc < 520:
        hb_map[resid_don-1, resid_acc-1] += 1

Mut_map = hb_map/ u.trajectory.n_frames
# %% ============================ Calculate test statistics =======================================================
#two proportion Z-test using the default naming 
#further information see wikipedia: https://en.wikipedia.org/wiki/Two-proportion_Z-test
n1 = 2001
n2 = 2001
p1 = WT_map
p2 = Mut_map
p = (p1*n1+p2*n2)/(n1+n2)
z = (p1-p2)/np.sqrt(p*(1-p)*(1/n1 + 1/n2))
SE = np.sqrt(p*(1-p)*(1/n1 + 1/n2))
z[np.isnan(z)] = 0

#%% =============================== Write TCL script ==========================================================
script = []
u = mda.Universe(plot_structure)
ca = u.select_atoms("name CA")

#in wt>mut
i_set, j_set = np.where(z>3+0.3/SE)
script.append("draw color 0\n")
for n in range(len(i_set)):
    i = i_set[n]
    j = j_set[n]
    if i>j:
        p_i = ca.atoms.positions[i,:]
        p_j = ca.atoms.positions[j,:]
        script.append(f"draw cylinder {{{p_i[0]:.3f} {p_i[1]:.3f} {p_i[2]:.3f}}} {{{p_j[0]:.3f} {p_j[1]:.3f} {p_j[2]:.3f}}} radius 0.2\n")

#in wt<mut
i_set, j_set = np.where(z< -(3 + 0.3/SE))
script.append("draw color 1\n")
for n in range(len(i_set)):
    i = i_set[n]
    j = j_set[n]
    if i>j:
        p_i = ca.atoms.positions[i,:]
        p_j = ca.atoms.positions[j,:]
        script.append(f"draw cylinder {{{p_i[0]:.3f} {p_i[1]:.3f} {p_i[2]:.3f}}} {{{p_j[0]:.3f} {p_j[1]:.3f} {p_j[2]:.3f}}} radius 0.2\n")

with open("draw_hbonds_z_test.tcl", "w") as file:
    file.writelines(script)


# %% ================================ Quality controll test ====================================================


plt.hist((np.abs(z))[np.where(np.abs(z) > 3)])
plt.savefig("z_test_histogram.png")
