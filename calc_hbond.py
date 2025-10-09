import MDAnalysis
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA
import numpy as np
import matplotlib.pyplot as plt

top_wt = "WT.tpr"
traj_wt = "WT_0_1_fit1.xtc"
u = MDAnalysis.Universe(top_wt, traj_wt)

hbonds = HBA(universe=u)#, between=["protein", "protein"]
hbonds.hydrogens_sel = hbonds.guess_hydrogens("protein")
hbonds.acceptors_sel = hbonds.guess_acceptors("protein")
hbonds.run(backend = "multiprocessing", n_workers = 11)
results = hbonds.results["hbonds"]

#%%
protein = u.select_atoms("protein")
hb_map = np.zeros([protein.atoms.n_residues, protein.atoms.n_residues])
for i in range(results.shape[0]):
    resid_don = u.atoms.resids[int(results[i,1])]
    resid_acc = u.atoms.resids[int(results[i,3])]
    #rint(resid_acc, ", ", resid_don)
    if resid_don < 520 and resid_acc < 520:
        hb_map[resid_don-1, resid_acc-1] += 1

WT_map = hb_map/ u.trajectory.n_frames
#%%
top_wt = "md_Mut_0_1.tpr"
traj_wt = "fitted_mutant.xtc"
u = MDAnalysis.Universe(top_wt, traj_wt)

hbonds = HBA(universe=u)#, between=["protein", "protein"]
hbonds.hydrogens_sel = hbonds.guess_hydrogens("protein")
hbonds.acceptors_sel = hbonds.guess_acceptors("protein")
hbonds.run(backend = "multiprocessing", n_workers = 11)
results = hbonds.results["hbonds"]

#%%
protein = u.select_atoms("protein")
hb_map = np.zeros([protein.atoms.n_residues, protein.atoms.n_residues])
for i in range(results.shape[0]):
    resid_don = u.atoms.resids[int(results[i,1])]
    resid_acc = u.atoms.resids[int(results[i,3])]
    #rint(resid_acc, ", ", resid_don)
    if resid_don < 520 and resid_acc < 520:
        hb_map[resid_don-1, resid_acc-1] += 1

Mut_map = hb_map/ u.trajectory.n_frames
# %%
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

script = []
new_wt_file = "/localscratch/projects/dennis/debug1/vis_structure.gro"
u = mda.Universe(new_wt_file)
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
        script.append(f"draw cylinder ö{p_i[0]:.3f} {p_i[1]:.3f} {p_i[2]:.3f}ü ö{p_j[0]:.3f} {p_j[1]:.3f} {p_j[2]:.3f}ü radius 0.2\n".replace("ö", "{").replace("ü", "}"))

#in wt<mut
i_set, j_set = np.where(z< -(3 + 0.3/SE))
script.append("draw color 1\n")
for n in range(len(i_set)):
    i = i_set[n]
    j = j_set[n]
    if i>j:
        p_i = ca.atoms.positions[i,:]
        p_j = ca.atoms.positions[j,:]
        script.append(f"draw cylinder ö{p_i[0]:.3f} {p_i[1]:.3f} {p_i[2]:.3f}ü ö{p_j[0]:.3f} {p_j[1]:.3f} {p_j[2]:.3f}ü radius 0.2\n".replace("ö", "{").replace("ü", "}"))

with open("draw_hbonds_z_test.tcl", "w") as file:
    file.writelines(script)


# %%


plt.hist((np.abs(z))[np.where(np.abs(z) > 3)])