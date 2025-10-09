import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda
import MDAnalysis.analysis.align as align
from MDAnalysis.lib import distances

#%%
#%% ---------- Einstellungen ----------
top_wt = "md_WT_0_1.gro"
traj_wt = "WT_0_1_fit1.xtc"
top_mut = "md_Mut_0_1.gro"
traj_mut = "fitted_mutant.xtc"
 
selection = "name CA"       # Residue-Vertreter
cutoff = 8.0                # Å
align_ref_sel = "backbone"  # Alignment
step = 1                    # jedes Frame (du kannst auf 5 oder 10 erhöhen)
 
# ---------- Hilfsfunktion ----------
def compute_residue_contacts(top, traj, selection, cutoff, align_ref_sel=None, step=1):
    u = mda.Universe(top, traj)
    sel = u.select_atoms(selection)
    residues = sel.residues
    n_res = len(residues)
 
    if align_ref_sel:
        ref = mda.Universe(top)
        R = align.AlignTraj(u, ref, select=align_ref_sel, in_memory=True)
        R.run()
 
    counts = np.zeros((n_res, n_res), dtype=np.uint16)
    n_frames = 0
 
    for ts in u.trajectory[::step]:
        coords = residues.center_of_mass(compound='residues')
        D = distances.distance_array(coords, coords, box=u.dimensions)
 
        # nur bool, keine große Matrix aufheben
        for i in range(n_res):
            for j in range(i+1, n_res):
                if D[i, j] < cutoff:
                    counts[i, j] += 1
                    counts[j, i] += 1
        n_frames += 1
 
    freq = counts.astype(np.float32) / float(n_frames)
    return residues, freq, n_frames
 
# ---------- Berechnungen ----------
res_wt, freq_wt, frames_wt = compute_residue_contacts(top_wt, traj_wt, selection, cutoff, align_ref_sel, step)
res_mut, freq_mut, frames_mut = compute_residue_contacts(top_mut, traj_mut, selection, cutoff, align_ref_sel, step)
 
assert len(res_wt) == len(res_mut)
 
diff = freq_wt - freq_mut


#%% generate TCL script


"""
script = []
new_wt_file = "/localscratch/projects/dennis/debug1/vis_structure.gro"
u = mda.Universe(new_wt_file)
ca = u.select_atoms("name CA")

#in WT but not in MUT
i_set, j_set = np.where(diff>0.8)
script.append("draw color 0\n")
for n in range(len(i_set)):
    i = i_set[n]
    j = j_set[n]
    if i>j:
        p_i = ca.atoms.positions[i,:]
        p_j = ca.atoms.positions[j,:]
        script.append(f"draw cylinder ö{p_i[0]:.3f} {p_i[1]:.3f} {p_i[2]:.3f}ü ö{p_j[0]:.3f} {p_j[1]:.3f} {p_j[2]:.3f}ü radius 0.2\n".replace("ö", "{").replace("ü", "}"))

#in mut but not in wt
i_set, j_set = np.where(diff<-0.8)
script.append("draw color 1\n")
for n in range(len(i_set)):
    i = i_set[n]
    j = j_set[n]
    if i>j:
        p_i = ca.atoms.positions[i,:]
        p_j = ca.atoms.positions[j,:]
        script.append(f"draw cylinder ö{p_i[0]:.3f} {p_i[1]:.3f} {p_i[2]:.3f}ü ö{p_j[0]:.3f} {p_j[1]:.3f} {p_j[2]:.3f}ü radius 0.2\n".replace("ö", "{").replace("ü", "}"))

with open("draw_bonds.tcl", "w") as file:
    file.writelines(script)

#in mut but not in wt
i_set, j_set = np.where(np.logical_and(freq_wt > 0.8, freq_mut >0.8))
script.append("draw color 16\n")
script.append("draw material Glass1\n")
for n in range(len(i_set)):
    i = i_set[n]
    j = j_set[n]
    if i>j+2:
        p_i = ca.atoms.positions[i,:]
        p_j = ca.atoms.positions[j,:]
        script.append(f"draw cylinder ö{p_i[0]:.3f} {p_i[1]:.3f} {p_i[2]:.3f}ü ö{p_j[0]:.3f} {p_j[1]:.3f} {p_j[2]:.3f}ü radius 0.2\n".replace("ö", "{").replace("ü", "}"))

with open("draw_bonds.tcl", "w") as file:
    file.writelines(script)
"""
# %%

n1 = frames_wt
n2 = frames_mut
p1 = freq_wt
p2 = freq_mut
p = (p1*n1+p2*n2)/(n1+n2)
z = (p1-p2)/np.sqrt(p*(1-p)*(1/n1 + 1/n2))
z[np.isnan(z)] = 0
SE = np.sqrt(p*(1-p)*(1/n1 + 1/n2))
script = []
new_wt_file = "/localscratch/projects/dennis/debug1/vis_structure.gro"
u = mda.Universe(new_wt_file)
ca = u.select_atoms("name CA")

#in wt>mut
i_set, j_set = np.where(z>3+0.6/SE)
script.append("draw color 0\n")
for n in range(len(i_set)):
    i = i_set[n]
    j = j_set[n]
    if i>j:
        p_i = ca.atoms.positions[i,:]
        p_j = ca.atoms.positions[j,:]
        script.append(f"draw cylinder ö{p_i[0]:.3f} {p_i[1]:.3f} {p_i[2]:.3f}ü ö{p_j[0]:.3f} {p_j[1]:.3f} {p_j[2]:.3f}ü radius 0.2\n".replace("ö", "{").replace("ü", "}"))

#in wt<mut
i_set, j_set = np.where(z< -(3 + 0.6/SE))
script.append("draw color 1\n")
for n in range(len(i_set)):
    i = i_set[n]
    j = j_set[n]
    if i>j:
        p_i = ca.atoms.positions[i,:]
        p_j = ca.atoms.positions[j,:]
        script.append(f"draw cylinder ö{p_i[0]:.3f} {p_i[1]:.3f} {p_i[2]:.3f}ü ö{p_j[0]:.3f} {p_j[1]:.3f} {p_j[2]:.3f}ü radius 0.2\n".replace("ö", "{").replace("ü", "}"))

with open("draw_bonds_z_test.tcl", "w") as file:
    file.writelines(script)


# %%


plt.hist((np.abs(z)*SE)[np.where(np.abs(z) > 3)])