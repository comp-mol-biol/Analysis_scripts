from multiprocessing import cpu_count
import MDAnalysis
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA
from MDAnalysis.lib import distances
import numpy as np

# Documentation for the contacts.py module
# This module provides functions for analyzing protein contacts and hydrogen bonds in molecular dynamics trajectories.
# The main functions are:
# 1. hbonds: Calculates hydrogen bonds between protein residues
# 2. compute_residue_contacts: Calculates residue-residue contacts based on a distance cutoff
def hbonds(top, traj, between = None):
    """
    Calculate hydrogen bonds between protein residues in a molecular dynamics trajectory.
    
    Parameters:
    top (str): Path to the topology file (e.g., .tpr, .gro, .pdb)
    traj (str): Path to the trajectory file (e.g., .xtc, .trr, .dcd)
    between (list, optional): List of selections for donor-acceptor pairs. If None, all protein atoms are considered.
    
    Returns:
    tuple: A tuple containing:
        - hb_map_normalized (numpy.ndarray): Normalized hydrogen bond map (n_residues x n_residues)
        - n_frames (int): Number of frames in the trajectory
    
    Note: The hydrogen bond map is symmetric, with each entry representing the frequency of hydrogen bonds
    between residue pairs across all frames.
    """

    u = MDAnalysis.Universe(top, traj)

    # Initialize hydrogen bond analysis
    hbonds = HBA(universe=u, between=between)
    # Guess hydrogen and acceptor atoms for protein
    hbonds.hydrogens_sel = hbonds.guess_hydrogens("protein")
    hbonds.acceptors_sel = hbonds.guess_acceptors("protein")
    # Run analysis using multiprocessing for better performance
    hbonds.run(backend="multiprocessing", n_workers=cpu_count()-1)
    results = hbonds.results["hbonds"]

    # Select protein atoms for residue-based analysis
    protein = u.select_atoms("protein")
    # Initialize hydrogen bond map (n_residues x n_residues)
    hb_map = np.zeros([protein.atoms.n_residues, protein.atoms.n_residues])
    
    # Process hydrogen bond results
    for i in range(results.shape[0]):
        # Get residue IDs for donor and acceptor
        resid_don = u.atoms.resids[int(results[i,1])]
        resid_acc = u.atoms.resids[int(results[i,3])]
        # Check if both residues are part of the protein
        max_protein_resid = protein.resids[-1]
        if resid_don < max_protein_resid and resid_acc < max_protein_resid:
            # Increment hydrogen bond count for this residue pair
            hb_map[resid_don-1, resid_acc-1] += 1
            # If not the same residue, increment the symmetric position (undirected graph)
            if resid_don != resid_acc:
                hb_map[resid_acc-1, resid_don-1] += 1

    # Normalize by number of frames to get frequency
    hb_map_normalized = hb_map / u.trajectory.n_frames
    return hb_map_normalized, u.trajectory.n_frames

def compute_residue_contacts(top, traj, selection="name CA", cutoff=8.0, align_ref_sel=None, step=1):
    """
    Calculate residue-residue contacts based on a distance cutoff in a molecular dynamics trajectory.
    
    Parameters:
    top (str): Path to the topology file (e.g., .tpr, .gro, .pdb)
    traj (str): Path to the trajectory file (e.g., .xtc, .trr, .dcd)
    selection (str, optional): Atom selection for calculating residue centers. Default is "name CA".
    cutoff (float, optional): Distance cutoff for contacts in Angstroms. Default is 8.0.
    align_ref_sel (str, optional): Selection for alignment reference. If provided, trajectory will be aligned.
    step (int, optional): Step size for trajectory frames. Default is 1 (every frame).
    
    Returns:
    tuple: A tuple containing:
        - residues (MDAnalysis.core.groups.ResidueGroup): Residue group for the selection
        - freq (numpy.ndarray): Contact frequency matrix (n_residues x n_residues)
        - n_frames (int): Number of frames processed
    
    Note: The contact matrix is symmetric, with each entry representing the frequency of contacts
    between residue pairs across all processed frames.
    """
    u = MDAnalysis.Universe(top, traj)
    sel = u.select_atoms(selection)
    residues = sel.residues
    n_res = len(residues)
 
    # Align trajectory if reference selection is provided
    if align_ref_sel:
        ref = mda.Universe(top)
        R = align.AlignTraj(u, ref, select=align_ref_sel, in_memory=True)
        R.run()
 
    # Initialize contact count matrix
    counts = np.zeros((n_res, n_res), dtype=np.uint16)
    n_frames = 0
 
    # Process trajectory frames
    for ts in u.trajectory[::step]:
        # Calculate center of mass for each residue
        coords = residues.center_of_mass(compound='residues')
        # Calculate distance matrix between all residue centers
        D = distances.distance_array(coords, coords, box=u.dimensions)
 
        # Count contacts for each residue pair
        for i in range(n_res):
            for j in range(i+1, n_res):
                if D[i, j] < cutoff:
                    counts[i, j] += 1
                    counts[j, i] += 1
        n_frames += 1
 
    # Calculate contact frequency by normalizing by number of frames
    freq = counts.astype(np.float32) / float(n_frames)
    return residues, freq, n_frames
