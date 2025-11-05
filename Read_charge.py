import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt

file = input()

u = mda.Universe(file)
protein = u.select_atoms("protein")
charges = protein.atoms.charges
plt.plot(charges)