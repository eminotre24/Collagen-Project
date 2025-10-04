import MDAnalysis as md
import matplotlib.pyplot as plt
import numpy as np

# Load data
t_set = 100

if t_set == 25:
    # 25 ns data
    file_dir = "/Users/alejandro/Desktop/mddocs/processing-analysis/fullcol-0308/results-0808"
elif t_set == 100:
    # 100 ns data
    file_dir = "/Users/alejandro/Desktop/mddocs/processing-analysis/fullcol-100ns-2008"

u = md.Universe(file_dir + "/md_1.tpr", file_dir + "/md_1_noPBC.xtc")
u.trajectory[-1]
protein = u.select_atoms("protein")

t_eq = 20000 # ps
# Time Data
ti = u.trajectory[0].time
tf = u.trajectory[-1].time
nt = len(u.trajectory)-1
deltat = (tf-ti)/nt
t = np.linspace(ti,tf,nt+1)
maskt = (t >= t_eq)

# --- Functions ---
def get_mean_coords_mol(molecule, mask, group):
    # Get the Mean Coordinates of the molecule
    mean_coords = None
    chains = protein.fragments[3*molecule:3*molecule + 3]
    colmol = chains[0] + chains[1] + chains[2]
    colmol = colmol.select_atoms(group)
    for ts in u.trajectory[mask]:
        coords = colmol.positions / 10
        if mean_coords is None:
            mean_coords = coords.astype(float)
            n = 1
        else:
            n += 1
            mean_coords += (coords - mean_coords) / n

    # Ordering the z for unwrapping
    order = np.argsort(mean_coords[:, 2])
    return mean_coords[order]

def mean_com(mask, group = "backbone"):
    fib_mean = None
    fib = protein.select_atoms(group)
    for ts in u.trajectory[mask]:
        com = fib.center_of_mass() / 10
        if fib_mean is None:
            fib_mean = com
            n = 1
        else:
            n += 1
            fib_mean += (com - fib_mean) / n
    
    return com[0], com[1], com[2]

# --- Center of the Fibril ---
x_cfp, y_cfp, _ = mean_com(maskt)

# Molecule Selections
mols = [14, 21, 28, 33, 27, 33]

fig = plt.figure(figsize=(6,6), dpi=150)
ax = fig.add_subplot(111, projection='3d')

cmap = plt.get_cmap("tab10")   # Coloring
colors = [cmap(j) for j in range(len(mols))]
for i, mol in enumerate(mols):
    molecule_mean = get_mean_coords_mol(mol, maskt, "name CA")

    x = molecule_mean[:, 0] - x_cfp
    y = molecule_mean[:, 1] - y_cfp
    z = molecule_mean[:, 2]

    ax.scatter(x, y, z, s=10, color=colors[i], edgecolor='k', alpha=0.1, label=f'Mol {mol}')


ax.set_xlim3d([-30, 30])
ax.set_ylim3d([-30, 30])
ax.set_zlim3d([0, 72])

# Centerline of fibril
z_line = np.linspace(0, 72, 100)
ax.plot(np.zeros_like(z_line), np.zeros_like(z_line), z_line, color='r', linestyle='dashed', linewidth=2)

plt.title("3D Visualization of Molecules")
plt.tight_layout()
plt.legend()
plt.show()

