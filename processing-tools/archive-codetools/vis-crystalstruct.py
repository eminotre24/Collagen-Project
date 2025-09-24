import numpy as np
import matplotlib.pyplot as plt
import gemmi as gm

path = "/Users/alejandro/Desktop/mddocs/processing-analysis/PDBs-Codes/fibrilGen/"
filepath = path + "collagen-molecule.pdb"

st = gm.read_structure(filepath)
cell = st.cell

# Fibril Structure
'''
lattice_points = [[0,0,0], [1,0,0], [0,1,0], [-1,0,0], [0,-1,0],
                  [1,1,0], [-1,1,0], [1,-1,0], [-1,-1,0],

                  [0,0,1], [1,0,1], [-1,0,1], [0,1,1], [0,-1,1],
                  [1,1,1], [-1,1,1], [1,-1,1], [-1,-1,1],

                  [0,0,-1], [1,0,-1], [-1,0,-1], [0,1,-1], [0,-1,-1],
                  [1,1,-1], [-1,1,-1], [1,-1,-1], [-1,-1,-1]]
'''
lattice_points = [[0,0,0], [1,0,0], [0,1,0], [0,0,1]]


x = np.zeros(len(lattice_points))
y = np.zeros(len(lattice_points))
z = np.zeros(len(lattice_points))

# Generation of Points
for i, point in enumerate(lattice_points):
    frac = gm.Fractional(point[0], point[1], point[2])
    cart = cell.orthogonalize(frac)
    x[i], y[i], z[i] = cart.tolist()

# Vis
fig = plt.figure(figsize=(6,6), dpi=200)
ax = fig.add_subplot(111, projection='3d')

ax.scatter(z, y, x, s=10, color='black', edgecolor='k', alpha=0.7)

# --- Add labels beside each point ---
'''
for i, (xi, yi, zi) in enumerate(zip(x, y, z)):
    ax.text(zi, yi, xi, f"{lattice_points[i]}", fontsize=8, color="b")
'''

# --- Visualization Primitive Vectors ---
#'''
basis = [[1,0,0], [0,1,0], [0,0,1]]
labels = ['a', 'b', 'c']
for vec, label in zip(basis, labels):
    fraci = gm.Fractional(vec[0], vec[1], vec[2])
    carti = cell.orthogonalize(fraci)
    xi, yi, zi = carti.tolist()
    dx, dy, dz = xi - x[0], yi - y[0], zi - z[0]
    ax.quiver(z[0], y[0], x[0], dz, dy, dx,
              arrow_length_ratio=0.1, color='r', linewidth=2)
    ax.text(zi/2, yi, xi, label, fontsize=12, color='r')
#'''


ax.set_xlim3d([-700, 700])
ax.set_ylim3d([-700, 700])
ax.set_zlim3d([-700, 700])

plt.title("Visualization of Lattice Points (Cartesian Coordinates)")
ax.set_xlabel(r'Z ($\AA$)')
ax.set_ylabel(r'Y ($\AA$)')
ax.set_zlabel(r'X ($\AA$)')
plt.tight_layout()
plt.show()

print(cell)