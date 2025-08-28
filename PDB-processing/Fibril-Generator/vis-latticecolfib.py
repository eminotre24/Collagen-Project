# We can see the lattice structure (what we care of, that is the shape in x,y), with a projection in that plane
import numpy as np
import matplotlib.pyplot as plt
import gemmi as gm

# Read collagen crystal structure
path = "./"
filepath = path + "col.pdb"

col_diam = 15 #Arms
st = gm.read_structure(filepath)
cell = st.cell

# Lattice descriptors
a = gm.Fractional(1, 0, 0)
b = gm.Fractional(0, 1, 0)
c = gm.Fractional(0, 0, 1)

axy = np.array(cell.orthogonalize(a).tolist()[0:2])
bxy = np.array(cell.orthogonalize(b).tolist()[0:2])
cxy = np.array(cell.orthogonalize(c).tolist()[0:2])

def proj(basis, point):
    ab, bb, cb = basis
    return ab*point[0] + bb*point[1] + cb*point[2]

# --- Visual of Lattice Construction: Primitive Vectors & Diameter of Molecule ---
scale_factor = 1

points = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]]
'''
points = [[0,0,0], [1,0,0], [0,1,0], [-1,0,0], [0,-1,0],
                  [1,1,0], [-1,1,0], [1,-1,0], [-1,-1,0],

                  [0,0,1], [1,0,1], [-1,0,1], [0,1,1], [0,-1,1],
                  [1,1,1], [-1,1,1], [1,-1,1], [-1,-1,1],

                  [0,0,-1], [1,0,-1], [-1,0,-1], [0,1,-1], [0,-1,-1],
                  [1,1,-1], [-1,1,-1], [1,-1,-1], [-1,-1,-1]]
'''

num_points = len(points)
coords = np.zeros((num_points, 2))

for i, point in enumerate(points):
    cart_point = proj((axy, bxy, cxy), point)
    coords[i, :] = scale_factor*cart_point

fig, ax = plt.subplots(figsize=(5,5), dpi=200)

colmol = plt.Circle((0, 0), col_diam/2, color='g', fill=False)
ax.add_patch(colmol)

ax.scatter(coords[:, 0], coords[:, 1], color='k', alpha=0.5, edgecolor = 'k', zorder=3)

# Add text
'''
for i in range(num_points):
    ax.text(coords[i, 0], coords[i, 1], f"[{points[i][0]}, {points[i][1]}, {points[i][2]}]", fontsize=6, ha='center', va='bottom')
'''

ax.axis('equal')
plt.xlabel(r'x ($\AA$)')
plt.ylabel(r'y ($\AA$)')
plt.title('XY Projection of Lattice Points - Collagen Fibril')
plt.grid(True)
plt.show()