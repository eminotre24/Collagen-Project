# Cylindrical Transformation Approach module
# Module to store the functions used - keep code cleaner
import numpy as np
import matplotlib.pyplot as plt

def mean_coords_generation(universe, protein, mask, center = True):
    # Whole System
    mean_coords = None
    box_mean = None
    # Mean
    for ts in universe.trajectory[mask]:
        if center:
            # Center routine
            com = protein.center_of_mass()
            coords = universe.select_atoms("all").positions - np.array([com[0], com[1], 0.0])
        else:
            coords = universe.select_atoms("all").positions
        if mean_coords is None:
            mean_coords = coords.astype(float)
            box_mean = ts.dimensions.astype(np.float64, copy=True)
            n = 1
        else:
            n += 1
            mean_coords += (coords - mean_coords) / n
            box_mean += (ts.dimensions - box_mean) / n

    new_u = universe.copy()
    new_u.atoms.positions = mean_coords
    new_u.dimensions = box_mean

    return new_u

def get_mean_coords_mol(universe, protein, molecule, mask, group):
    # Get the Mean Coordinates of the molecule
    mean_coords = None
    chains = protein.fragments[3*molecule:3*molecule + 3]
    colmol = chains[0] + chains[1] + chains[2]
    colmol = colmol.select_atoms(group)
    for ts in universe.trajectory[mask]:
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

def get_centered_mean_coords_mol(universe, protein, molecule, mask, group):
    # Get the Mean Coordinates of the molecule - COM included each step - already centralized
    mean_coords = None
    chains = protein.fragments[3*molecule:3*molecule + 3]
    colmol = chains[0] + chains[1] + chains[2]
    colmol = colmol.select_atoms(group)
    for ts in universe.trajectory[mask]:
        com = protein.center_of_mass() / 10 # check as its not the same as the group
        coords_nc = colmol.positions / 10 # not centered
        com_p = np.array([com[0], com[1], 0.0]) # leave z untouched
        coords = coords_nc - com_p # centralize
        if mean_coords is None:
            mean_coords = coords.astype(float)
            n = 1
        else:
            n += 1
            mean_coords += (coords - mean_coords) / n
    
    # Ordering the z for unwrapping
    order = np.argsort(mean_coords[:, 2]) 
    return mean_coords[order]

def mean_com(universe, protein, mask, group):
    fib_mean = None
    fib = protein.select_atoms(group)
    for ts in universe.trajectory[mask]:
        com = fib.center_of_mass() / 10
        if fib_mean is None:
            fib_mean = com
            n = 1
        else:
            n += 1
            fib_mean += (com - fib_mean) / n
    
    return com[0], com[1], com[2]

def clean_phi(y, x):
    # Get the phi effectively getting rid of the jumps
    phi_p = np.atan2(y, x)
    phi_uw = np.unwrap(phi_p)
    return phi_uw

def cyl_proj(molecule, y_c = 0, x_c = 0, r_min = False):
    x = molecule[:, 0] - x_c
    y = molecule[:, 1] - y_c
    z = molecule[:, 2]

    if r_min == False:
        r = np.sqrt(y**2 + x**2)
        phi = clean_phi(y, x)
    else: # Filter near 0 r
        r = np.sqrt(y**2 + x**2)
        mask_rmin = (r > r_min)
        r = r[mask_rmin]
        x_filt = x[mask_rmin] 
        y_filt = y[mask_rmin] 
        phi = clean_phi(y_filt, x_filt)
        z = z[mask_rmin]

    return r, phi, z

def linfitandr2(x, y):
    fit_res = np.polynomial.polynomial.polyfit(x, y, 1, full=True)
    b, m = fit_res[0]
    # R2
    SSE = fit_res[1][0] 
    SST = np.sum((y - np.mean(y)) ** 2) 
    R2 = 1 - SSE/SST
    return m, b, R2

def remove_ends(percentage, molecule):
    z = molecule[:, 2]
    z_cut = (z.max() - z.min())*percentage/2
    # Filter
    maskz = (z >= z.min() + z_cut) & (z <= z.max() - z_cut)
    return molecule[maskz]

# --- Deprecated Routines ---
# Bin cataloging and visualization
def bin_cat_vis(rad_mol, phi_mol, z_mol, mol, bins):
    cat_bins = np.linspace(rad_mol.min(), rad_mol.max(), bins+1)
    bin_ids = np.digitize(rad_mol, cat_bins) - 1
    bin_rad = 0.5 * (cat_bins[:-1] + cat_bins[1:])

    # --- Build figure with GridSpec: left & middle span both rows; right column is split ---
    fig = plt.figure(figsize=(15, 8), dpi=150)

    gs = fig.add_gridspec(2, 3, height_ratios=[1, 0.4])

    ax0 = fig.add_subplot(gs[:, 0])   # left: per-chain phi(z)
    ax1 = fig.add_subplot(gs[0, 1])   # middle: combined phi(z) + fit
    ax2 = fig.add_subplot(gs[1, 1], sharex=ax1)   # top-right: z vs r
    ax3 = fig.add_subplot(gs[:, 2])  # bottom-right: histogram of r

    # --- Left: Angle Projection of Molecule & Fit ---
    ax0.scatter(phi_mol, z_mol, s=10, edgecolor='k', alpha=0.3, label="Data")
    m, b, R2 = linfitandr2(phi_mol, z_mol)
    fit = np.poly1d((m, b))
    ax0.plot(phi_mol, fit(phi_mol), linestyle='dashed', color='r', label=rf"Fit $R^2$={R2:.2f}")
    ax0.set_title("All chains combined")
    ax0.set_xlabel(r"$\varphi$")
    ax0.legend()

    # --- Middle: Radius & Histogram ---
    # Radius
    ax1.scatter(rad_mol, z_mol, s=10, color='g', edgecolor='k', alpha=0.3)
    ax1.set_xlim(0, np.sqrt(xd**2 + yd**2) / 2)
    ax1.set_title(r"$z(r)$ - radius of molecule along fibril")
    ax1.tick_params(labelbottom=False)

    # Histogram
    ax2.hist(rad_mol, bins=bins, color='g', edgecolor='k', alpha=0.7, density=True)
    ax2.set_xlabel(r"$r$ (nm)")
    ax2.set_ylabel("Ratio")
    ax2.set_xlim(ax2.get_xlim())

    # --- Right: Separation and fit of each Group ---
    cmap = plt.get_cmap("tab10")   # Coloring
    colors = [cmap(j) for j in range(bins)]
    for i in range(bins):
        maskr = bin_ids == i
        sample_phi = phi_mol[maskr]
        sample_z = z_mol[maskr]

        color = colors[i]

        ax3.scatter(sample_phi, sample_z, edgecolor='k', color=color, alpha=0.3)

        # Fit Analysis
        m, b, R2 = linfitandr2(sample_phi, sample_z)
        fit = np.poly1d((m, b))
        ax3.plot(sample_phi, fit(sample_phi), linestyle='dashed', color=color, label=rf"Fit Rad {bin_rad[i]:.2f}, $R^2$={R2:.2f}")

    ax0.set_title("Separated Radius Bins")
    ax0.set_xlabel(r"$\varphi$")
    ax3.legend()

    # --- Shared y label & overall title ---
    fig.supylabel(r"$z$ (nm)")
    plt.suptitle(f"Analysis of Molecule {mol} in relation with Fibril", fontsize=15)
    plt.tight_layout()
    plt.show()