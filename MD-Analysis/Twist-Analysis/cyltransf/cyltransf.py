# Cylindrical Transformation Approach module
# Module to store the functions used - keep code cleaner

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

def clean_phi(y, x):
    # Get the phi effectively getting rid of the jumps
    phi_p = np.atan2(y, x)
    phi_uw = np.unwrap(phi_p)
    return phi_uw

def cyl_proj(molecule, y_c, x_c, r_min = False):
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
    fit_res = np.polyfit(x, y, 1, full=True) 
    m, b = fit_res[0]
    # R2
    SSE = fit_res[1][0] 
    SST = np.sum((y - np.mean(y)) ** 2) 
    R2 = 1 - SSE/SST
    return m, b, R2