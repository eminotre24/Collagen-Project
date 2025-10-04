# Cylindrical Transformation Approach module
# Module to store the functions used - keep code cleaner
import numpy as np
from scipy import odr

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

def lin_mdl(beta, x):
    m, c = beta
    return m * x + c

def ortreglinfit(x, y, beta0 = [0,0]):
    fun = odr.Model(lin_mdl)
    data = odr.Data(x, y)
    odr_run = odr.ODR(data, fun, beta0=beta0) # Guess Parameters
    res = odr_run.run()
    m_res, b_res = res.beta
    r_var = res.res_var
    return m_res, b_res, r_var

def remove_ends(percentage, molecule):
    z = molecule[:, 2]
    z_cut = (z.max() - z.min())*percentage/2
    # Filter
    maskz = (z >= z.min() + z_cut) & (z <= z.max() - z_cut)
    return molecule[maskz]