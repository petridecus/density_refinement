import numpy as np
import plotly.graph_objects as go
import math
import sys
import tqdm

# this will be different array of scatter factors for every atom type
# this appears to be arranged as:
# [c, a1, a2, a3, a4, b1, b2, b3, b4]
# http://lampx.tugraz.at/~hadley/ss1/crystaldiffraction/atomicformfactors/formfactors.php
h_scat = [0.001305, 0.489918, 0.262003, 
            0.196767, 0.049879, 20.6593,
            7.74039, 49.5519, 2.20159]
c_scat = [0.215600, 2.310000, 1.020000, 
            1.588600, 0.865000, 20.843899, 
            10.207500, 0.568700, 51.651199]
n_scat = [-11.529, 12.2126, 3.1322, 
            2.0125, 1.1663, 0.0057, 
            9.8933, 28.9975, 0.5826]
o_scat = [0.2508, 3.0485, 2.2868, 
            1.5463, 0.867, 13.2771, 
            5.7011, 0.3239, 32.9089]
s_scat = [0.8669, 6.9053, 5.2034,
            1.4379, 1.5863, 1.4679,
            22.2151, 0.2536, 56.172]

scatter_vals = dict(H = h_scat, C = c_scat, N = n_scat, O = o_scat, S = s_scat)

# globals 
lin_r = [[], [], []]
lin_f = [[], [], []]
cryst = [0, 0, 0]
coord_dims = [[], []]
reso = 5

# reads crystal size and resolution from cif file
def read_cif(filename: str):
    global cryst, reso
    lines = open(filename).readlines()

    num_datapoints = 0

    for ii in tqdm.tqdm(range(len(lines)), "Parsing diffraction data from cif file"):
        line = lines[ii]
        if '_cell.length_a' in line:
            cryst[0] = float(line[20:])
            continue
        if '_cell.length_b' in line:
            cryst[1] = float(line[20:])
            continue
        if '_cell.length_c' in line:
            cryst[2] = float(line[20:])
            continue
        if '_diffrn_reflns.pdbx_d_res_high' in line:
            reso = float(line[33:]) / 2
            continue
        if '_diffrn_reflns.number' in line:
            num_datapoints = int(line[33:])

    data_starting_line = len(lines) - (num_datapoints + 1)
    data_lines = lines[data_starting_line:]

    return data_lines


# reads in 'max_atoms' element types and coordinates from given file 'f' 
def read_pdb(filename: str, max_atoms: int = 30000):
    atoms = []
    element_types = []
    num_atoms = 0

    lines =  open(filename).readlines()
    for ii in tqdm.tqdm(range(len(lines)), "Parsing PDB file"): 
        line = lines[ii]
        if line[0:4] == 'ATOM':
            atoms.append([float(line[31:38]), float(line[39:46]), float(line[47:54])])
            element_types.append(line[13])
            num_atoms = num_atoms + 1
        if num_atoms > max_atoms:
            break

    global coord_dims
    coord_dims[0] = np.amin(atoms, axis=0)
    coord_dims[1] = np.amax(atoms, axis=0)

    return atoms, element_types

'''
https://www.xtal.iqfr.csic.es/Cristalografia/parte_07-en.html
equation to calculate density model from .cif file (hkl intensities aka Fhkl)
p(xyz) = 1/V * sum(|Fhkl| .* e^(-i * (2pihx + 2piky + 2pilz))) ==> e^-it = cost - isint
         1/V * sum(|Fhkl| .* (cos(2pihx + 2piky + 2pilz) - isin(2pihx + 2piky + 2pilz)))
need to determine volume and do this operation for every value of hkl and xyz

simplified when unit cell is "centrosymmetric"
2/V * sum(|Fhkl| .* cos2pi(hx + ky + lz)) - so basically it gets rid of the sin part of the equation?
'''
def reciprocal_space_density(cif_data):
    xv_f, yv_f, zv_f = np.meshgrid(lin_f[0], lin_f[1], lin_f[2])
    f_c = np.zeros((len(xv_f), len(xv_f[0]), len(xv_f[0][0])))
    
    # create linespaces to iterate thru xyz
    # for each xyz point calculate p value based on equation from site
    cryst_vol = cryst[0] * cryst[1] * cryst[2]

    for ii in tqdm.tqdm(range(len(cif_data) - 1), "Creating density map from raw crystal data"):
        data = cif_data[ii].split()

        if len(data) < 7:
            print("short line: " + str(data))
            continue

        if data[6] != 'o' and data[6] != 'f':
            continue
    
        # h k l are data[3-5]
        hkl_coef = 2 * np.pi * (float(data[3]) * xv_f + float(data[4]) * yv_f + float(data[5]) * zv_f)
        complex_val = np.cos(hkl_coef) - 1j* np.sin(hkl_coef)
        f_c += (abs(float(data[7])) * hkl_coef)
    
    f_c = f_c / cryst_vol

    f_realspace = np.fft.irfftn(f_c, f_c.shape) # , f_c.shape)
    fig1 = go.Figure(data=go.Volume(x=xv_f.flatten(), y=yv_f.flatten(), z=zv_f.flatten(), 
                                     value=f_realspace.flatten(), isomin=0.2, isomax=1.0, 
                                     surface_count=5, opacity=0.5))
    fig1.show()

    return f_c 

# atoms: coordinates of atoms read from pdb
# element_types: element types of atoms read from pdb
# crystal: dimensions in 3D of crystal
# b: b factor(s) #TODO
def real_space_density(atoms, element_types, b):
    global cryst

    max_vals = coord_dims[1] - coord_dims[0]
    xv_r, yv_r, zv_r = np.meshgrid(lin_r[0], lin_r[1], lin_r[2])
    xv_f, yv_f, zv_f = np.meshgrid(lin_f[0], lin_f[1], lin_f[2])

    # rho starts as all zeros and every atom iteration will add info to it 
    rho_3d = np.zeros((len(xv_r), len(xv_r[0]), len(xv_r[0][0]))) # each meshgrid array has correct dimensions so we can use xv
    atommask = np.full((len(xv_r), len(xv_r[0]), len(xv_r[0][0])), False, dtype=bool)

    rho_3d = np.zeros((len(xv_f), len(xv_f[0]), len(xv_f[0][0]))) # each meshgrid array has correct dimensions so we can use xv
    atommask = np.full((len(xv_f), len(xv_f[0]), len(xv_f[0][0])), False, dtype=bool)
    
    pi = np.pi

    # 3d equivalent to operations on rho, just add in Z dimension? 
    for ii in tqdm.tqdm(range(len(atoms)), "Calculating diffraction data from atom cooords"):
        atom = atoms[ii]
        scat = scatter_vals[element_types[ii]]

        # option to plot on real space axes based on atom coordinates- as of now this doesn't play nice w other axes
        # d2X = ((xv_r - atom[0] + max_vals[0]/2) % max_vals[0]) - (max_vals[0] / 2)
        # d2Y = ((yv_r - atom[1] + max_vals[1]/2) % max_vals[1]) - (max_vals[1] / 2)
        # d2Z = ((zv_r - atom[2] + max_vals[2]/2) % max_vals[2]) - (max_vals[2] / 2)

        # NOTE - swapping Y and Z coords is something specific to 3cad data
        d2X = ((xv_f - atom[0] + cryst[0]/2) % cryst[0]) - (cryst[0] / 2)
        d2Y = ((yv_f - atom[2] + cryst[1]/2) % cryst[1]) - (cryst[1] / 2)
        d2Z = ((zv_f - atom[1] + cryst[2]/2) % cryst[2]) - (cryst[2] / 2)
        
        d2 = np.power(d2X, 2) + np.power(d2Y, 2) + np.power(d2Z, 2)
        d2_coef = -d2*pi*pi*4 # makes it a little faster

        atommask = np.logical_or(atommask, d2 < 3.2 * 3.2 * 3.2) # don't think this does anything

        # NOTE what data can be thrown away?
        rho_3d += (scat[0] * np.sqrt(pi*4/b) * np.exp(d2_coef/b) + 
                    scat[1] * np.sqrt(pi*4/(b+scat[5])) * np.exp(d2_coef/(b+scat[5])) + 
                    scat[2] * np.sqrt(pi*4/(b+scat[6])) * np.exp(d2_coef/(b+scat[6])) + 
                    scat[3] * np.sqrt(pi*4/(b+scat[7])) * np.exp(d2_coef/(b+scat[7])) + 
                    scat[4] * np.sqrt(pi*4/(b+scat[8])) * np.exp(d2_coef/(b+scat[8])))
    
    fig1 = go.Figure(data=go.Volume(x=xv_f.flatten(), y=yv_f.flatten(), z=zv_f.flatten(), 
                                     value=rho_3d.flatten(), isomin=0.2, isomax=1.0, 
                                     surface_count=5, opacity=0.5))
    fig1.show()

    return rho_3d, atommask


def main():
    # args 
    # 1: pdb relative file location
    # 2: max number of atoms (optional)
    if len(sys.argv) < 2:
        print("Please pass in a pdb file location as the (only) argument for this script!")
        quit()
    elif len(sys.argv) > 3:
        print("Too many arguments!")
        quit()

    atoms3d, element_types = read_pdb(sys.argv[1])
    b = 1

    if len(sys.argv) == 3:
        cif_lines = read_cif(sys.argv[2])

    # create XYZ linespaces in real space based on atom coords
    global lin_r
    lin_r[0] = np.linspace(coord_dims[0][0], coord_dims[1][0], math.ceil((coord_dims[1][0] - coord_dims[0][0]) / reso))
    lin_r[1] = np.linspace(coord_dims[0][1], coord_dims[1][1], math.ceil((coord_dims[1][1] - coord_dims[0][1]) / reso))
    lin_r[2] = np.linspace(coord_dims[0][2], coord_dims[1][2], math.ceil((coord_dims[1][2] - coord_dims[0][2]) / reso))

    # create XYZ linespaces in fourier space based on crystal dims
    global lin_f
    lin_f[0] = np.linspace(0, cryst[0], math.ceil(cryst[0] / reso))
    lin_f[1] = np.linspace(0, cryst[1], math.ceil(cryst[1] / reso))
    lin_f[2] = np.linspace(0, cryst[2], math.ceil(cryst[2] / reso))
    
    f_c = reciprocal_space_density(cif_lines)
    f_realspace = np.fft.ifftn(f_c, f_c.shape)

    rho_3d, atommask = real_space_density(atoms3d, element_types, b)
    rho_3d = np.ma.array(rho_3d, mask = atommask) 
    
    # n-dimensionl fourier tranform to get into reciprocal space
    r_recip = np.fft.fftn(rho_3d, rho_3d.shape)
    f_recip = f_c # np.fft.fftn(f_c, rho_3d.shape)

    # isolate real and imaginary elements at every point
    r_recip_x = r_recip.real
    r_recip_y = r_recip.imag

    f_recip_x = f_recip.real
    f_recip_y = f_recip.imag

    # use sqrt of x^2 + y^2 to get intensities at every point
    intensities_exp = np.sqrt(np.square(r_recip_x) + np.square(r_recip_y))
    intensities_obs = np.sqrt(np.square(f_recip_x) + np.square(f_recip_y))

    # use intensities and x values to get the phases
    phases_exp = np.arccos(r_recip_x / intensities_exp)

    # now do the opposite of all of these operations to get back to the original model
    combined_x = intensities_obs * np.cos(phases_exp)
    combined_y = intensities_obs * np.sin(phases_exp)
    
    # reconverted_y specifically is what is breaking everything right now :/
    # the phase data must be ok if reconverted_x + recip_y * 1j works... what's the difference??
    combined_recip = combined_x + combined_y * 1j
    combined_model = np.fft.irfftn(combined_recip, rho_3d.shape)
    combined_model = np.abs(combined_model) / 10 # 000
    print(combined_model)

    xv_f, yv_f, zv_f = np.meshgrid(lin_f[0], lin_f[1], lin_f[2])
    '''
    real_space_aligned = np.zeros((len(xv_r), len(xv_r[0]), len(xv_r[0][0])))
    for xx in range(len(xv_r)):
        for yy in range(len(xv_r[0])):
            for zz in range(len(xv_r[0][0])):
                real_space_aligned[xx][yy][zz] += combined_model[xx][yy][zz]
    '''
    fig3 = go.Figure(data=go.Volume(x=xv_f.flatten(), y=yv_f.flatten(), z=zv_f.flatten(),
                                     value=combined_model.flatten(), isomin=0.2, isomax=1.0,
                                     surface_count=5, opacity=0.5))
    fig3.show()

main()
