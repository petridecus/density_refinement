import numpy as np
import plotly.graph_objects as go
import math
import sys
import tqdm

'''
https://www.xtal.iqfr.csic.es/Cristalografia/parte_07-en.html
equation to calculate density model from .cif file (hkl intensities aka Fhkl)
p(xyz) = 1/V * sum(|Fhkl| .* e^(-i * (2pihx + 2piky + 2pilz))) ==> e^-it = cost - isint
         1/V * sum(|Fhkl| .* (cos(2pihx + 2piky + 2pilz) - isin(2pihx + 2piky + 2pilz)))
need to determine volume and do this operation for every value of hkl and xyz

simplified when unit cell is "centrosymmetric"
2/V * sum(|Fhkl| .* cos2pi(hx + ky + lz)) - so basically it gets rid of the sin part of the equation? 

"Patterson function"
1/v * sum(|Fhkl|^2 .* cos2pi(hu + kv + lw)) -> this actually creates inaccurate atom coordinates unless processed
'''

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

xv = []
yv = []
zv = []

cryst = []
# reso = 0.05

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
        # TODO figure out if resolution is supposed to be this bad
        if '_diffrn_reflns.pdbx_d_res_high' in line:
            # reso = float(line[33:]) / 10
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

    return atoms, element_types

# atoms: coordinates of atoms read from pdb
# element_types: element types of atoms read from pdb
# crystal: dimensions in 3D of crystal
# b: b factor(s) #TODO
def real_space_density(atoms, element_types, cryst, b):
    # rho starts as all zeros and every atom iteration will add info to it 
    rho_3d = np.zeros((len(xv), len(xv[0]), len(xv[0][0])))
    atommask = np.full((len(xv), len(xv[0]), len(xv[0][0])), False, dtype=bool)
    pi = np.pi

    # 3d equivalent to operations on rho, just add in Z dimension? 
    for ii in tqdm.tqdm(range(len(atoms)), "Calculating diffraction data from atom cooords"):
        atom = atoms[ii]
        scat = scatter_vals[element_types[ii]]

        d2X = ((xv - atom[0] + cryst[0]/2) % cryst[0]) - (cryst[0] / 2)
        d2Y = ((yv - atom[1] + cryst[1]/2) % cryst[1]) - (cryst[1] / 2)
        d2Z = ((zv - atom[2] + cryst[2]/2) % cryst[2]) - (cryst[2] / 2)
        d2 = np.power(d2X, 2) + np.power(d2Y, 2) + np.power(d2Z, 2)
        d2_coef = -d2*pi*pi*4 # makes it a little faster

        atommask = np.logical_or(atommask, d2 < 3.2 * 3.2 * 3.2) # don't think this does anything

        # NOTE what data can be thrown away?
        rho_3d += (scat[0] * np.sqrt(pi*4/b) * np.exp(d2_coef/b) + 
                    scat[1] * np.sqrt(pi*4/(b+scat[5])) * np.exp(d2_coef/(b+scat[5])) + 
                    scat[2] * np.sqrt(pi*4/(b+scat[6])) * np.exp(d2_coef/(b+scat[6])) + 
                    scat[3] * np.sqrt(pi*4/(b+scat[7])) * np.exp(d2_coef/(b+scat[7])) + 
                    scat[4] * np.sqrt(pi*4/(b+scat[8])) * np.exp(d2_coef/(b+scat[8])))
    
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

    global cryst, reso

    # find min and max values of x, y, and z coords to create bounds in 3d space
    max_vals = np.amax(atoms3d, axis=0) 
    min_vals = np.amin(atoms3d, axis=0)
    cryst = [round(int((max_vals[0] - min_vals[0]) / 10)) * 10 + 10, 
             round(int((max_vals[1] - min_vals[1]) / 10)) * 10 + 10, 
             round(int((max_vals[2] - min_vals[2]) / 10)) * 10 + 10]

    reso = 0.05
    b = 1

    if len(sys.argv) == 3:
        cif_lines = read_cif(sys.argv[2])
    
    step = reso / 3

    global xv, yv, zv
    x_lin = np.linspace(0, int(cryst[0]), int(1 / step)) # int(min_vals[0]), int(max_vals[0]), int(1 / stepXYZ[0]))
    y_lin = np.linspace(0, int(cryst[1]), int(1 / step)) # int(min_vals[1]), int(max_vals[1]), int(1 / stepXYZ[1]))
    z_lin = np.linspace(0, int(cryst[2]), int(1 / step)) # int(min_vals[2]), int(max_vals[2]), int(1 / stepXYZ[2]))
    xv, yv, zv = np.meshgrid(x_lin, y_lin, z_lin)
    
    rho_3d, atommask = real_space_density(atoms3d, element_types, cryst, b)

    # doesn't seem to do anything...
    rho_3d = np.ma.array(rho_3d, mask = atommask) 
    fig1 = go.Figure(data=go.Volume(x=xv.flatten(), y=yv.flatten(), z=zv.flatten(), 
                                     value=rho_3d.flatten(), isomin=0.1, isomax=1.0, 
                                     surface_count=5, opacity=0.5))
    fig1.show()

    '''
    # n-dimensionl fourier tranform to get into reciprocal space
    recip = np.fft.rfftn(rho_3d)

    # isolate real and imaginary elements at every point
    recip_x = recip.real
    recip_y = recip.imag

    # use sqrt of x^2 + y^2 to get intensities at every point
    intensities = np.sqrt(np.square(recip_x) + np.square(recip_y))

    # use intensities and x values to get the phases
    phases = np.arccos(recip_x / intensities)

    # now do the opposite of all of these operations to get back to the original model
    reconverted_x = intensities * np.cos(phases)
    reconverted_y = intensities * np.sin(phases)
    
    # reconverted_y specifically is what is breaking everything right now :/
    # the phase data must be ok if reconverted_x + recip_y * 1j works... what's the difference??
    reconverted_recip = reconverted_x + recip_y * 1j # reconverted_y * 1j
    reconverted_model = np.fft.irfftn(reconverted_recip)

    fig2 = go.Figure(data=go.Volume(x=xv.flatten(), y=yv.flatten(), z=zv.flatten(), 
                                     value=reconverted_model.flatten(), isomin=0.1, isomax=1.0, 
                                     surface_count=5, opacity=0.5))
    fig2.show()
    '''

main()

# adjustments to step size (still don't totally understnd why this is necessary)
# ceil1 = math.ceil(cryst[0]/step) # 5 * 300 = 1500
# ceil2 = math.ceil(cryst[1]/step) # 5 * 300 = 1500
# ceil3 = math.ceil(cryst[2]/step)
    
# are there cases where this actually changes the step size?
# stepXYZ = [ cryst[0] / ceil1, cryst[1] / ceil2, cryst[2] / ceil3];

# create XYZ linespaces for matrices based on step sizes and calculated max/min values
