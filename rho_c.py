import numpy as np
import plotly.graph_objects as go
import math
import sys
import tqdm

# equation to calculate density model from .cif file (hkl intensities aka Fhkl)
# p(xyz) = 1/V * sum(|Fhkl| .* e^(-i * (2pihx + 2piky + 2pilz))) ==> e^-it = cost - isint
#          1/V * sum(|Fhkl| .* (cos(2pihx + 2piky + 2pilz) - isin(2pihx + 2piky + 2pilz)))
# need to determine volume and do this operation for every value of hkl and xyz

def main():
    max_atoms = 30000

    # args 
    # 1: pdb relative file location
    # 2: max number of atoms (optional)
    if len(sys.argv) < 2:
        print("Please pass in a pdb file location as the (only) argument for this script!")
        quit()
    elif len(sys.argv) > 3:
        print("Too many arguments!")
        quit()

    f =  open(sys.argv[1])
    if len(sys.argv) > 2:
        max_atoms = int(sys.argv[2])

    atoms3d = []
    element_types = []
    num_atoms = 0

    # this loop reads in element types and xyz coords from PDB file
    for line in f: 
        if line[0:4] == 'ATOM':
            atoms3d.append([float(line[31:38]), float(line[39:46]), float(line[47:54])])
            element_types.append(line[13])
            num_atoms = num_atoms + 1
        if num_atoms > max_atoms:
            break

    # find min and max values of x, y, and z coords to create bounds in 3d space
    max_vals = np.amax(atoms3d, axis=0) 
    min_vals = np.amin(atoms3d, axis=0)
    cryst = max_vals - min_vals;
    print("cryst: " + str(cryst))
    # cryst = [round(int((max_vals[0] - min_vals[0]) / 10)) * 10 + 10, 
    #         round(int((max_vals[1] - min_vals[1]) / 10)) * 10 + 10, 
    #         round(int((max_vals[2] - min_vals[2]) / 10)) * 10 + 10]

    # will eventually be inputs to the method
    reso = 0.1;
    b = 1;
    pi = np.pi;

    # step is called 'grid' in the matlab code
    step = reso / 3

    # adjustments to step size (still don't totally understnd why this is necessary)
    ceil1 = math.ceil(cryst[0]/step) # 5 * 300 = 1500
    ceil2 = math.ceil(cryst[1]/step) # 5 * 300 = 1500
    ceil3 = math.ceil(cryst[2]/step)
    
    # are there cases where this actually changes the step size?
    stepXYZ = [ cryst[0] / ceil1, cryst[1] / ceil2, cryst[2] / ceil3];
    print("step sizes: " + str(stepXYZ))

    # create XYZ linespaces for matrices based on step sizes and calculated max/min values
    x_lin = np.linspace(min_vals[0], max_vals[0], 1 / stepXYZ[0])
    y_lin = np.linspace(min_vals[1], max_vals[1], 1 / stepXYZ[1])
    z_lin = np.linspace(min_vals[2], max_vals[2], 1 / stepXYZ[2])
    xv_3d, yv_3d, zv_3d = np.meshgrid(x_lin, y_lin, z_lin)
    print("x meshgrid dims: " + str(xv_3d.shape))

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

    # rho starts as all zeros and every atom iteration will add info to it 
    rho_3d = np.zeros((len(xv_3d), len(xv_3d[0]), len(xv_3d[0][0])))

    atommask = np.full((len(xv_3d), len(xv_3d[0]), len(xv_3d[0][0])), False, dtype=bool)

    # 3d equivalent to operations on rho, just add in Z dimension? 
    for ii in tqdm.tqdm(range(len(atoms3d))):
        atom = atoms3d[ii]
        scat = scatter_vals[element_types[ii]]

        d2X = ((xv_3d - atom[0] + cryst[0]/2) % cryst[0]) - (cryst[0] / 2)
        d2Y = ((yv_3d - atom[1] + cryst[1]/2) % cryst[1]) - (cryst[1] / 2)
        d2Z = ((zv_3d - atom[2] + cryst[2]/2) % cryst[2]) - (cryst[2] / 2)
        d2 = np.power(d2X, 2) + np.power(d2Y, 2) + np.power(d2Z, 2)
        d2_coef = -d2*pi*pi*4 # makes it a little faster

        atommask = np.logical_or(atommask, d2 < 3.2 * 3.2 * 3.2) # don't think this does anything

        # NOTE what data can be thrown away?
        rho_3d += (scat[0] * np.sqrt(pi*4/b) * np.exp(d2_coef/b) + 
                    scat[1] * np.sqrt(pi*4/(b+scat[5])) * np.exp(d2_coef/(b+scat[5])) + 
                    scat[2] * np.sqrt(pi*4/(b+scat[6])) * np.exp(d2_coef/(b+scat[6])) + 
                    scat[3] * np.sqrt(pi*4/(b+scat[7])) * np.exp(d2_coef/(b+scat[7])) + 
                    scat[4] * np.sqrt(pi*4/(b+scat[8])) * np.exp(d2_coef/(b+scat[8])))

    # doesn't seem to do anything...
    rho_3d = np.ma.array(rho_3d, mask = atommask) 
    fig1 = go.Figure(data=go.Volume(x=xv_3d.flatten(), y=yv_3d.flatten(), z=zv_3d.flatten(), 
                                     value=rho_3d.flatten(), isomin=0.1, isomax=1.0, 
                                     surface_count=5, opacity=0.5))
    fig1.show()

    # n-dimensionl fourier tranform to get into reciprocal space
    recip = np.fft.rfftn(rho_3d)

    # isolate real and reciprocal elements at every point
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
    reconverted_recip = reconverted_x + reconverted_y * 1j
    reconverted_model = np.fft.irfftn(reconverted_recip)

    fig2 = go.Figure(data=go.Volume(x=xv_3d.flatten(), y=yv_3d.flatten(), z=zv_3d.flatten(), 
                                     value=reconverted_model.flatten(), isomin=0.1, isomax=1.0, 
                                     surface_count=5, opacity=0.5))
    fig2.show()

main()
