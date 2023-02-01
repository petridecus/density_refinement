Basic tool for refinining Electron Density model.
 1. Reciprocal space data:
    - Reads in reciprocal space data from a .cif file.
    - Calculates a real space density based on intensities recorded thru x-ray diffraction (which contains no phase data).
 2. Real space data:
    - Reads in atom coordinate and element type data from a .pdb file.
    - Calculates another real space density based on atom coordinates of model in PDB.
 3. Rephasing:
    - Isolates the magnitudes of Electron density from the x-ray datas model, by returning to fourier space and converting to polar coords.
    - Isolates the PHASES of Electron density model from the real space data, by converting to fourier space and taking angles from polar coords.
    - Combines magnitudes from x-ray data, angles from real-space data, converts back to complex data and does one final FFT to calculate rephased model!
