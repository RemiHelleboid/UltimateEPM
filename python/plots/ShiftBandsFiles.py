import numpy as np
import glob
import os, sys
import re

def get_bands_file(prefix: str, suffix: str = ".csv"):
    files = glob.glob(prefix + "*" + suffix)
    # Extract the band index from the file name
    # band_index = [int(re.findall(r'\d+', file)[-1]) for file in files]
    for i in range(len(files)):
        file = files[i]
        if re.findall(r'\d+', file):
            band_index = int(re.findall(r'\d+', file)[-1])
            files[i] = (band_index, file)
        else:
            continue

    # Sort the files by band index
    files.sort(key=lambda x: x[0])
    # Get the files names
    files = [file[1] for file in files]
    print("Files: ", files)
    return files

def shift_bands_to_zero():
    IndexFirst_CB_band = 4
    IndexFirst_VB_band = 3

    # Get the bands files
    bands_files = get_bands_file("band_")
    # Get the number of bands
    num_bands = len(bands_files)

    # Get the band gap
    file_first_cb = bands_files[IndexFirst_CB_band]
    file_first_vb = bands_files[IndexFirst_VB_band]
    print("First CB file: ", file_first_cb)
    print("First VB file: ", file_first_vb)
    CB_1 = np.loadtxt(bands_files[IndexFirst_CB_band], unpack=True)
    VB_3 = np.loadtxt(bands_files[IndexFirst_VB_band], unpack=True)
    band_gap = np.min(CB_1) - np.max(VB_3)
    print("Band gap: ", band_gap)

    # Check that the max of the VB is 0
    VB_max = np.max(VB_3)
    if VB_max != 0:
        print("Error: VB max is not 0")
        sys.exit(1)
    
    # Shift the CB bands by the band gap
    for i in range(IndexFirst_CB_band, num_bands):
        print("CB file: ", bands_files[i])
        CB_i = np.loadtxt(bands_files[i], unpack=True)
        CB_i = CB_i - band_gap

        new_filename = f"Si100Ge000_CB_{i-IndexFirst_CB_band+1}_EPM"
        np.savetxt(new_filename, CB_i, fmt="%.6f", delimiter=" ")
    
    # Export the VB bands
    for i in range(IndexFirst_VB_band+1):
        print("VB file: ", bands_files[i])
        VB_i = np.loadtxt(bands_files[i], unpack=True)
        new_filename = f"Si100Ge000_VB_{i+1}_EPM"
        np.savetxt(new_filename, VB_i, fmt="%.6f", delimiter=" ")

    return

def check_bands():
    bands_files_cb = get_bands_file("Si100Ge000_CB_*", suffix="_EPM")
    bands_files_vb = get_bands_file("Si100Ge000_VB_*", suffix="_EPM")
    print("CB files: ", bands_files_cb)
    print("VB files: ", bands_files_vb)
    for file in bands_files_cb:
        print("CB file: ", file)
        CB_i = np.loadtxt(file, unpack=True)
        print("Min CB: ", np.min(CB_i))
        print("Max CB: ", np.max(CB_i))
    for file in bands_files_vb:
        print("VB file: ", file)
        VB_i = np.loadtxt(file, unpack=True)
        print("Min VB: ", np.min(VB_i))
        print("Max VB: ", np.max(VB_i))
    return

if __name__ == "__main__":
    shift_bands_to_zero()
    check_bands()
    
