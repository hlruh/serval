#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
river.py

@author: Henrik Ruh
9 Sept. 2024

Script to merge serval result files into image cubes of dimensions
(number of orders, number of spectra, number of pixels). 
An image cube is thus structured as a series of rivermaps.

The serval result files can be obtained by using the -outfmt and -outchi options of serval.
 

To make executable :
    chmod u+x $PATH/river.py
    ln -s $PATH/river.py  ~/bin/river
    
In command line write:
    river tag
    
"""

### imports
import argparse
import numpy as np
from glob import glob
from astropy.io import fits

### functions

# Define a list of available extensions
AVAILABLE_EXTENSIONS = ["wave", "waverest", "err", "fmod", "res", "spec", "bpmap", "ratio", "chi2","chi2norm","diffrv"]

# Process the fits files
def fits_to_river(folder, ext=None, savename="rivermap"):
    if ext is None:
        ext = ["wave", "waverest", "err", "fmod", "res", "spec", "bpmap", "ratio","chi2","chi2norm","diffrv"]
    else:
        ext = np.ravel([ext])

    for ext_name in ext:
        print('Create rivermap: '+ext_name)
        try:
            if ext_name == "chi2":
                filenames = np.sort(glob(folder + "/*_chi2map.fits"))
                cube = []

                # Append chi2 maps
                for fn in filenames:
                    with fits.open(fn) as hdul:
                        dat = np.copy(hdul[0].data)

                    cube.append(dat)
                    
                # Swap axes to get the desired cube format
                cube = list(np.array(cube).swapaxes(0, 1))

                # Append vgrid
                vgrid = []
                for fn in filenames:
                    with fits.open(fn) as hdul:
                        hdr = hdul[0].header
                        dat = np.copy(hdul[0].data)
                        nspec,nsteps = len(filenames),dat.shape[1]
                        v_lo = hdr["CRVAL1"] 
                        v_step = hdr["CDELT1"]
                        berv = hdr['HIERARCH SERVAL BERV']
                        vgrid.append(v_lo + v_step * np.arange(nsteps) - berv)
                cube = np.concatenate(([vgrid],cube))
                
                # Add chi2 sum
                cube = np.concatenate((cube,[np.nansum(cube[1:], axis=0)]))
                
            elif ext_name == "chi2norm":
                filenames = np.sort(glob(folder + "/*_chi2map.fits"))
                cube = []
                chi2 = []
                
                # Append chi2 maps
                for fn in filenames:
                    with fits.open(fn) as hdul:
                        dat = np.copy(hdul[0].data)
                    chi2.append(dat)
                    
                    # Normalize
                    for d in dat:
                        d /= np.nanmax(d)
                    cube.append(dat)
                    
                # Add chi2 sum
                chi2sum = np.nansum(chi2, axis=1)
                # Normalize
                for d in chi2sum:
                    d /= np.nanmax(d)
                        
                # Swap axes to get the desired cube format
                cube = np.array(cube).swapaxes(0, 1)
                cube = np.concatenate((cube,[chi2sum]))

                # Append vgrid
                vgrid = []
                for fn in filenames:
                    with fits.open(fn) as hdul:
                        hdr = hdul[0].header
                        dat = np.copy(hdul[0].data)
                        nspec,nsteps = len(filenames),dat.shape[1]
                        v_lo = hdr["CRVAL1"] 
                        v_step = hdr["CDELT1"]
                        berv = hdr['HIERARCH SERVAL BERV']
                        vgrid.append(v_lo + v_step * np.arange(nsteps) - berv)
                cube = np.concatenate(([vgrid],cube))
                
            elif ext_name == 'diffrv':
                filenames = np.sort(glob(folder + "/*_mod.fits"))
                fmod = []
                res = []
                wave = []
                bpmap = []
                for fn in filenames:
                    with fits.open(fn) as hdul:
                        fmod.append(hdul['fmod'].data)
                        res.append(hdul['res'].data)
                        wave.append(hdul['wave'].data)
                        bpmap.append(hdul['bpmap'].data)
                                                
                # Swap axes to get the desired cube format
                fmod = np.array(fmod).swapaxes(0, 1)
                res = np.array(res).swapaxes(0, 1)
                wave = np.array(wave).swapaxes(0, 1)
                bpmap = np.array(bpmap).swapaxes(0, 1)
                                            
                # compute gradient
                nord = len(wave)
                grad = np.array([np.gradient(fmod[o],wave[o,0],axis=1) for o in range(nord)])

                c=3e8
                cube = res/grad*c # save differential RV
                cube[bpmap!=0] = np.nan         
                
            else:
                filenames = np.sort(glob(folder + "/*_mod.fits"))
                cube = []
                for fn in filenames:
                    with fits.open(fn) as hdul:
                        cube.append(hdul[ext_name].data)
                        
                # Swap axes to get the desired cube format
                cube = list(np.array(cube).swapaxes(0, 1))

            # Create a new FITS file with the processed data
            hdu = fits.PrimaryHDU(cube)
            hdu.writeto(folder + '/' + savename + "_" + ext_name + ".fits", overwrite=True)
            
        except Exception as e:
            print(f"An error occurred for extension {ext_name}: {str(e)}")
            
# Define the command-line arguments
def get_command_line_arguments():
    parser = argparse.ArgumentParser(description="Convert serval result fits files to rivermaps")
    parser.add_argument("obj", type=str, help="Tag, output directory and file prefix (e.g. Object name).")
    parser.add_argument("--ext", type=str, nargs="+", default=AVAILABLE_EXTENSIONS, choices=AVAILABLE_EXTENSIONS, help="List of extensions to process")

    return parser.parse_args()

### MAIN

if __name__ == "__main__":
    # get command line arguments
    args = get_command_line_arguments()
    
    obj = args.obj
    if obj[-1] == '/':
        obj = obj[:-1]
        
    # load observations
    folder = obj+'/res'
   
    # create river maps
    fits_to_river(folder, args.ext)
