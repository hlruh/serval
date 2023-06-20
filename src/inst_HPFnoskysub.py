#! /usr/bin/env python

# no sky flux , sky wavl present in pascals tac files so ignore this part

from read_spec import *
from astropy.time import Time
import numpy as np
import cspline as spl

# Instrument parameters
name = 'HPF'
obsname = "hpf" # for barycorrpy
       # wikipedia
obsloc = dict(lat=30.681444,    # 30d40'53" N
              lon=-104.014722,  # 104d00'53" W
              elevation=2026)
R = 53000 # resolving power

pat = '*.fits'
pmax = 2048 - 300
iomax = 28
oset = "[4,5,6,14,15,16,17,18,26]"
coset = "[4,5,6,14,15,16,17,18,26]"

maskfile = 'telluric_mask_carm_short.dat'
skyfile = 'sky_mask_5_sigma.dat'
blazefile = 'hpf_blaze_spectra.fits'  # https://github.com/grzeimann/Goldilocks_Documentation/blob/master/hpf_blaze_spectra.fits


# Instrument read functions
def scan(self, s, orders=None, pfits=True, verb=True):
   """
   Returns
   -------
   namedtuple('spectrum', 'w f berv bjd blaze drift timeid sn55 ')
           w    - wavelength
           f    - flux
           berv - Barycentric Earth Radial Velocity
           bjd  - Barycentric Julian Day
           blaze - Blaze filename
           drift - Used RV Drift
           sn55  - S_N order center55
   Example
   -------
   >>> read_carm_vis(filename)

   """
   HIERARCH = 'HIERARCH '
   hdulist = self.hdulist = pyfits.open(s) # slow 30 ms
   self.header = hdr = hdulist[0].header
   self.instname = hdr['INSTRUME']
   self.drsberv = hdr.get('BERV', np.nan)
   self.drsbjd = hdr.get('BJD', np.nan) + 2400000
   self.dateobs = hdr['DATE-OBS']
   self.mjd = Time(self.dateobs, format='isot', scale='utc').mjd
   self.sn55 = hdr.get('SNR', 50)
   self.fileid = hdr.get('DATE-OBS', 0) #fileid[fileid.index('(')+1:fileid.index(')')]
   self.timeid = self.fileid

   self.calmode = "%s,%s,%s" % (hdr.get('SCI-OBJ', ''), hdr.get('CAL-OBJ', ''), hdr.get('SKY-OBJ', ''))

   self.ra = hdr['RA']
   self.de = hdr['DEC']
   self.airmass = hdr.get('AIRMASS', np.nan)
   self.exptime = hdr['ITIME']

   # exclude files from files depending if the observation was charged
   # this can be internally decided via snr
   # a different flag might be useful
   #self.charged = hdr['CHARGED']
   #if not self.charged:
   #    self.flag |= sflag.user
   
   if 'Goldilocks' in self.filename:
      self.drift = hdr.get(HIERARCH+'LRVCORR', np.nan)
      # no error given (https://github.com/grzeimann/Goldilocks_Documentation)
      self.e_drift = 0
      
      self.tmmean = 0.5
   
   elif 1:
      # for HPF spectra the drift is already included in the wavelength solution
      self.drift = hdr.get(HIERARCH+'CARACAL DRIFT FP RV', hdr.get(HIERARCH+'CARACAL DRIFT RV', np.nan))
      self.e_drift = hdr.get(HIERARCH+'CARACAL DRIFT FP E_RV', hdr.get(HIERARCH+'CARACAL DRIFT RVERR', np.nan))

      self.tmmean = hdr.get(HIERARCH+'CARACAL TMEAN', 0.0)
      if self.exptime: self.tmmean /= self.exptime   # normalise
      if self.tmmean == 0: self.tmmean = 0.5

   if self.mjd > 59731:
       # add drift offset post downtime (estimate based on sample)
      self.drift = self.drift - 60

def data(self, orders, pfits=True):
   if 'Goldilocks' in self.filename:
      return data_goldilocks(self, orders, pfits=True)
   elif 1:
      return data_psu(self, orders, pfits=True)
      
def data_psu(self, orders, pfits=True):
   hdulist = self.hdulist
   if 1:  # read order data
      # science
      w,f,e = hdulist['Sci Wavl'].data, hdulist['Sci Flux'].data, hdulist['Sci Variance'].data
      w,f,e = w[orders],f[orders],e[orders]
      e = e**0.5        

      f, e = deblaze(f,e,orders) # deblaze spectra
      sky, esky = deblaze(sky,esky,orders,channel=0) # deblaze sky      
      f, e = skysub_orders([w,f,e,wsky,sky,esky]) # perform sky subtraction

   bpmap = np.isnan(f).astype(int)            # flag 1 for nan   

   with np.errstate(invalid='ignore'):
      bpmap[f < -3*e] |= flag.neg
      bpmap[e==0] |= flag.nan
      bpmap[np.isnan(e)] |= flag.nan

      return w, f, e, bpmap

def data_goldilocks(self, orders, pfits=True):
   hdulist = self.hdulist
   if 1:  # read order data
      # science
      w,f,e = hdulist['Sci Wavl'].data, hdulist['Sci Flux'].data, hdulist['Sci Error'].data
      w,f,e = w[orders],f[orders],e[orders]      

   bpmap = np.isnan(f).astype(int)            # flag 1 for nan

   with np.errstate(invalid='ignore'):
      bpmap[f < -3*e] |= flag.neg
      bpmap[e==0] |= flag.nan
      bpmap[np.isnan(e)] |= flag.nan
     
      return w, f, e, bpmap

