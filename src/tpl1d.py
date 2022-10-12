#! /usr/bin/env python
from __future__ import division

# merge orders of a template
# only for blaze-normalised spectra
# sampling currently hardcoded

from astropy.io import fits
from gplot import *
import numpy as np
import cspline

def tpl1d(filename, out=None, show=False):
    hdu = fits.open(filename)
    hdr = hdu[0].header

    w, f = hdu['wave'].data.astype(float), hdu['spec'].data.astype(float)
    ok = w != 1.    
    ok[w==0] = 0

    # parabolic window to downweight of edges to reduce egde effects
    nx = ok[0].size
    x = np.arange(nx)/nx - 0.5
    y = 1 - x**2/0.25
    ymap = ok * y
    K = int((w[ok].max()-w[ok].min())*3e5/0.2)
    spl = cspline.ucbspl_fit(w[ok], f[ok], w=ymap[ok], K=K, lam=0.1)

    if show:
        gplot( spl.osamp(0.1), 'w l lc 2 t "tpl",',
                                w[ok][::20], f[ok][::20], 'lc 1 pt 1 ps 0.3',)
        #gplot(np.exp(w[ok]), f[ok], 'w l t "2d", ', np.exp(spl.osamp(1)[0]), spl.osamp(1)[1], 'w l t "1d"')
        from pause import pause; pause()

    # get gaps
    w2 = w[~np.all(w==0,axis=1)]
    gaps = np.array([w2[:-1,-1],w2[1:,0]]).T
    gaps = gaps[gaps[:,1]-gaps[:,0]>0]
    
    tw,tf = spl.osamp(1)
    # set tpl in gaps to nan
    for g in gaps:
        m = np.logical_and(tw>g[0],tw<g[1])
        tf[m] = np.nan

    tpl = np.rec.array((tw,tf), names='lnwave,flux')

    if show:
        gplot( spl.osamp(0.1), 'w l lc 2 t "tpl",',
                                tw[::10], tf[::10], 'lc 1 pt 1 ps 0.3',)
        from pause import pause; pause()

    if not out:
        out = filename.replace('.fits','') + '_1d.fits'

    fits.writeto(out, tpl, hdr[20:], overwrite=True)

    print("-> "+out)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('filename')
    parser.add_argument('out', nargs='?')
    parser.add_argument('-show', action='store_true')
    args = parser.parse_args()
    tpl1d(**vars(args))

# ./tpl1d.py /home/astro115/carmenes/data/svn/serval/CARM_VIS/J12479+097/template.fits J12479+097_tpl.fits
