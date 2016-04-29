#!/usr/bin/env python

import fitsio
import numpy as np
import os
import suchyta_utils.balrog


if __name__=='__main__':
   
    datadir = '/gpfs01/astro/workarea/esuchyta/git-repos/BalrogDirs/2015-Nov/BalrogDB/'
    spto_file = os.path.join(datadir, 'y1a1_spto', 'y1a1_spto-truth.fits')
    p = suchyta_utils.balrog.Y1Processing(dir='/gpfs01/astro/workarea/esuchyta/Y1-files/', systematics='Y1A1NEW_COADD_SPT')
    truth = fitsio.read(spto_file, ext=1, columns=['ra','dec'])
    truth = suchyta_utils.balrog.Y1Dataset(truth, processing=p)

    name = 'y1a1_spto_s82'
    outdir = os.path.join('FITS', 'samples', name)
    if not os.path.exist(outdir):
        os.makedirs(outdir)

    truth.UsualMasking(ra='ra', dec='dec')
    outfile = os.path.join(outdir, 'uniform-pseudo.fits')
    if os.path.exists(outfile):
        os.remove(outfile)
    f = fitsio.FITS(outfile,'rw')
    f.write(truth.data)

    truth.ApplyDepth(cond='>', val=22, ra='ra', dec='dec')
    outfile = os.path.join(outdir, 'uniform-benchmark.fits')
    if os.path.exists(outfile):
        os.remove(outfile)
    f = fitsio.FITS(outfile,'rw')
    f.write(truth.data)
