#!/usr/bin/env python

import sys
import fitsio
import numpy as np
import os

import suchyta_utils.jk
import kmeans_radec



if __name__ == "__main__":
    file = 'FITS/samples/y1a1_spto_s82/des-benchmark.fits'

    njack = 64
    ra = 'alphawin_j2000_i'
    dec = 'deltawin_j2000_i'

    sample = fitsio.read(file, ext=1)
    f = 'y1a1_spto_s82_des-benchmark'
    jkfile = 'jk-regions/%s-%i.txt'%(f,njack)

    if not os.path.exists(jkfile):
        suchyta_utils.jk.GenerateJKRegions(sample[ra], sample[dec], njack, jkfile, maxiter=1000)

