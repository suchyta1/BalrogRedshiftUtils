#!/usr/bin/env python

import esutil
import fitsio
import numpy as np
import os
import numpy.lib.recfunctions as rec

import suchyta_utils.balrog
import suchyta_utils.mpi
from mpi4py import MPI



if __name__=='__main__':
    #name = 'y1a1_spto'
    #field = 0
    name = 'y1a1_s82'
    field = 1

    ddir = os.path.join('FITS', 'des')
    zfile = os.path.join(ddir, 'ANNz2_Y1A1_LSS_Mar16.fits')
    zzfile = os.path.join(ddir, 'ANNz2_Y1A1_LSS_Mar16_sorted.fits')

    inc = 500000
    bdir = '/gpfs01/astro/workarea/esuchyta/git-repos/BalrogDirs/2015-Nov/BalrogDB/%s'%(name)
    dfile = os.path.join(bdir, '%s-des.fits'%(name))
    outdir = os.path.join(ddir, name)

    if MPI.COMM_WORLD.Get_rank()==0:
        if not os.path.exists(outdir):
            os.makedirs(outdir)

    if not os.path.exists(zzfile):
        if MPI.COMM_WORLD.Get_rank()==0:
            z = fitsio.read(zfile, ext=-1)
            z = np.sort(z, order='COADD_OBJECTS_ID')
            names = {}
            for n in z.dtype.names:
                names[n] = n.lower()
            z = rec.rename_fields(z, names)
            f = fitsio.FITS(zzfile, 'rw')
            f.write(z)


    if MPI.COMM_WORLD.Get_rank()==0:
        size = fitsio.read_header(zzfile, ext=-1)['NAXIS2']
        z = fitsio.read(zzfile, ext=-1, columns=['COADD_OBJECTS_ID'])
        num = size/inc
        print num
        starts = np.arange(num)*inc
        ends = starts + inc
        mod = size%inc
        if mod != 0:
            starts = np.append(starts, ends[-1])
            ends = np.append(ends, ends[-1]+mod)
    else:
        starts = None
        ends = None
   

    starts, ends = suchyta_utils.mpi.Scatter(starts, ends)
    for i in range(len(starts)):
        z = fitsio.read(zzfile, ext=-1, rows=np.arange(starts[i],ends[i]))
        df = fitsio.FITS(dfile, 'r')
        print starts[i], ends[i]
        w = df[1].where('coadd_objects_id >= %i && coadd_objects_id <= %i'%(z['coadd_objects_id'][0],z['coadd_objects_id'][-1]))
        if len(w)==0:
            continue
        d = df[1][w]
        dz = rec.join_by('coadd_objects_id', d, z, usemask=False)
        dz = rec.append_fields(dz, 'field', [field]*len(dz))

        file = os.path.join(outdir, '%i-%i.fits'%(starts[i],ends[i]))
        if os.path.exists(file):
            os.remove(file)
        f = fitsio.FITS(file, 'rw')
        f.write(dz)

    MPI.COMM_WORLD.Barrier()
    if MPI.COMM_WORLD.Get_rank()==0:
        outfile = os.path.join(outdir, '%s.fits'%(name))
        if os.path.exists(outfile):
            os.remove(outfile)
        files = os.listdir(outdir)
        of = fitsio.FITS(outfile, 'rw')
        for i in range(len(files)):
            f = os.path.join(outdir, files[i])
            data = fitsio.read(f, ext=1)
            if i==0:
                of.write(data)
            else:
                of[1].append(data)
            os.remove(f)
