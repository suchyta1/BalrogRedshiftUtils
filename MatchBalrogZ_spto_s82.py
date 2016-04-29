#!/usr/bin/env python

import esutil
import fitsio
import numpy as np
import os
import sys
import copy
import numpy.lib.recfunctions as rec


def DoJoin(balrog, row, size, odir, zz, names, end=None, cols=False, field=None):
    if not os.path.exists(odir):
        os.makedirs(odir)

    if end is None:
        end = row + len(zz)

    if end > size:
        end = size
    ee = end - row

    b = balrog[-1].read(rows=np.arange(row,end))
    d = []
    for name in names:
        d.append(zz[name][:ee])

    n = list(names)
    n.append('field')
    d.append(np.array([field]*len(b)))

    c = rec.append_fields(b, n, d)
    if 'table' in c.dtype.names:
        c = rec.drop_fields(c, 'table')

    ofile = os.path.join(odir, '%i-%i.fits'%(row,end))
    esutil.io.write(ofile, c, clobber=True)
    
    if cols:
        return end, c.dtype.names
    else:
        return end


def GetAll(files, outdir):
    datadir = '/gpfs01/astro/workarea/esuchyta/git-repos/BalrogDirs/2015-Nov/BalrogDB/'
    spto_file = os.path.join(datadir, 'y1a1_spto', 'y1a1_spto-sim.fits')
    s82_file = os.path.join(datadir, 'y1a1_s82', 'y1a1_s82-sim.fits')

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    names = ('balrog_id_z', 'intrainflag', 'annz_best', 'annz_best_err')
    formats = ('<i8', 'f4', 'f4', 'f4')

    stack = None
    row = 0
    field = 0
    balrog = fitsio.FITS(spto_file, 'r')
    size = balrog[-1].read_header()['NAXIS2']
    odir = os.path.join(outdir,'y1a1_spto')

    for i in range(len(files)):
        zz = np.loadtxt(files[i], delimiter=',', usecols=(1,2,3,4), dtype={'names':names, 'formats':formats})

        if len(zz)==0:
            end = row + 500000
        else:
            if i==0:
                end, cols1 = DoJoin(balrog, row, size, odir, zz, names, cols=True, field=field)
            else:
                end = DoJoin(balrog, row, size, odir, zz, names, cols=False, field=field)

            if (i!=len(files)-1) and (end==size):
                done = end - row
                ee = (len(zz)-done)
                zz = zz[done:]
                row = 0
                field = 1
                balrog = fitsio.FITS(s82_file, 'r')
                size = balrog[-1].read_header()['NAXIS2']
                odir = os.path.join(outdir,'y1a1_s82')
                end, cols2 = DoJoin(balrog, row, size, odir, zz, names, end=ee, cols=True, field=field)

        row = end
        print files[i]
    return cols1, cols2


def Stack(name, keep, outdir):
    outdir = os.path.join(outdir, name)
    outfile = os.path.join(outdir, '%s.fits'%(name))
    if os.path.exists(outfile):
        os.remove(outfile)

    files = np.sort(os.listdir(outdir))
    of = fitsio.FITS(outfile, 'rw')
    for i in range(len(files)):
        f = os.path.join(outdir, files[i])
        data = fitsio.read(f, ext=1, columns=keep)
        if i==0:
            of.write(data)
        else:
            of[1].append(data)
        os.remove(f)
        print files[i]
            


if __name__ == "__main__":

    photoz_dir = 'TXT/balrog_ANNz2_cat'
    files = np.sort(os.listdir(photoz_dir))
    usefiles = []
    for file in files:
        if file.startswith('ANNZ_'):
            f = os.path.join(photoz_dir, file)
            usefiles.append(f)

    outdir = 'FITS/balrog'
    cols1, cols2 = GetAll(usefiles, outdir)
    keep = np.intersect1d(cols1, cols2)
    
    Stack('y1a1_spto', keep, outdir)
    Stack('y1a1_s82', keep, outdir)
    
