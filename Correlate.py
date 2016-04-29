#!/usr/bin/env python

import sys
import esutil
import fitsio
import numpy as np
import os

import suchyta_utils.balrog
import suchyta_utils.mpi

from mpi4py import MPI
import treecorr


def CorrelationFunction(source, source_random, corrconfig=None, source_ra='ra', source_dec='dec', source_random_ra=None, source_random_dec=None):

    if corrconfig is None:
        corrconfig = {'sep_units': 'degrees',
                      'min_sep': 0.1,
                      'max_sep': 6,
                      'nbins': 25,
                      'bin_slop': 0.25,
                      'num_threads': 1}

    if source_random_ra is None:
        source_random_ra = source_ra
    if source_random_dec is None:
        source_random_dec = source_dec

    SourceCat = treecorr.Catalog(ra=source[source_ra], dec=source[source_dec], ra_units='degrees', dec_units='degrees')
    SourceRand = treecorr.Catalog(ra=source_random[source_random_ra], dec=source_random[source_random_dec], ra_units='degrees', dec_units='degrees')
    dd = treecorr.NNCorrelation(**corrconfig)
    dr = treecorr.NNCorrelation(**corrconfig)
    rr = treecorr.NNCorrelation(**corrconfig)

    dd.process(SourceCat)
    dr.process(SourceCat, SourceRand)
    rr.process(SourceRand)

    xi, varxi = dd.calculateXi(rr, dr=dr)
    r = np.exp(dd.logr)
    
    return [xi, r]



if __name__ == "__main__":

    kind = 'pseudo'
    rand = 'balrog'


    njack = 64
    f = 'y1a1_spto_s82_des-benchmark'
    jkfile = 'jk-regions/%s-%i.txt'%(f,njack)

    z1 = np.arange(0.6, 1.01, 0.1)[-1:]
    z2 = np.arange(0.7, 1.11, 0.1)[-1:]

    regions = [jkfile]*2
    pos = ['alphawin_j2000_i', 'deltawin_j2000_i']
    upos = ['ra', 'dec']

    des = fitsio.read('FITS/samples/y1a1_spto_s82/des-%s.fits'%(kind), ext=1)
    sim = fitsio.read('FITS/samples/y1a1_spto_s82/sim-%s.fits'%(kind), ext=1)
    if rand=='balrog':
        jkby = [pos,pos]
        nojkkwargs = {'source_ra':pos[0], 'source_dec':pos[1], 'source_random_ra':pos[0], 'source_random_dec':pos[1]}
    elif rand=='uniform':
        jkby = [pos,upos]
        nojkkwargs = {'source_ra':pos[0], 'source_dec':pos[1], 'source_random_ra':upos[0], 'source_random_dec':upos[1]}


    for i in range(len(z1)):
        d = suchyta_utils.balrog.ApplyCut(des, key='annz2_z', cond='>', val=z1[i])
        d = suchyta_utils.balrog.ApplyCut(d, key='annz2_z', cond='<', val=z2[i])

        if rand=='balrog':
            s = suchyta_utils.balrog.ApplyCut(sim, key='annz_best', cond='>', val=z1[i])
            s = suchyta_utils.balrog.ApplyCut(s, key='annz_best', cond='<', val=z2[i])
        elif rand=='uniform':
            if MPI.COMM_WORLD.Get_rank()==0:
                print len(d)
                ukeep = np.random.choice(len(sim), size=10*len(d), replace=False) 
            else:
                ukeep = None
            ukeep = suchyta_utils.mpi.Broadcast(ukeep)
            s = sim[ukeep]

        jkargs = [d,s]
        Corr = suchyta_utils.jk.SphericalJK(target=CorrelationFunction, jkargs=jkargs, jkargsby=jkby, regions=regions, nojkkwargs=nojkkwargs)
        Corr.DoJK(mpi=True)

        if MPI.COMM_WORLD.Get_rank()==0:
            Results = Corr.GetResults()

            zstr = '%.2f<z<%.2f'%(z1[i],z2[i])
            outdir = os.path.join('CorrResults-%s_%s'%(kind, rand), zstr)
            if not os.path.exists(outdir):
                os.makedirs(outdir)

            tfile = os.path.join(outdir, 'separation.dat')
            jfile = os.path.join(outdir, 'jk.dat')
            ffile = os.path.join(outdir, 'full.dat')
            
            np.savetxt(tfile, Results['full'][1])
            np.savetxt(ffile, Results['full'][0])
            np.savetxt(jfile, Results['jk'][:, 0])

            print 'done with %s'%(zstr)
