#!/usr/bin/env python
#-*- coding: utf-8 -*-

## temp imports

from matplotlib import pyplot as p

import numpy as np
from astropy.io import fits
from astropy import units as u
from astropy import wcs
from pyraf import iraf

# Logger
from logutil import *

# Loading necessary IRAF packages and configurations

iraf.digiphot(_doprint=0)
iraf.daophot(_doprint=0)
iraf.apphot(_doprint=0)
iraf.reset(min_lenuserarea='200000')

__all__ = []

class Frame(mylogger):
    def __init__(self, fname, ext, logfn='last.log'):

        ## Initialize DAOPHOT
        self.base = './'
        import daophot
        daophot.set_default(self.base)

        ## Initiate logger utility
        self.sdict={}
        self.log = mylogger(self.sdict, logfn)
        self.fname = fname
        self.pfname = fname.split('/')[-1]
        self.iname = "%s[%i]" % (fname, ext)
        self.ext = ext
        self.read_fits(fname, ext)
        self.read_prim(fname)  # Some infos are only found in the primary HDU
        self.read_wcs()

        # Change with fwhm estimator routine. Bellow for processed DECAM data
        self.fwhm = self.hdup.header['FWHMAV']
        self.fwhmph = self.hdup.header['FWHMAV']*(3600*self.hdup.header['CDELT2'])

        ## Utility names
        self.daofindfn = '%s%d.coo.1' % (self.pfname, ext)
        self.photfn = '%s%d.mag.1' % (self.pfname, ext)
        self.fitpsffn = self.pfname.replace('.fits', '.fitpsf')

    def pix2sky(self, x, y):
        return self.wcs.wcs_pix2world(np.array([x,y]).T, 1)

    def read_fits(self, fname, ext):
        self.log(1, 'READ', 1, 'Reading %s[%i]' % (fname, ext))
        self.hdu = fits.open(fname, memmap=True)[ext]

    def read_prim(self, fname):
        self.log(1, 'READP', 1, 'Reading %s header' % (fname))
        self.hdup = fits.open(fname, memmap=True)[0]

    def read_wcs(self):
        self.log(1, 'WCS READ', 1, 'Reading WCS for %s[%i]' % (self.fname, self.ext))
        self.wcs = wcs.WCS(self.hdu.header)

    def findsky(self, scl, nwin, rerun=False):
        """
        Find sky statistics at random windows.
        Window size set by scl (square) and number of windows by nwin
        """

        a = self.hdu.data
        avoid = 100*scl
        s = scl*0.5
        xsize = a.shape[1]
        ysize = a.shape[1]
        b = np.empty((nwin,4))
        for i in xrange(nwin):
            x = y = -10
            while x < avoid or x > xsize-avoid or y < avoid or y > ysize-avoid:
                x = 1 + int(np.random.rand()*xsize)
                y = 1 + int(np.random.rand()*ysize)
            tmp = a[y-s:y+s, x-s:x+s]
            b[i,0] = np.sum(tmp)/(scl*scl)
            b[i,1] = np.mean(tmp)
            b[i,2] = np.median(tmp)
            b[i,3] = np.std(tmp)

        sigma = np.median(b[:,3])
        sky = np.median(b[:,0])
        self.sigma = sigma
        self.sky = sky
        self.log(1, 'SKY', sky, 'Sky value median: %lf' % sky)
        self.log(1, 'SKYSIGMA', sigma, 'Sky variance median: %lf' % sigma)
        iraf.datapars.setParam('sigma', sigma)
        iraf.fitskypars.setParam('skyvalu', sky)

    def run_daofind(self):
        iraf.daofind.setParam('image',self.iname)
        iraf.datapars.setParam('fwhmpsf',self.fwhm)
        iraf.daofind.setParam('output',self.base+'default')
        iraf.daofind(mode='h',Stdout=1)

    def run_phot(self):
        iraf.phot.setParam('image',self.iname)
        iraf.fitskypars.setParam('skyvalue',self.sky)
        iraf.fitskypars.setParam('annulus',4.*self.fwhm)
        iraf.fitskypars.setParam('dannulus',2.*self.fwhm)
        iraf.phot(mode='h',Stdout=1)

    def run_fitpsf(self):

        # select some guess stars for PSF building
        # Based on median magnitude from apperture phot
        daofind = np.loadtxt(self.daofindfn)
        i = (daofind[:,2] < np.median(daofind[:,2]) + 0.12)&(daofind[:,2] > np.median(daofind[:,2]) - 0.12)
        np.savetxt(self.base+self.pfname+'.guess.coo', daofind[i,0:2],
                   fmt=['%-10.3f','%-10.3f'])

        iraf.fitpsf.setParam('image', self.iname)
        iraf.fitpsf.setParam('output',self.pfname.replace('.fits', '.fitpsf')) # preliminary psf fit
        iraf.fitpsf.setParam('coords', self.base+self.pfname+'.guess.coo')
        iraf.fitpsf(mode='h',Stdout=1)

    def merge(self):
        iraf.txdump(textfile=self.photfn,
                    fields='xcenter,ycenter,mag,msky,merr,id',
                    expr='mag!=INDEF && merr!=INDEF',
                    Stdout=self.base+'phottmp')
        mags = np.loadtxt(self.base+'phottmp')

        iraf.txdump(textfile=self.daofindfn,
                    fields='sharpness,sround,ground,id',
                    expr='sharpness!=INDEF && sround!=INDEF && ground!=INDEF',
                    Stdout=self.base+'daofindtmp')
        daofind = np.loadtxt(self.base+'daofindtmp')

        iraf.txdump(textfile=self.fitpsffn,
                    fields='rsigma,id',
                    expr='rsigma!=INDEF && rsigma < 7.0',
                    Stdout=self.base+'fitpsftmp')
        fitpsf = np.loadtxt(self.base+'fitpsftmp')


        I = reduce(lambda l,r: np.intersect1d(l,r,True), (i[:,-1] for i in
                                                          (daofind, mags,
                                                           fitpsf)))


        oo = np.c_[daofind[np.searchsorted(daofind[:,-1], I)],
                  mags[np.searchsorted(mags[:,-1], I)],
                  fitpsf[np.searchsorted(fitpsf[:,-1], I)]]

        return oo



if __name__=='__main__':
    fname = "/scratch/gc_survey/raw_data/c4d_150715_013102_osi_g_v1.fits"
    tmp = Frame(fname, 2, 'bunda.log')
    t = tmp.pix2sky(np.random.rand(10)*1000, np.random.rand(10)*1000)
    tmp.findsky(10, 1000)
    #tmp.run_daofind()
    #tmp.run_phot()
    #tmp.run_fitpsf()

    f = tmp.merge()
    print f[0]

    p.plot(f[:,0], f[:,1], 'k.', ms=1)
    p.show()


