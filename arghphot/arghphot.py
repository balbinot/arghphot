#!/usr/bin/env python
#-*- coding: utf-8 -*-

## temp imports

from matplotlib import pyplot as p

import numpy as np
from astropy.io import fits, ascii
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy import wcs
import aplpy
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
        self.pstfile = '%s%d.pst.1' % (self.pfname, ext)
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

    def trim_phot(self):
        a = ascii.read(self.photfn)
        std = a['STDEV']
        sum = a['SUM']/a['AREA']
        sky = a['MSKY']
        sn = np.abs(sum-sky)/std
        i = sn > 3
        tmp = a[i]
        tmp.write(self.photfn+'trim', format='ascii')

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
        ## Use trimmed photometry to avoid variable sky background spurious
        ## detections
        # x,y, msky, stdev, sum, area, mag, merr, id
        mags = np.loadtxt(self.photfn+'trim', usecols=(6,7,14,15,26,27,29,30,3), skiprows=1)

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

        ## I have no idea how this works
        I = reduce(lambda l,r: np.intersect1d(l,r,True), (i[:,-1] for i in
                                                          (mags, daofind,
                                                           fitpsf)))
        oo = np.c_[mags[np.searchsorted(mags[:,-1], I)],
                   daofind[np.searchsorted(daofind[:,-1], I)],
                   fitpsf[np.searchsorted(fitpsf[:,-1], I)]]

        return oo

    def select_psf(self):
        f = self.merge()
        w = self.pix2sky(f[:,0], f[:,1])
        coo = SkyCoord(w[:,0]*u.deg, w[:,1]*u.deg)
        nid, nsep, _ = coo.match_to_catalog_sky(coo, nthneighbor=2)
        x = f[:,0]
        y = f[:,1]
        id = f[:,-1]
        sky = f[:,2]
        skystd = f[:,3]
        mag = f[:,6]
        merr = f[:,7]
        sharp = f[:,9]
        fwhm = f[:,13]

        ## Set of constrains for PSF stars
        i = (mag < np.mean(mag) - 1.2*np.std(mag))
        i *= (merr < 0.02)
        i *= (np.abs(sharp-np.median(sharp)) < 0.04)
        i *= (fwhm < 1.05*tmp.fwhm/2.355)
        i *= (fwhm > 0.8*tmp.fwhm/2.355)
        i *= (nsep.arcsec > 15*tmp.fwhmph)
        i *= (x > 60*tmp.fwhm)*(y > 60*tmp.fwhm)
        i *= (x < tmp.hdu.data.shape[1] - 60*tmp.fwhm)
        i *= (y < tmp.hdu.data.shape[0] - 60*tmp.fwhm)
        i *= (np.abs(skystd - np.mean(skystd)) < np.std(skystd))
        i *= (sky - np.mean(sky) < np.std(sky))

        self._parse_pst(id[i], x[i], y[i], mag[i], sky[i])
        return (id[i], x[i], y[i], mag[i], sky[i])

    def tvmark(self, ra, dec):
        gc = aplpy.FITSFigure(self.hdu)
        gc.show_grayscale(stretch='arcsinh')
        gc.set_tick_labels_font(size='small')
        gc.show_markers(ra,dec,layer='scatter',edgecolor='red',
                        facecolor='none',marker='o',s=10,alpha=0.5)

    def _parse_pst(self, id, x, y, mag, msky):
        pstfile = open(self.pstfile, 'w')
        pstfile.write("#N ID    XCENTER   YCENTER   MAG         MSKY                                  \\\n")
        pstfile.write("#U ##    pixels    pixels    magnitudes  counts                                \\\n")
        pstfile.write("#F %-9d  %-10.3f   %-10.3f   %-12.3f     %-15.7g                                 \n")
        pstfile.write("#\n")
        np.savetxt(pstfile, np.array([id, x,y, mag, msky]).T,
                   fmt=['%-9d','%-10.3f','%-10.3f','%-12.3f','%-15.7g'])
        pstfile.close()

    def grid_psf(self):
        from mpl_toolkits.axes_grid1 import ImageGrid
        id, x, y = np.loadtxt(self.pstfile, usecols=(0,1,2), unpack=True)


        grid1 = ImageGrid(fig, 111,
                nrows_ncols = (2, 4),
                axes_pad = 0.07,
                share_all=True,
                label_mode = "L",
                cbar_location = "right",
                cbar_mode="single",
                cbar_size="7%",
                cbar_pad="7%",
                aspect = True
                )


        xbox = int(f[i,0] - rad)
        Xbox = int(f[i,0] + rad)
        ybox = int(f[i,1] - rad)
        Ybox = int(f[i,1] + rad)
        block = img[xbox:Xbox,ybox:Ybox]





if __name__=='__main__':
    fname = "/scratch/gc_survey/raw_data/c4d_150715_013102_osi_g_v1.fits"
    tmp = Frame(fname, 2, 'bunda.log')
    tmp.findsky(10, 1000)
    #tmp.run_daofind()
    #tmp.run_phot()
    tmp.trim_phot()
    #tmp.run_fitpsf()

    f = tmp.select_psf()

    t = tmp.pix2sky(f[1], f[2])
    tmp.tvmark(t[:,0], t[:,1])

    p.show()


