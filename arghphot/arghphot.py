#!/usr/bin/env python
#-*- coding: utf-8 -*-

## temp imports

from matplotlib import pyplot as p
from matplotlib import cm

import tempfile

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

__all__ = ['Frame']

class Frame(mylogger):
    def __init__(self, fname, ext, mask, logfn='last.log'):

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

        # Read image and maks
        self.hdu = self.read_fits(fname, ext)
        self.mask = self.read_fits(mask, ext)
        self.read_prim(fname)  # Some infos are only found in the primary HDU
        self.read_wcs()

        # Change with fwhm estimator routine. Bellow for processed DECAM data
        self.fwhm = self.hdup.header['FWHMAV']
        self.fwhmph = self.hdup.header['FWHMAV']*(3600*self.hdup.header['CDELT2'])

        ## Utility names
        self.daofindfn = '%s%d.coo.1' % (self.pfname, ext)
        self.photfn = '%s%d.mag.1' % (self.pfname, ext)
        self.pstfile = '%s%d.pst.1' % (self.pfname, ext)
        self.fitpsffn = '%s%d.fitpsf.1' % (self.pfname, ext)
        self.guess = '%s%d.guess.1' % (self.pfname, ext)
        self.psfgridname = '%s%d.psfgrid.png' % (self.pfname, ext)
        self.psf = '%s%d.psf.1.fits' % (self.pfname, ext)
        self.psfimg = '%s%d.psf.1.img.fits' % (self.pfname, ext)
        self.psfselectplot = '%s%d.psfselect.png' % (self.pfname, ext)

    def pix2sky(self, x, y):
        return self.wcs.wcs_pix2world(np.array([x,y]).T, 1)

    def read_fits(self, fname, ext):
        self.log(1, 'READ', 1, 'Reading %s[%i]' % (fname, ext))
        return fits.open(fname, memmap=True)[ext]

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
        m = self.mask.data
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
            tmpm = m[y-s:y+s, x-s:x+s]
            if np.any(tmpm < 1):
                b[i,0] = np.nan
                b[i,1] = np.nan
                b[i,2] = np.nan
                b[i,3] = np.nan
            else:
                b[i,0] = np.sum(tmp)/(scl*scl)
                b[i,1] = np.mean(tmp)
                b[i,2] = np.median(tmp)
                b[i,3] = np.std(tmp)

        sigma = np.nanmedian(b[:,3])
        sky = np.nanmedian(b[:,0])
        self.sigma = sigma
        self.sky = sky
        self.log(1, 'SKY', sky, 'Sky value median: %lf' % sky)
        self.log(1, 'SKYSIGMA', sigma, 'Sky variance median: %lf' % sigma)
        iraf.datapars.setParam('sigma', sigma)
        iraf.fitskypars.setParam('skyvalu', sky)

    def run_daofind(self, coofn=False):
        iraf.daofind.setParam('image',self.iname)
        iraf.datapars.setParam('fwhmpsf',self.fwhm)
        if coofn==False:
            coofn = self.base+'default'
        iraf.daofind.setParam('output', coofn)
        iraf.daofind(mode='h',Stdout=1)

        return coofn

    def run_phot(self, coofn=False, outfn=False):
        if coofn==False:
            iraf.phot.setParam('coords',self.base+'default')
        else:
            iraf.phot.setParam('coords', coofn)

        if outfn==False:
            iraf.phot.setParam('output',self.base+'default')
        else:
            iraf.phot.setParam('output', outfn)

        iraf.phot.setParam('image',self.iname)
        iraf.fitskypars.setParam('skyvalue',self.sky)
        iraf.fitskypars.setParam('annulus',4.*self.fwhm)
        iraf.fitskypars.setParam('dannulus',2.*self.fwhm)
        iraf.photpars.setParam('zmag', self.hdup.header['MAGZERO']) # Use DECAM estimate of zeropoint
        iraf.phot(mode='h',Stdout=1)

    def trim_phot(self, infn=False, outfn=False):
        if infn==False:
            a = ascii.read(self.photfn)
        else:
            a = ascii.read(infn)
        std = a['STDEV']
        sum = a['SUM']/a['AREA']
        sky = a['MSKY']
        sn = np.abs(sum-sky)/std
        i = sn > 3
        tmp = a[i]
        if outfn==False:
           outfn = self.photfn+'trim'
        tmp.write(outfn, format='ascii')

        return outfn, np.where(i)[0]

    def run_fitpsf(self):

        # select some guess stars for PSF building
        # Based on median magnitude from apperture phot
        daofind = np.loadtxt(self.daofindfn, usecols=(0,1,2))
        i = (daofind[:,2] < np.median(daofind[:,2]) + 0.12)&(daofind[:,2] > np.median(daofind[:,2]) - 0.12)
        np.savetxt(self.base+self.guess, daofind[i,0:2],
                   fmt=['%-10.3f','%-10.3f'])

        iraf.fitpsf.setParam('image', self.iname)
        iraf.fitpsf.setParam('output', self.base+self.fitpsffn) # preliminary psf fit
        iraf.fitpsf.setParam('coords', self.base+self.guess)
        iraf.fitpsf(mode='h',Stdout=1)

    def merge(self):
        ## Use trimmed photometry to avoid variable sky background spurious
        ## detections
        # x,y, msky, stdev, sum, area, mag, merr, id
#        mags = np.loadtxt(self.photfn+'trim', usecols=(6,7,14,15,26,27,29,30,3), skiprows=1)
        mags = np.genfromtxt(self.photfn+'trim', usecols=(6,7,14,15,26,27,29,30,3), skip_header=1, dtype='|S5')
        j = True
        for i in np.arange(mags.shape[1]):
            j *= (mags[:,i] != '--')
        mags = mags[j]
        mags = mags.astype(np.float64)

        tf = tempfile.NamedTemporaryFile(dir=self.base)
        iraf.txdump(textfile=self.daofindfn,
                    fields='sharpness,sround,ground,id',
                    expr='sharpness!=INDEF && sround!=INDEF && ground!=INDEF',
                    Stdout=tf.name+'coo')
        daofind = np.loadtxt(tf.name+'coo')

        tf = tempfile.NamedTemporaryFile(dir=self.base)
        iraf.txdump(textfile=self.fitpsffn,
                    fields='rsigma,id',
                    expr='rsigma!=INDEF && rsigma < 7.0',
                    Stdout=tf.name+'psf')
        fitpsf = np.loadtxt(tf.name+'psf')

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

        self.maglim = np.mean(mag) -1.2
        self.merrlim = 0.08
        self.sharplim = 0.04
        self.fwhmlimup = 1.15*self.fwhm/2.355
        self.fwhmlimlow = 0.8*self.fwhm/2.355

        p.figure(figsize=(12,12))
        p.subplot(331)
        p.xlabel('mag')
        p.hist(mag, range=[12,32], bins=30, color='k', alpha=0.6)
        p.axvline(x=self.maglim, c='k')
        p.subplot(332)
        p.xlabel('merr')
        p.hist(merr, range=[0,0.5], bins=30, color='k', alpha=0.6)
        p.axvline(x=self.merrlim, c='k')
        p.subplot(333)
        p.xlabel('sharpness')
        p.hist(sharp, range=[-0.5,0.5], bins=30, color='k', alpha=0.6)
        p.axvline(x=np.median(sharp), c='k')
        p.axvline(x=np.median(sharp)+self.sharplim, ls='--', c='k')
        p.axvline(x=np.median(sharp)-self.sharplim, ls='--', c='k')
        p.subplot(334)
        p.xlabel('fwhm [px]')
        p.hist(fwhm, range=[1,10], bins=30, color='k', alpha=0.6)
        p.axvline(x=self.fwhmlimup, ls='--', c='k')
        p.axvline(x=self.fwhmlimlow, ls='--', c='k')
        p.subplot(335)
        p.xlabel('separation [arcsec]')
        p.hist(nsep.arcsec, bins=30, color='k', alpha=0.6)
        p.axvline(x=12*self.fwhmph, ls='--', c='k')
        p.subplot(336)
        p.xlabel('sky std [counts]')
        p.hist(skystd, bins=30, color='k', alpha=0.6)
        p.axvline(x=np.median(skystd) - np.std(skystd), ls='--', c='k')
        p.axvline(x=np.median(skystd) + np.std(skystd), ls='--', c='k')
        p.axvline(x=np.median(skystd), ls='-', c='k')
        p.subplot(337)
        p.xlabel('sky [counts]')
        p.hist(sky, range=[40,90], bins=30, color='k', alpha=0.6)
        p.axvline(x=np.mean(sky) - np.std(sky), ls='--', c='k')
        p.axvline(x=np.mean(sky) + np.std(sky), ls='--', c='k')
        p.axvline(x=np.mean(sky), ls='-', c='k')
        p.savefig(self.psfselectplot)

        ## Set of constrains for PSF stars
        i = (mag < self.maglim)  # works if zeropoint guess was ok
        i *= (merr < self.merrlim)
        i *= (np.abs(sharp-np.median(sharp)) < self.sharplim)
        i *= (fwhm < self.fwhmlimup)
        i *= (fwhm > self.fwhmlimlow)
        i *= (nsep.arcsec > 17*self.fwhmph)
        i *= (x > 60*self.fwhm)*(y > 60*self.fwhm)
        i *= (x < self.hdu.data.shape[1] - 60*self.fwhm)
        i *= (y < self.hdu.data.shape[0] - 60*self.fwhm)
        i *= (np.abs(skystd - np.median(skystd)) < np.std(skystd))
        i *= (np.abs(sky - np.mean(sky)) < np.std(sky))

        if len(id[i]) <= 2:
            self.log(3, 'NPSF', len(id[i]), 'Number of PSF stars less than 2')
        else:
            self.log(1, 'NPSF', len(id[i]), 'Number of PSF stars is %i' % len(id[i]))

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
        side = int(np.ceil(np.sqrt(len(x))))
        rad = int(6*self.fwhm)
        fig = p.figure(figsize=(12,12))
        grid =  ImageGrid(fig, 111,
                          nrows_ncols = (side, side),
                          axes_pad = 0.0,
                          share_all=True,
                          label_mode = "L",
                          cbar_location = "right",
                          cbar_mode=None,
#                          cbar_size="5%",
#                          cbar_pad="5%",
                          aspect = True
                )

        for i in np.arange(len(x)):
            xbox = int(x[i] - rad)
            Xbox = int(x[i] + rad)
            ybox = int(y[i] - rad)
            Ybox = int(y[i] + rad)
            block = self.hdu.data[ybox:Ybox,xbox:Xbox]
            grid[i].imshow(block.T, origin='lower', cmap=cm.gray,
                           vmin=self.sky-5*self.sigma, vmax=300,
                           interpolation='nearest')

        p.savefig(self.psfgridname)

    def run_psf(self):
        fwhm = self.fwhm
        iraf.daopars.setParam('matchra',fwhm)
        iraf.daopars.setParam('psfrad',4*fwhm+1)
        iraf.daopars.setParam('fitrad',fwhm)
        iraf.daopars.setParam('sannulu',2*fwhm)
        iraf.daopars.setParam('wsannul',4*fwhm)
        iraf.psf.setParam('image',self.iname)
        iraf.psf(mode='h')
        iraf.seepsf(psfimage=self.base+self.psf, image=self.base+self.psfimg,
                    magnitu='18.0')

    def run_allstar(self):
        fwhm = self.fwhm
        iraf.daopars.setParam('matchra',fwhm)
        iraf.daopars.setParam('psfrad',4*fwhm+1)
        iraf.daopars.setParam('fitrad',fwhm)
        iraf.daopars.setParam('sannulu',2*fwhm)
        iraf.daopars.setParam('wsannul',4*fwhm)
        iraf.allstar.setParam('image',self.iname)
        iraf.allstar(mode='h',verbose='no')


if __name__=='__main__':
    fname = "/scratch/gc_survey/raw_data/c4d_150715_013102_osi_g_v1.fits"
    tmp = Frame(fname, 2, 'bunda.log')
    tmp.findsky(10, 1000)
    #tmp.run_daofind()
    #tmp.run_phot()
    #tmp.trim_phot()
    #tmp.run_fitpsf()
    #f = tmp.select_psf()
    #tmp.grid_psf()
    #tmp.run_psf()
    tmp.run_allstar()

    #t = tmp.pix2sky(f[1], f[2])
    #tmp.tvmark(t[:,0], t[:,1])

    #p.show()


