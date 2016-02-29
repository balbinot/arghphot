#!/usr/bin/env python
#-*- coding: utf-8 -*-

import logging
from colorlog import ColoredFormatter

LOG_LEVEL = logging.DEBUG
LOGFORMAT = "[%(bold)s%(log_color)s%(levelname)-8s%(reset)s %(asctime)s] %(white)s%(message)s%(reset)s"
logging.root.setLevel(LOG_LEVEL)
formatter = ColoredFormatter(LOGFORMAT)
stream = logging.StreamHandler()
stream.setLevel(LOG_LEVEL)
stream.setFormatter(formatter)
fh = logging.FileHandler('spam.log')
fh.setLevel(logging.DEBUG)
fh.setFormatter(formatter)

import numpy as np
from astropy.io import fits
from astropy import units as u
from astropy import wcs
from pyraf import iraf

__all__ = []

# Loading necessary IRAF packages
iraf.digiphot(_doprint=0)
iraf.daophot(_doprint=0)
iraf.apphot(_doprint=0)
iraf.reset(min_lenuserarea='200000')


class Frame():
    def __init__(self, fname, ext):
        self.log = logging.getLogger(fname.replace('.fits','.log'))
        fh = logging.FileHandler('spam.log')
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(formatter)
        self.log.addHandler(stream)
        self.log.addHandler(fh)
        self.fname = fname
        self.ext = ext
        self.read_fits(fname, ext)
        self.read_wcs()

    def pix2sky(self, x, y):
        return self.wcs.wcs_pix2sky(np.array([x,y]).T, 1)

    def read_fits(self, fname, ext):
        self.log.info('Reading %s[%i]' % (fname, ext))
        self.hdu = fits.open(fname, memmap=True)[ext]

    def read_wcs(self):
        self.log.info('Reading WCS for %s[%i]' % (self.fname, self.ext))
        self.wcs = wcs.WCS(self.hdu.header)

    def get_params(self):
        pass


if __name__=='__main__':
    fname = "/scratch/gc_survey/raw_data/c4d_150715_013102_osi_g_v1.fits"
    tmp = Frame(fname, 2)

