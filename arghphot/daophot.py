#!/usr/bin/env python
#-*- coding: utf-8 -*-

from pyraf import iraf

def set_default(base):
    # Loading necessary IRAF packages
    iraf.digiphot(_doprint=0)
    iraf.daophot(_doprint=0)
    iraf.apphot(_doprint=0)
    iraf.reset(min_lenuserarea='200000')

    iraf.datapars.setParam('datamin','INDEF')
    iraf.datapars.setParam('datamax','50000')
    #iraf.datapars.setParam('readnoise','RDNOISE')
    #iraf.datapars.setParam('epadu','GAIN')
    iraf.datapars.setParam('exposure','EXPTIME')
    iraf.datapars.setParam('airmass', 'AIRMASS')
    iraf.datapars.setParam('filter', 'FILTER')

    iraf.findpars.setParam('threshold',4.0)
    iraf.findpars.setParam('sharphi',  1)
    iraf.findpars.setParam('roundhi', 1.3)
    iraf.findpars.setParam('roundlo', -1.3)

    iraf.daofind.setParam('verify','no')
    iraf.daofind.setParam('interactive','no')
    iraf.daofind.setParam('verify','no')
    iraf.daofind.setParam('interactive','no')

    iraf.photpars.setParam('zmag',25.0)
    iraf.photpars.setParam('weighti','constant')
    iraf.photpars.setParam('apertur',3.0)

    iraf.phot.setParam('output',base+'default')
    iraf.phot.setParam('coords',base+'default')
    iraf.phot.setParam('verify','no')
    iraf.phot.setParam('interactive','no')

    iraf.fitpsf.setParam('box',10.0)
    iraf.fitpsf.setParam('verify','no')
    iraf.fitpsf.setParam('interactive','no')

    iraf.centerpars.setParam('calgori','none')
    iraf.fitskypars.setParam('salgorithm','mode')

    iraf.daopars.setParam('functio','moffat15')
    iraf.daopars.setParam('varorde','0')
    iraf.daopars.setParam('nclean','0')
    iraf.daopars.setParam('saturat','no')
    iraf.daopars.setParam('fitsky','yes')
    iraf.daopars.setParam('recenter','yes')
    iraf.daopars.setParam('groupsk','yes')
    iraf.daopars.setParam('maxnsta','20000')

    iraf.psf.setParam('photfile',base+'default')
    iraf.psf.setParam('pstfile',base+'default')
    iraf.psf.setParam('psfimage',base+'default')
    iraf.psf.setParam('opstfile',base+'default')
    iraf.psf.setParam('groupfil',base+'default')
    iraf.psf.setParam('interac','no')
    iraf.psf.setParam('matchby','yes')
    iraf.psf.setParam('verify','no')
    iraf.psf.setParam('showplo','no')

    iraf.allstar.setParam('photfile',base+'default')
    iraf.allstar.setParam('psfimage',base+'default')
    iraf.allstar.setParam('allstarf',base+'default')
    iraf.allstar.setParam('rejfile',base+'default')
    iraf.allstar.setParam('subimage',base+'default')
    iraf.allstar.setParam('verify','no')
