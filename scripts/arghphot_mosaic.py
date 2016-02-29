#!/usr/bin/python

#===============================================================================#
# This is ARGHphot, a software to perform source photometry on each chip of     #
# MOSAIC2 images.                                                               #
#                                                                               #
# Eduardo Balbinot - June 2010                                                  #
#===============================================================================#

import sys
import os
from numpy import *
import numpy as np
from pyraf import iraf
import pyfits
import pywcs

# Loading necessary IRAF packages
iraf.digiphot(_doprint=0)
iraf.daophot(_doprint=0)
iraf.apphot(_doprint=0)
iraf.reset(min_lenuserarea='200000')

class bcolor:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m\033[1m'
    OKGREEN = '\033[92m\033[1m'
    WARNING = '\033[91m\033[1m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'

    def disable(self):
        self.HEADER = ''
        self.OKBLUE = ''
        self.OKGREEN = ''
        self.WARNING = ''
        self.FAIL = ''
        self.ENDC = ''

def addflags(filename,ext,img):
    f = open(filename,'r')
    g = open(filename+'.flag','w+')
    for line in f:
        a = ''.join(line[:])
        a = a.replace('\n','')
        a = a+'  '+str(ext)+'  '+str(img)+'\n'
        g.write(a)
    f.close()
    g.close()

def coconv(f,img,ext):
    hdu = pyfits.open(img)
    wcs = pywcs.WCS(hdu[0].header)
    g = np.empty((f.shape[0],f.shape[1]+2))
    g[:,0:3] = f[:,0:3]
    g[:,4:] = f[:,2:]
    for i in range(f.shape[0]):
        g[i,2],g[i,3] = wcs.wcs_pix2sky(f[i,0],f[i,1],1)
    import ipdb; ipdb.set_trace()  # XXX BREAKPOINT
    return(g)

def cutbad(f,img,low,clip,high,rad):

    img = np.transpose(pyfits.getdata(img))
    rad = rad
    out = np.empty((0,f.shape[1]))

    for i in range(f.shape[0]):
        xbox = int(f[i,0] - rad)
        Xbox = int(f[i,0] + rad)
        ybox = int(f[i,1] - rad)
        Ybox = int(f[i,1] + rad)
        block = img[xbox:Xbox,ybox:Ybox]
        side1 = np.mean(block[0:2,0])
        side2 = np.mean(block[int(rad)-1:int(rad)+1,0])
        side3 = np.mean(block[0,0:2])
        side4 = np.mean(block[0,int(rad)-1:int(rad)+1])
        sides = np.array([side1,side2,side3,side4])
        if np.any((block>=high)):
            print 'Star (',f[i,0],';',f[i,1],') eliminated. ' + bcolor.WARNING + 'Error: bad pixel neighbour' + bcolor.ENDC
        elif np.any((block<=low)):
            print 'Star (',f[i,0],';',f[i,1],') eliminated. ' + bcolor.WARNING + 'Error: low neighbour' + bcolor.ENDC
        elif np.any((sides>=clip)):
            print 'Star (',f[i,0],';',f[i,1],') eliminated. ' + bcolor.WARNING + 'Error: high neighbour' + bcolor.ENDC
        else:
            out = np.vstack((out,f[i]))
    return(out)

def rad(f,g,lim):

    f = f
    g = g
    lim = lim

    if f.shape[0] > g.shape[0]:
        print "Shape mismatch"
    else:

        out = np.empty((0,f.shape[1]))
        d = np.empty((g.shape[0],1))

        for i in range(f.shape[0]):
           for j in range(g.shape[0]):
               if f[i,-1] != g[j,-1]:
                   d[j] = np.sqrt((f[i,0] - g[j,0])**2 + (f[i,1] - g[j,1])**2)
               else:
                   d[j] = 100*lim
           if min(d) > lim:
               out = np.vstack((out,f[i]))

        return(out)

def readconfig(filename):
    """ A simple config file parser """
    options = {}
    COMMENT_CHAR = '#'
    OPTION_CHAR = '='
    f = open(filename)
    for line in f:
        if COMMENT_CHAR in line:
            line, comment = line.split(COMMENT_CHAR, 1)
        if OPTION_CHAR in line:
            option, value = line.split(OPTION_CHAR, 1)
            option = option.strip()
            value = value.strip()
            options[option] = float(value)
    f.close()
    return options

def writeconfig(filename,options,mode='w+'):
    """ Write the configuratin file """
    options = options
    OPTION_CHAR = '  =  '
    f = open(filename,mode)
    key = sorted(options)
    for i in key:
        f.write(str(i)+str(OPTION_CHAR)+str(options[i])+'\n')
    f.close()

def merge(a,b,c):
    """
    This is a propotype routine to merge files based on the last column of data
    """
    import numpy as np
    d = np.empty((0,a.shape[1]+b.shape[1]-1))
    for i in range(a.shape[0]):
        for j in range(b.shape[0]):
            if a[i,-1] == b[j,-1]:
                d = np.vstack((d,np.hstack((a[i,0:-1],b[j,:]))))

    e = np.empty((0,d.shape[1]+c.shape[1]-1))
    for i in range(c.shape[0]):
        for j in range(d.shape[0]):
            if c[i,-1] == d[j,-1]:
                e = np.vstack((e,np.hstack((d[j,0:-1],c[i,:]))))
    return(e)

def findsky(a,scl):
    """ Find the standard deviation of sky counts. Return median(stdv(sky))
        and the median(sum(sky)/nsky)
    """
    a = pyfits.getdata(a)
    scl = int(scl)
    avoid = 20
    xsize = a.shape[1]/scl
    ysize = a.shape[0]/scl

    b = empty((ysize*xsize,4))
    n = 0

    for i in arange(avoid,ysize-avoid,1):
        for j in arange(avoid,xsize-avoid,1):
            tmp = a[i*scl:(i+1)*scl,j*scl:(j+1)*scl]
            b[n,0] = sum(tmp)/(scl*scl)
            b[n,1] = mean(tmp)
            b[n,2] = median(tmp)
            b[n,3] = std(tmp)
            n = n+1

    return(median(b[:,3]), median(b[:,0]))

#Run tests and see how far does the script got in the last run
basedir = os.path.exists(sys.argv[1])
splits = os.path.exists(sys.argv[1]+'/im1.fits')
base = sys.argv[1]+'/'
infile = sys.argv[1]+'.fits'

if basedir==False:
    os.mkdir(base)

# This should be modified for non-MOSAIC2 data or non-conventional orientation
file = pyfits.open(infile)
header = pyfits.getheader(infile)
img = file[0].data
nx = shape(img)[1]
ny = shape(img)[0]
sy = ny/2
sx = nx/4

print bcolor.OKBLUE + "\nSpliting chips..." + bcolor.ENDC
if splits == False:

    for i in range(4):
        for j in range(2):
            im =  img[j*sy:(j+1)*sy,i*sx:(i+1)*sx]
            im = pyfits.HDUList([pyfits.PrimaryHDU(im,header=header)])
            fname = 'im'+str(j*4 + i+1)+'.fits'
            print "chip "+str(j*4+i+1)+" -> "+fname
            im.writeto(base+fname)
else:
    print bcolor.OKGREEN + "Chips already splited" + bcolor.ENDC

gain = header['GAIN']
rdnoise = header['RDNOISE']
print "Global parameters \n"
print "GAIN:    "+str(gain)
print "RDNOISE: "+str(rdnoise)

#Some parameters that NEVER CHANGE FOR THE WHOLE OBS.

iraf.datapars.setParam('datamin','INDEF')  #Min good data value ~Sky-3Sigma
iraf.datapars.setParam('datamax','50000')  #Max good data value
iraf.datapars.setParam('readnoise',rdnoise)#readnoise
iraf.datapars.setParam('epadu',gain)     #e- per ADU; calculated above
iraf.datapars.setParam('exposure','EXPTIME')    #Integration time
iraf.datapars.setParam('airmass', 'AIRMASS')  #Airmass at middle
iraf.datapars.setParam('filter', 'FILTER')    #Filter used
iraf.findpars.setParam('threshold',4.0)
iraf.findpars.setParam('sharphi',  1)
iraf.findpars.setParam('roundhi', 1.3)
iraf.findpars.setParam('roundlo', -1.3)
iraf.daofind.setParam('verify','no')
iraf.daofind.setParam('interactive','no')
iraf.fitpsf.setParam('box',10.0)
iraf.fitpsf.setParam('verify','no')
iraf.fitpsf.setParam('interactive','no')
iraf.daofind.setParam('verify','no')
iraf.daofind.setParam('interactive','no')
iraf.phot.setParam('output',base+'default')
iraf.phot.setParam('coords',base+'default')
iraf.phot.setParam('verify','no')
iraf.phot.setParam('interactive','no')
iraf.centerpars.setParam('calgori','none')
iraf.fitskypars.setParam('salgorithm','mode')
iraf.photpars.setParam('zmag',25.0)
iraf.photpars.setParam('weighti','constant')
iraf.photpars.setParam('apertur',3.0)
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

#begin chip loop
#change range for all chips
override = False

for imgnu in range(8):

    ext = str(imgnu + 1)
    print "\nRunning chip #"+ext
    FitsFileName = base+"im"+ext+".fits"
    print FitsFileName

    if os.path.exists(base+'im'+ext+'.param') and override == False:
        print 'Using existent param file: '+'im'+ext+'.param'
        pars = readconfig(base+'im'+ext+'.param')
        sigma = pars['sigma']
        sky = pars['sky']
        fwhm = pars['fwhm']
        print 'Skyvalue: ',sky
        print 'Stdv_sky: ',sigma
        print 'Best FWHM: ',fwhm
        iraf.datapars.setParam('sigma',sigma)
        iraf.fitskypars.setParam('skyvalu',sky)
    else:
        pars = {}
        print 'Calculating sky and sigma values ...'
        sigma, sky = findsky(FitsFileName,10)
        print 'Skyvalue: ',sky
        print 'Stdv_sky: ',sigma
        iraf.datapars.setParam('sigma',sigma)
        iraf.fitskypars.setParam('skyvalu',sky)
        pars['gain'] = gain
        pars['rdnoise'] = rdnoise
        pars['sigma'] = sigma
        pars['sky'] = sky

        iraf.datapars.setParam('fwhmpsf',4.0)
        iraf.daofind.setParam('image',FitsFileName)
        iraf.daofind.setParam('output',base+'im'+ext+'.guess.coo')
        iraf.daofind(mode='h',Stdout=1)

        iraf.txdump(textfile=base+'im'+ext+'.guess.coo',fields='xcenter,ycenter,sharpness',expr='sharpness!=INDEF',
        Stdout=base+'daofindtmp')
        daofind = np.loadtxt(base+'daofindtmp')
        if os.path.exists(base+'im'+ext+'.guess.coo') == True: os.remove(base+'im'+ext+'.guess.coo')
        i = (daofind[:,2] < np.median(daofind[:,2]) + 0.12) & (daofind[:,2] > np.median(daofind[:,2]) - 0.12)
        np.savetxt(base+'im'+ext+'.guess.coo',daofind[i,0:2],fmt=['%-10.3f','%-10.3f'])
        iraf.fitpsf.setParam('image',FitsFileName)
        iraf.fitpsf.setParam('coords',base+'im'+ext+'.guess.coo')

        val = linspace(2.,6.,11.)
        best = 0.0
        res = 99999.99
        for i in val:
            iraf.datapars.setParam('fwhmpsf',i)
            iraf.fitpsf.setParam('output',base+'fitpsf')
            iraf.fitpsf(mode='h',Stdout=1)
            asd = iraf.txdump(textfile=base+'fitpsf',fields='rsigma',expr='rsigma!=INDEF && rsigma < 7.0',Stdout=1)
            os.remove(base+'fitpsf')
            tmp = np.array(asd,dtype=float64)
            tmp = tmp*2.35482
#            tmp2 = np.empty((tmp.shape[0],2))
#            tmp2[:,0] = tmp
#            tmp2[:,1] = i
#            np.savetxt(base+'fitpsf'+str(i)+'.out',tmp2,fmt='%10.3f')
            if np.std(tmp) < res:
                best = np.median(tmp)
                res = np.std(tmp)
            else:
                best = best
                res = res

        pars['fwhm'] = best
        print 'Best FWHM: ',best
        fwhm = best
        os.remove(base+'im'+ext+'.guess.coo')

    iraf.daofind.setParam('image',FitsFileName)
    iraf.datapars.setParam('fwhmpsf',fwhm)
    iraf.daofind.setParam('output',base+'default')
    iraf.daofind(mode='h',Stdout=1)

    #config phot
    iraf.phot.setParam('image',FitsFileName)
    iraf.fitskypars.setParam('skyvalue',sky)
    iraf.fitskypars.setParam('annulus',4.*fwhm)
    iraf.fitskypars.setParam('dannulus',2*fwhm)
    iraf.phot(mode='h',Stdout=1)
    iraf.fitpsf.setParam('output',base+'im'+ext+'.fitpsf')
    iraf.fitpsf.setParam('coords',base+'im'+ext+'.fits.coo.1')
    iraf.fitpsf(mode='h',Stdout=1)
    writeconfig(base+'im'+ext+'.param',pars)

    #begin finding best PSF stars

    if os.path.exists(base+'phottmp') == True: os.remove(base+'phottmp')
    if os.path.exists(base+'fitpsftmp') == True: os.remove(base+'fitpsftmp')
    if os.path.exists(base+'daofindtmp') == True: os.remove(base+'daofindtmp')

    iraf.txdump(textfile=FitsFileName+'.mag.1',fields='xcenter,ycenter,mag,msky,merr,id',expr='mag\
    != INDEF && merr != INDEF', Stdout=base+'phottmp')
    mags = np.loadtxt(base+'phottmp')

    iraf.txdump(textfile=FitsFileName+'.coo.1',fields='sharpness,sround,ground,id',expr='sharpness!=INDEF\
    && sround!=INDEF && ground!=INDEF', Stdout=base+'daofindtmp')
    daofind = np.loadtxt(base+'daofindtmp')

    iraf.txdump(textfile=base+'im'+ext+'.fitpsf',fields='rsigma,id',expr='rsigma!=INDEF\
    && rsigma < 7.0', Stdout=base+'fitpsftmp')
    fitpsf = np.loadtxt(base+'fitpsftmp')

    if os.path.exists(base+'phottmp') == True: os.remove(base+'phottmp')
    if os.path.exists(base+'fitpsftmp') == True: os.remove(base+'fitpsftmp')
    if os.path.exists(base+'daofindtmp') == True: os.remove(base+'daofindtmp')

    f = merge(mags,daofind,fitpsf)

    i = (f[:,2] < mean(f[:,2])-0.6*std(f[:,2])) & (f[:,4] < 0.03) &\
    (f[:,5] < median(f[:,5]) + 0.12) & (f[:,5] > median(f[:,5]) - 0.12) & (f[:,3] < sky + 1.5*sigma) &\
    (f[:,8] < 1.05*fwhm/2.355) & (f[:,8] > 0.8*fwhm/2.355) & (f[:,0] > 40*fwhm) & (f[:,0] < 2056-40*fwhm) &\
    (f[:,1] > 40*fwhm) & (f[:,1] < 4060-40*fwhm)

    print bcolor.OKBLUE + '\nQUALITY CUTS ' + bcolor.ENDC
    print 'mag < ',mean(f[:,2])-0.8*std(f[:,2])
    print 'error < 0.03'
    print 's < ',median(f[:,5]) + 0.12
    print 's > ',median(f[:,5]) - 0.12
    print 'fwhm < ',1.05*fwhm
    print 'fwhm > ',0.8*fwhm

    #removes crowded psf stars
    g = rad(f[i],f[:],9*fwhm)
    #removes bad pixel or too low pixel stars
    g = cutbad(g,FitsFileName,sky-4*sigma,sky+3*sigma,48000.0,5*fwhm)

#insert new filter that excludes any source that has any bad pixel (50000) or any too low data value like
# data < sky - 2*stdv

#    np.savetxt(base+'qualityim'+ext+'.dat',f,fmt='%10.3f')
#    np.savetxt(base+'bestim'+ext+'.dat',g,fmt='%10.3f')

    pstfile = open(FitsFileName+'.pst.1', 'a+r')
    pstfile.write("#N ID    XCENTER   YCENTER   MAG         MSKY                                  \\\n")
    pstfile.write("#U ##    pixels    pixels    magnitudes  counts                                \\\n")
    pstfile.write("#F %-9d  %-10.3f   %-10.3f   %-12.3f     %-15.7g                                 \n")
    pstfile.write("#\n")
    np.savetxt(pstfile,g.take([-1,0,1,2,3],axis=1),fmt=['%-9d','%-10.3f','%-10.3f','%-12.3f','%-15.7g'])
    pstfile.close()

    iraf.daopars.setParam('matchra',fwhm)
    iraf.daopars.setParam('psfrad',4*fwhm+1)
    iraf.daopars.setParam('fitrad',fwhm)
    iraf.daopars.setParam('sannulu',2*fwhm)
    iraf.daopars.setParam('wsannul',4*fwhm)
    iraf.psf.setParam('image',FitsFileName)
    iraf.psf(mode='h')
    iraf.seepsf(psfimage=FitsFileName+'.psf.1.fits',image=base+'psfim'+ext+'.fits',magnitu='18.0')
    iraf.allstar.setParam('image',FitsFileName)
    iraf.allstar(mode='h',verbose='no')

    if os.path.exists(base+'allstartmp') == True: os.remove(base+'allstartmp')
    iraf.txdump(textfile=FitsFileName+'.als.1',fields='xcenter,ycenter,mag,merr,id',expr='mag\
    != INDEF && merr != INDEF', Stdout=base+'allstartmp')
    outmags = np.loadtxt(base+'allstartmp')
    if os.path.exists(base+'allstartmp') == True: os.remove(base+'allstartmp')

    if int(ext) <= 4:
        shifty = 0.0
        shiftx = (float(ext) - 1.)*sx
    else:
        shifty = sy
        shiftx = (float(ext) - 5.)*sx

    outmags[:,0] = outmags[:,0] + shiftx
    outmags[:,1] = outmags[:,1] + shifty
    outmags = coconv(outmags,infile,int(ext))
    np.savetxt(base+'im'+ext+'.PHOT',outmags,fmt=['%-10.3f','%-10.3f','%-15.8f','%-15.8f','%-10.3f','%-10.3f','%-6d'])
    addflags(base+'im'+ext+'.PHOT',str(ext),str(infile))
