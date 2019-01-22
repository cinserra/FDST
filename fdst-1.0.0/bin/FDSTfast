#!/usr/bin/env python

#import pyfits
import astropy.io.fits
import shutil
import time
import string
import re
import os
import fdst
from numpy import *
from scipy import *
from pylab import *

######################## for the help ########################
from optparse import OptionParser

description = " Fast reduction of LT+SPRAT spectra "
usage = "%prog "
if __name__ == "__main__":
    parser = OptionParser(usage=usage, description=description, version="%prog " + str(fdst.__version__))
    parser.add_option("-v", "--verbose",dest="verbose",\
                  action="store_true",default=False,
                  help='Print tasks description')
    option,args = parser.parse_args()

##############################################################
from pyraf import iraf

iraf.noao(_doprint=0)
iraf.imred(_doprint=0)
iraf.ccdred(_doprint=0)
iraf.twodspec(_doprint=0)
iraf.longslit(_doprint=0)
iraf.onedspec(_doprint=0)
iraf.specred(_doprint=0)
toforget = ['ccdproc', 'imcopy', 'specred.apall', 'longslit.identify', 'longslit.reidentify', 'specred.standard',
            'longslit.fitcoords', 'onedspec.wspectext']
for t in toforget:
    iraf.unlearn(t)

##################### opening setup ##########################

print ''
question = raw_input('Do you have a list of spectra ? ([yes]/no) ')
if not question:
    question = 'yes'

print ''
if question == 'yes' or question == 'y' or question == 'Y' or question == 'Yes' or question == 'YES':
	files = raw_input('List ? [e.g. list, list.txt, list.dat] ')
	lcf = open(files,'r')
	riga = lcf.readlines()
	lcf.close()
	snlist = []
	for line in riga:
		p = line.split()
		snlist.append(p[0])
else:
	files = raw_input('List the spectra (space separated list) to use: ')
	snlist = string.split(files)

######################### load metadata files ################
direc=fdst.__path__[0]+'/metadata/'
sensfile = direc+'sens/sensfile.fits'
atmofile = direc+'sens/atmofile.fits'
atmofile05 = direc+'sens/atmofile05.fits'
atmofile06 = direc+'sens/atmofile06.fits'
atmofile07 = direc+'sens/atmofile07.fits'
atmofile08 = direc+'sens/atmofile08.fits'
atmofile09 = direc+'sens/atmofile09.fits'
atmofile11 = direc+'sens/atmofile11.fits'
atmofile12 = direc+'sens/atmofile12.fits'
atmofile13 = direc+'sens/atmofile13.fits'
atmofile14 = direc+'sens/atmofile14.fits'
atmofile15 = direc+'sens/atmofile15.fits'
extinctionfile = direc+'extinction/lapalma_ext.dat'
##############################################################
os.system("rm -rf spectrum*.fits")
################ loop extranction and calibration ############
now = time.time()
ii = 0
while ii != len(snlist):
    _snname = snlist[ii]
    iraf.imcopy(_snname + "[3]",'spectrum.fits')
    spec = 'spectrum.fits'

    ################ load kheader keywords ###############
    hdulist = astropy.io.fits.open(_snname)
    prihdr = hdulist[0].header
    AIRMASS = prihdr['AIRMASS']
    EXPTIME = prihdr['EXPTIME']
    MJD = prihdr['MJD']
    MJDs = "%0.3f" % (MJD)
    DATE = prihdr['date-obs']
    DATEf = DATE.replace("-","")
    DATEs = DATEf[0:8]
    OBJECT = prihdr['OBJECT']
    #######################################################

    specf = re.sub('.fits', '_f.fits', spec)
    iraf.onedspec.calibrate(input=spec, output=specf, sensiti=sensfile, extinct='yes',flux='yes', ignorea='yes', extinction=extinctionfile, observatory='lapalma', airmass=AIRMASS, exptime=EXPTIME, fnu='no')

    spect10 = re.sub('.fits', '_e10.fits', specf)
    iraf.sarit(input1=specf,op='/',input2=atmofile, output=spect10)
    spect09 = re.sub('.fits', '_e09.fits', specf)
    iraf.sarit(input1=specf,op='/',input2=atmofile09, output=spect09)
    spect08 = re.sub('.fits', '_e08.fits', specf)
    iraf.sarit(input1=specf,op='/',input2=atmofile08, output=spect08)
    spect07 = re.sub('.fits', '_e07.fits', specf)
    iraf.sarit(input1=specf,op='/',input2=atmofile07, output=spect07)
    spect06 = re.sub('.fits', '_e06.fits', specf)
    iraf.sarit(input1=specf,op='/',input2=atmofile06, output=spect06)
    spect05 = re.sub('.fits', '_e05.fits', specf)
    iraf.sarit(input1=specf,op='/',input2=atmofile05, output=spect05)
    spect11 = re.sub('.fits', '_e11.fits', specf)
    iraf.sarit(input1=specf,op='/',input2=atmofile11, output=spect11)
    spect12 = re.sub('.fits', '_e12.fits', specf)
    iraf.sarit(input1=specf,op='/',input2=atmofile12, output=spect12)
    spect13 = re.sub('.fits', '_e13.fits', specf)
    iraf.sarit(input1=specf,op='/',input2=atmofile13, output=spect13)
    spect14 = re.sub('.fits', '_e14.fits', specf)
    iraf.sarit(input1=specf,op='/',input2=atmofile14, output=spect14)
    spect15 = re.sub('.fits', '_e15.fits', specf)
    iraf.sarit(input1=specf,op='/',input2=atmofile15, output=spect15)

    ax= subplot(111)
    xl = array([6500, 8000])

    listatmo =[spect05,spect06,spect07,spect08,spect09,spect10,spect11,spect12,spect13,spect14,spect15]
    i = 0
    while i != len(listatmo):
        iraf.imdel('tmp.fits',verify='no')
        sp = listatmo[i] + "[*,1,1]"
        iraf.imcopy(sp,'tmp.fits',verbose='no')
        xx = iraf.listpixel('tmp.fits',wcs='world',Stdout=1)
        lam,fl =[],[]
        for x in xx:
             lam.append(float(string.split(x)[0]))
             fl.append(float(string.split(x)[1]))

        fl2 = []
        for z in fl:
            fl2.append(float(z))

        cr =['#ee560b','orange','#e6d114','g','c','k','b','purple','m','pink','r']
        plot(lam,fl2,color=cr[i],ls='-',lw='1')

        i = i +1

    os.remove('tmp.fits')
    xlim(xl[0],xl[1])
    ax.set_xlabel(r'$\rm{Observed\hspace{0.3}Wavelength\hspace{0.3} (\AA)}$', fontsize=18)
    ax.set_ylabel(r'$\rm{F_\lambda\,(erg\, s^{-1}\, cm^{-2}\, \AA^{-1})}$', fontsize=18)
    ax.minorticks_on()
    ax.xaxis.set_minor_locator(MultipleLocator(100))
    show()

    print ''
    _questionat = raw_input('Which one does remove better the telluric lines (give the number) ? (1=dark orange, 2=orange, 3=yellow, 4=green, 5=cyan, [6]=black, 7=blue, 8=purple, 9=magenta, 10=pink, 11=red) ')
    if not _questionat:
        questionat = 6
    else:
        questionat = float(_questionat)

    if questionat == 1:
        spect = re.sub('.fits', '_t.fits', specf)
        iraf.sarit(input1=specf,op='/',input2=atmofile05, output=spect)
        os.system("rm -rf spec*_e*.fits")
    elif questionat == 2:
        spect = re.sub('.fits', '_t.fits', specf)
        iraf.sarit(input1=specf,op='/',input2=atmofile06, output=spect)
        os.system("rm -rf spec*_e*.fits")
    elif questionat == 3:
        spect = re.sub('.fits', '_t.fits', specf)
        iraf.sarit(input1=specf,op='/',input2=atmofile07, output=spect)
        os.system("rm -rf spec*_e*.fits")
    elif questionat == 4:
        spect = re.sub('.fits', '_t.fits', specf)
        iraf.sarit(input1=specf,op='/',input2=atmofile08, output=spect)
        os.system("rm -rf spec*_e*.fits")
    elif questionat == 5:
        spect = re.sub('.fits', '_t.fits', specf)
        iraf.sarit(input1=specf,op='/',input2=atmofile09, output=spect)
        os.system("rm -rf spec*_e*.fits")
    elif questionat == 6:
        spect = re.sub('.fits', '_t.fits', specf)
        iraf.sarit(input1=specf,op='/',input2=atmofile, output=spect)
        os.system("rm -rf spec*_e*.fits")
    elif questionat == 7:
        spect = re.sub('.fits', '_t.fits', specf)
        iraf.sarit(input1=specf,op='/',input2=atmofile11, output=spect)
        os.system("rm -rf spec*_e*.fits")
    elif questionat == 8:
        spect = re.sub('.fits', '_t.fits', specf)
        iraf.sarit(input1=specf,op='/',input2=atmofile12, output=spect)
        os.system("rm -rf spec*_e*.fits")
    elif questionat == 9:
        spect = re.sub('.fits', '_t.fits', specf)
        iraf.sarit(input1=specf,op='/',input2=atmofile13, output=spect)
        os.system("rm -rf spec*_e*.fits")
    elif questionat == 10:
        spect = re.sub('.fits', '_t.fits', specf)
        iraf.sarit(input1=specf,op='/',input2=atmofile14, output=spect)
        os.system("rm -rf spec*_e*.fits")
    elif questionat == 11:
        spect = re.sub('.fits', '_t.fits', specf)
        iraf.sarit(input1=specf,op='/',input2=atmofile15, output=spect)
        os.system("rm -rf spec*_e*.fits")

    fileName, fileExtension = os.path.splitext(spect)
    shutil.copyfile(fileName+fileExtension,str(OBJECT)+'_'+DATEs+'_'+str(MJDs)+'_e'+fileExtension)
    iraf.wspec(str(OBJECT)+'_'+DATEs+'_'+str(MJDs)+'_e'+fileExtension,str(OBJECT)+'_'+DATEs+'_'+str(MJDs)+'_e.asci', header='no')
    os.remove('spectrum_f_t.fits')
    os.remove('spectrum.fits')
    os.remove('spectrum_f.fits')

    ii = ii + 1

then = time.time()
time = then -now
print ''
print 'Extraction and calibration done in %.0is ' % (time)
print ''
