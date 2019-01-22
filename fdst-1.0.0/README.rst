###########################################################################

							 FDSTfast package
							 ----------------

					C. Inserra  v1.0.0 21/01/2019

###########################################################################
The package contains:
- FDSTfast

Usage:
> FDSTfast (in terminal wherever you want)

Help:
FDSTfast -h
###########################################################################

----------------------------- F.D.S.T. fast -------------------------------

The Fast and Dark Side of Transient experiment Fast extraction pipeline
(FDSTfast) is based on scipy, numpy, pylab and pyraf packages

The user can introduce the files/valuew or hit enter. In the last case the
values between brackets will be taken as default value.

The input is a Liverpool telescope spectrum.
The program will give a wavelength and flux calibrated spectrum both in
.fits and .asci format.

###########################################################################

					ACKNOWLEDGMENTS AND FAQ

Please report any problems to InserraC@cardiff.ac.uk
or open a ticket on https://github.com/cinserra (best way)

###########################################################################

						INSTALLATION

FDSTfast is written in python and requires the following package:

- Python 2.7
   these modules have to be installed:
        - numpy
        - scipy
	    	- pyraf
	    	- astropy
- Iraf

If you have ASTROCONDA installed is even better. Otherwise you can
retrieve it here: https://astroconda.readthedocs.io/en/latest/
The ideal distribution is that with the IRAF legacy value
https://astroconda.readthedocs.io/en/latest/installation.html#legacy-software-stack-with-iraf
(activate conda iraf27)

###########################################################################
1) extract the files from the tarball
> tar -xvf FDSTfast.tar
> cd FDST

2) install the programs as user

YOU DO NOT NEED TO BE ROOT!!!!
Between parentheses there are optionals commands

> python setup.py install  (--prefix=/home/user/fdst) (--record file.txt)

It is preferable to install it under your user directory (in a OSX
system it will be /Users/nameoftheuser/fdst)

You can also install the scripts in the same directory where you downloaded
the file avoiding the --prefix=/home/user/xxx

file.txt will be written, creating a log of the installation

###########################################################################

To uninstall a previous version

3) As a first attempt try to install it again, it should work,
if the version doesn't match proceed to next step

4)
- delete the fdst directory in your site-package path
- delete the fdst****.egg-info from the same directory
- delete all the executable:
FDSTfast

5) to do so type the following commands:

> which FDSTfast
(e.g. /Users/cosmo/anaconda2/envs/iraf27/bin/FDSTfast )
> rm -rf /Users/cosmo/anaconda2/envs/iraf27/bin/FDSTfast
(check if you do not have other bin with FDSTfast)

6) open python:

> python

>>> import fdst
>>> fdst.__path__
(e.g. ['/Users/spxci1/anaconda2/envs/iraf27/lib/python2.7/site-packages/fdst'])

7) then on the terminal (or X11):

> rm -rf /Users/cosmo/anaconda2/envs/iraf27/lib/python2.7/site-packages/fdst
> rm -rf /Users/cosmo/anaconda2/envs/iraf27/lib/python2.7/site-packages/fdst*.egg-info

ALTERNATIVE:
If you used the option " --record files.txt " during the installation
you can run the following command in the terminal:
> cat files.txt | xargs sudo rm -rf


8) Now you can repeat points 1 to 2

###########################################################################
