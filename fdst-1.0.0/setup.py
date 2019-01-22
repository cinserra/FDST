from distutils.core import setup
from distutils.command.install import INSTALL_SCHEMES
from os import sys, path
import os,shutil,re
from glob import glob
for scheme in INSTALL_SCHEMES.values():
    scheme['data'] = scheme['purelib']


from imp import find_module
try: find_module('numpy')
except: sys.exit('### Error: python module numpy not found')

try: find_module('astropy')
except: sys.exit('### Error: python module astropy not found')

try: find_module('pyraf')
except: sys.exit('### Error: python module pyraf not found')

try: find_module('matplotlib')
except: sys.exit('### Error: python module matplotlib not found')

try: find_module('scipy')
except: sys.exit('### Error: python module matplotlib not found')


setup(
    name='fdst',
    version='1.0.0',
    author='C.Inserra',
    author_email='InserraC@cardiff.ac.uk',
    classifiers=[
        # How mature is this project?
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 4 - Beta',
        'Intended Audience :: General users',
        'Topic :: Astronomy :: spectra extraction',
        'Programming Language :: Python :: 2.7',
    ],
    scripts=['bin/FDSTfast'],
    url='https://github.com/cinserra',
    license=open('LICENSE.rst').read(),
    description='FDSTfast is a package for Liverpool Telescope (SPRAT) spectra extraction',
    long_description=open('README.rst').read(),
    keywords='LT spectra SPRAT',
    install_requires = ['numpy','astropy','pyraf','matplotlib','scipy'],
    packages=['fdst'],
    package_dir={'':'src'},
    package_data = {'fdst' : ["metadata/*.txt","metadata/arc/*.dat","metadata/extinction/*.dat",\
                                    "metadata/sens/*.txt","metadata/sens/*.fits"]}
    #entry_points={
    #    'console_scripts': [
    #        'sample=sample:main',
    #    ],
    #}
)
