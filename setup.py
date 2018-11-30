from distutils.core import setup
from distutils.command.install import INSTALL_SCHEMES
from os import sys
import os
import re

for scheme in INSTALL_SCHEMES.values():
    scheme['data'] = scheme['purelib']


from imp import find_module
try:
    find_module('numpy')
except ImportError:
    sys.exit('### Error: python module numpy not found')
    
try:
    find_module('scipy')
except ImportError:
    sys.exit('### Error: python module scipy not found')
    
try:
    find_module('astropy')
except ImportError:
    sys.exit('### Error: Neither astropy nor pyfits found.')

try:
    find_module('matplotlib')
except ImportError:
    sys.exit('### Error: python module matplotlib not found')

try:
    find_module('CDFutils')
except ImportError:
    sys.exit('\n*** Error: python module CDFutils not found. ***\n '
             'Download and install: https://github.com/cdfassnacht/CDFutils\n')

try:
    find_module('specim')
except ImportError:
    sys.exit('\n*** Error: python module specim not found. ***\n'
             'Download and install: https://github.com/cdfassnacht/specim\n')

verstr = "unknown"
try:
    parentdir = os.getcwd()+'/'
    verstrline = open(parentdir+'/src/_version.py', "rt").read()
except EnvironmentError:
    pass # Okay, there is no version file.
else:
    VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
    mo = re.search(VSRE, verstrline, re.M)
    if mo:
        verstr = mo.group(1)
    else:
        errstr = 'unable to find version in %s+src/_version.py' % parentdir
        raise RuntimeError(errstr)


setup(
    name = 'keck_code',
    version = verstr,
    author = 'Chris Fassnacht',
    author_email = 'cdfassnacht@ucdavis.edu',
    scripts=[],
    #url = 'lcogt.net',
    license = 'LICENSE.txt',
    description = 'Code for analyzing data from various Keck instruments',
    long_description = open('README.txt').read(),
    requires = ['numpy', 'scipy', 'astropy', 'matplotlib'],
    packages = ['keck_code', 'keck_code.osiris'],
#    package_dir = {'':'src'},
#    package_data = {'esi_redux' : ['data/*']}
)
