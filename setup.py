######################################################################
# cgkit - setup script
#
# Copyright (C) 2003, Matthias Baas (baas@ira.uka.de)
#
# http://cgkit.sourceforge.net
#
# You may distribute under the terms of the BSD license, as
# specified in the file license.txt.
####################################################################

from setuptools import setup, Extension
from Cython.Build import cythonize

######################################################################

setup(
    ext_modules=cythonize([Extension("cgtypes", ["cgtypes.pyx"]),]),
)
