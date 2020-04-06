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

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from distutils.core import setup, Extension
import shutil, os, sys, os.path, time
from io import open
from Cython.Build import cythonize
from os import path

######################################################################

# read the contents of README.md
this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, "README.md"), encoding="utf-8") as f:
    long_description = f.read()

setup(
    name="cgkit1",
    version="1.2.4",  # also update cgkitinfo.py
    description="Python Computer Graphics Kit v1",
    long_description=long_description,
    long_description_content_type="text/markdown",
    maintainer="Andrew Straw",
    maintainer_email="strawman@astraw.com",
    url="https://github.com/strawlab/cgkit1",
    license="BSD license, see license.txt",
    py_modules=[
        "cgkitinfo",
        "pycgtypes.vec3",
        "pycgtypes.mat3",
        "pycgtypes.vec4",
        "pycgtypes.mat4",
        "pycgtypes.quat",
    ],
    ext_modules=cythonize([Extension("cgtypes", ["cgtypes.pyx"]),]),
)
