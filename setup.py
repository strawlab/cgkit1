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
from Cython.Build import cythonize

######################################################################

setup(name="cgkit1",
      version="1.2.1", # also update cgkitinfo.py
      description="Python Computer Graphics Kit",
      author="Matthias Baas",
      author_email="baas@ira.uka.de",
      url="http://cgkit.sourceforge.net",
      license="BSD license, see license.txt",
      py_modules=["cgkitinfo",
                  "pycgtypes.vec3","pycgtypes.mat3",
                  "pycgtypes.vec4","pycgtypes.mat4","pycgtypes.quat",
                  ],
      ext_modules=cythonize([Extension("cgtypes", ["cgtypes.pyx"]),
                   ]),
      )
