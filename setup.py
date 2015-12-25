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

from distutils.core import setup, Extension
import shutil, os, sys, os.path, time
from Cython.Build import cythonize

######################################################################

setup(name="cgkit",
      version="1.2.0",
      description="Python Computer Graphics Kit",
      author="Matthias Baas",
      author_email="baas@ira.uka.de",
      url="http://cgkit.sourceforge.net",
      license="BSD license, see license.txt",
      py_modules=["cgkitinfo", "ri","riutil",
                  "pycgtypes.vec3","pycgtypes.mat3",
                  "pycgtypes.vec4","pycgtypes.mat4","pycgtypes.quat",
                  "sltokenize","sl","slparams", "_slparser"],
      ext_modules=cythonize([Extension("cgtypes", ["cgtypes.pyx"]),
                   Extension("noise", ["src/noisemodule.cpp"]),
                   ]),
      )
