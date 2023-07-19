# cgkit1

[![build](https://github.com/strawlab/cgkit1/workflows/build-and-test/badge.svg?branch=master)](https://github.com/strawlab/cgkit1/actions?query=branch%3Amaster)
[![PyPI version](https://badge.fury.io/py/cgkit1.svg)](https://badge.fury.io/py/cgkit1)

This is a fork of Matthias Bass's python-cgkit v1.2.0.

It has been updated to be compatible with Python 3 and has the renderman and
noise generation code removed.

<hr>
Python Computer Graphics Kit v1.2.0

Copyright (C) 2002, Matthias Baas (see license.txt)
<hr>

The Python Computer Graphics Kit is a collection of Python modules
that contain the basic types and functions to be able to create 3D
computer graphics images. The kit mainly focuses on Pixar's RenderMan
interface, but some modules can also be used for OpenGL programs or
non-RenderMan compliant renderers like POV-Ray, for example.

## Installing the package

The package uses the Python
[build](https://pypa-build.readthedocs.io/en/stable/), so compiling and
installing looks the same on every platform, you simply have to call:

```
python -m build
```

This will compile the C-modules and install everything in the standard
location.

The cgtypes module uses Cython.
