####################################################################
# Utility functions
#
# Copyright (C) 2002, Matthias Baas (baas@ira.uka.de)
#
# http://cgkit.sourceforge.net
#
# You may distribute under the terms of the BSD license, as
# specified in the file license.txt.
####################################################################

import types, sys
from cgtypes import vec3
from ri import *
try:
    import Image
    _PIL_installed = 1
except ImportError:
    _PIL_installed = 0


def RiuDefaultHeader():
    """Outputs a default header into the RIB stream.

    This function can be called right after RiBegin() to write the
    following information into the RIB stream:

    ##RenderMan RIB-Structure 1.1
    ##Creator <Filename>
    ##CreationDate <Date>
    ##For <User>

    The "For" information is left out if the user name can't be determined.
    """
    RiArchiveRecord(RI_STRUCTURE,"RenderMan RIB-Structure 1.1")
    RiArchiveRecord(RI_STRUCTURE,"Creator %s",sys.argv[0])
    RiArchiveRecord(RI_STRUCTURE,"CreationDate %s",time.ctime())
    try:
        RiArchiveRecord(RI_STRUCTURE,"For %s",getpass.getuser())
    except:
        pass


def RiuArrow(height=1.0, thickness=0.1, headheight=0.2, headscale=1.7):
    """Outputs an arrow.

    The arrow starts in the origin and ends in (0,0,height). thickness is the
    overall thickness of the arrow, headheight and headscale specify the
    height and a factor for the radius of the arrow head.
    """

    # Radius for arrow head (paraboloid)
    rmax = headscale*thickness
    a = headheight/(rmax*rmax)
    offset = a*thickness*thickness
    RiTransformBegin()
    RiCylinder(thickness, 0,height-0.01-offset, 360)
    RiTranslate(0,0,height)
#    RiRotate(180,1,0,0)
    RiParaboloid(rmax, 0,-headheight,360)
    RiTransformEnd()

def RiuCoordSystem(thickness=0.06, shader="matte"):
    """Outputs a coordinate system.

    The X-,Y- and Z-axis are colored red, green and blue. thickness
    specifies the thickness of the arrows and shader is used as
    surface shader for the entire coordinate system.
    """

    RiAttributeBegin()
    RiSurface(shader)
    # X axis
    RiColor((1,0,0))
    RiTransformBegin()
    RiRotate(90, 0,1,0)
    RiuArrow(thickness=thickness)
    RiTransformEnd()

    # Y axis
    RiColor((0,1,0))
    RiTransformBegin()
    RiRotate(-90, 1,0,0)
    RiuArrow(thickness=thickness)
    RiTransformEnd()

    # Z axis
    RiColor((0,0,1))
    RiuArrow(thickness=thickness)

    RiAttributeEnd()


def RiuGrid(thickness=0.02, cells=6, shader="matte", color=(0.9,0.9,0.9)):
    """Output a grid primitive.

    thickness determines the thickness of the grid lines and cells the
    number of gridlines.  The grid lies on the XY plane and is
    centered at the origin. The grid spacing is 1 unit. The grid uses
    the given shader and color.
    """

    if cells%2==1: cells+=1

    RiAttributeBegin()
    RiSurface(shader)
    RiColor(color)

    RiTransformBegin()
    RiRotate(-90,1,0,0)
    RiTranslate(-cells/2,0,-cells/2)
    for i in range(cells+1):
        RiCylinder(thickness,0,cells,360)
        RiTranslate(1,0,0)
    RiTransformEnd()

    RiTransformBegin()
    RiRotate(90,0,1,0)
    RiTranslate(0,-cells/2,-cells/2)
    for i in range(cells+1):
        RiCylinder(thickness,0,cells,360)
        RiTranslate(0,1,0)
    RiTransformEnd()

    RiAttributeEnd()


def RiuHeightfield(image):
    if not _PIL_installed:
        raise ImportError, "the Python Imaging Library (PIL) is not installed"

    if type(image)==types.StringType or type(image)==types.UnicodeType:
        image = Image.open(image)

    points = []
    width, height = image.size
    for j in range(height):
        y = float(j)/(height-1)
        for i in range(width):
            x = float(i)/(height-1)
            col=image.getpixel((i,j))
            z=col[0]/255.0
            p=vec3(x,y,z)
            points.append(p)

    RiPatchMesh(RI_BILINEAR, width, RI_NONPERIODIC, height, RI_NONPERIODIC, P=points)


##########################################################

if __name__=='__main__':
    pass
