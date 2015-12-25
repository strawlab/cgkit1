#! /usr/bin/python
# Python version of the colorspheres example that comes with BMRT.

import sys
from ri import *

NFRAMES  = 100 
NSPHERES = 4
FRAMEROT = 15.0

def ColorSpheres(n,  s):
    if (n <= 0):
	return

    RiAttributeBegin()
  
    RiTranslate(-0.5, -0.5, -0.5)
    RiScale(1.0/n, 1.0/n, 1.0/n)

    for x in range(n):
        for y in range(n):
            for z in range(n):
		color = (float(x+1) / n, float(y+1) / n, float(z+1) / n)
	
		RiColor(color)
		RiTransformBegin()
		RiTranslate(x+.5, y+.5, z+.5)
		RiScale(s, s, s)
		RiSphere(0.5, -0.5, 0.5, 360.0)
		RiTransformEnd()

    RiAttributeEnd()


renderer = RI_NULL;

if (len(sys.argv) != 2):
    print >> sys.stderr, "USAGE: %s ribFile|rgl|rendrib\n" % sys.argv[0]
    sys.exit(-1)

renderer = sys.argv[1]

# if the variable renderer is "rgl" or "rendrib", RiBegin() will
# attempt to start up that renderer and pipe its output directly to
# that renderer.  Since RiDisplay is set to put its output to the
# framebuffer, that renderer will attempt to open the framebuffer
# and render directly to it.  If the variable renderer is set to
# some other value, RiBegin() will open that as a file and put the
# RIB commands in it.

RiBegin(renderer)

for frame in range(NFRAMES+1):
    filename = "colorSpheres.%03d.tif" % frame
    
    RiFrameBegin(frame)
    
    RiProjection("perspective")
    RiTranslate(0.0, 0.0, 1.5)
    RiRotate(40.0, -1.0, 1.0, 0.0)

    RiDisplay(filename, RI_FRAMEBUFFER, RI_RGBA)
    RiFormat(256, 192, -1.0)
    RiShadingRate(1.0)

    RiWorldBegin()
    
    RiLightSource("distantlight")
	  
    RiSides(1)
	  
    scale = (NFRAMES-(frame-1.0)) / NFRAMES
    RiRotate(FRAMEROT*frame, 0.0, 0.0, 1.0)
    RiSurface("plastic")
    ColorSpheres(NSPHERES, scale)
    RiWorldEnd()
    RiFrameEnd()

RiEnd()
