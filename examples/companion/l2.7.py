#! /usr/bin/env python
#
# Python version of listing 2.6 on p.30 from
#
# Steve Upstill
# The RenderMan Companion - A programmer's Guide to Realistic Computer Graphics
# Addison Wesley, 1992

from ri import *

# ColorCube(): create a unit color cube from smaller cubes
# Parameters:
#  n: the number of minicubes on a side
#  s: a scale factor for each minicube
def ColorCube(n,s):
    if (n<=0):
        return

    RiAttributeBegin()
    RiTranslate(-.5,-.5,-.5)
    RiScale(1.0/n, 1.0/n, 1.0/n)

    for x in range(n):
        for y in range(n):
            for z in range(n):
                color = [ (x+1.0)/n, (y+1.0)/n, (z+1.0)/n ]
                RiColor(color)
                RiTransformBegin()
                RiTranslate(x+.5, y+.5, z+.5)
                RiScale(s,s,s)
                UnitCube()
                RiTransformEnd()
                
    RiAttributeEnd()


# UnitCube(): Enter a unit cube into the scene
def UnitCube():
    square = [[.5,.5,.5],[-.5,.5,.5],[-.5,-.5,.5],[.5,-.5,.5]]

    RiTransformBegin()
    
    # far square
    RiPolygon(RI_P, square)

    # right face
    RiRotate(90.0, 0.0, 1.0, 0.0)
    RiPolygon(RI_P, square)

    # near face
    RiRotate(90.0, 0.0, 1.0, 0.0)
    RiPolygon(RI_P, square)

    # left face
    RiRotate(90.0, 0.0, 1.0, 0.0)
    RiPolygon(RI_P, square)

    RiTransformEnd()

    RiTransformBegin()

    # bottom face
    RiRotate(90.0, 1.0, 0.0, 0.0)
    RiPolygon(RI_P, square)

    RiTransformEnd()

    RiTransformBegin()

    # top face
    RiRotate(-90.0, 1.0, 0.0, 0.0)
    RiPolygon(RI_P, square)

    RiTransformEnd()
    
          

RiBegin(RI_NULL)
RiLightSource("distantlight")
RiProjection("perspective")
RiTranslate(0.0, 0.0, 1.5)
RiRotate(40.0, -1.0, 1.0, 0.0)

RiWorldBegin()
RiSurface("matte")
ColorCube(5,0.3)
RiWorldEnd()
RiEnd()
