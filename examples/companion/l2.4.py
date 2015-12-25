#! /usr/bin/env python
#
# Python version of listing 2.4 on p.26 from
#
# Steve Upstill
# The RenderMan Companion - A programmer's Guide to Realistic Computer Graphics
# Addison Wesley, 1992

from ri import *

Color = [.2,.4,.6]

# UnitCube(): Enter a unit cube into the scene
def UnitCube():
    square = [[.5,.5,.5],[-.5,.5,.5],[-.5,-.5,.5],[.5,-.5,.5]]

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

    # bottom face
    RiRotate(90.0, 1.0, 0.0, 0.0)
    RiPolygon(RI_P, square)

    # top face
    RiRotate(180.0, 1.0, 0.0, 0.0)
    RiPolygon(RI_P, square)
    
          

RiBegin(RI_NULL)
RiLightSource("distantlight")
RiProjection("perspective")
RiTranslate(0.0, 0.0, 1.5)
RiRotate(40.0, -1.0, 1.0, 0.0)

RiWorldBegin()
RiSurface("matte")
RiColor(Color)
UnitCube()
RiWorldEnd()
RiEnd()
