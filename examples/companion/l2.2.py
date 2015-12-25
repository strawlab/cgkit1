#! /usr/bin/env python
#
# Python version of listing 2.2 on p.20 from
#
# Steve Upstill
# The RenderMan Companion - A programmer's Guide to Realistic Computer Graphics
# Addison Wesley, 1992

from ri import *

Square = [[.5,.5,0], [.5,-.5,0],[-.5,-.5,0],[-.5,.5,0]]

Color = [.2,.4,.6]

RiBegin(RI_NULL)
RiLightSource("distantlight")

RiProjection("perspective")
RiTranslate(0.0, 0.0, 1.0)
RiRotate(40.0, -1.0, 1.0, 0.0)

RiWorldBegin()

RiSurface("matte")
RiColor(Color)
RiPolygon(RI_P, Square)

RiWorldEnd()
RiEnd()
