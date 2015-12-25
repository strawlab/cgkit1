#! /usr/bin/env python
#
# Python version of listing 4.1 on p.61 from
#
# Steve Upstill
# The RenderMan Companion - A programmer's Guide to Realistic Computer Graphics
# Addison Wesley, 1992

from ri import *


OFFSET = 1.2

def ShowQuads():
    RiRotate(-90.0, 1.0, 0.0, 0.0)

    RiTranslate(-OFFSET, 0.0, (OFFSET/2.0) )
    RiSphere(0.5, -0.5, 0.5, 360)

    RiTranslate(OFFSET, 0.0, 0.0)
    RiTranslate(0.0, 0.0, -0.5)
    RiCone(1.0, 0.5, 360)

    RiTranslate(0.0, 0.0, 0.5)
    RiTranslate(OFFSET, 0.0, 0.0)
    RiCylinder(0.5, -0.5, 0.5, 360)

    RiTranslate(-(OFFSET*2),0.0, -OFFSET)
    RiHyperboloid([0.4,-0.4,-0.4], [0.4,0.4,0.4], 360)

    RiTranslate(OFFSET, 0.0, -0.5)
    RiParaboloid(0.5, 0.0, 0.9, 360)

    RiTranslate(OFFSET, 0.0, 0.5)
    RiTorus(0.4, 0.15, 0.0, 360, 360)
    

# main()

RiBegin(RI_NULL)

RiLightSource("distantlight")

RiProjection("perspective", fov=32)
RiTranslate(0,0,5)

RiWorldBegin()
ShowQuads()
RiWorldEnd()

RiEnd()
