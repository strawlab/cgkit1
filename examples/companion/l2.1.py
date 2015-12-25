#! /usr/bin/env python
#
# Python version of listing 2.1 on p.18 from
#
# Steve Upstill
# The RenderMan Companion - A programmer's Guide to Realistic Computer Graphics
# Addison Wesley, 1992

from ri import *

Square = [[.5,.5,.5], [.5,-.5,.5],[-.5,-.5,.5],[-.5,.5,.5]]

RiBegin(RI_NULL)
RiWorldBegin()
RiSurface("constant")
RiPolygon(RI_P, Square)
RiWorldEnd()
RiEnd()
