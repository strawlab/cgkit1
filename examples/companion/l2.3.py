#! /usr/bin/env python
#
# Python version of listing 2.3 on p.25 from
#
# Steve Upstill
# The RenderMan Companion - A programmer's Guide to Realistic Computer Graphics
# Addison Wesley, 1992

from ri import *

Color = [.2,.4,.6]

# UnitCube(): define a cube in the graphics environment
def UnitCube():
    L = -.5    # For x: left side
    R = .5     # For x: right side
    D = -.5    # For y: down side
    U = .5     # For y: upper side
    F = .5     # For z: far side
    N = -.5    # For z: near side

    Cube = [
        [ [L,D,F], [R,D,F], [R,D,N], [L,D,N] ],   # Bottom face
        [ [L,D,F], [L,U,F], [L,U,N], [L,D,N] ],   # Left face
        [ [R,U,N], [L,U,N], [L,U,F], [R,U,F] ],   # Top face
        [ [R,U,N], [R,U,F], [R,D,F], [R,D,N] ],   # Right face
        [ [R,D,F], [R,U,F], [L,U,F], [L,D,F] ],   # Far face
        [ [L,U,N], [R,U,N], [R,D,N], [L,D,N] ]    # Near face
    ]

    # declare the cube
    for i in range(6):
        RiPolygon(RI_P, Cube[i])


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
