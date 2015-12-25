#! /usr/bin/env python
#
# Python version of listing 2.8 on p.34 from
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
          

# main()

NFRAMES  = 10    # number of frames in the animation
NCUBES   = 5     # # of minicubes on a side of the color cube
FRAMEROT = 5.0   # # of degrees to rotate cube between frames

RiBegin(RI_NULL)

RiLightSource("distantlight")

RiProjection("perspective")
RiTranslate(0.0, 0.0, 1.5)
RiRotate(40.0, -1.0, 1.0, 0.0)

for frame in range(1,NFRAMES+1):
    filename = "anim%03d.tif" % frame
    RiFrameBegin(frame)
    RiDisplay(filename, RI_FILE, RI_RGBA)
    RiWorldBegin()
    scale = (NFRAMES-(frame-1.0))/NFRAMES
    RiRotate(FRAMEROT*frame, 0.0, 0.0, 1.0)
    RiSurface("matte")
    ColorCube(NCUBES,scale)
    RiWorldEnd()
    RiFrameEnd()
    
RiEnd()
