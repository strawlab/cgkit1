/*
 Noise module

 The basic noise implementation is a Perlin noise as described in
 Texturing & Modeling.

 Copyright (C) 2002, Matthias Baas (baas@ira.uka.de)

 You may distribute under the terms of the BSD license, as
 specified in the file license.txt.
*/

#include <iostream.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <Python.h>

// Table mask which is used to create indexing values
// (must be length of the tables - 1)
#define TABMASK 0xff

#include "noisetabs.h"

// Vector class objects from the cgtypes module.
// They're initialized in the modules init function and are used to 
// create return values for the vector functions (v...)
static PyObject* vec3 = NULL;
static PyObject* vec4 = NULL;

// Period values for the pnoise() function family
// (this should actually be parameters to the pnoise() functions, but
// then the tabindex templates wouldn't work anymore)
static int xperiod = 1;
static int yperiod = 1;
static int zperiod = 1;
static int tperiod = 1;
static int poffset = 0;

/*----------------------------------------------------------------------
  Modulo function.
  Returns a mod b.
----------------------------------------------------------------------*/
inline int imod(int a, int b)
{
  int n = (int)(a/b);
  a-=n*b;
  if (a<0)
    a+=b;
  return a;
}

/*----------------------------------------------------------------------
  Tabindex functions.

  They create a "random" index between 0 and 255 which is dependend
  only on up to 4 integer input values.
----------------------------------------------------------------------*/

inline unsigned char tabindex2(int ix, int iy)
{
  return perm[(ix + perm[iy&TABMASK])&TABMASK];
}

inline unsigned char tabindex3(int ix, int iy, int iz)
{
  return perm[(ix + perm[(iy + perm[iz&TABMASK])&TABMASK])&TABMASK];
}

inline unsigned char tabindex4(int ix, int iy, int iz, int it)
{
  return perm[(it + perm[(ix + perm[(iy + perm[iz&TABMASK])&TABMASK])&TABMASK])&TABMASK];
}

/*----------------------------------------------------------------------
  Periodic tabindex functions.

  They create a "random" index between 0 and 255 which is dependend
  only on up to 4 integer input values.
  The return value is periodic with a perdiod stored in the global
  variables xperiod, yperiod, zperiod and tperiod.
  These period variables have to be set before calling the functions.
----------------------------------------------------------------------*/

inline unsigned char ptabindex2(int ix, int iy)
{
  ix=imod(ix,xperiod);
  iy=imod(iy,yperiod);
  return perm[(ix + perm[iy&TABMASK])&TABMASK];
}

inline unsigned char ptabindex3(int ix, int iy, int iz)
{
  ix=imod(ix,xperiod);
  iy=imod(iy,yperiod);
  iz=imod(iz,zperiod);
  return perm[(ix + perm[(iy + perm[iz&TABMASK])&TABMASK])&TABMASK];
}

inline unsigned char ptabindex4(int ix, int iy, int iz, int it)
{
  ix=imod(ix,xperiod);
  iy=imod(iy,yperiod);
  iz=imod(iz,zperiod);
  it=imod(it,tperiod);
  return perm[(it + perm[(ix + perm[(iy + perm[iz&TABMASK])&TABMASK])&TABMASK])&TABMASK];
}

/*----------------------------------------------------------------------
  Periodic tabindex functions with offset.

  Same as before, but they use an additional offset (global variable
  poffset) that is incorporated in the calculation.
  These versions are used for the vector version of the pnoise()
  functions.
----------------------------------------------------------------------*/

inline unsigned char ptabindex2offset(int ix, int iy)
{
  ix=imod(ix,xperiod);
  iy=imod(iy,yperiod)+poffset;
  return perm[(ix + perm[iy&TABMASK])&TABMASK];
}

inline unsigned char ptabindex3offset(int ix, int iy, int iz)
{
  ix=imod(ix,xperiod);
  iy=imod(iy,yperiod);
  iz=imod(iz,zperiod)+poffset;
  return perm[(ix + perm[(iy + perm[iz&TABMASK])&TABMASK])&TABMASK];
}

inline unsigned char ptabindex4offset(int ix, int iy, int iz, int it)
{
  ix=imod(ix,xperiod);
  iy=imod(iy,yperiod);
  iz=imod(iz,zperiod)+poffset;
  it=imod(it,tperiod);
  return perm[(it + perm[(ix + perm[(iy + perm[iz&TABMASK])&TABMASK])&TABMASK])&TABMASK];
}

/*----------------------------------------------------------------------
 lerp - Linear interpolation between a and b
 ----------------------------------------------------------------------*/
inline double lerp(double t, double a, double b)
{
  return a+(b-a)*t;
}

/*----------------------------------------------------------------------
  cellnoise

  Returns a value between 0 and 1 which is constant between integer
  lattice points.
----------------------------------------------------------------------*/
static double cellnoise(double x, double y, double z, double t)
{
  int ix = int(floor(x));
  int iy = int(floor(y));
  int iz = int(floor(z));
  int it = int(floor(t));
  return uniform[tabindex4(ix,iy,iz,it)];
};

/*----------------------------------------------------------------------
  vcellnoise (2D)

  Returns a value between 0 and 1 which is constant between integer
  lattice points.
----------------------------------------------------------------------*/
static void vcellnoise(double x, double y,
			 double& ox, double& oy)
{
  ox = cellnoise(x,y,0.0,0.0);
  x+=10.0;
  oy = cellnoise(x,y,0.0,0.0);
};

/*----------------------------------------------------------------------
  vcellnoise (3D)

  Returns a value between 0 and 1 which is constant between integer
  lattice points.
----------------------------------------------------------------------*/
static void vcellnoise(double x, double y, double z,
			 double& ox, double& oy, double& oz)
{
  ox = cellnoise(x,y,z,0.0);
  x+=10.0;
  oy = cellnoise(x,y,z,0.0);
  y+=10.0;
  oz = cellnoise(x,y,z,0.0);
};

/*----------------------------------------------------------------------
  vcellnoise (4D)

  Returns a value between 0 and 1 which is constant between integer
  lattice points.
----------------------------------------------------------------------*/
static void vcellnoise(double x, double y, double z, double t,
			 double& ox, double& oy, double& oz, double& ot)
{
  ox = cellnoise(x,y,z,t);
  x+=10.0;
  oy = cellnoise(x,y,z,t);
  y+=10.0;
  oz = cellnoise(x,y,z,t);
  z+=10.0;
  ot = cellnoise(x,y,z,t);
};


/*----------------------------------------------------------------------
  noise template (2D)

  Noise template for a 2D noise or pnoise function.
  The template parameter specifies the function that's used to create
  an index for the random values. If this function is periodic you'll
  also get a periodic noise.
  Returns a value between -1 and 1 
----------------------------------------------------------------------*/
template<class T_IndexFunc>
static inline double noise_template(const T_IndexFunc idxfunc, double x, double y)
{
  // Determine cell
  int ix = int(floor(x));
  int iy = int(floor(y));

  // Coordinate components relative to the cell vertices
  double rx0 = x-ix;
  double ry0 = y-iy;
  double rx1 = rx0-1.0;
  double ry1 = ry0-1.0;

  double sx = rx0*rx0*(3.0-2.0*rx0);
  double sy = ry0*ry0*(3.0-2.0*ry0);

  double *g;
  double u,v;
  double a,b;

  g = grads2[idxfunc(ix,iy)];
  u = g[0]*rx0 + g[1]*ry0;
  g = grads2[idxfunc(ix+1,iy)];
  v = g[0]*rx1 + g[1]*ry0;
  a = lerp(sx,u,v);

  g = grads2[idxfunc(ix,iy+1)];
  u = g[0]*rx0 + g[1]*ry1;
  g = grads2[idxfunc(ix+1,iy+1)];
  v = g[0]*rx1 + g[1]*ry1;
  b = lerp(sx,u,v);

  return lerp(sy,a,b);
};

/*---------------------------------------------------------------------- 
  noise template (3D)

  Noise template for a 3D noise or pnoise function.
  The template parameter specifies the function that's used to create
  an index for the random values. If this function is periodic you'll
  also get a periodic noise.
  Returns a value between -1 and 1 
----------------------------------------------------------------------*/
template<class T_IndexFunc>
static inline double noise_template(const T_IndexFunc idxfunc, double x, double y, double z)
{
  // Determine cell
  int ix = int(floor(x));
  int iy = int(floor(y));
  int iz = int(floor(z));

  // Coordinate components relative to the cell vertices
  double rx0 = x-ix;
  double ry0 = y-iy;
  double rz0 = z-iz;
  double rx1 = rx0-1.0;
  double ry1 = ry0-1.0;
  double rz1 = rz0-1.0;

  double sx = rx0*rx0*(3.0-2.0*rx0);
  double sy = ry0*ry0*(3.0-2.0*ry0);
  double sz = rz0*rz0*(3.0-2.0*rz0);

  double *g;
  double u,v;
  double a,b,c,d;

  g = grads3[idxfunc(ix,iy,iz)];
  u = g[0]*rx0 + g[1]*ry0 + g[2]*rz0;
  g = grads3[idxfunc(ix+1,iy,iz)];
  v = g[0]*rx1 + g[1]*ry0 + g[2]*rz0;
  a = lerp(sx,u,v);

  g = grads3[idxfunc(ix,iy+1,iz)];
  u = g[0]*rx0 + g[1]*ry1 + g[2]*rz0;
  g = grads3[idxfunc(ix+1,iy+1,iz)];
  v = g[0]*rx1 + g[1]*ry1 + g[2]*rz0;
  b = lerp(sx,u,v);

  c = lerp(sy,a,b);

  g = grads3[idxfunc(ix,iy,iz+1)];
  u = g[0]*rx0 + g[1]*ry0 + g[2]*rz1;
  g = grads3[idxfunc(ix+1,iy,iz+1)];
  v = g[0]*rx1 + g[1]*ry0 + g[2]*rz1;
  a = lerp(sx,u,v);

  g = grads3[idxfunc(ix,iy+1,iz+1)];
  u = g[0]*rx0 + g[1]*ry1 + g[2]*rz1;
  g = grads3[idxfunc(ix+1,iy+1,iz+1)];
  v = g[0]*rx1 + g[1]*ry1 + g[2]*rz1;
  b = lerp(sx,u,v);

  d = lerp(sy,a,b);
  
  return lerp(sz,c,d);
};

/*----------------------------------------------------------------------
  noise template (4D)

  Noise template for a 4D noise or pnoise function.
  The template parameter specifies the function that's used to create
  an index for the random values. If this function is periodic you'll
  also get a periodic noise.
  Returns a value between -1 and 1 
----------------------------------------------------------------------*/
template<class T_IndexFunc>
static inline double noise_template(const T_IndexFunc idxfunc, double x, double y, double z, double t)
{
  // Determine cell
  int ix = int(floor(x));
  int iy = int(floor(y));
  int iz = int(floor(z));
  int it = int(floor(t));

  // Coordinate components relative to the cell vertices
  double rx0 = x-ix;
  double ry0 = y-iy;
  double rz0 = z-iz;
  double rt0 = t-it;
  double rx1 = rx0-1.0;
  double ry1 = ry0-1.0;
  double rz1 = rz0-1.0;
  double rt1 = rt0-1.0;

  double sx = rx0*rx0*(3.0-2.0*rx0);
  double sy = ry0*ry0*(3.0-2.0*ry0);
  double sz = rz0*rz0*(3.0-2.0*rz0);
  double st = rt0*rt0*(3.0-2.0*rt0);

  double *g;
  double u,v;
  double a,b,c,d,e,f;

  // at it
  g = grads4[idxfunc(ix,iy,iz,it)];
  u = g[0]*rx0 + g[1]*ry0 + g[2]*rz0 + g[3]*rt0;
  g = grads4[idxfunc(ix+1,iy,iz,it)];
  v = g[0]*rx1 + g[1]*ry0 + g[2]*rz0 + g[3]*rt0;
  a = lerp(sx,u,v);

  g = grads4[idxfunc(ix,iy+1,iz,it)];
  u = g[0]*rx0 + g[1]*ry1 + g[2]*rz0 + g[3]*rt0;
  g = grads4[idxfunc(ix+1,iy+1,iz,it)];
  v = g[0]*rx1 + g[1]*ry1 + g[2]*rz0 + g[3]*rt0;
  b = lerp(sx,u,v);

  c = lerp(sy,a,b);

  g = grads4[idxfunc(ix,iy,iz+1,it)];
  u = g[0]*rx0 + g[1]*ry0 + g[2]*rz1 + g[3]*rt0;
  g = grads4[idxfunc(ix+1,iy,iz+1,it)];
  v = g[0]*rx1 + g[1]*ry0 + g[2]*rz1 + g[3]*rt0;
  a = lerp(sx,u,v);

  g = grads4[idxfunc(ix,iy+1,iz+1,it)];
  u = g[0]*rx0 + g[1]*ry1 + g[2]*rz1 + g[3]*rt0;
  g = grads4[idxfunc(ix+1,iy+1,iz+1,it)];
  v = g[0]*rx1 + g[1]*ry1 + g[2]*rz1 + g[3]*rt0;
  b = lerp(sx,u,v);

  d = lerp(sy,a,b);

  e = lerp(sz,c,d);

  // at it+1
  g = grads4[idxfunc(ix,iy,iz,it+1)];
  u = g[0]*rx0 + g[1]*ry0 + g[2]*rz0 + g[3]*rt1;
  g = grads4[idxfunc(ix+1,iy,iz,it+1)];
  v = g[0]*rx1 + g[1]*ry0 + g[2]*rz0 + g[3]*rt1;
  a = lerp(sx,u,v);

  g = grads4[idxfunc(ix,iy+1,iz,it+1)];
  u = g[0]*rx0 + g[1]*ry1 + g[2]*rz0 + g[3]*rt1;
  g = grads4[idxfunc(ix+1,iy+1,iz,it+1)];
  v = g[0]*rx1 + g[1]*ry1 + g[2]*rz0 + g[3]*rt1;
  b = lerp(sx,u,v);

  c = lerp(sy,a,b);

  g = grads4[idxfunc(ix,iy,iz+1,it+1)];
  u = g[0]*rx0 + g[1]*ry0 + g[2]*rz1 + g[3]*rt1;
  g = grads4[idxfunc(ix+1,iy,iz+1,it+1)];
  v = g[0]*rx1 + g[1]*ry0 + g[2]*rz1 + g[3]*rt1;
  a = lerp(sx,u,v);

  g = grads4[idxfunc(ix,iy+1,iz+1,it+1)];
  u = g[0]*rx0 + g[1]*ry1 + g[2]*rz1 + g[3]*rt1;
  g = grads4[idxfunc(ix+1,iy+1,iz+1,it+1)];
  v = g[0]*rx1 + g[1]*ry1 + g[2]*rz1 + g[3]*rt1;
  b = lerp(sx,u,v);

  d = lerp(sy,a,b);

  f = lerp(sz,c,d);
  
  return lerp(st,e,f);
};

/*----------------------------------------------------------------------
  noise (2D)
  Returns a value between -1 and 1 
----------------------------------------------------------------------*/
static double noise(double x, double y)
{
  return noise_template(tabindex2,x,y);
};

/*---------------------------------------------------------------------- 
  noise (3D)
  Returns a value between -1 and 1 
----------------------------------------------------------------------*/
static double noise(double x, double y, double z)
{
  return noise_template(tabindex3,x,y,z);
};

/*----------------------------------------------------------------------
  noise (4D)
  Returns a value between -1 and 1 
----------------------------------------------------------------------*/
static double noise(double x, double y, double z, double t)
{
  return noise_template(tabindex4,x,y,z,t);
};

/*---------------------------------------------------------------------- 
   vnoise (2D)
   Returns values between -1 and 1 
----------------------------------------------------------------------*/
static void vnoise(double x, double y, 
		   double& ox, double& oy)
{
  ox = noise_template(tabindex2,x,y);
  x += 10.0;
  oy = noise_template(tabindex2,x,y);
};

/*---------------------------------------------------------------------- 
   vnoise (3D)
   Returns values between -1 and 1 
----------------------------------------------------------------------*/
static void vnoise(double x, double y, double z, 
		   double& ox, double& oy, double& oz)
{
  ox = noise_template(tabindex3,x,y,z);
  x += 10.0;
  oy = noise_template(tabindex3,x,y,z);
  y += 10.0;
  oz = noise_template(tabindex3,x,y,z);
};

/*---------------------------------------------------------------------- 
   vnoise (4D)
   Returns values between -1 and 1 
----------------------------------------------------------------------*/
static void vnoise(double x, double y, double z, double t,
		   double& ox, double& oy, double& oz, double& ot)
{
  ox = noise_template(tabindex4,x,y,z,t);
  x += 10.0;
  oy = noise_template(tabindex4,x,y,z,t);
  y += 10.0;
  oz = noise_template(tabindex4,x,y,z,t);
  z += 10.0;
  ot = noise_template(tabindex4,x,y,z,t);
};


/*----------------------------------------------------------------------
   pnoise (2D)
   Returns a value between -1 and 1 
----------------------------------------------------------------------*/
static double pnoise(double x, double y, int px, int py)
{
  xperiod = px;
  yperiod = py;
  return noise_template(ptabindex2,x,y);
};

/*----------------------------------------------------------------------
   pnoise (3D)
   Returns a value between -1 and 1 
----------------------------------------------------------------------*/
static double pnoise(double x, double y, double z, int px, int py, int pz)
{
  xperiod = px;
  yperiod = py;
  zperiod = pz;
  return noise_template(ptabindex3,x,y,z);
};

/*----------------------------------------------------------------------
   pnoise (4D)
   Returns a value between -1 and 1 
----------------------------------------------------------------------*/
static double pnoise(double x, double y, double z, double t, int px, int py, int pz, int pt)
{
  xperiod = px;
  yperiod = py;
  zperiod = pz;
  tperiod = pt;
  return noise_template(ptabindex4,x,y,z,t);
};

/*---------------------------------------------------------------------- 
   vpnoise (2D)
   Returns values between -1 and 1 
----------------------------------------------------------------------*/
static void vpnoise(double x, double y, 
		    int px, int py,
		    double& ox, double& oy)
{
  xperiod = px;
  yperiod = py;
  poffset = 0;
  ox = noise_template(ptabindex2offset,x,y);
  x += 10.0;
  poffset = 37;
  oy = noise_template(ptabindex2offset,x,y);
};

/*---------------------------------------------------------------------- 
   vpnoise (3D)
   Returns values between -1 and 1 
----------------------------------------------------------------------*/
static void vpnoise(double x, double y, double z, 
		    int px, int py, int pz,
		    double& ox, double& oy, double& oz)
{
  xperiod = px;
  yperiod = py;
  zperiod = pz;
  poffset = 0;
  ox = noise_template(ptabindex3offset,x,y,z);
  x += 10.0;
  poffset = 37;
  oy = noise_template(ptabindex3offset,x,y,z);
  y += 10.0;
  poffset = 99;
  oz = noise_template(ptabindex3offset,x,y,z);
};

/*---------------------------------------------------------------------- 
   vpnoise (4D)
   Returns values between -1 and 1 
----------------------------------------------------------------------*/
static void vpnoise(double x, double y, double z, double t,
		    int px, int py, int pz, int pt,
		    double& ox, double& oy, double& oz, double& ot)
{
  xperiod = px;
  yperiod = py;
  zperiod = pz;
  tperiod = pt;
  poffset = 0;
  ox = noise_template(ptabindex4offset,x,y,z,t);
  x += 10.0;
  poffset = 37;
  oy = noise_template(ptabindex4offset,x,y,z,t);
  y += 10.0;
  poffset = 99;
  oz = noise_template(ptabindex4offset,x,y,z,t);
  z += 10.0;
  poffset = 105;
  ot = noise_template(ptabindex4offset,x,y,z,t);
};


/*----------------------------------------------------------------------
  Fractional Brownian Motion
  (as described in Advanced RenderMan)
----------------------------------------------------------------------*/
static double fBm(double x, double y, double z, int octaves, double lacunarity, double gain)
{
  double res=0.0;
  double amp=1.0;

  for(int i=0; i<octaves; i+=1)
  {
    res += amp*noise(x,y,z);
    amp *= gain;
    x *= lacunarity;
    y *= lacunarity;
    z *= lacunarity;
  }
  return res;
}

/*----------------------------------------------------------------------
  Vector Fractional Brownian Motion
----------------------------------------------------------------------*/
static void vfBm(double x, double y, double z, 
                 int octaves, double lacunarity, double gain,
		 double& ox, double& oy, double& oz)
{
  double amp=1.0;
  double nx,ny,nz;

  ox=0.0;
  oy=0.0;
  oz=0.0;

  for(int i=0; i<octaves; i+=1)
  {
    vnoise(x,y,z, nx,ny,nz);
    ox += amp*nx;
    oy += amp*ny;
    oz += amp*nz;
    amp *= gain;
    x *= lacunarity;
    y *= lacunarity;
    z *= lacunarity;
  }
}

/*----------------------------------------------------------------------
  Turbulence
  (as described in Advanced RenderMan)
----------------------------------------------------------------------*/
static double turbulence(double x, double y, double z, int octaves, double lacunarity, double gain)
{
  double res=0.0;
  double amp=1.0;

  for(int i=0; i<octaves; i+=1)
  {
    res += amp*fabs(noise(x,y,z));
    amp *= gain;
    x *= lacunarity;
    y *= lacunarity;
    z *= lacunarity;
  }
  return res;
}

/*----------------------------------------------------------------------
  Vector Turbulence
----------------------------------------------------------------------*/
static void vturbulence(double x, double y, double z, 
			int octaves, double lacunarity, double gain,
			double& ox, double& oy, double& oz)
{
  double amp=1.0;
  double nx,ny,nz;

  ox=0.0;
  oy=0.0;
  oz=0.0;

  for(int i=0; i<octaves; i+=1)
  {
    vnoise(x,y,z, nx,ny,nz);
    ox += amp*fabs(nx);
    oy += amp*fabs(ny);
    oz += amp*fabs(nz);
    amp *= gain;
    x *= lacunarity;
    y *= lacunarity;
    z *= lacunarity;
  }
}


/*######################################################################*/

// Point structure which is used for parsing points.
struct Point
{
  double x,y,z,t;
  int n;
}; 

/*----------------------------------------------------------------------
  Converter function to be used for PyArg_ParseTuple().
  The function converts "points" (=sequences with up to 4 floating
  point values) into a Point structure.
----------------------------------------------------------------------*/
static int point_converter(PyObject* seq, Point* point)
{
  PyObject* tup = PySequence_Fast(seq, "point argument must be a sequence");
  if (tup==NULL)
  {
    point->n = 0;
    return 0;
  }

  /* Number of element */
  int tuplen = PySequence_Size(tup);
  point->n = tuplen;

  switch(tuplen)
  {
  case 1:
    if (!PyArg_ParseTuple(tup, "d;invalid point argument", &point->x)) 
      point->n=0;
    break;
  case 2:
    if (!PyArg_ParseTuple(tup, "dd;invalid point argument", &point->x, &point->y)) 
      point->n=0;
    break;
  case 3:
    if (!PyArg_ParseTuple(tup, "ddd;invalid point argument", &point->x, &point->y, &point->z)) 
      point->n=0;
    break;
  case 4:
    if (!PyArg_ParseTuple(tup, "dddd;invalid point argument", &point->x, &point->y, &point->z, &point->t)) 
      point->n=0;
    break;
  default:
    PyErr_SetString(PyExc_ValueError,"point argument can only have between 1 and 4 values");
    point->n=0;
  }
   
  Py_DECREF(tup);
  return (point->n == 0)? 0 : 1;
}

/*----------------------------------------------------------------------
  parse_seqpoint 

  Parses a point in an argument tuple. The point can be given as
  a tuple or as a tuple + float. The first argument of args has to be
  a tuple!

  args: Argument tuple (len(args) must be greater than 0)
  x,y,z,t: Output values
  Return value: Number of elements in seq or 0 (=Error)
----------------------------------------------------------------------*/
static int parse_seqpoint(PyObject* args, double& x, double& y, double& z, double& t)
{
  Point pnt;
  double f;
  int arglen;

  /* Number of arguments */
  arglen = PySequence_Size(args);

  if (arglen==1)
  {
    if (!PyArg_ParseTuple(args, "O&", point_converter, &pnt)) return 0;
    switch(pnt.n)
    {
    case 1: x=pnt.x; break;
    case 2: x=pnt.x; y=pnt.y; break;
    case 3: x=pnt.x; y=pnt.y; z=pnt.z; break;     
    case 4: x=pnt.x; y=pnt.y; z=pnt.z; t=pnt.t; break;
    }
    return pnt.n;
  }
  else
  {
    if (!PyArg_ParseTuple(args, "O&d", point_converter, &pnt, &f)) return 0;
    switch(pnt.n)
    {
    case 1: x=pnt.x; y=f; break;
    case 2: x=pnt.x; y=pnt.y; z=f; break;
    case 3: x=pnt.x; y=pnt.y; z=pnt.z; t=f; break;     
    default:
      PyErr_SetString(PyExc_TypeError,"4D vectors are not allowed when time is specified separately");
      return 0;
    }
    return pnt.n+1;
  }
}

/*---------------------------------------------------------
  Parsing the arguments to the noise functions.  

  Accepts:  - (float)
            - (float, float)
            - (float, float, float)
            - (float, float, float, float)
            - ((float))
            - ((float, float))
            - ((float, float, float))
            - ((float, float, float, float))

  The last four don't necessarily have to be tuples, but can be any
  type that supports the sequence protocol (e.g. vec3).

  The return value is the number of parsed arguments.
---------------------------------------------------------*/
static int parse_args(PyObject* args, double& x, double& y, double& z, double& t, bool& external_t)
{
  //  Point pnt;
  int arglen;
  int res=0;
  PyObject* val;
  PyObject* val2;
  PyObject* s;
  PyObject* errmsg;

  external_t = false;

  /* Number of arguments */
  arglen = PySequence_Size(args);

  switch(arglen)
  {
  /* --- 1 argument given --- */
  case 1:
    /* val: First (and only) argument */
    val = PySequence_GetItem(args,0);
    /* Is argument a sequence or single value? */
    if (PySequence_Check(val))
    {
      res = parse_seqpoint(args, x,y,z,t);
    }
    else  /* ...a single value instead of a sequence */
    {
      if (PyArg_ParseTuple(args, "d", &x)) res=1;
    }
    Py_DECREF(val);
    break;

  /* --- 2 arguments given --- */
  case 2:
    /* val: First argument */

    val = PySequence_GetItem(args,0);
    /* Is argument a sequence or single value? */
    if (PySequence_Check(val))
    {
      external_t = true;
      res = parse_seqpoint(args, x,y,z,t);
    }
    else  /* ...two floats */
    {
      if (PyArg_ParseTuple(args, "dd", &x,&y)) res=2;
    }
    Py_DECREF(val);
    break;

  /* --- 3 arguments given --- */
  case 3:
    if (PyArg_ParseTuple(args, "ddd", &x,&y,&z)) res=3;
    break;

  /* --- 4 arguments given --- */
  case 4:
    if (PyArg_ParseTuple(args, "dddd", &x,&y,&z,&t)) res=4;
    break;

  /* --- the wrong number of arguments given --- */
  default:
    val2=Py_BuildValue("i",arglen);
    s=PyString_FromString("the function takes between 1 and 4 arguments (%i given)");
    errmsg=PyString_Format(s,val2);
    PyErr_SetObject(PyExc_TypeError,errmsg);
    Py_DECREF(errmsg);
    Py_DECREF(s);
    Py_DECREF(val2);
  }

  return res;
}

// same as above, but without external_t Flag.
static inline int parse_args(PyObject* args, double& x, double& y, double& z, double& t)
{
  bool unused;
  return parse_args(args,x,y,z,t,unused);
}

/*----------------------------------------------------------------------
  Parse arguments for pnoise
----------------------------------------------------------------------*/
static int parse_args_pnoise(PyObject* args, 
			     double& x, double &y, double& z, double &t, 
			     int& ipx, int& ipy, int& ipz, int &ipt,
			     bool& external_t)
{
  //  Point pnt1;
  //  Point pnt2;
  double px=1.0,py=1.0,pz=1.0,pt=1.0;
  PyObject* seq1;
  PyObject* seq2;
  int arglen;
  int n=0;
  int m=0;

  external_t = false;

  arglen = PySequence_Size(args);
  switch(arglen)
  {
   // two arguments (seq,seq)
   case 2:
     seq1 = PySequence_GetSlice(args,0,1);
     seq2 = PySequence_GetSlice(args,1,2);
     n = parse_seqpoint(seq1,x,y,z,t);
     if (n>0) m = parse_seqpoint(seq2,px,py,pz,pt);
     Py_DECREF(seq1);
     Py_DECREF(seq2);
     break;
   // four arguments (seq,float,seq,float)
   case 4:
     external_t = true;
     seq1 = PySequence_GetSlice(args,0,2);
     seq2 = PySequence_GetSlice(args,2,4);
     n = parse_seqpoint(seq1,x,y,z,t);
     if (n>0) m = parse_seqpoint(seq2,px,py,pz,pt);
     Py_DECREF(seq1);
     Py_DECREF(seq2);
     break;
   // wrong number of arguments
   default:
     PyErr_SetString(PyExc_TypeError,"only 2 or 4 arguments allowed");
     return 0;
  }

  if ((n==0) || (m==0)) return 0;

  if (n!=m)
  {
    PyErr_SetString(PyExc_TypeError,"the period must have as many values as the point");
    return 0;
  }

  ipx=int(px);
  ipy=int(py);
  ipz=int(pz);
  ipt=int(pt);

  if ((ipx==0) || (ipy==0) || (ipz==0) || (ipt==0))
  {
    PyErr_SetString(PyExc_ValueError,"a period value is zero");
    return 0;
  }

  return n;
}

// same as above, but without external_t Flag.
static inline int parse_args_pnoise(PyObject* args,
				   double& x, double &y, double& z, double &t, 
				   int& ipx, int& ipy, int& ipz, int &ipt)
{
  bool unused;
  return parse_args_pnoise(args,x,y,z,t,ipx,ipy,ipz,ipt,unused);
}

/*######################################################################*/
/*######################################################################*/
/*######################################################################*/

/*--------------------------------------------------
  Noise-function (Perlin noise).

  Returns a value between 0 and 1.
---------------------------------------------------*/
static PyObject* noise_noise(PyObject* self, PyObject* args)
{
  double x=0.0, y=0.0, z=0.0, t=0.0;
  int n;

  n = parse_args(args, x,y,z,t);
  if (n>0)
  {
    switch(n)
    {
    case 1: 
    case 2: return Py_BuildValue("d",0.5*(noise(x,y)+1.0));
    case 3: return Py_BuildValue("d",0.5*(noise(x,y,z)+1.0));
    case 4: return Py_BuildValue("d",0.5*(noise(x,y,z,t)+1.0));
    }
    // Should never be reached (but to get rid of the warning)
    return Py_BuildValue("");
  }
  else
  {
    return NULL;
  }
}

/*--------------------------------------------------
  Signed Noise-function (Perlin noise).

  Returns a value between -1 and 1.
---------------------------------------------------*/
static PyObject* noise_snoise(PyObject* self, PyObject* args)
{
  double x=0.0, y=0.0, z=0.0, t=0.0;
  int n;

  n = parse_args(args, x,y,z,t);
  if (n>0)
  {
    return Py_BuildValue("d", noise(x,y,z));
  }
  else
  {
    return NULL;
  }
}

/*--------------------------------------------------
  Vector Noise-function (Perlin noise).
  Returns values between 0 and 1.
---------------------------------------------------*/
static PyObject* noise_vnoise(PyObject* self, PyObject* args)
{
  double x=0.0, y=0.0, z=0.0, t=0.0;
  double resx,resy,resz,rest;
  int n;
  bool external_t;

  n = parse_args(args, x,y,z,t, external_t);

  if (n>0)
  {
    switch(n)
    {
    case 1: 
      return Py_BuildValue("d",0.5*(noise(x,y)+1.0));
    case 2: 
      vnoise(x,y, resx,resy);
      resx = 0.5*(resx+1.0);
      resy = 0.5*(resy+1.0);
      break;
    case 3:
      vnoise(x,y,z, resx,resy,resz);
      resx = 0.5*(resx+1.0);
      resy = 0.5*(resy+1.0);
      resz = 0.5*(resz+1.0);
      break;
    case 4:
      vnoise(x,y,z,t, resx,resy,resz,rest);
      resx = 0.5*(resx+1.0);
      resy = 0.5*(resy+1.0);
      resz = 0.5*(resz+1.0);
      rest = 0.5*(rest+1.0);
      break;
    }

    // Was t provided separately? Then decrease output length
    if (external_t) n-=1;

    switch(n)
    {
    case 1: 
      return Py_BuildValue("d", resx);
    case 2: 
      return Py_BuildValue("(dd)", resx,resy);
    case 3:
      return PyObject_CallFunction(vec3, "ddd", resx,resy,resz);
    case 4:
      return PyObject_CallFunction(vec4, "dddd", resx,resy,resz,rest);
    // should never occur (but there'd be a warning without)
    default:
      return Py_BuildValue("");
    }
  }
  else
  {
    return NULL;
  }
}


/*--------------------------------------------------
  Unsigned vector Noise-function (Perlin noise).
  Returns values between -1 and 1.
---------------------------------------------------*/
static PyObject* noise_vsnoise(PyObject* self, PyObject* args)
{
  double x=0.0, y=0.0, z=0.0, t=0.0;
  double resx,resy,resz,rest;
  int n;  
  bool external_t;

  n = parse_args(args, x,y,z,t, external_t);

  if (n>0)
  {
    switch(n)
    {
    case 1: 
      return Py_BuildValue("d",noise(x,y));
    case 2: 
      vnoise(x,y, resx,resy);
      break;
    case 3:
      vnoise(x,y,z, resx,resy,resz);
      break;
    case 4:
      vnoise(x,y,z,t, resx,resy,resz,rest);
      break;
    }

    // Was t provided separately? Then decrease output length
    if (external_t) n-=1;
    switch(n)
    {
    case 1: 
      return Py_BuildValue("d",resx);
    case 2: 
      return Py_BuildValue("(dd)", resx,resy);
    case 3:
      return PyObject_CallFunction(vec3, "ddd", resx,resy,resz);
    case 4:
      return PyObject_CallFunction(vec4, "dddd", resx,resy,resz,rest);
    // should never occur (but there'd be a warning without)
    default:
      return Py_BuildValue("");
    }
  }
  else
  {
    return NULL;
  }
}

/*--------------------------------------------------
  Periodic Noise.
  Returns a value between 0 and 1.
---------------------------------------------------*/
static PyObject* noise_pnoise(PyObject* self, PyObject* args)
{
  double x=0.0, y=0.0, z=0.0, t=0.0;
  int ipx,ipy,ipz,ipt;
  int n=0;

  /*  Point p;

  if (!PyArg_ParseTuple(args, "O&", point_converter, &p)) return NULL;

  printf("p.x = %f\n",p.x);
  printf("p.y = %f\n",p.y);
  printf("p.z = %f\n",p.z);
  printf("p.t = %f\n",p.t);
  printf("p.n = %d\n",p.n);
  return Py_BuildValue("");*/
  

  n = parse_args_pnoise(args,x,y,z,t,ipx,ipy,ipz,ipt);

  if (n==0) return NULL;
  
  if ((n==1) || (n==2))
    return Py_BuildValue("d",0.5*(pnoise(x,y,ipx,ipy)+1.0));
  if (n==3)
    return Py_BuildValue("d",0.5*(pnoise(x,y,z,ipx,ipy,ipz)+1.0));
  if (n==4)
    return Py_BuildValue("d",0.5*(pnoise(x,y,z,t,ipx,ipy,ipz,ipt)+1.0));

  // Should never occur
  return Py_BuildValue("");
}

/*--------------------------------------------------
  Signed Periodic Noise.
  Returns a value between -1 and 1.
---------------------------------------------------*/
static PyObject* noise_spnoise(PyObject* self, PyObject* args)
{
  double x=0.0, y=0.0, z=0.0, t=0.0;
  int ipx,ipy,ipz,ipt;
  int n=0;

  n = parse_args_pnoise(args,x,y,z,t,ipx,ipy,ipz,ipt);

  if (n==0) return NULL;
  
  if ((n==1) || (n==2))
    return Py_BuildValue("d",pnoise(x,y,ipx,ipy));
  if (n==3)
    return Py_BuildValue("d",pnoise(x,y,z,ipx,ipy,ipz));
  if (n==4)
    return Py_BuildValue("d",pnoise(x,y,z,t,ipx,ipy,ipz,ipt));

  // Should never occur
  return Py_BuildValue("");
}

/*--------------------------------------------------
  Vector Periodic Noise.
  Returns a value between 0 and 1.
---------------------------------------------------*/
static PyObject* noise_vpnoise(PyObject* self, PyObject* args)
{
  double x=0.0, y=0.0, z=0.0, t=0.0;
  double resx,resy,resz,rest;
  int ipx,ipy,ipz,ipt;
  int n=0;
  bool external_t;

  n = parse_args_pnoise(args,x,y,z,t,ipx,ipy,ipz,ipt,external_t);

  if (n==0) return NULL;

  switch(n)
  {
  case 1: 
    return Py_BuildValue("d",0.5*(pnoise(x,y,ipx,ipy)+1.0));
  case 2: 
    vpnoise(x,y, ipx,ipy, resx,resy);
    resx = 0.5*(resx+1.0);
    resy = 0.5*(resy+1.0);
    break;
  case 3:
    vpnoise(x,y,z, ipx,ipy,ipz, resx,resy,resz);
    resx = 0.5*(resx+1.0);
    resy = 0.5*(resy+1.0);
    resz = 0.5*(resz+1.0);
    break;
  case 4:
    vpnoise(x,y,z,t, ipx,ipy,ipz,ipt, resx,resy,resz,rest);
    resx = 0.5*(resx+1.0);
    resy = 0.5*(resy+1.0);
    resz = 0.5*(resz+1.0);
    rest = 0.5*(rest+1.0);
    break;
  }

  // Was t provided separately? Then decrease output length
  if (external_t) n-=1;

  switch(n)
  {
  case 1: 
    return Py_BuildValue("d",resx);
  case 2: 
    return Py_BuildValue("(dd)", resx,resy);
  case 3:
    return PyObject_CallFunction(vec3, "ddd", resx,resy,resz);
  case 4:
    return PyObject_CallFunction(vec4, "dddd", resx,resy,resz,rest);
  }

  // Should never occur
  return Py_BuildValue("");
}

/*--------------------------------------------------
  Vector Signed Periodic Noise.
  Returns a value between -1 and 1.
---------------------------------------------------*/
static PyObject* noise_vspnoise(PyObject* self, PyObject* args)
{
  double x=0.0, y=0.0, z=0.0, t=0.0;
  double resx,resy,resz,rest;
  int ipx,ipy,ipz,ipt;
  int n=0;
  bool external_t;

  n = parse_args_pnoise(args,x,y,z,t,ipx,ipy,ipz,ipt,external_t);

  if (n==0) return NULL;

  switch(n)
  {
  case 1: 
    return Py_BuildValue("d",0.5*(pnoise(x,y,ipx,ipy)+1.0));
  case 2: 
    vpnoise(x,y, ipx,ipy, resx,resy);
    break;
  case 3:
    vpnoise(x,y,z, ipx,ipy,ipz, resx,resy,resz);
    break;
  case 4:
    vpnoise(x,y,z,t, ipx,ipy,ipz,ipt, resx,resy,resz,rest);
    break;
  }

  // Was t provided separately? Then decrease output length
  if (external_t) n-=1;

  switch(n)
  {
  case 1: 
    return Py_BuildValue("d",resx);
  case 2: 
    return Py_BuildValue("(dd)", resx,resy);
  case 3:
    return PyObject_CallFunction(vec3, "ddd", resx,resy,resz);
  case 4:
    return PyObject_CallFunction(vec4, "dddd", resx,resy,resz,rest);
  }

  // Should never occur
  return Py_BuildValue("");
}


/*--------------------------------------------------
  Cellnoise
  Returns a value between 0 and 1.
---------------------------------------------------*/
static PyObject* noise_cellnoise(PyObject* self, PyObject* args)
{
  double x=0.0, y=0.0, z=0.0, t=0.0;
  int n;

  n = parse_args(args, x,y,z,t);
  if (n>0)
  {
    return Py_BuildValue("d",cellnoise(x,y,z,t));  
  }
  else
  {
    return NULL;
  }
}

/*--------------------------------------------------
  Signed cellnoise
  Returns a value between -1 and 1.
---------------------------------------------------*/
static PyObject* noise_scellnoise(PyObject* self, PyObject* args)
{
  double x=0.0, y=0.0, z=0.0, t=0.0;
  int n;

  n = parse_args(args, x,y,z,t);
  if (n>0)
  {
    return Py_BuildValue("d",2.0*cellnoise(x,y,z,t)-1.0);
  }
  else
  {
    return NULL;
  }
}

/*--------------------------------------------------
  Vector cellnoise.
  Returns values between 0 and 1.
---------------------------------------------------*/
static PyObject* noise_vcellnoise(PyObject* self, PyObject* args)
{
  double x=0.0, y=0.0, z=0.0, t=0.0;
  double resx,resy,resz,rest;
  int n;
  bool external_t;

  n = parse_args(args, x,y,z,t, external_t);

  if (n>0)
  {
    switch(n)
    {
    case 1: 
      return Py_BuildValue("d",cellnoise(x,y,z,t));
    case 2: 
      vcellnoise(x,y, resx,resy);
      break;
    case 3:
      vcellnoise(x,y,z, resx,resy,resz);
      break;
    case 4:
      vcellnoise(x,y,z,t, resx,resy,resz,rest);
      break;
    }

    // Was t provided separately? Then decrease output length
    if (external_t) n-=1;

    switch(n)
    {
    case 1: 
      return Py_BuildValue("d", resx);
    case 2: 
      return Py_BuildValue("(dd)", resx,resy);
    case 3:
      return PyObject_CallFunction(vec3, "ddd", resx,resy,resz);
    case 4:
      return PyObject_CallFunction(vec4, "dddd", resx,resy,resz,rest);
    // should never occur (but there'd be a warning without)
    default:
      return Py_BuildValue("");
    }
  }
  else
  {
    return NULL;
  }
}

/*--------------------------------------------------
  Vector scellnoise.
  Returns values between -1 and 1.
---------------------------------------------------*/
static PyObject* noise_vscellnoise(PyObject* self, PyObject* args)
{
  double x=0.0, y=0.0, z=0.0, t=0.0;
  double resx,resy,resz,rest;
  int n;
  bool external_t;

  n = parse_args(args, x,y,z,t, external_t);

  if (n>0)
  {
    switch(n)
    {
    case 1: 
      return Py_BuildValue("d",2.0*cellnoise(x,y,z,t)-1.0);
    case 2: 
      vcellnoise(x,y, resx,resy);
      resx = 2.0*resx-1.0;
      resy = 2.0*resy-1.0;
      break;
    case 3:
      vcellnoise(x,y,z, resx,resy,resz);
      resx = 2.0*resx-1.0;
      resy = 2.0*resy-1.0;
      resz = 2.0*resz-1.0;
      break;
    case 4:
      vcellnoise(x,y,z,t, resx,resy,resz,rest);
      resx = 2.0*resx-1.0;
      resy = 2.0*resy-1.0;
      resz = 2.0*resz-1.0;
      rest = 2.0*rest-1.0;
      break;
    }

    // Was t provided separately? Then decrease output length
    if (external_t) n-=1;

    switch(n)
    {
    case 1: 
      return Py_BuildValue("d", resx);
    case 2: 
      return Py_BuildValue("(dd)", resx,resy);
    case 3:
      return PyObject_CallFunction(vec3, "ddd", resx,resy,resz);
    case 4:
      return PyObject_CallFunction(vec4, "dddd", resx,resy,resz,rest);
    // should never occur (but there'd be a warning without)
    default:
      return Py_BuildValue("");
    }
  }
  else
  {
    return NULL;
  }
}

/*--------------------------------------------------
  fBm
---------------------------------------------------*/
static PyObject* noise_fBm(PyObject* self, PyObject* args)
{
  double x=0.0, y=0.0, z=0.0, lacunarity=2.0, gain=0.5;
  int octaves;

  if (!PyArg_ParseTuple(args,"(dd)i|dd", &x,&y,&octaves,&lacunarity,&gain))
  {
    PyErr_Clear();
    if (!PyArg_ParseTuple(args,"(ddd)i|dd", &x,&y,&z,&octaves,&lacunarity,&gain))
	return NULL;
  }

  return Py_BuildValue("d",0.5*(fBm(x,y,z,octaves,lacunarity,gain)+1.0));  
}

/*--------------------------------------------------
  vfBm
---------------------------------------------------*/
static PyObject* noise_vfBm(PyObject* self, PyObject* args)
{
  double x=0.0, y=0.0, z=0.0, lacunarity=2.0, gain=0.5;
  int octaves;
  double resx,resy,resz;

  if (PyArg_ParseTuple(args,"(dd)i|dd", &x,&y,&octaves,&lacunarity,&gain))
  {
    vfBm(x,y,z,octaves,lacunarity,gain,resx,resy,resz);
    resx=0.5*(resx+1.0);
    resy=0.5*(resy+1.0);
    return Py_BuildValue("(dd)",resx,resy);  
  }
  PyErr_Clear();
  if (PyArg_ParseTuple(args,"(ddd)i|dd", &x,&y,&z,&octaves,&lacunarity,&gain))
  {
    vfBm(x,y,z,octaves,lacunarity,gain,resx,resy,resz);
    resx=0.5*(resx+1.0);
    resy=0.5*(resy+1.0);
    resz=0.5*(resz+1.0);
    return PyObject_CallFunction(vec3, "ddd", resx,resy,resz);
  }

  return NULL;
}

/*--------------------------------------------------
  turbulence
---------------------------------------------------*/
static PyObject* noise_turbulence(PyObject* self, PyObject* args)
{
  double x=0.0, y=0.0, z=0.0, lacunarity=2.0, gain=0.5;
  int octaves;

  if (!PyArg_ParseTuple(args,"(dd)i|dd", &x,&y,&octaves,&lacunarity,&gain))
  {
    PyErr_Clear();
    if (!PyArg_ParseTuple(args,"(ddd)i|dd", &x,&y,&z,&octaves,&lacunarity,&gain))
	return NULL;
  }

  return Py_BuildValue("d",turbulence(x,y,z,octaves,lacunarity,gain));  
}

/*--------------------------------------------------
  Vector turbulence
---------------------------------------------------*/
static PyObject* noise_vturbulence(PyObject* self, PyObject* args)
{
  double x=0.0, y=0.0, z=0.0, lacunarity=2.0, gain=0.5;
  int octaves;
  double resx,resy,resz;

  if (PyArg_ParseTuple(args,"(dd)i|dd", &x,&y,&octaves,&lacunarity,&gain))
  {
    vturbulence(x,y,z,octaves,lacunarity,gain,resx,resy,resz);
    return Py_BuildValue("(dd)",resx,resy);
  }
  PyErr_Clear();
  if (PyArg_ParseTuple(args,"(ddd)i|dd", &x,&y,&z,&octaves,&lacunarity,&gain))
  {
    vturbulence(x,y,z,octaves,lacunarity,gain,resx,resy,resz);
    return PyObject_CallFunction(vec3, "ddd", resx,resy,resz);
  }

  return NULL;
}


/*######################################################################*/

// Function table
static PyMethodDef NoiseMethods[] = {
  {"noise", noise_noise, METH_VARARGS, 
   "Returns a noise value (Perlin) in the range from 0 to 1.\n\nThe arguments can be up to 4 floating point values or a sequence\nwith up to 4 floating point values. Time can be specified separately."},
  {"snoise", noise_snoise, METH_VARARGS,
   "Returns a signed noise value (Perlin) in the range from -1 to 1.\n\nThe arguments can be up to 4 floating point values or a sequence\nwith up to 4 floating point values. Time can be specified separately."},
  {"vnoise", noise_vnoise, METH_VARARGS,
   "Returns a noise vector with components in the range from 0 to 1.\n\nThe arguments can be up to 4 floating point values or a sequence\nwith up to 4 floating point values. Time can be specified separately.\nThe output has the same dimension than the input."},
  {"vsnoise", noise_vsnoise, METH_VARARGS,
   "Returns a signed noise vector with components in the range from -1 to 1.\n\nThe arguments can be up to 4 floating point values or a sequence\nwith up to 4 floating point values. Time can be specified separately.\nThe output has the same dimension than the input."},
  {"pnoise", noise_pnoise, METH_VARARGS,
   "Returns a periodic noise value in the range from 0 to 1.\n\nAs arguments you have to provide the point and the period, each can be\nup to 4 floating point values or a sequence with up to 4 floating point\nvalues. Time can be specified separately."},
  {"spnoise", noise_spnoise, METH_VARARGS,
   "Returns a signed periodic noise value in the range from -1 to 1.\n\nAs arguments you have to provide the point and the period, each can be\nup to 4 floating point values or a sequence with up to 4 floating point\nvalues. Time can be specified separately."},
  {"vpnoise", noise_vpnoise, METH_VARARGS,
   "Returns a periodic noise vector with components in the range from 0 to 1.\n\nAs arguments you have to provide the point and the period, each can be\nup to 4 floating point values or a sequence with up to 4 floating point\nvalues. Time can be specified separately."},
  {"vspnoise", noise_vspnoise, METH_VARARGS,
   "Returns a signed periodic noise vector with components in the range from -1 to 1.\n\nAs arguments you have to provide the point and the period, each can be\nup to 4 floating point values or a sequence with up to 4 floating point\nvalues. Time can be specified separately."},
  {"cellnoise", noise_cellnoise, METH_VARARGS,
   "Returns a cellnoise value in the range from 0 to 1.\n\nThe return value is constant between integer lattice points.\nThe arguments can be up to 4 floating point values or a sequence\nwith up to 4 floating point values. Time can be specified separately."},
  {"scellnoise", noise_scellnoise, METH_VARARGS,
   "Returns a signed cellnoise value in the range from -1 to 1.\n\nThe return value is constant between integer lattice points.\nThe arguments can be up to 4 floating point values or a sequence\nwith up to 4 floating point values. Time can be specified separately."},
  {"vcellnoise", noise_vcellnoise, METH_VARARGS,
   "Returns a cellnoise vector with components in the range from 0 to 1.\n\nThe return value is constant between integer lattice points.\nThe arguments can be up to 4 floating point values or a sequence\nwith up to 4 floating point values. Time can be specified separately."},
  {"vscellnoise", noise_vscellnoise, METH_VARARGS,
   "Returns a signed cellnoise vector with components in the range from -1 to 1.\n\nThe return value is constant between integer lattice points.\nThe arguments can be up to 4 floating point values or a sequence\nwith up to 4 floating point values. Time can be specified separately."},
  {"fBm", noise_fBm, METH_VARARGS,
   "Fractional Brownian motion.\n\nfBm(point, octaves [,lacunarity [,gain]]). The default values for\nlacunarity and gain are 2.0 and 0.5."},
  {"vfBm", noise_vfBm, METH_VARARGS,
   "Vector version of fractional Brownian motion.\n\nvfBm(point, octaves [,lacunarity [,gain]]). The default values for\nlacunarity and gain are 2.0 and 0.5."},
  {"turbulence", noise_turbulence, METH_VARARGS,
   "Turbulence.\n\nturbulence(point, octaves [,lacunarity [,gain]]). The default values for\nlacunarity and gain are 2.0 and 0.5."},
  {"vturbulence", noise_vturbulence, METH_VARARGS,
   "Vector version of turbulence.\n\nvturbulence(point, octaves [,lacunarity [,gain]]). The default values for\nlacunarity and gain are 2.0 and 0.5."},
  {NULL, NULL}
};

/*---------------------------------------------------------------------- 
  Module initialization 
----------------------------------------------------------------------*/
extern "C" void initnoise()
{
  Py_InitModule("noise", NoiseMethods);

  // Import cgtypes and get the vec3 and vec4 class
  PyObject* mod = PyImport_ImportModule("cgtypes");
  if (mod!=NULL)
  {
    vec3 = PyObject_GetAttrString(mod,"vec3");
    if (vec3==NULL)
    {
      Py_DECREF(mod);
      return;
    }
    vec4 = PyObject_GetAttrString(mod,"vec4");
    if (vec4==NULL)
    {
      Py_DECREF(vec3);
      Py_DECREF(mod);
      return;
    }
    Py_DECREF(mod);
  }
}
