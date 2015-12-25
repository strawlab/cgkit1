# -*- coding: utf-8 -*-
######################################################################
# cgtypes - vec3, vec4, mat3, mat4, quat
#
# Copyright (C) 2002, Matthias Baas (baas@ira.uka.de)
#
# http://cgkit.sourceforge.net
#
# You may distribute under the terms of the BSD license, as
# specified in the file license.txt.
####################################################################

# Some math functions
cdef extern from "math.h":
    double sqrt(double x)
    double fabs(double x)
    double sin(double x)
    double cos(double x)
    double tan(double x)
    double acos(double x)
    double exp(double x)
    double log(double x)
    double atan2(double y, double x)

# Threshold for comparisons
cdef double eps
eps = 1E-16

# Forward declarations
cdef class vec3
cdef class vec4
cdef class mat3
cdef class mat4

######################################################################

# _fmod: A modulo function that also works with negative a
cdef double _fmod(double a, double b):
    cdef int n
    n = <int>(a/b)
    a = a-n*b
    if a<0:
        a=a+b
    return a

######################################################################

# setEpsilon
def setEpsilon(e):
    """Sets a new epsilon threshold value.

    Sets a new epsilon threshold value and returns the previously set
    value. Two values are considered to be equal if their absolute
    difference is less than or equal to epsilon."""
    global eps
    res = eps
    eps = e
    return res

# getEpsilon
def getEpsilon():
    """Return the epsilon threshold which is used for doing comparisons."""
    return eps


######################################################################


# vec3iter
cdef class vec3iter:
    """vec3 iterator."""
    cdef int idx
    cdef object v

    def __cinit__(self, v):
        self.idx = 0
        self.v   = v

    def __iter__(self):
        return self

    def __next__(self):
        if self.idx==0:
            self.idx=1
            return self.v.x
        elif self.idx==1:
            self.idx=2
            return self.v.y
        elif self.idx==2:
            self.idx=3
            return self.v.z
        else:
            raise StopIteration


# vec3
cdef class vec3:
    """Three-dimensional vector.

    This class can be used to represent points, vectors, normals
    or even colors. The usual vector operations are available.

    Comparisons are done within an epsilon environment.
    """

    cdef double x,y,z

    def __cinit__(self, *args):
        """Constructor.

        There are several possibilities how to initialize a vector:

        v = vec3()       -> v = <0,0,0>
        v = vec3(a)      -> v = <a,a,a>
        v = vec3(x,y)    -> v = <x,y,0>
        v = vec3(x,y,z)  -> v = <x,y,z>
        w = vec3(v)      -> w = <x,y,z>

        Note that specifying just one value sets all three components to
        that value (except when that single value is a another vec3, then
        that vector is copied).

        Additionally you can wrap those values in a list or a tuple or
        specify them as a string:

        v = vec3([1,2,3]) -> v = <1,2,3>
        v = vec3("4,5")   -> v = <4,5,0>
        """
        cdef int arglen, seqlen
        cdef vec3 v
        arglen = len(args)

        if arglen==0:
            self.x = 0.0
            self.y = 0.0
            self.z = 0.0
        elif arglen==1:
            T = type(args[0])
            # scalar
            if T==float or T==int or T==long:
                self.x = args[0]
                self.y = self.x
                self.z = self.x
            # vec3
            elif T==vec3:
                v = args[0]
                self.x = v.x
                self.y = v.y
                self.z = v.z
            # String
            elif T==str or T==unicode:
                s=args[0].replace(","," ").replace("  "," ").strip().split(" ")
                if s==[""]:
                    s=[]
                for i in range(len(s)):
                    try:
                        s[i] = float(s[i])
                    except:
                        raise ValueError,"vec3() arg is a malformed string"
                v = vec3(s)
                self.x = v.x
                self.y = v.y
                self.z = v.z
            # Sequence of floats
            else:
                seq = args[0]
                try:
                    seqlen = len(seq)
                except:
                    raise TypeError,"vec3() arg can't be converted to vec3"

                if seqlen==3:
                    self.x = seq[0]
                    self.y = seq[1]
                    self.z = seq[2]
                elif seqlen==2:
                    self.x = seq[0]
                    self.y = seq[1]
                    self.z = 0.0
                elif seqlen==1:
                    self.x = seq[0]
                    self.y = self.x
                    self.z = self.x
                elif seqlen==0:
                    self.x = 0.0
                    self.y = 0.0
                    self.z = 0.0
                else:
                    raise TypeError, "vec3() takes at most 3 arguments"

        elif arglen==2:
            self.x = args[0]
            self.y = args[1]
            self.z = 0.0

        elif arglen==3:
            self.x = args[0]
            self.y = args[1]
            self.z = args[2]

        else:
            raise TypeError, "vec3() takes at most 3 arguments"

    def __reduce__(self):
        return (_vec3construct, (self.x, self.y, self.z))

    def __repr__(self):
        return 'vec3(%s, %s, %s)'%(repr(self.x),repr(self.y),repr(self.z))

    def __str__(self):
        return "(%1.4f, %1.4f, %1.4f)"%(self.x, self.y, self.z)

    def __iter__(self):
        return vec3iter(self)

    def __richcmp__(a,b,op):
        cdef vec3 va, vb
        cdef int cop

        cop = op

        ta = type(a)
        tb = type(b)
        if (ta!=vec3 or tb!=vec3):
            if cop==3:
                return 1
            else:
                return 0

        va = a
        vb = b

        # <
        if cop==0:
            return (va.x<vb.x and va.y<vb.y and va.z<vb.z)
        # <=
        elif cop==1:
            return (va.x-eps<=vb.x and va.y-eps<=vb.y and va.z-eps<=vb.z)
        # ==
        elif cop==2:
#            return (va.x==vb.x and va.y==vb.y and va.z==vb.z)
            return (fabs(va.x-vb.x)<=eps and fabs(va.y-vb.y)<=eps and fabs(va.z-vb.z)<=eps)
        # !=
        elif cop==3:
#            return (va.x!=vb.x or va.y!=vb.y or va.z!=vb.z)
            return (fabs(va.x-vb.x)>eps or fabs(va.y-vb.y)>eps and fabs(va.z-vb.z)>eps)
        # >
        elif cop==4:
            return (va.x>vb.x and va.y>vb.y and va.z>vb.z)
        # >=
        elif cop==5:
            return (va.x+eps>=vb.x and va.y+eps>=vb.y and va.z+eps>=vb.z)
        # sonst (Fehler)
        else:
            raise ValueError,"internal error: illegal rich comparison number"


    def __add__(vec3 a, vec3 b):
        """Vector addition.

        >>> a=vec3(1.0, 0.5, -1.8)
        >>> b=vec3(-0.3, 0.75, 0.5)
        >>> print a+b
        (0.7000, 1.2500, -1.3000)
        """
        cdef vec3 res
        res = vec3()
        res.x = a.x+b.x
        res.y = a.y+b.y
        res.z = a.z+b.z
        return res

    def __sub__(vec3 a, vec3 b):
        """Vector subtraction.

        >>> a=vec3(1.0, 0.5, -1.8)
        >>> b=vec3(-0.3, 0.75, 0.5)
        >>> print a-b
        (1.3000, -0.2500, -2.3000)
        """
        cdef vec3 res
        res = vec3()
        res.x = a.x-b.x
        res.y = a.y-b.y
        res.z = a.z-b.z
        return res

    def __mul__(a, b):
        """Multiplication with a scalar or dot product.

        >>> a=vec3(1.0, 0.5, -1.8)
        >>> b=vec3(-0.3, 0.75, 0.5)
        >>> print a*2.0
        (2.0000, 1.0000, -3.6000)
        >>> print 2.0*a
        (2.0000, 1.0000, -3.6000)
        >>> print a*b
        -0.825
        """
        cdef vec3 res, va, vb
        cdef double r

        ta = type(a)
        tb = type(b)

        if ta==vec3:
            va = a
            if tb==vec3:
                # vec3*vec3   (dot product a*b)
                vb = b
                return va.x*vb.x + va.y*vb.y + va.z*vb.z
            elif tb==float or tb==int or tb==long:
                # vec3*scalar
                res = vec3()
                r   = b
                res.x = va.x*r
                res.y = va.y*r
                res.z = va.z*r
                return res
            elif tb==mat3:
                # vec3*mat3
                return b.__rmul__(a)
            elif tb==mat4:
                # vec3*mat4
                return b.__rmul__(a)

        elif ta==float or ta==int or ta==long:
            if tb==vec3:
                # scalar*vec3
                res = vec3()
                vb  = b
                r   = a
                res.x = r*vb.x
                res.y = r*vb.y
                res.z = r*vb.z
                return res

        raise TypeError, "unsupported operand type for *"

    def __div__(vec3 a, b):
        """Division by scalar

        >>> a=vec3(1.0, 0.5, -1.8)
        >>> print a/2.0
        (0.5000, 0.2500, -0.9000)
        """
        cdef vec3 res
        cdef double r

        tb = type(b)

        if tb==float or tb==int or tb==long:
            # vec3/scalar
            res = vec3()
            r = b
            if fabs(r)<=eps:
                raise ZeroDivisionError,"vec3 division"
            res.x = a.x/r
            res.y = a.y/r
            res.z = a.z/r
            return res
        else:
            raise TypeError, "unsupported operand type for /"

    def __mod__(vec3 self, b):
        """Modulo (component wise).

        >>> a=vec3(3.0, 2.5, -1.8)
        >>> print a%2.0
        (1.0000, 0.5000, 0.2000)
        """
        cdef vec3 res, vb
        cdef double r

        tb = type(b)

        if tb==float or tb==int or tb==long:
            # vec3%scalar
            res = vec3()
            r = b
            if fabs(r)<=eps:
                raise ZeroDivisionError,"vec3 modulo"
            res.x = _fmod(self.x,r)
            res.y = _fmod(self.y,r)
            res.z = _fmod(self.z,r)
            return res
        elif tb==vec3:
            # vec3%vec3
            res = vec3()
            vb = b
            if fabs(vb.x)<=eps or fabs(vb.y)<=eps or fabs(vb.z)<=eps:
                raise ZeroDivisionError,"vec3 modulo"
            res.x = _fmod(self.x, vb.x)
            res.y = _fmod(self.y, vb.y)
            res.z = _fmod(self.z, vb.z)
            return res
        else:
            raise TypeError, "unsupported operand type for %"

    def __neg__(self):
        """Negation

        >>> a=vec3(3.0, 2.5, -1.8)
        >>> print -a
        (-3.0000, -2.5000, 1.8000)
        """
        cdef vec3 res
        res = vec3()
        res.x = -self.x
        res.y = -self.y
        res.z = -self.z
        return res

    def __pos__(self):
        """
        >>> a=vec3(3.0, 2.5, -1.8)
        >>> print +a
        (3.0000, 2.5000, -1.8000)
        """
        cdef vec3 res
        res = vec3()
        res.x = self.x
        res.y = self.y
        res.z = self.z
        return res


    def __iadd__(self, vec3 other):
        """Inline vector addition.

        >>> a=vec3(1.0, 0.5, -1.8)
        >>> b=vec3(-0.3, 0.75, 0.5)
        >>> a+=b
        >>> print a
        (0.7000, 1.2500, -1.3000)
        """
        self.x=self.x+other.x
        self.y=self.y+other.y
        self.z=self.z+other.z
        return self

    def __isub__(self, vec3 other):
        """Inline vector subtraction.

        >>> a=vec3(1.0, 0.5, -1.8)
        >>> b=vec3(-0.3, 0.75, 0.5)
        >>> a-=b
        >>> print a
        (1.3000, -0.2500, -2.3000)
        """
        self.x=self.x-other.x
        self.y=self.y-other.y
        self.z=self.z-other.z
        return self

    def __imul__(self, other):
        """Inline multiplication (only with scalar)

        >>> a=vec3(1.0, 0.5, -1.8)
        >>> a*=2.0
        >>> print a
        (2.0000, 1.0000, -3.6000)
        """
        cdef double r
        tb = type(other)
        if tb==float or tb==int or tb==long:
            r=other
            self.x=self.x*r
            self.y=self.y*r
            self.z=self.z*r
            return self
        else:
            raise TypeError, "unsupported operand type for *="

    def __idiv__(self, other):
        """Inline division with scalar

        >>> a=vec3(1.0, 0.5, -1.8)
        >>> a/=2.0
        >>> print a
        (0.5000, 0.2500, -0.9000)
        """
        cdef double r
        tb = type(other)
        if tb==float or tb==int or tb==long:
            r=other
            if fabs(r)<=eps:
                raise ZeroDivisionError,"vec3 division"
            self.x=self.x/r
            self.y=self.y/r
            self.z=self.z/r
            return self
        else:
            raise TypeError, "unsupported operand type for /="

    def __imod__(self, b):
        """Inline modulo.

        >>> a=vec3(3.0, 2.5, -1.8)
        >>> a%=2.0
        >>> print a
        (1.0000, 0.5000, 0.2000)
        """
        cdef vec3 vb
        cdef double r

        tb = type(b)

        if tb==float or tb==int or tb==long:
            # vec3%=scalar
            r = b
            if fabs(r)<=eps:
                raise ZeroDivisionError,"vec3 modulo"
            self.x = _fmod(self.x,r)
            self.y = _fmod(self.y,r)
            self.z = _fmod(self.z,r)
            return self
        elif tb==vec3:
            # vec3%=vec3
            vb = b
            if fabs(vb.x)<=eps or fabs(vb.y)<=eps or fabs(vb.z)<=eps:
                raise ZeroDivisionError,"vec3 modulo"
            self.x = _fmod(self.x, vb.x)
            self.y = _fmod(self.y, vb.y)
            self.z = _fmod(self.z, vb.z)
            return self
        else:
            raise TypeError, "unsupported operand type for %="


    def __abs__(self):
        """Return the length of the vector.

        abs(v) is equivalent to v.length().

        >>> a=vec3(1.0, 0.5, -1.8)
        >>> print abs(a)
        2.11896201004
        """
        cdef double a
        a = self.x*self.x + self.y*self.y + self.z*self.z
        return sqrt(a)

    def __getattr__(self, name):
        if name=="x":
            return self.x
        elif name=="y":
            return self.y
        elif name=="z":
            return self.z
        else:
            raise AttributeError,"vec3 has no attribute '"+name+"'"

    def __setattr__(self, name, val):
        if name=="x":
            self.x = val
        elif name=="y":
            self.y = val
        elif name=="z":
            self.z = val
        else:
            raise AttributeError,"vec3 has no writable attribute '"+name+"'"

    def __len__(self):
        """Length of the sequence (always 3)"""
        return 3

    def __getitem__(vec3 self, key):
        """Return a component by index (0-based)

        >>> a=vec3(1.0, 0.5, -1.8)
        >>> print a[0]
        1.0
        >>> print a[1]
        0.5
        >>> print a[2]
        -1.8
        """
        cdef int k

        T=type(key)
        if T!=int and T!=long:
            raise TypeError, "index must be integer"

        k = key
        if   k==0: return self.x
        elif k==1: return self.y
        elif k==2: return self.z
        else:
            raise IndexError,"index out of range"

    def __setitem__(self, key, value):
        """Set a component by index (0-based)

        >>> a=vec3()
        >>> a[0]=1.5; a[1]=0.7; a[2]=-0.3
        >>> print a
        (1.5000, 0.7000, -0.3000)
        """
        cdef int k

        T=type(key)
        if T!=int and T!=long:
            raise TypeError, "index must be integer"

        k = key
        if   k==0: self.x = value
        elif k==1: self.y = value
        elif k==2: self.z = value
        else:
            raise IndexError,"index out of range"

    def cross(self, vec3 other):
        """Cross product.

        >>> a=vec3(1.0, 0.5, -1.8)
        >>> b=vec3(-0.3, 0.75, 0.5)
        >>> c=a.cross(b)
        >>> print c
        (1.6000, 0.0400, 0.9000)
        """
        cdef vec3 res

        if type(other)==vec3:
            res = vec3()
            res.x = self.y*other.z-self.z*other.y
            res.y = self.z*other.x-self.x*other.z
            res.z = self.x*other.y-self.y*other.x
            return res
        else:
            raise TypeError, "unsupported operand type for cross()"

    def length(self):
        """Return the length of the vector.

        v.length() is equivalent to abs(v).

        >>> a=vec3(1.0, 0.5, -1.8)
        >>> print a.length()
        2.11896201004
        """
        cdef double a
        a = self.x*self.x + self.y*self.y + self.z*self.z
        return sqrt(a)

    def normalize(self):
        """Return normalized vector.

        >>> a=vec3(1.0, 0.5, -1.8)
        >>> print a.normalize()
        (0.4719, 0.2360, -0.8495)
        """
        cdef vec3 res
        cdef double nlen

        nlen = sqrt(self.x*self.x + self.y*self.y + self.z*self.z)
        if nlen<=eps:
            raise ZeroDivisionError,"null vector can't be normalized"
        nlen = 1.0/nlen
        res = vec3()
        res.x = self.x*nlen
        res.y = self.y*nlen
        res.z = self.z*nlen
        return res

    def angle(self, vec3 other):
        """Return angle (in radians) between self and other.

        >>> a=vec3(1.0, 0.5, -1.8)
        >>> b=vec3(-0.3, 0.75, 0.5)
        >>> print a.angle(b)
        1.99306755584
        """

        cdef double a, ls, lo

        ls = sqrt(self.x*self.x + self.y*self.y + self.z*self.z)
        lo = sqrt(other.x*other.x + other.y*other.y + other.z*other.z)
        a  = ls*lo
        if a<=eps:
            raise ZeroDivisionError,"angle(): null vector can't be normalized"
        a  = (self.x*other.x + self.y*other.y + self.z*other.z)/a
        if a>1.0:
            a=1.0
        if a<-1.0:
            a=-1.0
        return acos(a)
#        if isinstance(other, vec3):
#            return math.acos((self*other) / (abs(self)*abs(other)))
#        else:
#            raise TypeError, "unsupported operand type for angle()"

    def reflect(self, vec3 N):
        """Return the reflection vector.

        N is the surface normal which has to be of unit length.

        >>> a=vec3(1.0, -1.0, 0.0)
        >>> print a.reflect(vec3(0,1,0)
        (1.0000, 1.0000, 0.0000)
        """
        cdef vec3 res
        cdef double a

        res = vec3()
        a = 2.0*(self.x*N.x + self.y*N.y + self.z*N.z)
        res.x = self.x - a*N.x
        res.y = self.y - a*N.y
        res.z = self.z - a*N.z
        return res

    def refract(self, vec3 N, eta):
        """Return the transmitted vector.

        N is the surface normal which has to be of unit length.
        eta is the relative index of refraction. If the returned
        vector is zero then there is no transmitted light because
        of total internal reflection.

        >>> a=vec3(1.0, -1.5, 0.8)
        >>> print a.refract(vec3(0,1,0), 1.33)
        (1.3300, -1.7920, 1.0640)
        """
        cdef vec3 res
        cdef double dot, k, ceta

        res = vec3()

        ceta = eta
        dot = self.x*N.x + self.y*N.y + self.z*N.z
        k   = 1.0 - ceta*ceta*(1.0 - dot*dot)
        if k<=eps:
            res.x = 0.0
            res.y = 0.0
            res.z = 0.0
        else:
            k = ceta*dot + sqrt(k)
            res.x = ceta*self.x - k*N.x
            res.y = ceta*self.y - k*N.y
            res.z = ceta*self.z - k*N.z

        return res

    def ortho(self):
        """Returns an orthogonal vector.

        Returns a vector that is orthogonal to self (where
        self*self.ortho()==0).

        >>> a=vec3(1.0, -1.5, 0.8)
        >>> print round(a*a.ortho(),8)
        0.0
        """
        cdef vec3 res
        cdef double x,y,z
        res = vec3()

        x=fabs(self.x)
        y=fabs(self.y)
        z=fabs(self.z)
        # Is z the smallest element? Then use x and y
        if z<=x and z<=y:
            res.x = -self.y
            res.y = self.x
            res.z = 0.0
        # Is y smallest element? Then use x and z
        elif y<=x and y<=z:
            res.x = -self.z
            res.y = 0.0
            res.z = self.x
        # x is smallest
        else:
            res.x = 0.0
            res.y = -self.z
            res.z = self.y

        return res


    def fadd(self, vec3 a, vec3 b):
        self.x = a.x+b.x
        self.y = a.y+b.y
        self.z = a.z+b.z
        return self

    def fnormalize(self, vec3 a):
        """Normalize vector a and store result in self.
        """
        cdef double nlen

        nlen = sqrt(a.x*a.x + a.y*a.y + a.z*a.z)
        if nlen<=eps:
            raise ZeroDivisionError,"null vector can't be normalized"
        nlen = 1.0/nlen
        self.x = a.x*nlen
        self.y = a.y*nlen
        self.z = a.z*nlen
        return self

######################################################################

# vec4iter
cdef class vec4iter:
    """vec4 iterator."""
    cdef int idx
    cdef object v

    def __cinit__(self, v):
        self.idx = 0
        self.v   = v

    def __iter__(self):
        return self

    def __next__(self):
        if self.idx==0:
            self.idx=1
            return self.v.x
        elif self.idx==1:
            self.idx=2
            return self.v.y
        elif self.idx==2:
            self.idx=3
            return self.v.z
        elif self.idx==3:
            self.idx=4
            return self.v.w
        else:
            raise StopIteration

# vec4
cdef class vec4:
    """Four-dimensional vector.

    This class represents a 4D vector.
    """

    cdef double x,y,z,w

    def __cinit__(self, *args):
        """Constructor.

        There are several possibilities how to initialize a vector:

        v = vec4()        -> v = <0,0,0,0>
        v = vec4(a)       -> v = <a,a,a,a>
        v = vec4(x,y)     -> v = <x,y,0,0>
        v = vec4(x,y,z)   -> v = <x,y,z,0>
        v = vec4(x,y,z,w) -> v = <x,y,z,w>

        Note that specifying just one value sets all four components to
        that value.

        Additionally you can wrap those values in a list or a tuple or
        specify them as a string:

        v = vec4([1,2,3]) -> v = <1,2,3,0>
        v = vec4("4,5")   -> v = <4,5,0,0>
        """

        cdef int arglen, seqlen
        cdef vec4 v
        arglen = len(args)

        if arglen==0:
            self.x = 0.0
            self.y = 0.0
            self.z = 0.0
            self.w = 0.0
        elif arglen==1:
            T = type(args[0])
            # scalar
            if T==float or T==int or T==long:
                self.x = args[0]
                self.y = self.x
                self.z = self.x
                self.w = self.x
            # vec4
            elif T==vec4:
                v = args[0]
                self.x = v.x
                self.y = v.y
                self.z = v.z
                self.w = v.w
            # String
            elif T==str or T==unicode:
                s=args[0].replace(","," ").replace("  "," ").strip().split(" ")
                if s==[""]:
                    s=[]
                for i in range(len(s)):
                    try:
                        s[i] = float(s[i])
                    except:
                        raise ValueError,"vec4() arg is a malformed string"
                v = vec4(s)
                self.x = v.x
                self.y = v.y
                self.z = v.z
                self.w = v.w
            # Sequence of floats
            else:
                seq = args[0]
                try:
                    seqlen = len(seq)
                except:
                    raise TypeError,"vec4() arg can't be converted to vec4"

                if seqlen==4:
                    self.x = seq[0]
                    self.y = seq[1]
                    self.z = seq[2]
                    self.w = seq[3]
                elif seqlen==3:
                    self.x = seq[0]
                    self.y = seq[1]
                    self.z = seq[2]
                    self.w = 0.0
                elif seqlen==2:
                    self.x = seq[0]
                    self.y = seq[1]
                    self.z = 0.0
                    self.w = 0.0
                elif seqlen==1:
                    self.x = seq[0]
                    self.y = self.x
                    self.z = self.x
                    self.w = self.x
                elif seqlen==0:
                    self.x = 0.0
                    self.y = 0.0
                    self.z = 0.0
                    self.w = 0.0
                else:
                    raise TypeError, "vec4() takes at most 4 arguments"

        elif arglen==2:
            self.x = args[0]
            self.y = args[1]
            self.z = 0.0
            self.w = 0.0

        elif arglen==3:
            self.x = args[0]
            self.y = args[1]
            self.z = args[2]
            self.w = 0.0

        elif arglen==4:
            self.x = args[0]
            self.y = args[1]
            self.z = args[2]
            self.w = args[3]

        else:
            raise TypeError, "vec4() takes at most 4 arguments"

    def __reduce__(self):
        return (_vec4construct, (self.x, self.y, self.z, self.w))

    def __repr__(self):
        return 'vec4(%s, %s, %s, %s)'%(repr(self.x),repr(self.y),repr(self.z),repr(self.w))

    def __str__(self):
        return "(%1.4f, %1.4f, %1.4f, %1.4f)"%(self.x, self.y, self.z, self.w)

    def __iter__(self):
        return vec4iter(self)

    def __richcmp__(a,b,op):
        cdef vec4 va, vb
        cdef int cop

        cop = op

        ta = type(a)
        tb = type(b)
        if (ta!=vec4 or tb!=vec4):
            if cop==3:
                return 1
            else:
                return 0

        va = a
        vb = b

        # <
        if cop==0:
            return (va.x<vb.x and va.y<vb.y and va.z<vb.z and va.w<vb.w)
        # <=
        elif cop==1:
            return (va.x-eps<=vb.x and va.y-eps<=vb.y and va.z-eps<=vb.z and va.w-eps<=vb.w)
        # ==
        elif cop==2:
            return (fabs(va.x-vb.x)<=eps and fabs(va.y-vb.y)<=eps and fabs(va.z-vb.z)<=eps and fabs(va.w-vb.w)<=eps)
        # !=
        elif cop==3:
            return (fabs(va.x-vb.x)>eps or fabs(va.y-vb.y)>eps and fabs(va.z-vb.z)>eps and fabs(va.w-vb.w)>eps)
        # >
        elif cop==4:
            return (va.x>vb.x and va.y>vb.y and va.z>vb.z and va.w>vb.w)
        # >=
        elif cop==5:
            return (va.x+eps>=vb.x and va.y+eps>=vb.y and va.z+eps>=vb.z and va.w+eps>=vb.w)
        # sonst (Fehler)
        else:
            raise ValueError,"internal error: illegal rich comparison number"

    def __add__(vec4 a, vec4 b):
        """Vector addition.

        >>> a=vec4(1.0, 0.5, -1.8, 0.2)
        >>> b=vec4(-0.3, 0.75, 0.5, 0.3)
        >>> print a+b
        (0.7000, 1.2500, -1.3000, 0.5000)
        """
        cdef vec4 res
        res = vec4()
        res.x = a.x+b.x
        res.y = a.y+b.y
        res.z = a.z+b.z
        res.w = a.w+b.w
        return res

    def __sub__(vec4 a, vec4 b):
        """Vector subtraction.

        >>> a=vec4(1.0, 0.5, -1.8, 0.2)
        >>> b=vec4(-0.3, 0.75, 0.5, 0.3)
        >>> print a-b
        (1.3000, -0.2500, -2.3000, -0.1000)
        """
        cdef vec4 res
        res = vec4()
        res.x = a.x-b.x
        res.y = a.y-b.y
        res.z = a.z-b.z
        res.w = a.w-b.w
        return res

    def __mul__(a, b):
        """Multiplication with a scalar or dot product.

        >>> a=vec4(1.0, 0.5, -1.8, 0.2)
        >>> b=vec4(-0.3, 0.75, 0.5, 0.3)
        >>> print a*2.0
        (2.0000, 1.0000, -3.6000, 0.4000)
        >>> print 2.0*a
        (2.0000, 1.0000, -3.6000, 0.4000)
        >>> print a*b
        -0.765
        """

        cdef vec4 res, va, vb
        cdef double r

        ta = type(a)
        tb = type(b)

        if ta==vec4:
            va = a
            if tb==vec4:
                # vec4*vec4   (dot product a*b)
                vb = b
                return va.x*vb.x + va.y*vb.y + va.z*vb.z + va.w*vb.w
            elif tb==float or tb==int or tb==long:
                # vec4*scalar
                res = vec4()
                r   = b
                res.x = va.x*r
                res.y = va.y*r
                res.z = va.z*r
                res.w = va.w*r
                return res
            elif tb==mat4:
                # vec4*mat4
                return b.__rmul__(a)
        elif ta==float or ta==int or ta==long:
            if tb==vec4:
                # scalar*vec4
                res = vec4()
                vb  = b
                r   = a
                res.x = r*vb.x
                res.y = r*vb.y
                res.z = r*vb.z
                res.w = r*vb.w
                return res

        raise TypeError, "unsupported operand type for *"

    def __div__(vec4 a, b):
        """Division by scalar.

        >>> a=vec4(1.0, 0.5, -1.8, 0.2)
        >>> print a/2.0
        (0.5000, 0.2500, -0.9000, 0.1000)
        """

        cdef vec4 res
        cdef double r

        tb = type(b)

        if tb==float or tb==int or tb==long:
            # vec4/scalar
            res = vec4()
            r = b
            if fabs(r)<=eps:
                raise ZeroDivisionError,"vec4 division"
            res.x = a.x/r
            res.y = a.y/r
            res.z = a.z/r
            res.w = a.w/r
            return res
        else:
            raise TypeError, "unsupported operand type for /"

    def __mod__(vec4 self, b):
        """Modulo (component wise)

        >>> a=vec4(3.0, 2.5, -1.8, 0.2)
        >>> print a%2.0
        (1.0000, 0.5000, 0.2000, 0.2000)
        """
        cdef vec4 res, vb
        cdef double r

        tb = type(b)

        if tb==float or tb==int or tb==long:
            # vec4%scalar
            res = vec4()
            r = b
            if fabs(r)<=eps:
                raise ZeroDivisionError,"vec4 modulo"
            res.x = _fmod(self.x,r)
            res.y = _fmod(self.y,r)
            res.z = _fmod(self.z,r)
            res.w = _fmod(self.w,r)
            return res
        elif tb==vec4:
            # vec4%vec4
            res = vec4()
            vb = b
            if (fabs(vb.x)<=eps or fabs(vb.y)<=eps or
                fabs(vb.z)<=eps or fabs(vb.w)<=eps):
                raise ZeroDivisionError,"vec4 modulo"
            res.x = _fmod(self.x, vb.x)
            res.y = _fmod(self.y, vb.y)
            res.z = _fmod(self.z, vb.z)
            res.w = _fmod(self.w, vb.w)
            return res
        else:
            raise TypeError, "unsupported operand type for %"

    def __neg__(self):
        """Negation.

        >>> a=vec4(3.0, 2.5, -1.8, 0.2)
        >>> print -a
        (-3.0000, -2.5000, 1.8000, -0.2000)
        """
        cdef vec4 res
        res = vec4()
        res.x = -self.x
        res.y = -self.y
        res.z = -self.z
        res.w = -self.w
        return res

    def __pos__(self):
        """
        >>> a=vec4(3.0, 2.5, -1.8, 0.2)
        >>> print +a
        (3.0000, 2.5000, -1.8000, 0.2000)
        """
        cdef vec4 res
        res = vec4()
        res.x = self.x
        res.y = self.y
        res.z = self.z
        res.w = self.w
        return res

    def __iadd__(self, vec4 other):
        """Inline vector addition.

        >>> a=vec4(1.0, 0.5, -1.8, 0.2)
        >>> b=vec4(-0.3, 0.75, 0.5, 0.3)
        >>> a+=b
        >>> print a
        (0.7000, 1.2500, -1.3000, 0.5000)
        """
        self.x=self.x+other.x
        self.y=self.y+other.y
        self.z=self.z+other.z
        self.w=self.w+other.w
        return self

    def __isub__(self, vec4 other):
        """Inline vector subtraction.

        >>> a=vec4(1.0, 0.5, -1.8, 0.2)
        >>> b=vec4(-0.3, 0.75, 0.5, 0.3)
        >>> a-=b
        >>> print a
        (1.3000, -0.2500, -2.3000, -0.1000)
        """
        self.x=self.x-other.x
        self.y=self.y-other.y
        self.z=self.z-other.z
        self.w=self.w-other.w
        return self

    def __imul__(self, other):
        """Inline multiplication (only with scalar)

        >>> a=vec4(1.0, 0.5, -1.8, 0.2)
        >>> a*=2.0
        >>> print a
        (2.0000, 1.0000, -3.6000, 0.4000)
        """
        cdef double r
        tb = type(other)
        if tb==float or tb==int or tb==long:
            r=other
            self.x=self.x*r
            self.y=self.y*r
            self.z=self.z*r
            self.w=self.w*r
            return self
        else:
            raise TypeError, "unsupported operand type for *="

    def __idiv__(self, other):
        """Inline division with scalar

        >>> a=vec4(1.0, 0.5, -1.8, 0.2)
        >>> a/=2.0
        >>> print a
        (0.5000, 0.2500, -0.9000, 0.1000)
        """
        cdef double r
        tb = type(other)
        if tb==float or tb==int or tb==long:
            r=other
            if fabs(r)<=eps:
                raise ZeroDivisionError,"vec4 division"
            self.x=self.x/r
            self.y=self.y/r
            self.z=self.z/r
            self.w=self.w/r
            return self
        else:
            raise TypeError, "unsupported operand type for /="

    def __imod__(self, b):
        """Inline modulo.

        >>> a=vec4(3.0, 2.5, -1.8, 0.2)
        >>> a%=2.0
        >>> print a
        (1.0000, 0.5000, 0.2000, 0.2000)
        """
        cdef vec4 vb
        cdef double r

        tb = type(b)

        if tb==float or tb==int or tb==long:
            # vec4%=scalar
            r = b
            if fabs(r)<=eps:
                raise ZeroDivisionError,"vec4 modulo"
            self.x = _fmod(self.x,r)
            self.y = _fmod(self.y,r)
            self.z = _fmod(self.z,r)
            self.w = _fmod(self.w,r)
            return self
        elif tb==vec4:
            # vec4%=vec4
            vb = b
            if (fabs(vb.x)<=eps or fabs(vb.y)<=eps or
                fabs(vb.z)<=eps or fabs(vb.w)<=eps):
                raise ZeroDivisionError,"vec4 modulo"
            self.x = _fmod(self.x, vb.x)
            self.y = _fmod(self.y, vb.y)
            self.z = _fmod(self.z, vb.z)
            self.w = _fmod(self.w, vb.w)
            return self
        else:
            raise TypeError, "unsupported operand type for %="


    def __abs__(self):
        """Return the length of the vector.

        abs(v) is equivalent to v.length().

        >>> a=vec4(1.0, 0.5, -1.8, 0.2)
        >>> print abs(a)
        2.12837966538
        """
        cdef double a
        a = self.x*self.x + self.y*self.y + self.z*self.z + self.w*self.w
        return sqrt(a)

    def __getattr__(self, name):
        if name=="x":
            return self.x
        elif name=="y":
            return self.y
        elif name=="z":
            return self.z
        elif name=="w" or name=="t":
            return self.w
        else:
            raise AttributeError,"vec4 has no attribute '"+name+"'"

    def __setattr__(self, name, val):
        if name=="x":
            self.x = val
        elif name=="y":
            self.y = val
        elif name=="z":
            self.z = val
        elif name=="w" or name=="t":
            self.w = val
        else:
            raise AttributeError,"vec4 has no writable attribute '"+name+"'"

    def __len__(self):
        """Length of the sequence (always 4)."""
        return 4

    def __getitem__(self, key):
        """Return a component by index (0-based).

        >>> a=vec4(1.0, 0.5, -1.8, 0.2)
        >>> print a[0]
        1.0
        >>> print a[1]
        0.5
        >>> print a[2]
        -1.8
        >>> print a[3]
        0.2
        """
        cdef int k

        T=type(key)
        if T!=int and T!=long:
            raise TypeError, "index must be integer"

        k = key
        if   k==0: return self.x
        elif k==1: return self.y
        elif k==2: return self.z
        elif k==3: return self.w
        else:
            raise IndexError,"index out of range"

    def __setitem__(self, key, value):
        """Set a component by index (0-based).

        >>> a=vec4()
        >>> a[0]=1.5; a[1]=0.7; a[2]=-0.3; a[3]=0.2
        >>> print a
        (1.5000, 0.7000, -0.3000, 0.2000)
        """
        cdef int k

        T=type(key)
        if T!=int and T!=long:
            raise TypeError, "index must be integer"

        k = key
        if   k==0: self.x = value
        elif k==1: self.y = value
        elif k==2: self.z = value
        elif k==3: self.w = value
        else:
            raise IndexError,"index out of range"

    def length(self):
        """Return the length of the vector.

        v.length() is equivalent to abs(v).

        >>> a=vec4(1.0, 0.5, -1.8, 0.2)
        >>> print a.length()
        2.12837966538
        """
        cdef double a
        a = self.x*self.x + self.y*self.y + self.z*self.z + self.w*self.w
        return sqrt(a)

    def normalize(self):
        """Return normalized vector.

        >>> a=vec4(1.0, 0.5, -1.8, 1.2)
        >>> print a.normalize()
        (0.4107, 0.2053, -0.7392, 0.4928)
        """
        cdef vec4 res
        cdef double nlen

        nlen = sqrt(self.x*self.x + self.y*self.y + self.z*self.z + self.w*self.w)
        if nlen<=eps:
            raise ZeroDivisionError,"null vector can't be normalized"
        nlen = 1.0/nlen
        res = vec4()
        res.x = self.x*nlen
        res.y = self.y*nlen
        res.z = self.z*nlen
        res.w = self.w*nlen
        return res

######################################################################

# mat3iter
cdef class mat3iter:
    """mat iterator."""
    cdef int idx
    cdef object m

    def __cinit__(self, m):
        self.idx = 0
        self.m   = m

    def __iter__(self):
        return self

    def __next__(self):
        if self.idx==0:
            self.idx=1
            return self.m[0]
        elif self.idx==1:
            self.idx=2
            return self.m[1]
        elif self.idx==2:
            self.idx=3
            return self.m[2]
        else:
            raise StopIteration

# mat3
cdef class mat3:

    cdef double m11, m12, m13
    cdef double m21, m22, m23
    cdef double m31, m32, m33

    def __cinit__(self, *args):
        """Constructor.

        There are several possibilities how to initialize a matrix,
        depending on the number of arguments you provide to the constructor.

        - 0 arguments: Every component is zero.
        - 1 number argument: The diagonal is initialized to that number,
          all the other elements are zero.
        - 1 sequence argument: The elements are initialized with the numbers
          in the sequence (the sequence must contain 9 numbers).
        - 1 mat3 argument: The matrix is copied.
        - 3 sequence arguments: The columns are initialized with the
          respective sequence (each sequence must contain 3 numbers).
        - 9 number arguments: The matrix is initialized with those values
          (row-major order).
        """
        cdef int arglen
        cdef mat3 B
        arglen = len(args)

        # No arguments
        if arglen==0:
            m11=m12=m13=0.0
            m21=m22=m23=0.0
            m31=m32=m33=0.0

        # 1 argument (list, scalar or mat3)
        elif arglen==1:
            T = type(args[0])
            # scalar
            if T==float or T==int or T==long:
                self.m11 = args[0]
                self.m12 = 0.0
                self.m13 = 0.0
                self.m21 = 0.0
                self.m22 = self.m11
                self.m23 = 0.0
                self.m31 = 0.0
                self.m32 = 0.0
                self.m33 = self.m11

            # mat3
            elif T==mat3:
                B = args[0]
                self.m11 = B.m11
                self.m12 = B.m12
                self.m13 = B.m13
                self.m21 = B.m21
                self.m22 = B.m22
                self.m23 = B.m23
                self.m31 = B.m31
                self.m32 = B.m32
                self.m33 = B.m33
            # String
            elif T==str or T==unicode:
                s=args[0].replace(","," ").replace("  "," ").strip().split(" ")
                if s==[""]:
                    s=[]
                for i in range(len(s)):
                    try:
                        s[i] = float(s[i])
                    except:
                        raise ValueError,"mat3() arg is a malformed string"
                B = mat3(s)
                self.m11 = B.m11
                self.m12 = B.m12
                self.m13 = B.m13
                self.m21 = B.m21
                self.m22 = B.m22
                self.m23 = B.m23
                self.m31 = B.m31
                self.m32 = B.m32
                self.m33 = B.m33
            # Sequence of 9 floats
            else:
                seq = args[0]
                try:
                    seqlen = len(seq)
                except:
                    raise TypeError,"mat3() arg can't be converted to mat3"

                if seqlen==9:
                    self.m11 = seq[0]
                    self.m12 = seq[1]
                    self.m13 = seq[2]
                    self.m21 = seq[3]
                    self.m22 = seq[4]
                    self.m23 = seq[5]
                    self.m31 = seq[6]
                    self.m32 = seq[7]
                    self.m33 = seq[8]
                elif seqlen==3:
                    b = seq[0]
                    self.m11 = b[0]
                    self.m21 = b[1]
                    self.m31 = b[2]
                    b = seq[1]
                    self.m12 = b[0]
                    self.m22 = b[1]
                    self.m32 = b[2]
                    b = seq[2]
                    self.m13 = b[0]
                    self.m23 = b[1]
                    self.m33 = b[2]
                else:
                    raise TypeError, "mat3() takes at most 9 arguments"

        # 3 arguments (columns as sequences)
        elif arglen==3:
            a,b,c=args
            self.m11 = a[0]
            self.m12 = b[0]
            self.m13 = c[0]
            self.m21 = a[1]
            self.m22 = b[1]
            self.m23 = c[1]
            self.m31 = a[2]
            self.m32 = b[2]
            self.m33 = c[2]

        # 9 arguments
        elif arglen==9:
            self.m11 = args[0]
            self.m12 = args[1]
            self.m13 = args[2]
            self.m21 = args[3]
            self.m22 = args[4]
            self.m23 = args[5]
            self.m31 = args[6]
            self.m32 = args[7]
            self.m33 = args[8]

        else:
            raise TypeError,"mat3() arg can't be converted to mat3"

    def __reduce__(self):
        return (_mat3construct, (self[0], self[1], self[2]))

    def __repr__(self):
        return 'mat3(%s, %s, %s, %s, %s, %s, %s, %s, %s)'%(
            repr(self.m11),repr(self.m12),repr(self.m13),
            repr(self.m21),repr(self.m22),repr(self.m23),
            repr(self.m31),repr(self.m32),repr(self.m33))

    def __str__(self):
        return "[%9.4f, %9.4f, %9.4f]\n[%9.4f, %9.4f, %9.4f]\n[%9.4f, %9.4f, %9.4f]"%(
            self.m11, self.m12, self.m13,
            self.m21, self.m22, self.m23,
            self.m31, self.m32, self.m33)

    def __iter__(self):
        return mat3iter(self)

    def __richcmp__(a,b,op):
        cdef mat3 ma, mb
        cdef int cop

        cop = op

        ta = type(a)
        tb = type(b)
        if (ta!=mat3 or tb!=mat3):
            if cop==3:
                return 1
            else:
                return 0

        ma = a
        mb = b

        # <
        if cop==0:
            return 0
        # <=
        elif cop==1:
            return 0
        # ==
        elif cop==2:
            return (fabs(ma.m11-mb.m11)<=eps and fabs(ma.m12-mb.m12)<=eps and
                    fabs(ma.m13-mb.m13)<=eps and
                    fabs(ma.m21-mb.m21)<=eps and fabs(ma.m22-mb.m22)<=eps and
                    fabs(ma.m23-mb.m23)<=eps and
                    fabs(ma.m31-mb.m31)<=eps and fabs(ma.m32-mb.m32)<=eps and
                    fabs(ma.m33-mb.m33)<=eps)
        # !=
        elif cop==3:
            return (fabs(ma.m11-mb.m11)>eps or fabs(ma.m12-mb.m12)>eps or
                    fabs(ma.m13-mb.m13)>eps or
                    fabs(ma.m21-mb.m21)>eps or fabs(ma.m22-mb.m22)>eps or
                    fabs(ma.m23-mb.m23)>eps or
                    fabs(ma.m31-mb.m31)>eps or fabs(ma.m32-mb.m32)>eps or
                    fabs(ma.m33-mb.m33)>eps)
        # >
        elif cop==4:
            return 0
        # >=
        elif cop==5:
            return 0
        # sonst (Fehler)
        else:
            raise ValueError,"internal error: illegal rich comparison number"


    def __add__(mat3 a, mat3 b):
        """Matrix addition.
        """
        cdef mat3 res
        res = mat3()
        res.m11 = a.m11 + b.m11
        res.m12 = a.m12 + b.m12
        res.m13 = a.m13 + b.m13
        res.m21 = a.m21 + b.m21
        res.m22 = a.m22 + b.m22
        res.m23 = a.m23 + b.m23
        res.m31 = a.m31 + b.m31
        res.m32 = a.m32 + b.m32
        res.m33 = a.m33 + b.m33
        return res

    def __sub__(mat3 a, mat3 b):
        """Matrix subtraction.
        """
        cdef mat3 res
        res = mat3()
        res.m11 = a.m11 - b.m11
        res.m12 = a.m12 - b.m12
        res.m13 = a.m13 - b.m13
        res.m21 = a.m21 - b.m21
        res.m22 = a.m22 - b.m22
        res.m23 = a.m23 - b.m23
        res.m31 = a.m31 - b.m31
        res.m32 = a.m32 - b.m32
        res.m33 = a.m33 - b.m33
        return res

    def __mul__(a, b):
        """Multiplication.
        """
        cdef mat3 ma, mb, res
        cdef vec3 va, vb, vres
        cdef double r

        ta = type(a)
        tb = type(b)

        # mat3 on the left?
        if ta==mat3:
            ma = a
            # mat3 on the right?
            if tb==mat3:
                # mat3*mat3
                res = mat3()
                mb = b
                res.m11 = ma.m11*mb.m11+ma.m12*mb.m21+ma.m13*mb.m31
                res.m12 = ma.m11*mb.m12+ma.m12*mb.m22+ma.m13*mb.m32
                res.m13 = ma.m11*mb.m13+ma.m12*mb.m23+ma.m13*mb.m33
                res.m21 = ma.m21*mb.m11+ma.m22*mb.m21+ma.m23*mb.m31
                res.m22 = ma.m21*mb.m12+ma.m22*mb.m22+ma.m23*mb.m32
                res.m23 = ma.m21*mb.m13+ma.m22*mb.m23+ma.m23*mb.m33
                res.m31 = ma.m31*mb.m11+ma.m32*mb.m21+ma.m33*mb.m31
                res.m32 = ma.m31*mb.m12+ma.m32*mb.m22+ma.m33*mb.m32
                res.m33 = ma.m31*mb.m13+ma.m32*mb.m23+ma.m33*mb.m33
                return res
            # vec3 on the right?
            elif tb==vec3:
                # mat3*vec3
                wres = vec3()
                wb   = b
                wres.x = ma.m11*wb.x + ma.m12*wb.y + ma.m13*wb.z
                wres.y = ma.m21*wb.x + ma.m22*wb.y + ma.m23*wb.z
                wres.z = ma.m31*wb.x + ma.m32*wb.y + ma.m33*wb.z
                return wres
            # Scalar on the right?
            elif tb==float or tb==int or tb==long:
                res = mat3()
                r = b
                res.m11 = ma.m11*r
                res.m12 = ma.m12*r
                res.m13 = ma.m13*r
                res.m21 = ma.m21*r
                res.m22 = ma.m22*r
                res.m23 = ma.m23*r
                res.m31 = ma.m31*r
                res.m32 = ma.m32*r
                res.m33 = ma.m33*r
                return res

        # vec3 on the left? (then there must be a mat3 on the right)
        elif ta==vec3:
            # vec3*mat3
            wres = vec3()
            wa   = a
            mb   = b
            wres.x = wa.x*mb.m11 + wa.y*mb.m21 + wa.z*mb.m31
            wres.y = wa.x*mb.m12 + wa.y*mb.m22 + wa.z*mb.m32
            wres.z = wa.x*mb.m13 + wa.y*mb.m23 + wa.z*mb.m33
            return wres

        # Scalar on the left? (then there must be a mat3 on the right)
        elif ta==float or ta==int or ta==long:
            res = mat3()
            r  = a
            mb = b
            res.m11 = r*mb.m11
            res.m12 = r*mb.m12
            res.m13 = r*mb.m13
            res.m21 = r*mb.m21
            res.m22 = r*mb.m22
            res.m23 = r*mb.m23
            res.m31 = r*mb.m31
            res.m32 = r*mb.m32
            res.m33 = r*mb.m33
            return res

        # No valid combination found
        raise TypeError, "unsupported operand type for *"

    def __div__(mat3 a, b):
        """Division.
        """
        cdef mat3 res
        cdef double r

        tb = type(b)

        if tb==float or tb==int or tb==long:
            # mat3/scalar
            res = mat3()
            r = b
            if fabs(r)<=eps:
                raise ZeroDivisionError,"mat3 division"
            res.m11 = a.m11/r
            res.m12 = a.m12/r
            res.m13 = a.m13/r
            res.m21 = a.m21/r
            res.m22 = a.m22/r
            res.m23 = a.m23/r
            res.m31 = a.m31/r
            res.m32 = a.m32/r
            res.m33 = a.m33/r
            return res
        else:
            raise TypeError, "unsupported operand type for /"

    def __mod__(mat3 self, other):
        """Modulo.
        """
        cdef mat3 mb, res

        tb = type(other)

        if tb==float or tb==int or tb==long:
            r = other
            if fabs(r)<=eps:
                raise ZeroDivisionError,"mat3 modulo"
            res = mat3()
            res.m11 = _fmod(self.m11, r)
            res.m12 = _fmod(self.m12, r)
            res.m13 = _fmod(self.m13, r)
            res.m21 = _fmod(self.m21, r)
            res.m22 = _fmod(self.m22, r)
            res.m23 = _fmod(self.m23, r)
            res.m31 = _fmod(self.m31, r)
            res.m32 = _fmod(self.m32, r)
            res.m33 = _fmod(self.m33, r)
            return res

        elif tb==mat3:
            mb = other

            if (fabs(mb.m11)<=eps or fabs(mb.m12)<=eps or fabs(mb.m13)<=eps or
                fabs(mb.m21)<=eps or fabs(mb.m22)<=eps or fabs(mb.m23)<=eps or
                fabs(mb.m31)<=eps or fabs(mb.m32)<=eps or fabs(mb.m33)<=eps):
                raise ZeroDivisionError,"mat3 modulo"

            res=mat3()
            res.m11 = _fmod(self.m11, mb.m11)
            res.m12 = _fmod(self.m12, mb.m12)
            res.m13 = _fmod(self.m13, mb.m13)
            res.m21 = _fmod(self.m21, mb.m21)
            res.m22 = _fmod(self.m22, mb.m22)
            res.m23 = _fmod(self.m23, mb.m23)
            res.m31 = _fmod(self.m31, mb.m31)
            res.m32 = _fmod(self.m32, mb.m32)
            res.m33 = _fmod(self.m33, mb.m33)
            return res

        else:
            raise TypeError, "unsupported operand type for %"


    def __neg__(self):
        """Negation.
        """
        cdef mat3 res
        res = mat3()
        res.m11 = -self.m11
        res.m12 = -self.m12
        res.m13 = -self.m13
        res.m21 = -self.m21
        res.m22 = -self.m22
        res.m23 = -self.m23
        res.m31 = -self.m31
        res.m32 = -self.m32
        res.m33 = -self.m33
        return res

    def __pos__(self):
        cdef mat3 res
        res = mat3()
        res.m11 = self.m11
        res.m12 = self.m12
        res.m13 = self.m13
        res.m21 = self.m21
        res.m22 = self.m22
        res.m23 = self.m23
        res.m31 = self.m31
        res.m32 = self.m32
        res.m33 = self.m33
        return res

    def __iadd__(self, mat3 other):
        """Inline matrix addition.
        """
        self.m11 = self.m11 + other.m11
        self.m12 = self.m12 + other.m12
        self.m13 = self.m13 + other.m13
        self.m21 = self.m21 + other.m21
        self.m22 = self.m22 + other.m22
        self.m23 = self.m23 + other.m23
        self.m31 = self.m31 + other.m31
        self.m32 = self.m32 + other.m32
        self.m33 = self.m33 + other.m33
        return self

    def __isub__(self, mat3 other):
        """Inline matrix subtraction.
        """
        self.m11 = self.m11 - other.m11
        self.m12 = self.m12 - other.m12
        self.m13 = self.m13 - other.m13
        self.m21 = self.m21 - other.m21
        self.m22 = self.m22 - other.m22
        self.m23 = self.m23 - other.m23
        self.m31 = self.m31 - other.m31
        self.m32 = self.m32 - other.m32
        self.m33 = self.m33 - other.m33
        return self

    def __imul__(self, other):
        """Inline multiplication (scalar or matrix).
        """
        cdef double r, a,b,c
        cdef mat3 mb

        tb = type(other)
        if tb==float or tb==int or tb==long:
            r=other
            self.m11 = self.m11*r
            self.m12 = self.m12*r
            self.m13 = self.m13*r
            self.m21 = self.m21*r
            self.m22 = self.m22*r
            self.m23 = self.m23*r
            self.m31 = self.m31*r
            self.m32 = self.m32*r
            self.m33 = self.m33*r
            return self
        elif tb==mat3:
            mb = other
            a = self.m11
            b = self.m12
            c = self.m13
            self.m11 = a*mb.m11+b*mb.m21+c*mb.m31
            self.m12 = a*mb.m12+b*mb.m22+c*mb.m32
            self.m13 = a*mb.m13+b*mb.m23+c*mb.m33
            a = self.m21
            b = self.m22
            c = self.m23
            self.m21 = a*mb.m11+b*mb.m21+c*mb.m31
            self.m22 = a*mb.m12+b*mb.m22+c*mb.m32
            self.m23 = a*mb.m13+b*mb.m23+c*mb.m33
            a = self.m31
            b = self.m32
            c = self.m33
            self.m31 = a*mb.m11+b*mb.m21+c*mb.m31
            self.m32 = a*mb.m12+b*mb.m22+c*mb.m32
            self.m33 = a*mb.m13+b*mb.m23+c*mb.m33
            return self
        else:
            raise TypeError, "unsupported operand type for *="

    def __idiv__(self, other):
        """Inline division with scalar.
        """
        cdef double r
        tb = type(other)
        if tb==float or tb==int or tb==long:
            r=other
            if fabs(r)<=eps:
                raise ZeroDivisionError,"mat3 division"
            self.m11 = self.m11/r
            self.m12 = self.m12/r
            self.m13 = self.m13/r
            self.m21 = self.m21/r
            self.m22 = self.m22/r
            self.m23 = self.m23/r
            self.m31 = self.m31/r
            self.m32 = self.m32/r
            self.m33 = self.m33/r
            return self
        else:
            raise TypeError, "unsupported operand type for /="

    def __imod__(self, b):
        """Inline modulo.
        """
        cdef mat3 mb
        cdef double r

        tb = type(b)

        if tb==float or tb==int or tb==long:
            r = b
            if fabs(r)<=eps:
                raise ZeroDivisionError,"mat3 modulo"
            self.m11 = _fmod(self.m11, r)
            self.m12 = _fmod(self.m12, r)
            self.m13 = _fmod(self.m13, r)
            self.m21 = _fmod(self.m21, r)
            self.m22 = _fmod(self.m22, r)
            self.m23 = _fmod(self.m23, r)
            self.m31 = _fmod(self.m31, r)
            self.m32 = _fmod(self.m32, r)
            self.m33 = _fmod(self.m33, r)
            return self
        elif tb==mat3:
            mb = b
            if (fabs(mb.m11)<=eps or fabs(mb.m12)<=eps or fabs(mb.m13)<=eps or
                fabs(mb.m21)<=eps or fabs(mb.m22)<=eps or fabs(mb.m23)<=eps or
                fabs(mb.m31)<=eps or fabs(mb.m32)<=eps or fabs(mb.m33)<=eps):
                raise ZeroDivisionError,"mat3 modulo"
            self.m11 = _fmod(self.m11, mb.m11)
            self.m12 = _fmod(self.m12, mb.m12)
            self.m13 = _fmod(self.m13, mb.m13)
            self.m21 = _fmod(self.m21, mb.m21)
            self.m22 = _fmod(self.m22, mb.m22)
            self.m23 = _fmod(self.m23, mb.m23)
            self.m31 = _fmod(self.m31, mb.m31)
            self.m32 = _fmod(self.m32, mb.m32)
            self.m33 = _fmod(self.m33, mb.m33)
            return self
        else:
            raise TypeError, "unsupported operand type for %="


    def __len__(self):
        """Length of the sequence (always 3)."""
        return 3

    def __getitem__(mat3 self, key):
        """Return a column or a single matrix element.
        """
        cdef int i,j

        T=type(key)
        if T==int or T==long:
            i = key
            if i==0:
                return vec3(self.m11, self.m21, self.m31)
            elif i==1:
                return vec3(self.m12, self.m22, self.m32)
            elif i==2:
                return vec3(self.m13, self.m23, self.m33)
            else:
                raise IndexError, "index out of range"
        elif T==tuple:
            if len(key)!=2:
                raise ValueError, "index tuple must be a 2-tuple"
            i = key[0]
            j = key[1]
            if i<0 or i>2 or j<0 or j>2:
                raise IndexError, "index out of range"

            if i==0:
                if j==0:
                    return self.m11
                elif j==1:
                    return self.m12
                elif j==2:
                    return self.m13
            elif i==1:
                if j==0:
                    return self.m21
                elif j==1:
                    return self.m22
                elif j==2:
                    return self.m23
            elif i==2:
                if j==0:
                    return self.m31
                elif j==1:
                    return self.m32
                elif j==2:
                    return self.m33
        else:
            raise TypeError,"index must be integer or 2-tuple"

    def __setitem__(self, key, value):
        """Set a column or a single matrix element.

        The column must be given as a 3-sequence (this includes vec3 as well).
        """
        cdef int i,j
        cdef double v

        T=type(key)
        if T==int or T==long:
            if len(value)!=3:
                raise ValueError, "3-sequence expected"
            i = key
            if i==0:
                self.m11 = value[0]
                self.m21 = value[1]
                self.m31 = value[2]
            elif i==1:
                self.m12 = value[0]
                self.m22 = value[1]
                self.m32 = value[2]
            elif i==2:
                self.m13 = value[0]
                self.m23 = value[1]
                self.m33 = value[2]
            else:
                raise IndexError, "index out of range"
        elif T==tuple:
            if len(key)!=2:
                raise ValueError, "index tuple must be a 2-tuple"
            i = key[0]
            j = key[1]
            v = value
            if i<0 or i>2 or j<0 or j>2:
                raise IndexError, "index out of range"

            if i==0:
                if j==0:
                    self.m11 = v
                elif j==1:
                    self.m12 = v
                elif j==2:
                    self.m13 = v
            elif i==1:
                if j==0:
                    self.m21 = v
                elif j==1:
                    self.m22 = v
                elif j==2:
                    self.m23 = v
            elif i==2:
                if j==0:
                    self.m31 = v
                elif j==1:
                    self.m32 = v
                elif j==2:
                    self.m33 = v
        else:
            raise TypeError,"index must be integer or 2-tuple"

    def getRow(self, idx):
        """Return row (as vec3)."""
        cdef int i

        T=type(idx)
        if T!=int and T!=long:
            raise TypeError,"index must be integer"

        i = idx
        if i==0:
            return vec3(self.m11, self.m12, self.m13)
        elif i==1:
            return vec3(self.m21, self.m22, self.m23)
        elif i==2:
            return vec3(self.m31, self.m32, self.m33)
        else:
            raise IndexError, "index out of range"

    def setRow(self, idx, value):
        """Set row."""
        cdef int i

        T=type(idx)
        if T!=int and T!=long:
            raise TypeError,"index must be integer"

        if len(value)!=3:
            raise ValueError, "3-sequence expected"

        i = idx
        if i==0:
            self.m11 = value[0]
            self.m12 = value[1]
            self.m13 = value[2]
        elif i==1:
            self.m21 = value[0]
            self.m22 = value[1]
            self.m23 = value[2]
        elif i==2:
            self.m31 = value[0]
            self.m32 = value[1]
            self.m33 = value[2]
        else:
            raise IndexError, "index out of range"


    def getColumn(self, idx):
        """Return column (as vec3)."""
        cdef int i

        T=type(idx)
        if T!=int and T!=long:
            raise TypeError,"index must be integer"

        i = idx
        if i==0:
            return vec3(self.m11, self.m21, self.m31)
        elif i==1:
            return vec3(self.m12, self.m22, self.m32)
        elif i==2:
            return vec3(self.m13, self.m23, self.m33)
        else:
            raise IndexError, "index out of range"

    def setColumn(self, idx, value):
        """Set column."""
        cdef int i

        T=type(idx)
        if T!=int and T!=long:
            raise TypeError,"index must be integer"

        if len(value)!=3:
            raise ValueError, "3-sequence expected"

        i = idx
        if i==0:
            self.m11 = value[0]
            self.m21 = value[1]
            self.m31 = value[2]
        elif i==1:
            self.m12 = value[0]
            self.m22 = value[1]
            self.m32 = value[2]
        elif i==2:
            self.m13 = value[0]
            self.m23 = value[1]
            self.m33 = value[2]
        else:
            raise IndexError, "index out of range"


    def toList(self, rowmajor=0):
        """Return a list containing the matrix elements.

        By default the list is in column-major order. If you set the
        optional argument rowmajor to 1, you'll get the list in row-major
        order.
        """
        if rowmajor:
            return [self.m11, self.m12, self.m13,
                    self.m21, self.m22, self.m23,
                    self.m31, self.m32, self.m33]
        else:
            return [self.m11, self.m21, self.m31,
                    self.m12, self.m22, self.m32,
                    self.m13, self.m23, self.m33]


    def identity(self):
        """Return identity matrix.
        """
        cdef mat3 res
        res = mat3()
        res.m11 = 1.0
        res.m12 = 0.0
        res.m13 = 0.0
        res.m21 = 0.0
        res.m22 = 1.0
        res.m23 = 0.0
        res.m31 = 0.0
        res.m32 = 0.0
        res.m33 = 1.0
        return res

    def transpose(self):
        """Return transpose matrix.
        """
        cdef mat3 res
        res = mat3()
        res.m11 = self.m11
        res.m12 = self.m21
        res.m13 = self.m31
        res.m21 = self.m12
        res.m22 = self.m22
        res.m23 = self.m32
        res.m31 = self.m13
        res.m32 = self.m23
        res.m33 = self.m33
        return res

    def determinant(self):
        """Return determinant.
        """
        return self.m11*self.m22*self.m33+ \
               self.m12*self.m23*self.m31+ \
               self.m13*self.m21*self.m32- \
               self.m31*self.m22*self.m13- \
               self.m32*self.m23*self.m11- \
               self.m33*self.m21*self.m12

    def inverse(self):
        """Return inverse matrix."""
        cdef mat3 res
        cdef double d

        d = self.determinant()
        if fabs(d)<=eps:
            raise ZeroDivisionError, "matrix is not invertible"

        d = 1.0/d

        res=mat3()
        res.m11 = d*(self.m22*self.m33-self.m23*self.m32)
        res.m12 = d*(self.m32*self.m13-self.m12*self.m33)
        res.m13 = d*(self.m12*self.m23-self.m22*self.m13)
        res.m21 = d*(self.m23*self.m31-self.m21*self.m33)
        res.m22 = d*(self.m11*self.m33-self.m31*self.m13)
        res.m23 = d*(self.m21*self.m13-self.m11*self.m23)
        res.m31 = d*(self.m21*self.m32-self.m31*self.m22)
        res.m32 = d*(self.m31*self.m12-self.m11*self.m32)
        res.m33 = d*(self.m11*self.m22-self.m12*self.m21)
        return res

    def scaling(self, s):
        """Return scaling matrix."""
        cdef mat3 res
        res = mat3()
        res.m11 = s[0]
        res.m12 = 0.0
        res.m13 = 0.0
        res.m21 = 0.0
        res.m22 = s[1]
        res.m23 = 0.0
        res.m31 = 0.0
        res.m32 = 0.0
        res.m33 = s[2]
        return res

    def rotation(self, angle, axis):
        """Return rotation matrix.

        angle must be given in radians. axis must be a 3-sequence.
        """
        cdef double cangle, ax, ay, az
        cdef double sqr_a, sqr_b, sqr_c, len2
        cdef double k1, k2, k3, k1ab, k1ac, k1bc, k3a, k3b, k3c
        cdef mat3 res

        cangle   = angle
        ax    = axis[0]
        ay    = axis[1]
        az    = axis[2]
        sqr_a = ax*ax
        sqr_b = ay*ay
        sqr_c = az*az
        len2  = sqr_a+sqr_b+sqr_c

        k2    = cos(cangle)
        k1    = (1.0-k2)/len2
        k3    = sin(cangle)/sqrt(len2)
        k1ab  = k1*ax*ay
        k1ac  = k1*ax*az
        k1bc  = k1*ay*az
        k3a   = k3*ax
        k3b   = k3*ay
        k3c   = k3*az

        res = mat3()
        res.m11 = k1*sqr_a+k2
        res.m12 = k1ab-k3c
        res.m13 = k1ac+k3b
        res.m21 = k1ab+k3c
        res.m22 = k1*sqr_b+k2
        res.m23 = k1bc-k3a
        res.m31 = k1ac-k3b
        res.m32 = k1bc+k3a
        res.m33 = k1*sqr_c+k2
        return res


    def scale(self, s):
        """Concatenate a scaling."""
        cdef double sx, sy, sz
        sx = s[0]
        sy = s[1]
        sz = s[2]
        self.m11 = self.m11*sx
        self.m12 = self.m12*sy
        self.m13 = self.m13*sz
        self.m21 = self.m21*sx
        self.m22 = self.m22*sy
        self.m23 = self.m23*sz
        self.m31 = self.m31*sx
        self.m32 = self.m32*sy
        self.m33 = self.m33*sz
        return self

    def rotate(self, angle, axis):
        """Concatenate a rotation.

        angle must be given in radians. axis must be a 3-sequence.
        """
        cdef mat3 R
        cdef double c1, c2, c3
        R=self.rotation(angle, axis)

        # self*R  (in-place)
        c1 = self.m11
        c2 = self.m12
        c3 = self.m13
        self.m11 = c1*R.m11+c2*R.m21+c3*R.m31
        self.m12 = c1*R.m12+c2*R.m22+c3*R.m32
        self.m13 = c1*R.m13+c2*R.m23+c3*R.m33

        c1 = self.m21
        c2 = self.m22
        c3 = self.m23
        self.m21 = c1*R.m11+c2*R.m21+c3*R.m31
        self.m22 = c1*R.m12+c2*R.m22+c3*R.m32
        self.m23 = c1*R.m13+c2*R.m23+c3*R.m33

        c1 = self.m31
        c2 = self.m32
        c3 = self.m33
        self.m31 = c1*R.m11+c2*R.m21+c3*R.m31
        self.m32 = c1*R.m12+c2*R.m22+c3*R.m32
        self.m33 = c1*R.m13+c2*R.m23+c3*R.m33

        return self

    def ortho(self):
        """Return a matrix with orthogonal base vectors.

        Makes the x-, y- and z-axis orthogonal.
        """
        cdef mat3 res
        cdef vec3 x,y,z

        x = vec3(self.m11, self.m21, self.m31)
        y = vec3(self.m12, self.m22, self.m32)
        z = vec3(self.m13, self.m23, self.m33)

        xl = x.length()
        xl = xl*xl
        y = y - ((x*y)/xl)*x
        z = z - ((x*z)/xl)*x

        yl = y.length()
        yl = yl*yl
        z = z - ((y*z)/yl)*y

        res = mat3()
        res.m11 = x.x
        res.m12 = y.x
        res.m13 = z.x
        res.m21 = x.y
        res.m22 = y.y
        res.m23 = z.y
        res.m31 = x.z
        res.m32 = y.z
        res.m33 = z.z
        return res

    def decompose(self):
        """Decomposes the matrix into a rotation and scaling part.

        Returns a tuple (rotation, scaling). The scaling part is given
        as a vec3, the rotation is still a mat3.
        """
        cdef mat3 m
        cdef vec3 a,b,c
        cdef double al, bl, cl

        m = self.ortho()

        # a,b,c = Column 0,1,2 of m
        a = vec3()
        b = vec3()
        c = vec3()
        a.x = m.m11
        a.y = m.m21
        a.z = m.m31
        b.x = m.m12
        b.y = m.m22
        b.z = m.m32
        c.x = m.m13
        c.y = m.m23
        c.z = m.m33

        # al,bl,cl = length of a,b,c (= scaling factors)
        al = sqrt(a.x*a.x + a.y*a.y + a.z*a.z)
        bl = sqrt(b.x*b.x + b.y*b.y + b.z*b.z)
        cl = sqrt(c.x*c.x + c.y*c.y + c.z*c.z)
        if al<=eps or bl<=eps or cl<=eps:
            raise ZeroDivisionError,"transformation contains a 0 scaling"

        scale = vec3(al,bl,cl)

        # normalizing a,b,c
#        a/=al
#        b/=bl
#        c/=cl
        a.x = a.x/al
        a.y = a.y/al
        a.z = a.z/al
        b.x = b.x/bl
        b.y = b.y/bl
        b.z = b.z/bl
        c.x = c.x/cl
        c.y = c.y/cl
        c.z = c.z/cl

        m.m11 = a.x
        m.m21 = a.y
        m.m31 = a.z
        m.m12 = b.x
        m.m22 = b.y
        m.m32 = b.z
        m.m13 = c.x
        m.m23 = c.y
        m.m33 = c.z
        if m.determinant()<0.0:
            m.m11 = -m.m11
            m.m21 = -m.m21
            m.m31 = -m.m31
            scale.x = -scale.x

        return (m, scale)


######################################################################

# mat4iter
cdef class mat4iter:
    """mat iterator."""
    cdef int idx
    cdef object m

    def __cinit__(self, m):
        self.idx = 0
        self.m   = m

    def __iter__(self):
        return self

    def __next__(self):
        if self.idx==0:
            self.idx=1
            return self.m[0]
        elif self.idx==1:
            self.idx=2
            return self.m[1]
        elif self.idx==2:
            self.idx=3
            return self.m[2]
        elif self.idx==3:
            self.idx=4
            return self.m[3]
        else:
            raise StopIteration

# Calculates m3.determinant() where m3 is a mat3 that is based on m4
# where row i and column j are left out.
cdef double _subdet(double m4[16], short i, short j):
    cdef double m3[9]
    cdef double* p
    cdef short k, l

    p = m3
    k=0
    while 1:
        if k==i:
            k=k+1
        if k==4:
            break
        l=0
        while 1:
            if l==j:
                l=l+1
            if l==4:
                break
            p[0]=m4[k*4+l]
            p=p+1
            l=l+1
        k=k+1

    return m3[0]*m3[4]*m3[8]+ \
           m3[1]*m3[5]*m3[6]+ \
           m3[2]*m3[3]*m3[7]- \
           m3[6]*m3[4]*m3[2]- \
           m3[7]*m3[5]*m3[0]- \
           m3[8]*m3[3]*m3[1]



# mat4
cdef class mat4:

    cdef double m11, m12, m13, m14
    cdef double m21, m22, m23, m24
    cdef double m31, m32, m33, m34
    cdef double m41, m42, m43, m44

    def __cinit__(self, *args):
        cdef int arglen
        cdef mat4 B
        arglen = len(args)

        # No arguments
        if arglen==0:
            m11=m12=m13=m14=0.0
            m21=m22=m23=m24=0.0
            m31=m32=m33=m34=0.0
            m41=m42=m43=m44=0.0

        # 1 argument (list, scalar or mat4)
        elif arglen==1:
            T = type(args[0])
            # scalar
            if T==float or T==int or T==long:
                self.m11 = args[0]
                self.m12 = 0.0
                self.m13 = 0.0
                self.m14 = 0.0
                self.m21 = 0.0
                self.m22 = self.m11
                self.m23 = 0.0
                self.m24 = 0.0
                self.m31 = 0.0
                self.m32 = 0.0
                self.m33 = self.m11
                self.m34 = 0.0
                self.m41 = 0.0
                self.m42 = 0.0
                self.m43 = 0.0
                self.m44 = self.m11
            # mat4
            elif T==mat4:
                B = args[0]
                self.m11 = B.m11
                self.m12 = B.m12
                self.m13 = B.m13
                self.m14 = B.m14
                self.m21 = B.m21
                self.m22 = B.m22
                self.m23 = B.m23
                self.m24 = B.m24
                self.m31 = B.m31
                self.m32 = B.m32
                self.m33 = B.m33
                self.m34 = B.m34
                self.m41 = B.m41
                self.m42 = B.m42
                self.m43 = B.m43
                self.m44 = B.m44
            # String
            elif T==str or T==unicode:
                s=args[0].replace(","," ").replace("  "," ").strip().split(" ")
                if s==[""]:
                    s=[]
                for i in range(len(s)):
                    try:
                        s[i] = float(s[i])
                    except:
                        raise ValueError,"mat4() arg is a malformed string"
                B = mat4(s)
                self.m11 = B.m11
                self.m12 = B.m12
                self.m13 = B.m13
                self.m14 = B.m14
                self.m21 = B.m21
                self.m22 = B.m22
                self.m23 = B.m23
                self.m24 = B.m24
                self.m31 = B.m31
                self.m32 = B.m32
                self.m33 = B.m33
                self.m34 = B.m34
                self.m41 = B.m41
                self.m42 = B.m42
                self.m43 = B.m43
                self.m44 = B.m44
            # Sequence of 16 floats or 4 4-sequences
            else:
                seq = args[0]
                try:
                    seqlen = len(seq)
                except:
                    raise TypeError,"mat4() arg can't be converted to mat4"

                if seqlen==16:
                    self.m11 = seq[0]
                    self.m12 = seq[1]
                    self.m13 = seq[2]
                    self.m14 = seq[3]
                    self.m21 = seq[4]
                    self.m22 = seq[5]
                    self.m23 = seq[6]
                    self.m24 = seq[7]
                    self.m31 = seq[8]
                    self.m32 = seq[9]
                    self.m33 = seq[10]
                    self.m34 = seq[11]
                    self.m41 = seq[12]
                    self.m42 = seq[13]
                    self.m43 = seq[14]
                    self.m44 = seq[15]
                elif seqlen==4:
                    b = seq[0]
                    self.m11 = b[0]
                    self.m21 = b[1]
                    self.m31 = b[2]
                    self.m41 = b[3]
                    b = seq[1]
                    self.m12 = b[0]
                    self.m22 = b[1]
                    self.m32 = b[2]
                    self.m42 = b[3]
                    b = seq[2]
                    self.m13 = b[0]
                    self.m23 = b[1]
                    self.m33 = b[2]
                    self.m43 = b[3]
                    b = seq[3]
                    self.m14 = b[0]
                    self.m24 = b[1]
                    self.m34 = b[2]
                    self.m44 = b[3]
                else:
                    raise TypeError, "mat4() takes at most 16 arguments"

        # 4 arguments (columns as sequences)
        elif arglen==4:
            a,b,c,d=args
            self.m11 = a[0]
            self.m12 = b[0]
            self.m13 = c[0]
            self.m14 = d[0]
            self.m21 = a[1]
            self.m22 = b[1]
            self.m23 = c[1]
            self.m24 = d[1]
            self.m31 = a[2]
            self.m32 = b[2]
            self.m33 = c[2]
            self.m34 = d[2]
            self.m41 = a[3]
            self.m42 = b[3]
            self.m43 = c[3]
            self.m44 = d[3]

        # 16 arguments
        elif arglen==16:
            self.m11 = args[0]
            self.m12 = args[1]
            self.m13 = args[2]
            self.m14 = args[3]
            self.m21 = args[4]
            self.m22 = args[5]
            self.m23 = args[6]
            self.m24 = args[7]
            self.m31 = args[8]
            self.m32 = args[9]
            self.m33 = args[10]
            self.m34 = args[11]
            self.m41 = args[12]
            self.m42 = args[13]
            self.m43 = args[14]
            self.m44 = args[15]

        else:
            raise TypeError,"mat4() arg can't be converted to mat4"

    def __reduce__(self):
        return (_mat4construct, (self[0], self[1], self[2], self[3]))

    def __repr__(self):
        return 'mat4(%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)'%(
            repr(self.m11),repr(self.m12),repr(self.m13),repr(self.m14),
            repr(self.m21),repr(self.m22),repr(self.m23),repr(self.m24),
            repr(self.m31),repr(self.m32),repr(self.m33),repr(self.m34),
            repr(self.m41),repr(self.m42),repr(self.m43),repr(self.m44))

    def __str__(self):
        return "[%9.4f, %9.4f, %9.4f, %9.4f]\n[%9.4f, %9.4f, %9.4f, %9.4f]\n[%9.4f, %9.4f, %9.4f, %9.4f]\n[%9.4f, %9.4f, %9.4f, %9.4f]"%(
            self.m11, self.m12, self.m13, self.m14,
            self.m21, self.m22, self.m23, self.m24,
            self.m31, self.m32, self.m33, self.m34,
            self.m41, self.m42, self.m43, self.m44)

    def __iter__(self):
        return mat4iter(self)

    def __richcmp__(a,b,op):
        cdef mat4 ma, mb
        cdef int cop

        cop = op

        ta = type(a)
        tb = type(b)
        if (ta!=mat4 or tb!=mat4):
            if cop==3:
                return 1
            else:
                return 0

        ma = a
        mb = b

        # <
        if cop==0:
            return 0
        # <=
        elif cop==1:
            return 0
        # ==
        elif cop==2:
            return (fabs(ma.m11-mb.m11)<=eps and fabs(ma.m12-mb.m12)<=eps and
                    fabs(ma.m13-mb.m13)<=eps and fabs(ma.m14-mb.m14)<=eps and
                    fabs(ma.m21-mb.m21)<=eps and fabs(ma.m22-mb.m22)<=eps and
                    fabs(ma.m23-mb.m23)<=eps and fabs(ma.m24-mb.m24)<=eps and
                    fabs(ma.m31-mb.m31)<=eps and fabs(ma.m32-mb.m32)<=eps and
                    fabs(ma.m33-mb.m33)<=eps and fabs(ma.m34-mb.m34)<=eps and
                    fabs(ma.m41-mb.m41)<=eps and fabs(ma.m42-mb.m42)<=eps and
                    fabs(ma.m43-mb.m43)<=eps and fabs(ma.m44-mb.m44)<=eps)
        # !=
        elif cop==3:
            return (fabs(ma.m11-mb.m11)>eps or fabs(ma.m12-mb.m12)>eps or
                    fabs(ma.m13-mb.m13)>eps or fabs(ma.m14-mb.m14)>eps or
                    fabs(ma.m21-mb.m21)>eps or fabs(ma.m22-mb.m22)>eps or
                    fabs(ma.m23-mb.m23)>eps or fabs(ma.m24-mb.m24)>eps or
                    fabs(ma.m31-mb.m31)>eps or fabs(ma.m32-mb.m32)>eps or
                    fabs(ma.m33-mb.m33)>eps or fabs(ma.m34-mb.m34)>eps or
                    fabs(ma.m41-mb.m41)>eps or fabs(ma.m42-mb.m42)>eps or
                    fabs(ma.m43-mb.m43)>eps or fabs(ma.m44-mb.m44)>eps)
        # >
        elif cop==4:
            return 0
        # >=
        elif cop==5:
            return 0
        # sonst (Fehler)
        else:
            raise ValueError,"internal error: illegal rich comparison number"


    def __add__(mat4 a, mat4 b):
        """Matrix addition.

        >>> M=mat4(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)
        >>> print M+M
        [   2.0000,    4.0000,    6.0000,    8.0000]
        [  10.0000,   12.0000,   14.0000,   16.0000]
        [  18.0000,   20.0000,   22.0000,   24.0000]
        [  26.0000,   28.0000,   30.0000,   32.0000]
        """
        cdef mat4 res
        res = mat4()
        res.m11 = a.m11 + b.m11
        res.m12 = a.m12 + b.m12
        res.m13 = a.m13 + b.m13
        res.m14 = a.m14 + b.m14
        res.m21 = a.m21 + b.m21
        res.m22 = a.m22 + b.m22
        res.m23 = a.m23 + b.m23
        res.m24 = a.m24 + b.m24
        res.m31 = a.m31 + b.m31
        res.m32 = a.m32 + b.m32
        res.m33 = a.m33 + b.m33
        res.m34 = a.m34 + b.m34
        res.m41 = a.m41 + b.m41
        res.m42 = a.m42 + b.m42
        res.m43 = a.m43 + b.m43
        res.m44 = a.m44 + b.m44
        return res

    def __sub__(mat4 a, mat4 b):
        """Matrix subtraction.

        >>> M=mat4(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)
        >>> print M-M
        [   0.0000,    0.0000,    0.0000,    0.0000]
        [   0.0000,    0.0000,    0.0000,    0.0000]
        [   0.0000,    0.0000,    0.0000,    0.0000]
        [   0.0000,    0.0000,    0.0000,    0.0000]
        """
        cdef mat4 res
        res = mat4()
        res.m11 = a.m11 - b.m11
        res.m12 = a.m12 - b.m12
        res.m13 = a.m13 - b.m13
        res.m14 = a.m14 - b.m14
        res.m21 = a.m21 - b.m21
        res.m22 = a.m22 - b.m22
        res.m23 = a.m23 - b.m23
        res.m24 = a.m24 - b.m24
        res.m31 = a.m31 - b.m31
        res.m32 = a.m32 - b.m32
        res.m33 = a.m33 - b.m33
        res.m34 = a.m34 - b.m34
        res.m41 = a.m41 - b.m41
        res.m42 = a.m42 - b.m42
        res.m43 = a.m43 - b.m43
        res.m44 = a.m44 - b.m44
        return res

    def __mul__(a, b):
        """Multiplication.

        >>> M=mat4(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)
        >>> print M*2.0
        [   2.0000,    4.0000,    6.0000,    8.0000]
        [  10.0000,   12.0000,   14.0000,   16.0000]
        [  18.0000,   20.0000,   22.0000,   24.0000]
        [  26.0000,   28.0000,   30.0000,   32.0000]
        >>> print 2.0*M
        [   2.0000,    4.0000,    6.0000,    8.0000]
        [  10.0000,   12.0000,   14.0000,   16.0000]
        [  18.0000,   20.0000,   22.0000,   24.0000]
        [  26.0000,   28.0000,   30.0000,   32.0000]
        >>> print M*M
        [  90.0000,  100.0000,  110.0000,  120.0000]
        [ 202.0000,  228.0000,  254.0000,  280.0000]
        [ 314.0000,  356.0000,  398.0000,  440.0000]
        [ 426.0000,  484.0000,  542.0000,  600.0000]
        >>> print M*_vec3(1,2,3)
        (0.1765, 0.4510, 0.7255)
        >>> print _vec3(1,2,3)*M
        (0.7083, 0.8056, 0.9028)
        """
        cdef mat4 ma, mb, res
        cdef vec3 va, vb, vres
        cdef vec4 wa, wb, wres
        cdef double r, w

        ta = type(a)
        tb = type(b)

        # mat4 on the left?
        if ta==mat4:
            ma = a
            # mat4 on the right?
            if tb==mat4:
                # mat4*mat4
                res = mat4()
                mb = b
                res.m11 = ma.m11*mb.m11+ma.m12*mb.m21+ma.m13*mb.m31+ma.m14*mb.m41
                res.m12 = ma.m11*mb.m12+ma.m12*mb.m22+ma.m13*mb.m32+ma.m14*mb.m42
                res.m13 = ma.m11*mb.m13+ma.m12*mb.m23+ma.m13*mb.m33+ma.m14*mb.m43
                res.m14 = ma.m11*mb.m14+ma.m12*mb.m24+ma.m13*mb.m34+ma.m14*mb.m44
                res.m21 = ma.m21*mb.m11+ma.m22*mb.m21+ma.m23*mb.m31+ma.m24*mb.m41
                res.m22 = ma.m21*mb.m12+ma.m22*mb.m22+ma.m23*mb.m32+ma.m24*mb.m42
                res.m23 = ma.m21*mb.m13+ma.m22*mb.m23+ma.m23*mb.m33+ma.m24*mb.m43
                res.m24 = ma.m21*mb.m14+ma.m22*mb.m24+ma.m23*mb.m34+ma.m24*mb.m44
                res.m31 = ma.m31*mb.m11+ma.m32*mb.m21+ma.m33*mb.m31+ma.m34*mb.m41
                res.m32 = ma.m31*mb.m12+ma.m32*mb.m22+ma.m33*mb.m32+ma.m34*mb.m42
                res.m33 = ma.m31*mb.m13+ma.m32*mb.m23+ma.m33*mb.m33+ma.m34*mb.m43
                res.m34 = ma.m31*mb.m14+ma.m32*mb.m24+ma.m33*mb.m34+ma.m34*mb.m44
                res.m41 = ma.m41*mb.m11+ma.m42*mb.m21+ma.m43*mb.m31+ma.m44*mb.m41
                res.m42 = ma.m41*mb.m12+ma.m42*mb.m22+ma.m43*mb.m32+ma.m44*mb.m42
                res.m43 = ma.m41*mb.m13+ma.m42*mb.m23+ma.m43*mb.m33+ma.m44*mb.m43
                res.m44 = ma.m41*mb.m14+ma.m42*mb.m24+ma.m43*mb.m34+ma.m44*mb.m44
                return res
            # vec4 on the right?
            elif tb==vec4:
                # mat4*vec4
                wres = vec4()
                wb   = b
                wres.x = ma.m11*wb.x + ma.m12*wb.y + ma.m13*wb.z + ma.m14*wb.w
                wres.y = ma.m21*wb.x + ma.m22*wb.y + ma.m23*wb.z + ma.m24*wb.w
                wres.z = ma.m31*wb.x + ma.m32*wb.y + ma.m33*wb.z + ma.m34*wb.w
                wres.w = ma.m41*wb.x + ma.m42*wb.y + ma.m43*wb.z + ma.m44*wb.w
                return wres
            # vec3 on the right?
            elif tb==vec3:
                # mat4*vec3
                vres = vec3()
                vb   = b
                vres.x = ma.m11*vb.x + ma.m12*vb.y + ma.m13*vb.z + ma.m14
                vres.y = ma.m21*vb.x + ma.m22*vb.y + ma.m23*vb.z + ma.m24
                vres.z = ma.m31*vb.x + ma.m32*vb.y + ma.m33*vb.z + ma.m34
                w    = ma.m41*vb.x + ma.m42*vb.y + ma.m43*vb.z + ma.m44
                if fabs(w)>eps:
                    vres.x = vres.x/w
                    vres.y = vres.y/w
                    vres.z = vres.z/w
                return vres
            # Scalar on the right?
            elif tb==float or tb==int or tb==long:
                res = mat4()
                r = b
                res.m11 = ma.m11*r
                res.m12 = ma.m12*r
                res.m13 = ma.m13*r
                res.m14 = ma.m14*r
                res.m21 = ma.m21*r
                res.m22 = ma.m22*r
                res.m23 = ma.m23*r
                res.m24 = ma.m24*r
                res.m31 = ma.m31*r
                res.m32 = ma.m32*r
                res.m33 = ma.m33*r
                res.m34 = ma.m34*r
                res.m41 = ma.m41*r
                res.m42 = ma.m42*r
                res.m43 = ma.m43*r
                res.m44 = ma.m44*r
                return res

        # vec4 on the left? (then there must be a mat4 on the right)
        elif ta==vec4:
            # vec4*mat4
            wres = vec4()
            wa   = a
            mb   = b
            wres.x = wa.x*mb.m11 + wa.y*mb.m21 + wa.z*mb.m31 + wa.w*mb.m41
            wres.y = wa.x*mb.m12 + wa.y*mb.m22 + wa.z*mb.m32 + wa.w*mb.m42
            wres.z = wa.x*mb.m13 + wa.y*mb.m23 + wa.z*mb.m33 + wa.w*mb.m43
            wres.w = wa.x*mb.m14 + wa.y*mb.m24 + wa.z*mb.m34 + wa.w*mb.m44
            return wres

        # vec3 on the left? (then there must be a mat4 on the right)
        elif ta==vec3:
            # vec3*mat4
            vres = vec3()
            va   = a
            mb   = b
            vres.x = va.x*mb.m11 + va.y*mb.m21 + va.z*mb.m31 + mb.m41
            vres.y = va.x*mb.m12 + va.y*mb.m22 + va.z*mb.m32 + mb.m42
            vres.z = va.x*mb.m13 + va.y*mb.m23 + va.z*mb.m33 + mb.m43
            w = va.x*mb.m14 + va.y*mb.m24 + va.z*mb.m34 + mb.m44
            if fabs(w)>eps:
                vres.x = vres.x/w
                vres.y = vres.y/w
                vres.z = vres.z/w
            return vres

        # Scalar on the left? (then there must be a mat4 on the right)
        elif ta==float or ta==int or ta==long:
            res = mat4()
            r  = a
            mb = b
            res.m11 = r*mb.m11
            res.m12 = r*mb.m12
            res.m13 = r*mb.m13
            res.m14 = r*mb.m14
            res.m21 = r*mb.m21
            res.m22 = r*mb.m22
            res.m23 = r*mb.m23
            res.m24 = r*mb.m24
            res.m31 = r*mb.m31
            res.m32 = r*mb.m32
            res.m33 = r*mb.m33
            res.m34 = r*mb.m34
            res.m41 = r*mb.m41
            res.m42 = r*mb.m42
            res.m43 = r*mb.m43
            res.m44 = r*mb.m44
            return res

        # No valid combination found
        raise TypeError, "unsupported operand type for *"

    def __div__(mat4 a, b):
        """Division.

        >>> M=mat4(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)
        >>> print M/2.0
        [   0.5000,    1.0000,    1.5000,    2.0000]
        [   2.5000,    3.0000,    3.5000,    4.0000]
        [   4.5000,    5.0000,    5.5000,    6.0000]
        [   6.5000,    7.0000,    7.5000,    8.0000]
        """
        cdef mat4 res
        cdef double r

        tb = type(b)

        if tb==float or tb==int or tb==long:
            # mat4/scalar
            res = mat4()
            r = b
            if fabs(r)<=eps:
                raise ZeroDivisionError,"mat4 division"
            res.m11 = a.m11/r
            res.m12 = a.m12/r
            res.m13 = a.m13/r
            res.m14 = a.m14/r
            res.m21 = a.m21/r
            res.m22 = a.m22/r
            res.m23 = a.m23/r
            res.m24 = a.m24/r
            res.m31 = a.m31/r
            res.m32 = a.m32/r
            res.m33 = a.m33/r
            res.m34 = a.m34/r
            res.m41 = a.m41/r
            res.m42 = a.m42/r
            res.m43 = a.m43/r
            res.m44 = a.m44/r
            return res
        else:
            raise TypeError, "unsupported operand type for /"

    def __mod__(mat4 self, other):
        """Modulo.

        >>> M=mat4(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)
        >>> print M%5.0
        [   1.0000,    2.0000,    3.0000,    4.0000]
        [   0.0000,    1.0000,    2.0000,    3.0000]
        [   4.0000,    0.0000,    1.0000,    2.0000]
        [   3.0000,    4.0000,    0.0000,    1.0000]
        """
        cdef mat4 mb, res

        tb = type(other)

        if tb==float or tb==int or tb==long:
            r = other
            if fabs(r)<=eps:
                raise ZeroDivisionError,"mat4 modulo"
            res = mat4()
            res.m11 = _fmod(self.m11, r)
            res.m12 = _fmod(self.m12, r)
            res.m13 = _fmod(self.m13, r)
            res.m14 = _fmod(self.m14, r)
            res.m21 = _fmod(self.m21, r)
            res.m22 = _fmod(self.m22, r)
            res.m23 = _fmod(self.m23, r)
            res.m24 = _fmod(self.m24, r)
            res.m31 = _fmod(self.m31, r)
            res.m32 = _fmod(self.m32, r)
            res.m33 = _fmod(self.m33, r)
            res.m34 = _fmod(self.m34, r)
            res.m41 = _fmod(self.m41, r)
            res.m42 = _fmod(self.m42, r)
            res.m43 = _fmod(self.m43, r)
            res.m44 = _fmod(self.m44, r)
            return res

        elif tb==mat4:
            mb = other

            if (fabs(mb.m11)<=eps or fabs(mb.m12)<=eps or fabs(mb.m13)<=eps or
                fabs(mb.m14)<=eps or
                fabs(mb.m21)<=eps or fabs(mb.m22)<=eps or fabs(mb.m23)<=eps or
                fabs(mb.m24)<=eps or
                fabs(mb.m31)<=eps or fabs(mb.m32)<=eps or fabs(mb.m33)<=eps or
                fabs(mb.m34)<=eps or
                fabs(mb.m41)<=eps or fabs(mb.m42)<=eps or fabs(mb.m43)<=eps or
                fabs(mb.m44)<=eps):
                raise ZeroDivisionError,"mat4 modulo"

            res=mat4()
            res.m11 = _fmod(self.m11, mb.m11)
            res.m12 = _fmod(self.m12, mb.m12)
            res.m13 = _fmod(self.m13, mb.m13)
            res.m14 = _fmod(self.m14, mb.m14)
            res.m21 = _fmod(self.m21, mb.m21)
            res.m22 = _fmod(self.m22, mb.m22)
            res.m23 = _fmod(self.m23, mb.m23)
            res.m24 = _fmod(self.m24, mb.m24)
            res.m31 = _fmod(self.m31, mb.m31)
            res.m32 = _fmod(self.m32, mb.m32)
            res.m33 = _fmod(self.m33, mb.m33)
            res.m34 = _fmod(self.m34, mb.m34)
            res.m41 = _fmod(self.m41, mb.m41)
            res.m42 = _fmod(self.m42, mb.m42)
            res.m43 = _fmod(self.m43, mb.m43)
            res.m44 = _fmod(self.m44, mb.m44)
            return res

        else:
            raise TypeError, "unsupported operand type for %"


    def __neg__(self):
        """Negation.

        >>> M=mat4(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)
        >>> print -M
        [  -1.0000,   -2.0000,   -3.0000,   -4.0000]
        [  -5.0000,   -6.0000,   -7.0000,   -8.0000]
        [  -9.0000,  -10.0000,  -11.0000,  -12.0000]
        [ -13.0000,  -14.0000,  -15.0000,  -16.0000]
        """
        cdef mat4 res
        res = mat4()
        res.m11 = -self.m11
        res.m12 = -self.m12
        res.m13 = -self.m13
        res.m14 = -self.m14
        res.m21 = -self.m21
        res.m22 = -self.m22
        res.m23 = -self.m23
        res.m24 = -self.m24
        res.m31 = -self.m31
        res.m32 = -self.m32
        res.m33 = -self.m33
        res.m34 = -self.m34
        res.m41 = -self.m41
        res.m42 = -self.m42
        res.m43 = -self.m43
        res.m44 = -self.m44
        return res

    def __pos__(self):
        """
        >>> M=mat4(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)
        >>> print +M
        [   1.0000,    2.0000,    3.0000,    4.0000]
        [   5.0000,    6.0000,    7.0000,    8.0000]
        [   9.0000,   10.0000,   11.0000,   12.0000]
        [  13.0000,   14.0000,   15.0000,   16.0000]
        """
        cdef mat4 res
        res = mat4()
        res.m11 = self.m11
        res.m12 = self.m12
        res.m13 = self.m13
        res.m14 = self.m14
        res.m21 = self.m21
        res.m22 = self.m22
        res.m23 = self.m23
        res.m24 = self.m24
        res.m31 = self.m31
        res.m32 = self.m32
        res.m33 = self.m33
        res.m34 = self.m34
        res.m41 = self.m41
        res.m42 = self.m42
        res.m43 = self.m43
        res.m44 = self.m44
        return res

    def __iadd__(self, mat4 other):
        """Inline matrix addition.
        """
        self.m11 = self.m11 + other.m11
        self.m12 = self.m12 + other.m12
        self.m13 = self.m13 + other.m13
        self.m14 = self.m14 + other.m14
        self.m21 = self.m21 + other.m21
        self.m22 = self.m22 + other.m22
        self.m23 = self.m23 + other.m23
        self.m24 = self.m24 + other.m24
        self.m31 = self.m31 + other.m31
        self.m32 = self.m32 + other.m32
        self.m33 = self.m33 + other.m33
        self.m34 = self.m34 + other.m34
        self.m41 = self.m41 + other.m41
        self.m42 = self.m42 + other.m42
        self.m43 = self.m43 + other.m43
        self.m44 = self.m44 + other.m44
        return self

    def __isub__(self, mat4 other):
        """Inline matrix subtraction.
        """
        self.m11 = self.m11 - other.m11
        self.m12 = self.m12 - other.m12
        self.m13 = self.m13 - other.m13
        self.m14 = self.m14 - other.m14
        self.m21 = self.m21 - other.m21
        self.m22 = self.m22 - other.m22
        self.m23 = self.m23 - other.m23
        self.m24 = self.m24 - other.m24
        self.m31 = self.m31 - other.m31
        self.m32 = self.m32 - other.m32
        self.m33 = self.m33 - other.m33
        self.m34 = self.m34 - other.m34
        self.m41 = self.m41 - other.m41
        self.m42 = self.m42 - other.m42
        self.m43 = self.m43 - other.m43
        self.m44 = self.m44 - other.m44
        return self

    def __imul__(self, other):
        """Inline multiplication (scalar or matrix).
        """
        cdef double r, a,b,c,d
        cdef mat4 mb

        tb = type(other)
        if tb==float or tb==int or tb==long:
            r=other
            self.m11 = self.m11*r
            self.m12 = self.m12*r
            self.m13 = self.m13*r
            self.m14 = self.m14*r
            self.m21 = self.m21*r
            self.m22 = self.m22*r
            self.m23 = self.m23*r
            self.m24 = self.m24*r
            self.m31 = self.m31*r
            self.m32 = self.m32*r
            self.m33 = self.m33*r
            self.m34 = self.m34*r
            self.m41 = self.m41*r
            self.m42 = self.m42*r
            self.m43 = self.m43*r
            self.m44 = self.m44*r
            return self
        elif tb==mat4:
            mb = other
            a = self.m11
            b = self.m12
            c = self.m13
            d = self.m14
            self.m11 = a*mb.m11+b*mb.m21+c*mb.m31+d*mb.m41
            self.m12 = a*mb.m12+b*mb.m22+c*mb.m32+d*mb.m42
            self.m13 = a*mb.m13+b*mb.m23+c*mb.m33+d*mb.m43
            self.m14 = a*mb.m14+b*mb.m24+c*mb.m34+d*mb.m44
            a = self.m21
            b = self.m22
            c = self.m23
            d = self.m24
            self.m21 = a*mb.m11+b*mb.m21+c*mb.m31+d*mb.m41
            self.m22 = a*mb.m12+b*mb.m22+c*mb.m32+d*mb.m42
            self.m23 = a*mb.m13+b*mb.m23+c*mb.m33+d*mb.m43
            self.m24 = a*mb.m14+b*mb.m24+c*mb.m34+d*mb.m44
            a = self.m31
            b = self.m32
            c = self.m33
            d = self.m34
            self.m31 = a*mb.m11+b*mb.m21+c*mb.m31+d*mb.m41
            self.m32 = a*mb.m12+b*mb.m22+c*mb.m32+d*mb.m42
            self.m33 = a*mb.m13+b*mb.m23+c*mb.m33+d*mb.m43
            self.m34 = a*mb.m14+b*mb.m24+c*mb.m34+d*mb.m44
            a = self.m41
            b = self.m42
            c = self.m43
            d = self.m44
            self.m41 = a*mb.m11+b*mb.m21+c*mb.m31+d*mb.m41
            self.m42 = a*mb.m12+b*mb.m22+c*mb.m32+d*mb.m42
            self.m43 = a*mb.m13+b*mb.m23+c*mb.m33+d*mb.m43
            self.m44 = a*mb.m14+b*mb.m24+c*mb.m34+d*mb.m44
            return self
        else:
            raise TypeError, "unsupported operand type for *="

    def __idiv__(self, other):
        """Inline division with scalar.
        """
        cdef double r
        tb = type(other)
        if tb==float or tb==int or tb==long:
            r=other
            if fabs(r)<=eps:
                raise ZeroDivisionError,"mat4 division"
            self.m11 = self.m11/r
            self.m12 = self.m12/r
            self.m13 = self.m13/r
            self.m14 = self.m14/r
            self.m21 = self.m21/r
            self.m22 = self.m22/r
            self.m23 = self.m23/r
            self.m24 = self.m24/r
            self.m31 = self.m31/r
            self.m32 = self.m32/r
            self.m33 = self.m33/r
            self.m34 = self.m34/r
            self.m41 = self.m41/r
            self.m42 = self.m42/r
            self.m43 = self.m43/r
            self.m44 = self.m44/r
            return self
        else:
            raise TypeError, "unsupported operand type for /="

    def __imod__(self, b):
        """Inline modulo.
        """
        cdef mat4 mb
        cdef double r

        tb = type(b)

        if tb==float or tb==int or tb==long:
            r = b
            if fabs(r)<=eps:
                raise ZeroDivisionError,"mat4 modulo"
            self.m11 = _fmod(self.m11, r)
            self.m12 = _fmod(self.m12, r)
            self.m13 = _fmod(self.m13, r)
            self.m14 = _fmod(self.m14, r)
            self.m21 = _fmod(self.m21, r)
            self.m22 = _fmod(self.m22, r)
            self.m23 = _fmod(self.m23, r)
            self.m24 = _fmod(self.m24, r)
            self.m31 = _fmod(self.m31, r)
            self.m32 = _fmod(self.m32, r)
            self.m33 = _fmod(self.m33, r)
            self.m34 = _fmod(self.m34, r)
            self.m41 = _fmod(self.m41, r)
            self.m42 = _fmod(self.m42, r)
            self.m43 = _fmod(self.m43, r)
            self.m44 = _fmod(self.m44, r)
            return self
        elif tb==mat4:
            mb = b
            if (fabs(mb.m11)<=eps or fabs(mb.m12)<=eps or fabs(mb.m13)<=eps or
                fabs(mb.m14)<=eps or
                fabs(mb.m21)<=eps or fabs(mb.m22)<=eps or fabs(mb.m23)<=eps or
                fabs(mb.m24)<=eps or
                fabs(mb.m31)<=eps or fabs(mb.m32)<=eps or fabs(mb.m33)<=eps or
                fabs(mb.m34)<=eps or
                fabs(mb.m41)<=eps or fabs(mb.m42)<=eps or fabs(mb.m43)<=eps or
                fabs(mb.m44)<=eps):
                raise ZeroDivisionError,"mat4 modulo"
            self.m11 = _fmod(self.m11, mb.m11)
            self.m12 = _fmod(self.m12, mb.m12)
            self.m13 = _fmod(self.m13, mb.m13)
            self.m14 = _fmod(self.m14, mb.m13)
            self.m21 = _fmod(self.m21, mb.m21)
            self.m22 = _fmod(self.m22, mb.m22)
            self.m23 = _fmod(self.m23, mb.m23)
            self.m24 = _fmod(self.m24, mb.m23)
            self.m31 = _fmod(self.m31, mb.m31)
            self.m32 = _fmod(self.m32, mb.m32)
            self.m33 = _fmod(self.m33, mb.m33)
            self.m34 = _fmod(self.m34, mb.m33)
            self.m41 = _fmod(self.m41, mb.m41)
            self.m42 = _fmod(self.m42, mb.m42)
            self.m43 = _fmod(self.m43, mb.m43)
            self.m44 = _fmod(self.m44, mb.m43)
            return self
        else:
            raise TypeError, "unsupported operand type for %="

    def __len__(self):
        """Length of the sequence (always 4)."""
        return 4

    def __getitem__(mat4 self, key):
        """Return a column or a single matrix element.
        """
        cdef int i,j

        T=type(key)
        if T==int or T==long:
            i = key
            if i==0:
                return vec4(self.m11, self.m21, self.m31, self.m41)
            elif i==1:
                return vec4(self.m12, self.m22, self.m32, self.m42)
            elif i==2:
                return vec4(self.m13, self.m23, self.m33, self.m43)
            elif i==3:
                return vec4(self.m14, self.m24, self.m34, self.m44)
            else:
                raise IndexError, "index out of range"
        elif T==tuple:
            if len(key)!=2:
                raise ValueError, "index tuple must be a 2-tuple"
            i = key[0]
            j = key[1]
            if i<0 or i>3 or j<0 or j>3:
                raise IndexError, "index out of range"

            if i==0:
                if j==0:
                    return self.m11
                elif j==1:
                    return self.m12
                elif j==2:
                    return self.m13
                elif j==3:
                    return self.m14
            elif i==1:
                if j==0:
                    return self.m21
                elif j==1:
                    return self.m22
                elif j==2:
                    return self.m23
                elif j==3:
                    return self.m24
            elif i==2:
                if j==0:
                    return self.m31
                elif j==1:
                    return self.m32
                elif j==2:
                    return self.m33
                elif j==3:
                    return self.m34
            elif i==3:
                if j==0:
                    return self.m41
                elif j==1:
                    return self.m42
                elif j==2:
                    return self.m43
                elif j==3:
                    return self.m44
        else:
            raise TypeError,"index must be integer or 2-tuple"

    def __setitem__(self, key, value):
        """Set a column or a single matrix element.

        The column must be given as a 4-sequence (this includes vec4 as well).
        """
        cdef int i,j
        cdef double v

        T=type(key)
        if T==int or T==long:
            if len(value)!=4:
                raise ValueError, "4-sequence expected"
            i = key
            if i==0:
                self.m11 = value[0]
                self.m21 = value[1]
                self.m31 = value[2]
                self.m41 = value[3]
            elif i==1:
                self.m12 = value[0]
                self.m22 = value[1]
                self.m32 = value[2]
                self.m42 = value[3]
            elif i==2:
                self.m13 = value[0]
                self.m23 = value[1]
                self.m33 = value[2]
                self.m43 = value[3]
            elif i==3:
                self.m14 = value[0]
                self.m24 = value[1]
                self.m34 = value[2]
                self.m44 = value[3]
            else:
                raise IndexError, "index out of range"
        elif T==tuple:
            if len(key)!=2:
                raise ValueError, "index tuple must be a 2-tuple"
            i = key[0]
            j = key[1]
            v = value
            if i<0 or i>3 or j<0 or j>3:
                raise IndexError, "index out of range"

            if i==0:
                if j==0:
                    self.m11 = v
                elif j==1:
                    self.m12 = v
                elif j==2:
                    self.m13 = v
                elif j==3:
                    self.m14 = v
            elif i==1:
                if j==0:
                    self.m21 = v
                elif j==1:
                    self.m22 = v
                elif j==2:
                    self.m23 = v
                elif j==3:
                    self.m24 = v
            elif i==2:
                if j==0:
                    self.m31 = v
                elif j==1:
                    self.m32 = v
                elif j==2:
                    self.m33 = v
                elif j==3:
                    self.m34 = v
            elif i==3:
                if j==0:
                    self.m41 = v
                elif j==1:
                    self.m42 = v
                elif j==2:
                    self.m43 = v
                elif j==3:
                    self.m44 = v
        else:
            raise TypeError,"index must be integer or 2-tuple"

    def getRow(self, idx):
        """Return row (as vec4)."""
        cdef int i

        T=type(idx)
        if T!=int and T!=long:
            raise TypeError,"index must be integer"

        i = idx
        if i==0:
            return vec4(self.m11, self.m12, self.m13, self.m14)
        elif i==1:
            return vec4(self.m21, self.m22, self.m23, self.m24)
        elif i==2:
            return vec4(self.m31, self.m32, self.m33, self.m34)
        elif i==3:
            return vec4(self.m41, self.m42, self.m43, self.m44)
        else:
            raise IndexError, "index out of range"

    def setRow(self, idx, value):
        """Set row."""
        cdef int i

        T=type(idx)
        if T!=int and T!=long:
            raise TypeError,"index must be integer"

        if len(value)!=4:
            raise ValueError, "4-sequence expected"

        i = idx
        if i==0:
            self.m11 = value[0]
            self.m12 = value[1]
            self.m13 = value[2]
            self.m14 = value[3]
        elif i==1:
            self.m21 = value[0]
            self.m22 = value[1]
            self.m23 = value[2]
            self.m24 = value[3]
        elif i==2:
            self.m31 = value[0]
            self.m32 = value[1]
            self.m33 = value[2]
            self.m34 = value[3]
        elif i==3:
            self.m41 = value[0]
            self.m42 = value[1]
            self.m43 = value[2]
            self.m44 = value[3]
        else:
            raise IndexError, "index out of range"

    def getColumn(self, idx):
        """Return column (as vec4)."""
        cdef int i

        T=type(idx)
        if T!=int and T!=long:
            raise TypeError,"index must be integer"

        i = idx
        if i==0:
            return vec4(self.m11, self.m21, self.m31, self.m41)
        elif i==1:
            return vec4(self.m12, self.m22, self.m32, self.m42)
        elif i==2:
            return vec4(self.m13, self.m23, self.m33, self.m43)
        elif i==3:
            return vec4(self.m14, self.m24, self.m34, self.m44)
        else:
            raise IndexError, "index out of range"

    def setColumn(self, idx, value):
        """Set column."""
        cdef int i

        T=type(idx)
        if T!=int and T!=long:
            raise TypeError,"index must be integer"

        if len(value)!=4:
            raise ValueError, "4-sequence expected"

        i = idx
        if i==0:
            self.m11 = value[0]
            self.m21 = value[1]
            self.m31 = value[2]
            self.m41 = value[3]
        elif i==1:
            self.m12 = value[0]
            self.m22 = value[1]
            self.m32 = value[2]
            self.m42 = value[3]
        elif i==2:
            self.m13 = value[0]
            self.m23 = value[1]
            self.m33 = value[2]
            self.m43 = value[3]
        elif i==3:
            self.m14 = value[0]
            self.m24 = value[1]
            self.m34 = value[2]
            self.m44 = value[3]
        else:
            raise IndexError, "index out of range"


    def toList(self, rowmajor=0):
        """Return a list containing the matrix elements.

        By default the list is in column-major order (which can directly be
        used in OpenGL or RenderMan). If you set the optional argument
        rowmajor to 1, you'll get the list in row-major order.

        >>> M=mat4(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)
        >>> print M.toList()
        [1, 5, 9, 13, 2, 6, 10, 14, 3, 7, 11, 15, 4, 8, 12, 16]
        >>> print M.toList(rowmajor=1)
        [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
        """
        if rowmajor:
            return [self.m11, self.m12, self.m13, self.m14,
                    self.m21, self.m22, self.m23, self.m24,
                    self.m31, self.m32, self.m33, self.m34,
                    self.m41, self.m42, self.m43, self.m44]
        else:
            return [self.m11, self.m21, self.m31, self.m41,
                    self.m12, self.m22, self.m32, self.m42,
                    self.m13, self.m23, self.m33, self.m43,
                    self.m14, self.m24, self.m34, self.m44]


    def identity(self):
        """Return identity matrix.

        >>> print mat4().identity()
        [   1.0000,    0.0000,    0.0000,    0.0000]
        [   0.0000,    1.0000,    0.0000,    0.0000]
        [   0.0000,    0.0000,    1.0000,    0.0000]
        [   0.0000,    0.0000,    0.0000,    1.0000]
        """
        cdef mat4 res
        res = mat4()
        res.m11 = 1.0
        res.m12 = 0.0
        res.m13 = 0.0
        res.m14 = 0.0
        res.m21 = 0.0
        res.m22 = 1.0
        res.m23 = 0.0
        res.m24 = 0.0
        res.m31 = 0.0
        res.m32 = 0.0
        res.m33 = 1.0
        res.m34 = 0.0
        res.m41 = 0.0
        res.m42 = 0.0
        res.m43 = 0.0
        res.m44 = 1.0
        return res

    def transpose(self):
        """Return transpose matrix.

        >>> M=mat4(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)
        >>> print M.transpose()
        [   1.0000,    5.0000,    9.0000,   13.0000]
        [   2.0000,    6.0000,   10.0000,   14.0000]
        [   3.0000,    7.0000,   11.0000,   15.0000]
        [   4.0000,    8.0000,   12.0000,   16.0000]
        """
        cdef mat4 res
        res = mat4()
        res.m11 = self.m11
        res.m12 = self.m21
        res.m13 = self.m31
        res.m14 = self.m41
        res.m21 = self.m12
        res.m22 = self.m22
        res.m23 = self.m32
        res.m24 = self.m42
        res.m31 = self.m13
        res.m32 = self.m23
        res.m33 = self.m33
        res.m34 = self.m43
        res.m41 = self.m14
        res.m42 = self.m24
        res.m43 = self.m34
        res.m44 = self.m44
        return res

    def determinant(self):
        """Return determinant.

        >>> M=mat4(2.0,0,0,0, 0,2.0,0,0, 0,0,2.0,0, 0,0,0,2.0)
        >>> print M.determinant()
        16.0
        """
        return self.m11*self.m22*self.m33*self.m44 \
               -self.m11*self.m22*self.m34*self.m43 \
               +self.m11*self.m23*self.m34*self.m42 \
               -self.m11*self.m23*self.m32*self.m44 \
               +self.m11*self.m24*self.m32*self.m43 \
               -self.m11*self.m24*self.m33*self.m42 \
               -self.m12*self.m23*self.m34*self.m41 \
               +self.m12*self.m23*self.m31*self.m44 \
               -self.m12*self.m24*self.m31*self.m43 \
               +self.m12*self.m24*self.m33*self.m41 \
               -self.m12*self.m21*self.m33*self.m44 \
               +self.m12*self.m21*self.m34*self.m43 \
               +self.m13*self.m24*self.m31*self.m42 \
               -self.m13*self.m24*self.m32*self.m41 \
               +self.m13*self.m21*self.m32*self.m44 \
               -self.m13*self.m21*self.m34*self.m42 \
               +self.m13*self.m22*self.m34*self.m41 \
               -self.m13*self.m22*self.m31*self.m44 \
               -self.m14*self.m21*self.m32*self.m43 \
               +self.m14*self.m21*self.m33*self.m42 \
               -self.m14*self.m22*self.m33*self.m41 \
               +self.m14*self.m22*self.m31*self.m43 \
               -self.m14*self.m23*self.m31*self.m42 \
               +self.m14*self.m23*self.m32*self.m41

    def _submat(self, mat3 M3, i,j):
        cdef int k,l
        cdef double v

#        M=mat3()
        k=0
        while k<3:
            l=0
            while l<3:
                t=(k,l)
                if k>=i:
                    t=(k+1,t[1])
                if l>=j:
                    t=(t[0],l+1)

#                M3[k,l] = self[t]
                v = self[t]
                if k==0 and l==0:
                    M3.m11 = v
                if k==0 and l==1:
                    M3.m12 = v
                if k==0 and l==2:
                    M3.m13 = v
                if k==1 and l==0:
                    M3.m21 = v
                if k==1 and l==1:
                    M3.m22 = v
                if k==1 and l==2:
                    M3.m23 = v
                if k==2 and l==0:
                    M3.m31 = v
                if k==2 and l==1:
                    M3.m32 = v
                if k==2 and l==2:
                    M3.m33 = v

                l=l+1
            k=k+1
        return M3

    def inverse(self):
        """Return inverse matrix.

        >>> M=mat4(0,-2.0,0,0, 2.0,0,0,0, 0,0,2,0, 0,0,0,2)
        >>> print M.inverse()
        [   0.0000,    0.5000,    0.0000,    0.0000]
        [  -0.5000,    0.0000,    0.0000,    0.0000]
        [   0.0000,    0.0000,    0.5000,    0.0000]
        [   0.0000,    0.0000,    0.0000,    0.5000]
        """
        cdef mat4 res
        cdef short i,j, sign
        cdef double m4[16]
        cdef double m4_res[16]
        cdef double det, d

        m4[0] = self.m11
        m4[1] = self.m12
        m4[2] = self.m13
        m4[3] = self.m14
        m4[4] = self.m21
        m4[5] = self.m22
        m4[6] = self.m23
        m4[7] = self.m24
        m4[8] = self.m31
        m4[9] = self.m32
        m4[10] = self.m33
        m4[11] = self.m34
        m4[12] = self.m41
        m4[13] = self.m42
        m4[14] = self.m43
        m4[15] = self.m44

        det=self.determinant()
        if fabs(det)<=eps:
            raise ZeroDivisionError,"matrix not invertible"

        i=0
        while i<4:
            j=0
            while j<4:
                sign=1-((i+j)%2)*2
                d = _subdet(m4,i,j)
                m4_res[j*4+i]=sign*d/det
                j=j+1
            i=i+1

        res = mat4()
        res.m11 = m4_res[0]
        res.m12 = m4_res[1]
        res.m13 = m4_res[2]
        res.m14 = m4_res[3]
        res.m21 = m4_res[4]
        res.m22 = m4_res[5]
        res.m23 = m4_res[6]
        res.m24 = m4_res[7]
        res.m31 = m4_res[8]
        res.m32 = m4_res[9]
        res.m33 = m4_res[10]
        res.m34 = m4_res[11]
        res.m41 = m4_res[12]
        res.m42 = m4_res[13]
        res.m43 = m4_res[14]
        res.m44 = m4_res[15]
        return res


    def translation(self, t):
        """Return translation matrix.

        t must be a 3-sequence (e.g. a vec3).
        """
        cdef mat4 res
        res = mat4()
        res.m11 = 1.0
        res.m12 = 0.0
        res.m13 = 0.0
        res.m14 = t[0]
        res.m21 = 0.0
        res.m22 = 1.0
        res.m23 = 0.0
        res.m24 = t[1]
        res.m31 = 0.0
        res.m32 = 0.0
        res.m33 = 1.0
        res.m34 = t[2]
        res.m41 = 0.0
        res.m42 = 0.0
        res.m43 = 0.0
        res.m44 = 1.0
        return res

    def scaling(self, s):
        """Return scaling matrix.

        s must be a 3-sequence (e.g. a vec3).
        """
        cdef mat4 res
        res = mat4()
        res.m11 = s[0]
        res.m12 = 0.0
        res.m13 = 0.0
        res.m14 = 0.0
        res.m21 = 0.0
        res.m22 = s[1]
        res.m23 = 0.0
        res.m24 = 0.0
        res.m31 = 0.0
        res.m32 = 0.0
        res.m33 = s[2]
        res.m34 = 0.0
        res.m41 = 0.0
        res.m42 = 0.0
        res.m43 = 0.0
        res.m44 = 1.0
        return res

    def rotation(self, angle, axis):
        """Return rotation matrix.

        angle must be given in radians. axis must be a 3-sequence.
        """
        cdef double cangle, ax, ay, az
        cdef double sqr_a, sqr_b, sqr_c, len2
        cdef double k1, k2, k3, k1ab, k1ac, k1bc, k3a, k3b, k3c
        cdef mat4 res

        cangle   = angle
        ax    = axis[0]
        ay    = axis[1]
        az    = axis[2]
        sqr_a = ax*ax
        sqr_b = ay*ay
        sqr_c = az*az
        len2  = sqr_a+sqr_b+sqr_c

        k2    = cos(cangle)
        k1    = (1.0-k2)/len2
        k3    = sin(cangle)/sqrt(len2)
        k1ab  = k1*ax*ay
        k1ac  = k1*ax*az
        k1bc  = k1*ay*az
        k3a   = k3*ax
        k3b   = k3*ay
        k3c   = k3*az

        res = mat4()
        res.m11 = k1*sqr_a+k2
        res.m12 = k1ab-k3c
        res.m13 = k1ac+k3b
        res.m14 = 0.0
        res.m21 = k1ab+k3c
        res.m22 = k1*sqr_b+k2
        res.m23 = k1bc-k3a
        res.m24 = 0.0
        res.m31 = k1ac-k3b
        res.m32 = k1bc+k3a
        res.m33 = k1*sqr_c+k2
        res.m34 = 0.0
        res.m41 = 0.0
        res.m42 = 0.0
        res.m43 = 0.0
        res.m44 = 1.0
        return res


    def translate(self, t):
        """Concatenate a translation.

        t must be a 3-sequence (e.g. a vec3).
        """
        cdef double tx, ty, tz
        tx = t[0]
        ty = t[1]
        tz = t[2]

        self.m14 = self.m11*tx + self.m12*ty + self.m13*tz + self.m14
        self.m24 = self.m21*tx + self.m22*ty + self.m23*tz + self.m24
        self.m34 = self.m31*tx + self.m32*ty + self.m33*tz + self.m34
        self.m44 = self.m41*tx + self.m42*ty + self.m43*tz + self.m44
        return self

    def scale(self, s):
        """Concatenate a scaling.

        s must be a 3-sequence (e.g. a vec3).
        """
        cdef double sx, sy, sz
        sx = s[0]
        sy = s[1]
        sz = s[2]
        self.m11 = self.m11*sx
        self.m12 = self.m12*sy
        self.m13 = self.m13*sz
        self.m21 = self.m21*sx
        self.m22 = self.m22*sy
        self.m23 = self.m23*sz
        self.m31 = self.m31*sx
        self.m32 = self.m32*sy
        self.m33 = self.m33*sz
        self.m41 = self.m41*sx
        self.m42 = self.m42*sy
        self.m43 = self.m43*sz
        return self

    def rotate(self, angle, axis):
        """Concatenate a rotation.

        angle must be given in radians. axis must be a 3-sequence.
        """
        cdef mat4 R
        cdef double c1, c2, c3, c4
        R=self.rotation(angle, axis)

        # self*R  (in-place)
        c1 = self.m11
        c2 = self.m12
        c3 = self.m13
        c4 = self.m14
        self.m11 = c1*R.m11+c2*R.m21+c3*R.m31+c4*R.m41
        self.m12 = c1*R.m12+c2*R.m22+c3*R.m32+c4*R.m42
        self.m13 = c1*R.m13+c2*R.m23+c3*R.m33+c4*R.m43
        self.m14 = c1*R.m14+c2*R.m24+c3*R.m34+c4*R.m44

        c1 = self.m21
        c2 = self.m22
        c3 = self.m23
        c4 = self.m24
        self.m21 = c1*R.m11+c2*R.m21+c3*R.m31+c4*R.m41
        self.m22 = c1*R.m12+c2*R.m22+c3*R.m32+c4*R.m42
        self.m23 = c1*R.m13+c2*R.m23+c3*R.m33+c4*R.m43
        self.m24 = c1*R.m14+c2*R.m24+c3*R.m34+c4*R.m44

        c1 = self.m31
        c2 = self.m32
        c3 = self.m33
        c4 = self.m34
        self.m31 = c1*R.m11+c2*R.m21+c3*R.m31+c4*R.m41
        self.m32 = c1*R.m12+c2*R.m22+c3*R.m32+c4*R.m42
        self.m33 = c1*R.m13+c2*R.m23+c3*R.m33+c4*R.m43
        self.m34 = c1*R.m14+c2*R.m24+c3*R.m34+c4*R.m44

        c1 = self.m41
        c2 = self.m42
        c3 = self.m43
        c4 = self.m44
        self.m41 = c1*R.m11+c2*R.m21+c3*R.m31+c4*R.m41
        self.m42 = c1*R.m12+c2*R.m22+c3*R.m32+c4*R.m42
        self.m43 = c1*R.m13+c2*R.m23+c3*R.m33+c4*R.m43
        self.m44 = c1*R.m14+c2*R.m24+c3*R.m34+c4*R.m44

        return self

    def orthographic(self, left, right, bottom, top, near, far):
        """equivalent to the OpenGL command glOrtho()"""
        cdef mat4 res
        cdef double l,r,b,t,n,f, r_l, t_b, f_n

        l = left
        r = right
        b = bottom
        t = top
        n = near
        f = far

        r_l = r-l
        t_b = t-b
        f_n = f-n

        if r_l<=eps:
            raise ValueError, "right-value must be greater than left-value"
        if t_b<=eps:
            raise ValueError, "top-value must be greater than bottom-value"
        if f_n<=eps:
            raise ValueError, "far-value must be greater than near-value"

        res = mat4()

        res.m11 = 2.0/r_l
        res.m12 = 0.0
        res.m13 = 0.0
        res.m14 = -(r+l)/r_l
        res.m21 = 0.0
        res.m22 = 2.0/t_b
        res.m23 = 0.0
        res.m24 = -(t+b)/t_b
        res.m31 = 0.0
        res.m32 = 0.0
        res.m33 = -2.0/f_n
        res.m34 = -(f+n)/f_n
        res.m41 = 0.0
        res.m42 = 0.0
        res.m43 = 0.0
        res.m44 = 1.0
        return res

    def frustum(self, left, right, bottom, top, near, far):
        """equivalent to the OpenGL command glFrustum()"""
        cdef mat4 res
        cdef double l,r,b,t,n,f, r_l, t_b, f_n

        l = left
        r = right
        b = bottom
        t = top
        n = near
        f = far

        r_l = r-l
        t_b = t-b
        f_n = f-n

        if r_l<=eps:
            raise ValueError, "right-value must be greater than left-value"
        if t_b<=eps:
            raise ValueError, "top-value must be greater than bottom-value"
        if f_n<=eps:
            raise ValueError, "far-value must be greater than near-value"

        res = mat4()

        res.m11 = (2.0*n)/r_l
        res.m12 = 0.0
        res.m13 = (r+l)/r_l
        res.m14 = 0.0
        res.m21 = 0.0
        res.m22 = (2.0*n)/t_b
        res.m23 = (t+b)/t_b
        res.m24 = 0.0
        res.m31 = 0.0
        res.m32 = 0.0
        res.m33 = -(f+n)/f_n
        res.m34 = (-2.0*f*n)/f_n
        res.m41 = 0.0
        res.m42 = 0.0
        res.m43 = -1.0
        res.m44 = 0.0
        return res

    def perspective(self, fovy, aspect, near, far):
        """von Mesa übernommen (glu.c)"""
        cdef double f, a
        cdef double top, bottom, left, right

        f = fovy
        a = aspect

        top    = near * tan(f * 3.1415926535897931 / 360.0)
        bottom = -top
        left   = bottom * a
        right  = top * a

        return self.frustum(left, right, bottom, top, near, far)

    def lookAt(self, pos, target, up=(0,0,1)):
        """Look from pos to target.

        The resulting transformation moves the origin to pos and
        rotates so that The z-axis points to target. The y-axis is
        as close as possible to the up vector.
        """
        cdef mat4 res
        cdef vec3 right,vup,dir

        pos    = vec3(pos)
        target = vec3(target)
        vup    = vec3(up)

        dir = (target - pos).normalize()
        vup = vup.normalize()
        vup = vup - (vup * dir) * dir
        try:
            vup  = vup.normalize()
        except:
            # We're looking along the up direction, so choose
            # an arbitrary direction that is perpendicular to dir
            # as new up.
            vup = dir.ortho()

        right = vup.cross(dir).normalize()

        res = mat4()
        res.m11 = right.x
        res.m21 = right.y
        res.m31 = right.z
        res.m41 = 0.0
        res.m12 = vup.x
        res.m22 = vup.y
        res.m32 = vup.z
        res.m42 = 0.0
        res.m13 = dir.x
        res.m23 = dir.y
        res.m33 = dir.z
        res.m43 = 0.0
        res.m14 = pos.x
        res.m24 = pos.y
        res.m34 = pos.z
        res.m44 = 1.0
        return res

    def ortho(self):
        """Return a matrix with orthogonal base vectors.

        Makes the x-, y- and z-axis orthogonal.
        The fourth column and row remain untouched.
        """
        cdef mat4 res
        cdef vec3 x,y,z

        x = vec3(self.m11, self.m21, self.m31)
        y = vec3(self.m12, self.m22, self.m32)
        z = vec3(self.m13, self.m23, self.m33)

        xl = x.length()
        xl = xl*xl
        y = y - ((x*y)/xl)*x
        z = z - ((x*z)/xl)*x

        yl = y.length()
        yl = yl*yl
        z = z - ((y*z)/yl)*y

        res = mat4()
        res.m11 = x.x
        res.m12 = y.x
        res.m13 = z.x
        res.m14 = self.m14
        res.m21 = x.y
        res.m22 = y.y
        res.m23 = z.y
        res.m24 = self.m24
        res.m31 = x.z
        res.m32 = y.z
        res.m33 = z.z
        res.m34 = self.m34
        res.m41 = self.m41
        res.m42 = self.m42
        res.m43 = self.m43
        res.m44 = self.m44
        return res

    def decompose(self):
        """Decomposes the matrix into a translation, rotation and scaling part.

        Returns a tuple (translation, rotation, scaling). The
        translation and scaling parts are given as vec3's, the rotation
        is still given as a mat4.
        """
        cdef mat4 m
        cdef vec3 a,b,c
        cdef double al, bl, cl

        m = self.ortho()
        m.m14 = 0.0
        m.m24 = 0.0
        m.m34 = 0.0
        m.m44 = 1.0
        m.m41 = 0.0
        m.m42 = 0.0
        m.m43 = 0.0

        # a,b,c = Column 0,1,2 of m
        a = vec3()
        b = vec3()
        c = vec3()
        a.x = m.m11
        a.y = m.m21
        a.z = m.m31
        b.x = m.m12
        b.y = m.m22
        b.z = m.m32
        c.x = m.m13
        c.y = m.m23
        c.z = m.m33

        # al,bl,cl = length of a,b,c (= scaling factors)
        al = sqrt(a.x*a.x + a.y*a.y + a.z*a.z)
        bl = sqrt(b.x*b.x + b.y*b.y + b.z*b.z)
        cl = sqrt(c.x*c.x + c.y*c.y + c.z*c.z)
        if al<=eps or bl<=eps or cl<=eps:
            raise ZeroDivisionError,"transformation contains a 0 scaling"

        scale = vec3(al,bl,cl)

        # normalizing a,b,c
#        a/=al
#        b/=bl
#        c/=cl
        a.x = a.x/al
        a.y = a.y/al
        a.z = a.z/al
        b.x = b.x/bl
        b.y = b.y/bl
        b.z = b.z/bl
        c.x = c.x/cl
        c.y = c.y/cl
        c.z = c.z/cl

        m.m11 = a.x
        m.m21 = a.y
        m.m31 = a.z
        m.m12 = b.x
        m.m22 = b.y
        m.m32 = b.z
        m.m13 = c.x
        m.m23 = c.y
        m.m33 = c.z
        if m.determinant()<0.0:
            m.m11 = -m.m11
            m.m21 = -m.m21
            m.m31 = -m.m31
            scale.x = -scale.x

        return (vec3(self.m14, self.m24, self.m34), m, scale)

    def getMat3(self):
        """Convert to mat3 by discarding 4th row and column.
        """
        cdef mat3 res
        res = mat3()
        res.m11 = self.m11
        res.m12 = self.m12
        res.m13 = self.m13
        res.m21 = self.m21
        res.m22 = self.m22
        res.m23 = self.m23
        res.m31 = self.m31
        res.m32 = self.m32
        res.m33 = self.m33
        return res

    def setMat3(self, mat3 m):
        """Set the mat3 part of self.
        """
        self.m11 = m.m11
        self.m12 = m.m12
        self.m13 = m.m13
        self.m21 = m.m21
        self.m22 = m.m22
        self.m23 = m.m23
        self.m31 = m.m31
        self.m32 = m.m32
        self.m33 = m.m33

######################################################################

# quat
cdef class quat:
    """Quaternion class.

    Quaternions are an extension to complex numbers and can be used
    to store rotations. They are composed of four floats which can be
    seen as an angle and an axis of rotation.
    """

    cdef double w,x,y,z

    def __cinit__(self, *args):
        """Constructor.

        0 arguments: zeroes
        1 float argument:  w component, x,y,z = (0,0,0)
        1 quat argument: Make a copy
        1 mat3 argument: Initialize by rotation matrix
        1 mat4 argument: Initialize by rotation matrix
        2 arguments: angle & axis (doesn't have to be of unit length)
        4 arguments: components w,x,y,z
        """

        cdef int arglen, seqlen
        cdef quat q
        arglen = len(args)

        # 0 arguments
        if arglen==0:
            self.w = 0.0
            self.x = 0.0
            self.y = 0.0
            self.z = 0.0

        # 1 argument
        elif arglen==1:
            T = type(args[0])
            # scalar
            if T==float or T==int or T==long:
                self.w = args[0]
                self.x = 0.0
                self.y = 0.0
                self.z = 0.0
            # quat
            elif T==quat:
                q = args[0]
                self.w = q.w
                self.x = q.x
                self.y = q.y
                self.z = q.z
            # String
            elif T==str or T==unicode:
                s=args[0].replace(","," ").replace("  "," ").strip().split(" ")
                if s==[""]:
                    s=[]
                for i in range(len(s)):
                    try:
                        s[i] = float(s[i])
                    except:
                        raise ValueError,"quat() arg is a malformed string"
                q = quat(s)
                self.w = q.w
                self.x = q.x
                self.y = q.y
                self.z = q.z
            # Sequence of floats
            else:
                seq = args[0]
                try:
                    seqlen = len(seq)
                except:
                    raise TypeError,"quat() arg can't be converted to quat"

                if seqlen==0:
                    self.w = 0.0
                    self.x = 0.0
                    self.y = 0.0
                    self.z = 0.0
                elif seqlen==1:
                    self.w = seq[0]
                    self.x = 0.0
                    self.y = 0.0
                    self.z = 0.0
                elif seqlen==4:
                    self.w = seq[0]
                    self.x = seq[1]
                    self.y = seq[2]
                    self.z = seq[3]
                else:
                    raise TypeError, "quat() arg can't be converted to quat"

        # 2 arguments (angle & axis)
        elif arglen==2:
            angle, axis = args
            self.fromAngleAxis(angle,axis)

        # 4 arguments
        elif arglen==4:
            self.w = args[0]
            self.x = args[1]
            self.y = args[2]
            self.z = args[3]

        else:
            raise TypeError, "quat() arg can't be converted to quat"

    def __reduce__(self):
        return (_quatconstruct, (self.w, self.x, self.y, self.z))

    def __repr__(self):
        return 'quat(%s, %s, %s, %s)'%(repr(self.w),repr(self.x),repr(self.y),repr(self.z))

    def __str__(self):
        return "(%1.4f, %1.4f, %1.4f, %1.4f)"%(self.w, self.x, self.y, self.z)


    def __richcmp__(a,b,op):
        cdef quat va, vb
        cdef int cop

        cop = op

        ta = type(a)
        tb = type(b)
        if (ta!=quat or tb!=quat):
            if cop==3:
                return 1
            else:
                return 0

        va = a
        vb = b

        # <
        if cop==0:
            return 0
        # <=
        elif cop==1:
            return 0
        # ==
        elif cop==2:
            return (fabs(va.w-vb.w)<=eps and fabs(va.x-vb.x)<=eps and fabs(va.y-vb.y)<=eps and fabs(va.z-vb.z)<=eps)
        # !=
        elif cop==3:
            return (fabs(va.w-vb.w)>eps or fabs(va.x-vb.x)>eps or fabs(va.y-vb.y)>eps and fabs(va.z-vb.z)>eps)
        # >
        elif cop==4:
            return 0
        # >=
        elif cop==5:
            return 0
        # sonst (Fehler)
        else:
            raise ValueError,"internal error: illegal rich comparison number"

    def __add__(quat a, quat b):
        """Addition.

        >>> q=quat(0.9689, 0.2160, 0.1080, 0.0540)
        >>> print q+q
        (1.9378, 0.4320, 0.2160, 0.1080)
        """
        cdef quat res
        res = quat()
        res.w = a.w+b.w
        res.x = a.x+b.x
        res.y = a.y+b.y
        res.z = a.z+b.z
        return res

    def __sub__(quat a, quat b):
        """Subtraction.

        >>> q=quat(0.9689, 0.2160, 0.1080, 0.0540)
        >>> print q-q
        (0.0000, 0.0000, 0.0000, 0.0000)
        """
        cdef quat res
        res = quat()
        res.w = a.w-b.w
        res.x = a.x-b.x
        res.y = a.y-b.y
        res.z = a.z-b.z
        return res

    def __mul__(a, b):
        """Multiplication.

        >>> q=quat(0.9689, 0.2160, 0.1080, 0.0540)
        >>> print q*2.0
        (1.9378, 0.4320, 0.2160, 0.1080)
        >>> print 2.0*q
        (1.9378, 0.4320, 0.2160, 0.1080)
        >>> print q*q
        (0.8775, 0.4186, 0.2093, 0.1046)
        """
        cdef quat res, va, vb
        cdef double r

        ta = type(a)
        tb = type(b)

        if ta==quat:
            va = a
            if tb==quat:
                # quat*quat
                vb = b
                res = quat()
                res.w = va.w*vb.w - va.x*vb.x - va.y*vb.y - va.z*vb.z
                res.x = va.w*vb.x + va.x*vb.w + va.y*vb.z - va.z*vb.y
                res.y = va.w*vb.y + va.y*vb.w - va.x*vb.z + va.z*vb.x
                res.z = va.w*vb.z + va.z*vb.w + va.x*vb.y - va.y*vb.x
                return res
            elif tb==float or tb==int or tb==long:
                # quat*scalar
                res = quat()
                r   = b
                res.w = va.w*r
                res.x = va.x*r
                res.y = va.y*r
                res.z = va.z*r
                return res
        elif ta==float or ta==int or ta==long:
            if tb==quat:
                # scalar*quat
                res = quat()
                vb  = b
                r   = a
                res.w = r*vb.w
                res.x = r*vb.x
                res.y = r*vb.y
                res.z = r*vb.z
                return res

        raise TypeError, "unsupported operand type for *"

    def __div__(quat a, b):
        """Division.

        >>> q=quat(0.9689, 0.2160, 0.1080, 0.0540)
        >>> print q/2.0
        (0.4844, 0.1080, 0.0540, 0.0270)
        """
        cdef quat res
        cdef double r

        tb = type(b)

        if tb==float or tb==int or tb==long:
            # quat/scalar
            res = quat()
            r = b
            if fabs(r)<=eps:
                raise ZeroDivisionError,"quat division"
            res.w = a.w/r
            res.x = a.x/r
            res.y = a.y/r
            res.z = a.z/r
            return res
        else:
            raise TypeError, "unsupported operand type for /"

    def __pow__(quat self, other, modulo):
        """Return self**q."""
        cdef quat q
        if modulo!=None:
            raise TypeError, "unsupported operation"
        q = quat(other)
        return (q*self.log()).exp()

    def __neg__(self):
        """Negation.

        >>> q=quat(0.9689, 0.2160, 0.1080, 0.0540)
        >>> print -q
        (-0.9689, -0.2160, -0.1080, -0.0540)
        """
        cdef quat res
        res = quat()
        res.w = -self.w
        res.x = -self.x
        res.y = -self.y
        res.z = -self.z
        return res

    def __pos__(self):
        """
        >>> q=quat(0.9689, 0.2160, 0.1080, 0.0540)
        >>> print +q
        (0.9689, 0.2160, 0.1080, 0.0540)
        """
        cdef quat res
        res = quat()
        res.w = self.w
        res.x = self.x
        res.y = self.y
        res.z = self.z
        return res

    def __iadd__(self, quat other):
        """Inline quaternion addition.
        """
        self.w=self.w+other.w
        self.x=self.x+other.x
        self.y=self.y+other.y
        self.z=self.z+other.z
        return self

    def __isub__(self, quat other):
        """Inline quaternion subtraction.
        """
        self.w=self.w-other.w
        self.x=self.x-other.x
        self.y=self.y-other.y
        self.z=self.z-other.z
        return self

    def __imul__(self, other):
        """Inline multiplication.
        """
        cdef quat vb
        cdef double r
        cdef double w,x,y,z

        T = type(other)

        if T==quat:
            # quat*=quat
            vb = other
            w = self.w*vb.w - self.x*vb.x - self.y*vb.y - self.z*vb.z
            x = self.w*vb.x + self.x*vb.w + self.y*vb.z - self.z*vb.y
            y = self.w*vb.y + self.y*vb.w - self.x*vb.z + self.z*vb.x
            z = self.w*vb.z + self.z*vb.w + self.x*vb.y - self.y*vb.x
            self.w = w
            self.x = x
            self.y = y
            self.z = z
            return self

        elif T==float or T==int or T==long:
            # quat*=scalar
            r   = other
            self.w = self.w*r
            self.x = self.x*r
            self.y = self.y*r
            self.z = self.z*r
            return self

        else:
            raise TypeError, "unsupported operand type for *="

    def __idiv__(self, other):
        """Inline division with scalar.
        """
        cdef double r
        tb = type(other)
        if tb==float or tb==int or tb==long:
            r=other
            if fabs(r)<=eps:
                raise ZeroDivisionError,"quat division"
            self.w=self.w/r
            self.x=self.x/r
            self.y=self.y/r
            self.z=self.z/r
            return self
        else:
            raise TypeError, "unsupported operand type for /="

    def __getattr__(self, name):
        if name=="w":
            return self.w
        elif name=="x":
            return self.x
        elif name=="y":
            return self.y
        elif name=="z":
            return self.z
        else:
            raise AttributeError,"quat has no attribute '"+name+"'"

    def __setattr__(self, name, val):
        if name=="w":
            self.w = val
        elif name=="x":
            self.x = val
        elif name=="y":
            self.y = val
        elif name=="z":
            self.z = val
        else:
            raise AttributeError,"quat has no writable attribute '"+name+"'"


    def __abs__(self):
        """Return magnitude.

        >>> q=quat(0.9689, 0.2160, 0.1080, 0.0540)
        >>> print round(abs(q),5)
        1.0
        """
        cdef double a
        a = self.w*self.w + self.x*self.x + self.y*self.y + self.z*self.z
        return sqrt(a)

    def conjugate(self):
        """Return conjugate.

        >>> q=quat(0.9689, 0.2160, 0.1080, 0.0540)
        >>> print q.conjugate()
        (0.9689, -0.2160, -0.1080, -0.0540)
        """
        cdef quat res
        res=quat()
        res.w = self.w
        res.x = -self.x
        res.y = -self.y
        res.z = -self.z
        return res

    def normalize(self):
        """Return normalized quaternion.

        >>> q=quat(0.9, 0.5, 0.2, 0.3)
        >>> q=q.normalize()
        >>> print q
        (0.8250, 0.4583, 0.1833, 0.2750)
        >>> print abs(q)
        1.0
        """
        cdef quat res
        cdef double nlen

        nlen = sqrt(self.w*self.w+self.x*self.x+self.y*self.y+self.z*self.z)
        if nlen<=eps:
            raise ZeroDivisionError,"quat division"

        nlen = 1.0/nlen
        res = quat()
        res.w = self.w*nlen
        res.x = self.x*nlen
        res.y = self.y*nlen
        res.z = self.z*nlen
        return res

    def inverse(self):
        """Return inverse.

        >>> q=quat(0.9, 0.5, 0.2, 0.3)
        >>> print q.inverse()
        (0.7563, -0.4202, -0.1681, -0.2521)
        """
        cdef quat res

        len_2 = self.w*self.w+self.x*self.x+self.y*self.y+self.z*self.z
        if fabs(len_2)<=eps:
            raise ZeroDivisionError,"quat is not invertible"

        len_2 = 1.0/len_2
        res = quat()
        res.w =  self.w*len_2
        res.x = -self.x*len_2
        res.y = -self.y*len_2
        res.z = -self.z*len_2
        return res

    def toAngleAxis(self):
        """Return angle (in radians) and rotation axis.

        >>> q=quat(0.9, 0.5, 0.2, 0.3)
        >>> angle, axis = q.toAngleAxis()
        >>> print round(angle,4)
        1.2011
        >>> print axis
        (0.8111, 0.3244, 0.4867)
        """
        cdef quat q
        cdef double s,w

        q = self.normalize()

        # Clamp nself.w (since the quat has to be normalized it should
        # be between -1 and 1 anyway, but it might be slightly off due
        # to numerical inaccuracies)
        if q.w<-1.0:
            w = -1.0
        elif q.w>1.0:
            w = 1.0
        else:
            w = q.w

        w = acos(w)
        s = sin(w)
        if s<=eps:
            return (0.0, vec3(0.0,0.0,0.0))

        return (2.0*w, vec3(q.x/s, q.y/s, q.z/s))

    def fromAngleAxis(self, angle, axis):
        """Initialize self from an angle (in radians) and an axis and returns self."""
        cdef double a,n, x,y,z
        cdef quat q

        a = angle/2.0
        self.w = cos(a)

        x = axis[0]
        y = axis[1]
        z = axis[2]
        n = x*x+y*y+z*z
        if n<=eps:
            raise ValueError,"axis mustn't be the null vector"

        s = sin(a)/n
        self.x = x*s
        self.y = y*s
        self.z = z*s

        # Normalize the quat
        q = self.normalize()
        self.w = q.w
        self.x = q.x
        self.y = q.y
        self.z = q.z
        return self

    def toMat3(self):
        """Return rotation matrix as mat3."""
        cdef mat3 res
        cdef xx,yy,zz,xy,zw,xz,yw,yz,xw

        xx = 2.0*self.x*self.x
        yy = 2.0*self.y*self.y
        zz = 2.0*self.z*self.z
        xy = 2.0*self.x*self.y
        zw = 2.0*self.z*self.w
        xz = 2.0*self.x*self.z
        yw = 2.0*self.y*self.w
        yz = 2.0*self.y*self.z
        xw = 2.0*self.x*self.w

        res = mat3()
        res.m11 = 1.0-yy-zz
        res.m12 = xy-zw
        res.m13 = xz+yw
        res.m21 = xy+zw
        res.m22 = 1.0-xx-zz
        res.m23 = yz-xw
        res.m31 = xz-yw
        res.m32 = yz+xw
        res.m33 = 1.0-xx-yy
        return res

    def toMat4(self):
        """Return rotation matrix as mat4."""
        cdef mat4 res
        cdef xx,yy,zz,xy,zw,xz,yw,yz,xw

        xx = 2.0*self.x*self.x
        yy = 2.0*self.y*self.y
        zz = 2.0*self.z*self.z
        xy = 2.0*self.x*self.y
        zw = 2.0*self.z*self.w
        xz = 2.0*self.x*self.z
        yw = 2.0*self.y*self.w
        yz = 2.0*self.y*self.z
        xw = 2.0*self.x*self.w

        res = mat4()

        res.m11 = 1.0-yy-zz
        res.m12 = xy-zw
        res.m13 = xz+yw
        res.m14 = 0.0
        res.m21 = xy+zw
        res.m22 = 1.0-xx-zz
        res.m23 = yz-xw
        res.m24 = 0.0
        res.m31 = xz-yw
        res.m32 = yz+xw
        res.m33 = 1.0-xx-yy
        res.m34 = 0.0
        res.m41 = 0.0
        res.m42 = 0.0
        res.m43 = 0.0
        res.m44 = 1.0
        return res

    def fromMat(self, m):
        """Initialize self from either a mat3 or mat4 and returns self."""
        cdef double d1,d2,d3,ad1,ad2,ad3,t,s

        d1 = m[0,0]
        d2 = m[1,1]
        d3 = m[2,2]
        t = d1+d2+d3+1.0
        if t>0.0:
            s = 0.5/sqrt(t)
            self.w = 0.25/s
            self.x = (m[2,1]-m[1,2])*s
            self.y = (m[0,2]-m[2,0])*s
            self.z = (m[1,0]-m[0,1])*s
        else:
            ad1 = fabs(d1)
            ad2 = fabs(d2)
            ad3 = fabs(d3)
            if ad1>=ad2 and ad1>=ad3:
                s = sqrt(1.0+d1-d2-d3)*2.0
                self.x = 0.5/s
                self.y = (m[0,1]+m[1,0])/s
                self.z = (m[0,2]+m[2,0])/s
                self.w = (m[1,2]+m[2,1])/s
            elif ad2>=ad1 and ad2>=ad3:
                s = sqrt(1.0+d2-d1-d3)*2.0
                self.x = (m[0,1]+m[1,0])/s
                self.y = 0.5/s
                self.z = (m[1,2]+m[2,1])/s
                self.w = (m[0,2]+m[2,0])/s
            else:
                s = sqrt(1.0+d3-d1-d2)*2.0
                self.x = (m[0,2]+m[2,0])/s
                self.y = (m[1,2]+m[2,1])/s
                self.z = 0.5/s
                self.w = (m[0,1]+m[1,0])/s

        return self

    def dot(self, quat b):
        """Return the dot product of self and b."""
        return self.w*b.w + self.x*b.x + self.y*b.y + self.z*b.z

    def log(self):
        """Return the natural logarithm of self."""

        cdef double b, t, f, ct, r
        cdef quat res

        b = sqrt(self.x*self.x + self.y*self.y + self.z*self.z)
        res = quat()
        if fabs(b)<=eps:
            res.x = 0.0
            res.y = 0.0
            res.z = 0.0
            if self.w<=eps:
                raise ValueError, "math domain error"
            res.w = log(self.w)
        else:
            t = atan2(b, self.w)
            f = t/b
            res.x = f*self.x
            res.y = f*self.y
            res.z = f*self.z
            ct = cos(t)
            if fabs(ct)<=eps:
                raise ValueError, "math domain error"
            r = self.w/ct
            if r<=eps:
                raise ValueError, "math domain error"
            res.w = log(r)

        return res

    def exp(self):
        """Return the exponential of self."""

        cdef double b, f
        cdef quat res

        b = sqrt(self.x*self.x + self.y*self.y + self.z*self.z)
        res = quat()
        if fabs(b)<=eps:
            res.x = 0.0
            res.y = 0.0
            res.z = 0.0
            res.w = exp(self.w)
        else:
            f = sin(b)/b
            res.x = f*self.x
            res.y = f*self.y
            res.z = f*self.z
            res.w = exp(self.w)*cos(b)

        return res


def slerp(double t, quat q0, quat q1):
    """Spherical linear interpolation between two quaternions.

    The return value is an interpolation between q0 and q1. For t=0.0
    the return value equals q0, for t=1.0 it equals q1.
    q0 and q1 must be unit quaternions.
    """
    cdef double o,so,a,b

    o = acos(q0.dot(q1))
    so = sin(o)

    if (fabs(so)<eps):
        return quat(q0)

    a = sin(o*(1.0-t)) / so
    b = sin(o*t) / so
    return q0*a + q1*b

def squad(double t, quat a, quat b, quat c, quat d):
    """Spherical cubic interpolation."""
    return slerp(2*t*(1.0-t), slerp(t,a,d), slerp(t,b,c))

######################################################################

import sys
if sys.version_info[0] >= 3:
    import copyreg as copy_reg
else:
    import copy_reg

def _vec3construct(x,y,z):
    return vec3(x,y,z)

def _vec4construct(x,y,z,w):
    return vec4(x,y,z,w)

def _mat3construct(b1,b2,b3):
    return mat3(b1,b2,b3)

def _mat4construct(b1,b2,b3,b4):
    return mat4(b1,b2,b3,b4)

def _quatconstruct(w,x,y,z):
    return quat(w,x,y,z)

copy_reg.constructor(_vec3construct)
copy_reg.constructor(_vec4construct)
copy_reg.constructor(_mat3construct)
copy_reg.constructor(_mat4construct)
copy_reg.constructor(_quatconstruct)
