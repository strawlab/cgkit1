######################################################################
# cgkit - setup script
#
# Copyright (C) 2003, Matthias Baas (baas@ira.uka.de)
#
# http://cgkit.sourceforge.net
#
# You may distribute under the terms of the BSD license, as
# specified in the file license.txt.
####################################################################

from distutils.core import setup, Extension
import shutil, os, sys, os.path, time
#from Pyrex.Distutils import build_ext

# Check for an old version of cgtypes
exitflag = 0
for mp in sys.path:
    if mp==os.getcwd():
        continue
    p = os.path.join(mp,"cgtypes")
    if os.path.isdir(p):
        print "There is an old version of cgkit installed."
        print "Please uninstall the old version or simply remove the following directory:"
        print p
        print "(otherwise the new cgtypes module won't be in effect)"
        exitflag=1
        break

if exitflag:
    sys.exit(1)

######################################################################

# updateInfoModule
def updateInfoModule():
    """Updates the module cgkitinfo.

    Updates the info module so that the version string contains the
    current date and time.
    """

    infomod = "cgkitinfo.py"
    tmpname = "__cgkitinfo.py"

    # Read the cgkitinfo module
    try:
        f = file(infomod)
        lines = f.readlines()
        f.close()
    except:
        print 'Could not read file "%s"'%infomod
        print 'Make sure that "%s" is renamed to "%s" if it exists'%(tmpname,infomod)
        sys.exit(1)

    # Replace the version string...
    version_info = None
    for i in range(len(lines)):
        s = lines[i]
        f = s.split(" ")
        if len(f)==0:
            continue
        if f[0]=="version_info":
            try:
                exec s.strip()
                if type(version_info)!=tuple:
                    raise Exception
                if len(version_info)!=4:
                    raise Exception
                major,minor,micro,rlevel = version_info
                if type(major)!=int or type(minor)!=int or type(micro)!=int:
                    raise Exception
                if type(rlevel)!=str:
                    raise Exception
            except:
                print 'Invalid version tuple in file "%s": %s'%(infomod,s.strip())
                sys.exit(1)
        if f[0]=="version":
            if version_info==None:
                print '%s: version_info must occur before version'%infomod
                sys.exit(1)
            major,minor,micro,rlevel = version_info
            v = "%d.%d.%d"%(major,minor,micro)
            if rlevel!="final":
                v+=rlevel
            v+=" (%s)"%time.strftime("%b %d %Y, %H:%M")
            lines[i] = 'version = "%s"\n'%v

    # Save the new content to a temporary file...
    try:
        f = file(tmpname, "wt")
        f.writelines(lines)
        f.close()
    except:
        print 'Could not write temporary file "%s"'%tmpname
        sys.exit(1)

    # Remove the old module...
    try:
        os.remove(infomod)
    except:
        print 'Replacing module failed: Could not remove file "%s"'%infomod
        sys.exit(1)

    # ...and rename the new file
    try:
        os.rename(tmpname, infomod)
    except:
        print 'Replacing module failed: Could not rename file "%s"'%tmpname
        print 'Please rename the file "%s" into "%s" manually.'%(tmpname,infomod)
        sys.exit(1)


######################################################################

updateInfoModule()

pyrex_in  = "cgtypes.pyx"
pyrex_out = "cgtypes.c"
pyrex_cmd = "cython "+pyrex_in

print "Updating",pyrex_out
os.system(pyrex_cmd)

setup(name="cgkit",
      version="1.2.0",
      description="Python Computer Graphics Kit",
      author="Matthias Baas",
      author_email="baas@ira.uka.de",
      url="http://cgkit.sourceforge.net",
      license="BSD license, see license.txt",
      py_modules=["cgkitinfo", "ri","riutil",
                  "pycgtypes.vec3","pycgtypes.mat3",
                  "pycgtypes.vec4","pycgtypes.mat4","pycgtypes.quat",
                  "sltokenize","sl","slparams", "_slparser"],
      ext_modules=[Extension("cgtypes", ["cgtypes.c"]),
                   Extension("noise", ["src/noisemodule.cpp"])
                   ]
      )
