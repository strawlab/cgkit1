####################################################################
# slparams - Extract shader parameters from a shader source file
#
# Copyright (C) 2003, Matthias Baas (baas@ira.uka.de)
#
# http://cgkit.sourceforge.net
#
# You may distribute under the terms of the BSD license, as
# specified in the file license.txt.
####################################################################

import os, sys, string, StringIO, sltokenize
import cgtypes, math, sl
import _slparser

class SLParamsError(Exception):
    pass

class PreprocessorNotFound(SLParamsError):
    pass

class SyntaxError(SLParamsError):
    pass

class NoMoreTokens(SLParamsError):
    pass


# Parser class (subclassed from the generated Yapps parser)
class _SLParser(_slparser._SLParserBase):
    """SL parser class.

    This class is derived from the generated parser and implements
    some missing methods.
    """
    def __init__(self, scanner):
        _slparser._SLParserBase.__init__(self, scanner)

        # Current filename
        self.filename = "?"
        # Offset which has to be subtracted from the line number to
        # get the true line number within the file.
        self.linenroffset = 0

        # Parameter attributes...
        self.output   = ""
        self.detail   = ""
        self.type     = ""
        self.arraylen = None
        self.name     = ""
        self.space    = None
        self.default  = ""
        self.spaces   = []

        # Shader parameters
        self.params = []

    def newParams(self):
        """Start a new shader.

        The parameter list is cleared.
        """
        self.params = []

    def newType(self):
        """Clear all type parameters.

        This is called when a new type is declared (which is not equivalent
        to a new parameter, e.g. "float a,b,c")
        """
        self.output   = ""
        self.detail   = "uniform"
        self.type     = ""
        self.arraylen = None
        self.name     = ""
        self.space    = None
        self.default  = ""
        self.spaces   = []

    def defaultSpace(self, typ):
        """Return the default space for typ."""
        if typ in ["point", "vector", "normal", "matrix"]:
            return "current"
        elif typ=="color":
            return "rgb"
        else:
            return None
        
    def appendSpace(self):
        """Append self.space to self.spaces."""
        if self.space!=None:
            self.spaces.append(self.space)
        else:
            self.spaces.append(self.defaultSpace(self.type))
                
        self.space = None

    def storeParam(self):
        """Store the current set of attributes as a new parameter.

        The attributes are reset so that a new parameter can begin.
        """
        if self.arraylen==None:
            if self.space==None:
                self.space = self.defaultSpace(self.type)
            self.params.append((self.output, self.detail, self.type, None,
                            self.name, self.space, self.default))
        else:
            spaces = self.spaces
            if self.defaultSpace(self.type)==None:
                spaces = None
            self.params.append((self.output, self.detail, self.type,
                               self.arraylen, self.name, spaces, self.default))
        self.arraylen = None
        self.name = ""
        self.space = None
        self.default = ""
        self.spaces   = []
        

    def switchFile(self, cppline):
        """Switch to another file.

        This method is called when a preprocessor line (starting with #)
        is read. This line contains the current file name and the line number.
        """
        f = cppline.strip().split(" ")
        linenr = int(f[1])
        filename = f[2][1:-1]
        self.filename = filename
        self.linenroffset = self._scanner.get_line_number()-linenr+1

        
# _SLfilter
class _SLfilter:
    """Used by the sltokenizer to filter the Shading Language source.

    Only the shader and function definitions remain, the bodies will
    be dropped. The filtered result will be used as input for the
    actual parser.
    """
    def __init__(self):
        # Current {}-depth (0 = outside any {..})
        self.curly_depth = 0
        # Current ()-depth
        self.bracket_depth = 0
        # Receives the source code that'll get passed to the parser
        self.SLsource = ""
        self.stream_enabled = True

        self.current_filename = None
        
    def eater(self, type, tok, start, end, line, filename):
        """Record only tokens with depth 0."""
        if filename!=self.current_filename:
            self.current_filename = filename
            if self.SLsource!="" and self.SLsource[-1]!="\n":
                self.SLsource+="\n"
            self.SLsource+='# %d "%s"\n'%(start[0],filename)
        
        if tok=="}":
            self.curly_depth-=1
            if self.curly_depth==0 and self.bracket_depth==0:
                self.stream_enabled = True
        elif tok==")":
            self.bracket_depth-=1

        # Always record newlines so that line numbers won't get messed up
        if self.stream_enabled or tok=="\n":
            self.SLsource+=tok
            
        if tok=="{":
            if self.curly_depth==0 and self.bracket_depth==0:
                self.stream_enabled = False
            self.curly_depth+=1
        elif tok=="(":
            self.bracket_depth+=1


# slparams
def slparams(slname, cpp="cpp", cpperrstream=sys.stderr):
    """Extracts the shader parameters from a Shading Language source file.

    The argument slname is the name of the shader source file (*.sl).
    cpp is the name of the preprocessor to use which has to be
    installed separately. If the preprocessor doesn't produce any data
    a PreprocessorNotFound exception is thrown. If cpp is None, then
    it's assumed that the input file is already preprocessed.  The
    error stream of the preprocessor is written to the object that's
    specified by cpperrstream which must have a write() method.  If
    cpperrstream is None, the error stream is ignored.  The function
    returns a list of 3-tuples, one for each shader found in the
    file. The tuple contains the type, the name and the parameters of
    the shader. The parameters are given as a list of 7-tuples
    describing each parameter. The tuple contains the following
    information (in the given order):

    - The output specifier (either "output" or an empty string)
    - The storage class ("uniform" or "varying")
    - The parameter type
    - The array length or None if the parameter is not an array
    - The name of the parameter
    - The space in which a point-like type was defined
    - The default value (always given as a string)
    """

    # Check if the file exists and can be accessed (by trying to open it)
    # If the file doesn't exist, an exception is thrown.
    dummy = file(slname)
    dummy.close()
    
    # Run the preprocessor on the input file.
    # The preprocessed data will be in slsrc.
    if cpp!=None:
        cmd = "%s %s"%(cpp, slname)
        stdin, stdout, stderr = os.popen3(cmd)
        slsrc = stdout.read()
        stdout.close()
        if cpperrstream!=None:
            errs = stderr.read()
            cpperrstream.write(errs)
        # No data? Then it's assumed that the preprocessor couldn't be found
        if len(slsrc)==0:
            raise PreprocessorNotFound("Calling '%s' didn't produce any data."%cmd)
        f = StringIO.StringIO(slsrc)
    else:
        f = file(slname)

    # ...and filter it, so that only the shader and function
    # definitions remain...
    filter = _SLfilter()
    sltokenize.tokenize(f.readline, filter.eater)
    f.close()

#    print filter.SLsource

    # Parse the filtered source code...
    scanner = _slparser._SLParserBaseScanner(filter.SLsource)
    parser = _SLParser(scanner)
    
#    return wrap_error_reporter(parser, "definitions")

    try:
        return getattr(parser, "definitions")()
    except _slparser.NoMoreTokens, err:
        raise NoMoreTokens, "No more tokens"
    except _slparser.SyntaxError, err:
        scanner = parser._scanner
        input = scanner.input
        cpplineno = scanner.get_line_number()
        lineno = cpplineno-parser.linenroffset
        # '+1' at the end to exclude the newline. If no newline is found
        # rfind() returns -1, so by adding +1 we get 0 which is what we want.
        start = input.rfind("\n", 0, err.charpos)+1
        end   = input.find("\n", err.charpos, -1)
        origline = input[start:end].replace("\t", " ")
        line = " "+origline.lstrip()
        errpos = err.charpos-start-(len(origline)-len(line))
#        print "%s^"%((errpos)*" ")
        msg = 'Syntax error in "%s", line %d:\n'%(parser.filename, lineno)
        msg += '\n%s\n%s^'%(line, errpos*" ")
        exc = SyntaxError(msg)
        exc.filename = parser.filename
        exc.lineno = lineno
        exc.line = line
        exc.errpos = errpos
        raise exc



# Setup local namespace for convertdefault()
_local_namespace = {}
exec "from sl import *" in _local_namespace
   
def convertdefault(paramtuple):
    """Converts the default value of a shader parameter into a Python type.

    paramtuple must be a 7-tuple as returned by slparams(). The
    function returns a Python object that corresponds to the default
    value of the parameter. If the default value can't be converted
    then None is returned. Only the functions that are present in the
    sl module are evaluated. If a default value calls a user defined
    function then None is returned.

    The SL types will be converted into the following Python types:

    - float  -> float
    - string -> string
    - color  -> vec3
    - point  -> vec3
    - vector -> vec3
    - normal -> vec3
    - matrix -> mat4

    Arrays will be converted into lists of the corresponding type.    
    """
    global _local_namespace
    
    typ = paramtuple[2]
    arraylen = paramtuple[3]
    defstr = paramtuple[6]

    # Replace {} with [] so that SL arrays look like Python lists
    defstr = defstr.replace("{","[").replace("}","]")
    # If the parameter is not an array, then create an array with one
    # element (to unify further processing). It will be unwrapped in the end
    if arraylen==None:
        defstr = "[%s]"%defstr
    # Evaluate the string to create "raw" Python types (lists and tuples)
    try:
        rawres = eval(defstr, globals(), _local_namespace)
    except:
        return None

    # Convert into the appropriate type...
    if typ=="float":
        try:
            res = map(lambda x: float(x), rawres)
        except:
            return None
    elif typ=="color" or typ=="point" or typ=="vector" or typ=="normal":
        try:
            res = map(lambda x: cgtypes.vec3(x), rawres)
        except:
            return None
    elif typ=="matrix":
        try:
            res = map(lambda x: cgtypes.mat4(x), rawres)
        except:
            return None
    elif typ=="string":
        try:
            res = map(lambda x: str(x), rawres)
        except:
            return None

    if arraylen==None:
        res = res[0]

    return res

######################################################################

if __name__=="__main__":
    pass
