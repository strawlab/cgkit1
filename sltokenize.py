####################################################################
# RenderMan Shading Language Tokenizer
#
# Copyright (C) 2002, Matthias Baas (baas@ira.uka.de)
#
# http://cgkit.sourceforge.net
#
# You may distribute under the terms of the BSD license, as
# specified in the file license.txt.
####################################################################

import re

WHITESPACE = 0
NAME       = 1
NUMBER     = 2
STRING     = 3
NEWLINE    = 4
OPERATOR   = 5
CHARACTER  = 6
TYPE       = 7

# tokenize
def tokenize(readline, tokeater):
    """Reads a Shading Language input stream and creates tokens.

    The first parameter, readline, must be a callable object which
    provides the same interface as the readline() method of built-in
    file objects. Each call to the function should return one line of
    input as a string.

    The second parameter, tokeneater, must also be a callable object.
    It is called with six parameters: the token type, the token
    string, a tuple (srow, scol) specifying the row and column where
    the token begins in the source, a tuple (erow, ecol) giving the
    ending position of the token, the line on which the token was
    found and the filename of the current file.

    By default the filename argument is an empty string. It will only
    be the actual filename if you provide a preprocessed file stream
    as input (so you should first run cpp on any shader). The
    tokenizer actually expects preprocessed data as it doesn't handle
    comments.
    """

    regs =  ( (WHITESPACE, re.compile(r"[ \t]+")),
              (TYPE,       re.compile(r"float|point|vector|normal|matrix|color")),
              (NAME,       re.compile(r"[A-Za-z_][A-Za-z_0-9]*")),
              (NUMBER,     re.compile(r"[0-9]+(\.[0-9]+)?(E(\+|-)?[0-9]+)?")),
              (STRING,     re.compile(r"\"[^\"]*\"")),
              (OPERATOR,   re.compile(r"\+|-|!|\.|\*|/|\^|<|>|<=|>=|==|!=|&&|\|\||\?|:|=|\(|\)")),
              (NEWLINE,    re.compile(r"\n"))
            )

    linenr   = 0
    filename = ""
    while 1:
        # Read next line
        line = readline()
        # No more lines? then finish
        if line=="":
            break

        linenr+=1
        # Base for starting column...
        scolbase = 0

        # Process preprocessor lines...
        if line[0]=="#":
            try:
                f = line.strip().split(" ")
                linenr = int(f[1])-1
                filename = f[2][1:-1]
            except:
                pass
            continue

        s = line

        # Create tokens...
        while s!="":
            unmatched=1
            # Check all regular expressions...
            for r in regs:
                m=r[1].match(s)
                # Does it match? then the token is found
                if m!=None:
                    scol = m.start()
                    ecol = m.end()
                    tok = s[scol:ecol]
                    s   = s[ecol:]
                    tokeater(r[0], tok, (linenr, scolbase+scol), (linenr, scolbase+ecol), line, filename)
                    scolbase += ecol
                    unmatched=0
                    continue

            # No match? then report a single character...
            if unmatched:
                tok = s[0]
                tokeater(CHARACTER, tok, (linenr, scolbase), (linenr, scolbase+1), line, filename)
                s = s[1:]
                scolbase += 1



def _tokeater(type, s, start, end, line, filename):
    if type==WHITESPACE or type==NEWLINE:
        return
    print "Token:",type,s, start,end,'\t"%s"'%line.replace("\n",""),filename

######################################################################

if __name__=="__main__":

    f=open("test.sl")
    tokenize(f.readline, _tokeater)

