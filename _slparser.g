####################################################################
# Shading Language parser
#
# Copyright (C) 2003, Matthias Baas (baas@ira.uka.de)
#
# http://cgkit.sourceforge.net
#
# You may distribute under the terms of the BSD license, as
# specified in the file license.txt.
####################################################################
# Create parser with: yapps2 _slparser.g   ->  _slparser.py

parser _SLParserBase:
    ignore: "[ \n\t]+"
    token END: "$"
    token number: "[0-9]+(\.[0-9]*)?([eE][0-9]+)?|\.[0-9]+([eE][0-9]+)?"
    token stringconstant: "\"[^\"]*\""
    token shader_type: "(light|surface|volume|displacement|imager)" 
    token outputspec: "output"
    token type: "(float|string|color|point|vector|normal|matrix|void)"
    token detail: "(varying|uniform)"
    token singletypes: "(float|string)"
    token spacetypes: "(color|point|vector|normal|matrix)"
    token identifier: "[_a-zA-Z][_a-zA-Z0-9]*"
    token binop: "[+\-/^*.]"
    token preprocessorline: "#.*"

    ######################################################################

    # The main rule for the whole shader source
    rule definitions:                     {{ shaders = [] }}
                      (preprocessorline   {{ self.switchFile(preprocessorline) }}
                       | shader_def       {{ shaders.append(shader_def) }}
                       | function_def)* END {{ return shaders }}

    ######################################################################

    # A shader definition
    rule shader_def: shader_type  identifier  {{ self.newParams() }}
                     "\(" [formals] "\)"
                     "{" "}"       {{ return (shader_type,identifier,self.params) }}

    # A function definition
    rule function_def: [type] identifier "\(" [formals] "\)" "{" "}"

                     
    ######################################################################

    # Parameter list of a shader or function...
    # (the definition from the spec is slightly changed so that
    # a trailing ';' is also allowed)
    rule formals: formal_var_defs            
                  (";" [formal_var_defs])*

    rule formal_var_defs:                     {{ self.newType() }}
                          [outputspec  {{ self.output="output" }} ]
                          typespec
                          def_expressions
   
    rule typespec: [detail  {{ self.detail = detail }}]
                   type     {{ self.type = type }}

    rule def_expressions: def_expression ("," def_expression)*

    # (added array support)
    rule def_expression: identifier   {{ self.name = identifier }}
                         ["\[" number "\]"   {{ self.arraylen=int(number)}}]
                         [def_init  {{ self.default = def_init }}]
                                      {{ self.storeParam() }}

    rule def_init: "=" 
                   expression   {{ return expression }}


    ######################################################################


    # An expression (for the default values of the parameters)...
    rule expression:  {{ expr="" }}
                     (
                      primary        {{ expr+=primary }}
                      | "-" expression  {{ expr+="-"+expression }}
                      | typecast expression   {{ expr+=expression }}
                     )
                     (binop expression  {{ expr+=binop+expression}} )*
                                        {{ return expr }}
                     

    rule primary:   number          {{ return number }}
                  | stringconstant  {{ return stringconstant }}
#                  | identifier      {{ return identifier }}
                  | tuple           {{ return tuple }}
                  | array           {{ return array }}
                  | identifier_or_procedurecall   {{ return identifier_or_procedurecall }}

    # Rule tuple generalizes triple and sixteentuple...
    # (it's assumed that the shader source is syntactically correct,
    # so we don't have to check if the tuple really only contains 3 or 16
    # elements)
    rule tuple:  {{ tup = "" }}
                  "\(" expression   {{ tup+=expression }}
                     ("," expression {{ tup+=", "+expression }} )*
                  "\)"   {{ return "("+tup+")" }}

    rule typecast:   singletypes            
                   | spacetypes   
                   [spacename  {{ self.space = spacename[1:-1] }} ]  

    rule spacename: stringconstant   {{ return stringconstant }}

    # This rule is a combination if "identifier" and "procedurecall"
    # (because both start with the same token)
    rule identifier_or_procedurecall: 
                     identifier   {{ res = identifier; proc_args=""  }}
                     ["\("  [proc_args] "\)"  {{ res+="("+proc_args+")" }} ]
                                  {{ return res }}

    rule proc_args:  expression  {{ res = expression }}
                     ("," expression  {{ res+=", "+expression }} )*
                                 {{ return res }}


    rule array:  {{ arr = "" }}
                  "{" expression   {{ arr+=expression; self.appendSpace() }}
                    ("," expression {{ arr+=", "+expression; self.appendSpace() }} )*
                  "}"   {{ return "{"+arr+"}" }}
