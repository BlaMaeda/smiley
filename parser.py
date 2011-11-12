import re

STATEMENT_FORMAT = '<S'
VARIABLE_FORMAT = '(\<[^\>|^\s]+\>)'
LEFT_DEL = '<'
RIGHT_DEL = '>'

bnf = """
<S>                 ::= <expr>
<expr>              ::= <expr> <biop> <expr> | x | <uop> <expr> | <real> | <int-const>
                        math.log(abs(x)) | <pow> | math.sin(<exprf>) | math.cos(<exprf> )
                        (<exprf>) | math.exp(<exprf>) 
<exprf>             ::= <real> | math.log(abs(x)) | <powf> | math.sin(x) | x 
                        <int-const> | math.exp(x) | math.cos(x)
<biop>              ::= + | - | * | /
<uop>               ::= + | -
<pow>               ::= pow(abs(<exprf>), <real>) | pow(abs(<exprf>), <int-const>)
<powf>              ::= pow(abs(x), <real>)
<real>              ::= <int-const>.<int-const>
<int-const>         ::= <int-const> | 1 | 2 | 3 | 4 | 5
                        6 | 7 | 8 | 9 | 0
"""

def strip_spaces(key, values):
    """
    This removes white space unless it is a statement
    """
    if key.startswith(STATEMENT_FORMAT):
        values = [value.rstrip()
            for value in values.split('|') if value]
    else:
        values = [value.strip()
            for value in values.split('|') if value]

    return values

def parse_bnf(bnf):

    bnf_dict = {}
    for item in bnf.split('\n'):
        if item.find('::=') >= 0:
            key, values = item.split('::=')
            key = key.strip()
            bnf_dict[key] = strip_spaces(key, values)
        elif item:
            values = bnf_dict[key]
            values.extend(strip_spaces(key, item))
            if key.startswith(STATEMENT_FORMAT):
                #   Convert statements back to string
                values = ['\n'.join(values)]
            bnf_dict[key] = values
        else:
            #   blank line
            pass
    return bnf_dict

def parse_program(bnf, indiv):
    if not isinstance(bnf, dict):
        msg = "bnf must be a dict"
        raise ValueError(msg)

    complete = False

    program = bnf['<S>'][0]
    prg_list = re.split(VARIABLE_FORMAT, program)
    i = 0
    j = 0
    while True:
        if i == len(prg_list):
            complete = True
            break
        item = prg_list[i]
        if item != '' and item[0] == LEFT_DEL and item[-1] == RIGHT_DEL:
            if j == len(indiv):
                break
            ind = indiv[j] % len(bnf[item])
            replacement = bnf[item][ind]
            replacement = re.split(VARIABLE_FORMAT, replacement)
            prg_list = prg_list[0:i] + replacement + prg_list[i+1:]
            j += 1
        else:
            i += 1
    
    program = ''.join(prg_list)
    
    return complete, j, program

bnf = parse_bnf(bnf)
