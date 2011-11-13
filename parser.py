import re
from grammar import Grammar

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
#    if isinstance(bnf, str): # filename
#        bnf = ''.join(open(bnf, "r").readlines())
#        print bnf
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

def parse_program(grammar, int_list):
    """ Given a grammar and a list with integers, returns 
        the program obtained from using the grammar
        with the list.
    """
    assert isinstance(grammar, Grammar),\
            "'grammar' must be an instance of Grammar."
    assert all(isinstance(i, int) for i in int_list),\
            "'int_list' must be a list of integers."

    complete = False
    extended_cromosom = []

    program = grammar['<S>', 0]
    prg_list = re.split(VARIABLE_FORMAT, program)
    i = 0
    j = 0
    while True:
        if i == len(prg_list):
            complete = True
            break

        item = prg_list[i]
        if item.strip() == '':
            i += 1
        elif item[0] == LEFT_DEL and item[-1] == RIGHT_DEL:
            if j == len(int_list):
                break
            ind = int_list[j]
            replacement = grammar[item, ind]
            replacement = re.split(VARIABLE_FORMAT, replacement)
            prg_list = prg_list[0:i] + replacement + prg_list[i+1:]
            extended_cromosom.append("T" if len(replacement) == 1 else "S")
            print extended_cromosom[j], len([1 for k in replacement if re.match(VARIABLE_FORMAT, k)])
            j += 1
        else:
            i += 1
    
    program = ''.join(prg_list)
    
    return complete, extended_cromosom, program

#bnf = ''.join(open("bnf_paper.txt", "r").readlines())
#bnf = parse_bnf(bnf)
#grammar = Grammar(bnf)
#indiv = [7, 1, 3, 3, 1, 1, 1, 7]
