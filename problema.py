from parser import parse_bnf, parse_program
from copy import deepcopy
from pylab import subplot, plot, show
import math, scipy, re, sys

FORMA_ECUACION = "(_y.*?_)"
SEP_EC = "&"
ARG_COND = ".*\((.*)\).*"

class Problema:
    def __init__(self, bnf=None, bnf_filename=None, ec=None, ff=1e4, li=0, ls=5, step=0.1, pa=1.0, ps=1.0):
        if bnf == None and bnf_filename == None:
            self._bnf = """
            <S>                 ::= <expr>
            <expr>              ::= <expr> <biop> <expr> | <exprf> <biop> <exprf>
                                    <uop> <expr> | <uop> <exprf>
                                    math.log(abs(<expr>)) | math.exp(<expr>)
                                    math.log(abs(<exprf>)) | math.exp(<exprf>)
                                    <pow> | <powf> 
                                    math.sqrt(abs(<expr>)) | math.sqrt(abs(<exprf>))
                                    math.sin(<expr> ) | math.cos(<expr> ) 
                                    math.sin(<exprf> ) | math.cos(<exprf> ) 
                                    (<expr>) 
                                    (<exprf>) 
                                    <real> | <int-const>
                                    x | -x 
            <exprf>              ::= math.log(abs(x)) | math.exp(x)
                                    <powf>
                                    math.sqrt(abs(x))
                                    math.sin(x) | math.cos(x) 
                                    <real> | <int-const>
                                    x | -x 
            <biop>              ::= + | - | * | /
            <uop>               ::= + | -
            <pow>               ::= pow(abs(<expr> ), <real>) | pow(abs(<expr> ), <int-const>)
                                    pow(abs(<exprf> ), <real>) | pow(abs(<exprf> ), <int-const>)
            <powf>              ::= pow(abs(x), <real>)
            <real>              ::= <int-const>.<int-const>
            <int-const>         ::= <int-const> | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 0"""
        elif bnf != None:
            self._bnf = bnf
        elif bnf_filename != None:
            fd = open(bnf_filename, "r")
            self._bnf = ''.join(fd.readlines())
            fd.close()

        self._bnf = parse_bnf(self._bnf)

        if ec == None:
            self._ecuacion = "_y_  + x**2"
        else:
            partes_ec = re.split(SEP_EC, ec)
            self._ecuacion, self._condiciones = partes_ec[0], partes_ec[1:]
        self._generar_ecuacion()
        self._generar_condiciones()

        self._fitness_fail = ff
        self._lim_inf = li
        self._lim_sup = ls
        self._step = step
        self._peso_ajuste = pa
        self._peso_satisfaccion = ps

    def _generar_ecuacion(self):
        ecuacion = "abs("
        partes_ecuacion = re.split(FORMA_ECUACION, self._ecuacion)
        for parte in partes_ecuacion:
            if not parte:
                continue
            if parte[0:2] == "_y":
                orden = parte.count("'")
                if orden == 0:
                    ecuacion += "f(x)"
                else:
                    npuntos = orden + (1 if orden%2==0 else 2)
                    ecuacion += "scipy.derivative(f, x, dx=1e-3, n=%i, order=%i)" % (orden, npuntos)
            else:
                ecuacion += parte
        ecuacion += ")"
        self._ecuacion = ecuacion 

    def _generar_condiciones(self):
        condiciones = []
        for cond in self._condiciones:
            condicion = "abs("
            partes_condicion = re.split(FORMA_ECUACION, cond)
            for parte in partes_condicion:
                if not parte:
                    continue
                if parte[0:2] == "_y":
                    argumento = re.match(ARG_COND, parte).groups()[0]
                    orden = parte.count("'")
                    if orden == 0:
                        condicion += "f(%s)" % argumento
                    else:
                        npuntos = orden + (1 if orden%2==0 else 2)
                        condicion += "scipy.derivative(f, %s, dx=1e-3, n=%i, order=%i)" % (argumento, orden, npuntos)
                else:
                    condicion += parte
            condicion += ")"
            condiciones.append(condicion)
        self._condiciones = deepcopy(condiciones)


    def eval_fitness(self, fenotipo):
        assert(isinstance(fenotipo, list))
        x = self._lim_inf
        ajuste = 0.0
        parseo_completo, fin, program = parse_program(self._bnf, fenotipo)
        f = lambda x: eval(program)

        if parseo_completo:
            while x <= self._lim_sup:
                try:
                    fitness = eval(self._ecuacion)
                except:
                    return self._fitness_fail
                ajuste += fitness
                x += self._step
        else:
            return self._fitness_fail
        
        satisfaccion = 0.0
        mult_satisfaccion = round((self._lim_sup - self._lim_inf) / self._step)
        for condicion in self._condiciones:
            try:
                fitness = eval(condicion)
            except:
                return self._fitness_fail
            satisfaccion += fitness
                

        fitness = self._peso_ajuste*ajuste + self._peso_satisfaccion * mult_satisfaccion * satisfaccion
        return fin, fitness

    def compute_program(self, lista):
        _, __, program = parse_program(self._bnf, lista)
        return program

    def plotear(self, program):
        if (isinstance(program, str)):
            f = lambda x: eval(program)
            x = []
            i = self._lim_inf
            while i <= self._lim_sup:
                x.append(i)
                i += self._step
            plot(x, [f(i) for i in x])
            show()
        elif (isinstance(program, list)):
            for prog in program:
                f = lambda x: eval(prog)
                x = []
                i = self._lim_inf
                while i <= self._lim_sup:
                    x.append(i)
                    i += self._step
                plot(x, [f(i) for i in x])
            show()

