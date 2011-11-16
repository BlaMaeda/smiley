from parser import parse_bnf, parse_program
from copy import deepcopy
from pylab import subplot, plot, show
from grammar import Grammar
import math, scipy, re, sys

FORMA_ECUACION = "(_y.*?_)"
SEP_EC = "&"
ARG_COND = ".*\((.*)\).*"

class Problem:
    def __init__(self, ec, ff=1e4, li=0, ls=5, step=0.1, pa=1.0, ps=1.0):
        partes_ec = re.split(SEP_EC, ec)
        self._arg_ec = ec
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
                    ecuacion += "scipy.derivative(f, x, dx=0.01, n=%i, order=%i)" % \
                                 (orden, max(5, npuntos)) # TODO optimizar dx/order?
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
                        condicion += "scipy.derivative(f, %s, dx=0.01, n=%i, order=%i)" % \
                                      (argumento, orden, max(5, npuntos))
                else:
                    condicion += parte
            condicion += ")"
            condiciones.append(condicion)
        self._condiciones = deepcopy(condiciones)


    def eval_fitness(self, program):
        x = self._lim_inf
        ajuste = 0.0
        f = lambda x: eval(program)

        while x <= self._lim_sup:
            try:
                fitness = eval(self._ecuacion)
            except:
                return self._fitness_fail
            ajuste += fitness
            x += self._step
        
        satisfaccion = 0.0
        mult_satisfaccion = round((self._lim_sup - self._lim_inf) / self._step)
        for condicion in self._condiciones:
            try:
                fitness = eval(condicion)
            except:
                return self._fitness_fail
            satisfaccion += fitness
                

        fitness = self._peso_ajuste*ajuste + self._peso_satisfaccion * mult_satisfaccion * satisfaccion
        return fitness

    def get_fitness_fail(self):
        return self._fitness_fail

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

    def __str__(self):
        return self._arg_ec
