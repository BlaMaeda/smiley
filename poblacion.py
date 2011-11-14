import random as random_module
import math
import pylab
from random import randint, random, choice
from copy import copy
from pprint import pprint
from problem import Problem
from grammar import Grammar
from crom import Crom
from parser import parse_bnf #XXX

class Pair:
    def __init__(self, i, f):
        self.indice, self.fitness = i, f
    def __cmp__(self, p2):
        return -1 if self.fitness < p2.fitness else 1
    def __str__(self):
        return "%i-%f" % (self.indice, self.fitness)

class Poblacion:
    def __init__(self, n, l, problem, grammar, ml=None, pc=0.8, pm=0.05, bg=0.1, elit=True, cm="homologous"):
        random_module.seed()
        
        self._n = n
        self._l = l
        if ml == None:
            self._max_length = 2*l
        else:
            self._max_length = ml
        self._problem = problem
        self._grammar = grammar
        self._prob_cruza = pc
        self._prob_mutacion = pm
        self._brecha_gen = bg
        self._elitismo = elit
        self._next_generation = []
        self._crossover_method = cm
        
        self._individuos = []
        for i in xrange(n):
            self._individuos.append(Crom(l, self._max_length, self._problem, self._grammar,
                                         cross_meth = self._crossover_method))

####

    def _compute_fitness_list(self, individuos):
        fitness_list = []
        for i, indiv in enumerate(individuos):
            fitness_list.append(Pair(i, indiv.eval_fitness()))
        return fitness_list

####

    def _compute_fitness(self, ind_indiv): # XXX borrable
        fin, fitness = self._problem.eval_fitness(self._individuos[ind_indiv])
        if fin != -1: # individuo valido
            self._contract_indiv(ind_indiv, fin)
        return fitness

####

    def get_best_fitness(self):
        return min(self._fitness_list).fitness

####

    def get_best_index(self):
        return min(self._fitness_list).indice

####

    def get_best_member(self):
        best_indiv = self._individuos[self.get_best_index()]
        return best_indiv._program

####

    def ev_and_print(self, maxit, tol):
        assert(isinstance(maxit, int))
        
        self._medians = []
        i = 0
        while i < maxit:
            self._fitness_list = sorted(self._compute_fitness_list(self._individuos))
            mejor_fitness = self.get_best_fitness()
            self._medians.append(self._get_median())
            print "Generacion", i
            print "mejor fitness:", mejor_fitness, \
                  "fitness mediana:", self._fitness_list[self._n/2].fitness,\
                  "invalidos:", sum(1 for indiv in self._individuos if not indiv._valid),\
                  "promedio longitud", self._average_length()
            #self.show_all_indivs(); xxx = raw_input()
            #print str([len(j) for j in self._individuos])
            #print "fitness medio: ", sum(j.fitness for j in self._fitness_list)/float(self._n)
            #print [str(j) for j in sorted(self._fitness_list)]
            if mejor_fitness <= tol:
                break
            self.evolucionar()
            i += 1

        if i == maxit:
            print "No hubo convergencia"
        else:
            print "Convergencia en %i generaciones" % i

        self.plot_medians()

####

    def evolucionar(self):
        n = self._n
        cant_padres = round(n * self._brecha_gen)
        aptitudes = sorted(self._fitness_list)

        assert(len(self._next_generation) == 0)

        if (self._elitismo):
            self._next_generation.append(self._individuos[aptitudes[0].indice])
        
        i = 0
        while i < cant_padres:
            fin_ventana = n-i-1
            a = randint(0, fin_ventana)
            b = randint(0, fin_ventana)
            self._cruza(a, b)
            self._next_generation.append(self._individuos[a if randint(0,1) else b])
            i += 1
        
        while len(self._next_generation) < n:
            fin_ventana = n-i-1
            a = randint(0, fin_ventana)
            b = randint(0, fin_ventana)
            self._cruza(a, b)
            i += 1
        
        if len(self._next_generation) > n:
            ind_eliminar = randint(1 if self._elitismo else 0, n-1)
            self._next_generation = self._next_generation[:ind_eliminar] + self._next_generation[ind_eliminar+1:]
        
        assert(len(self._next_generation) == n)
        
        self._individuos = [copy(i) for i in self._next_generation]
        try:
            assert(self._individuos == self._next_generation)
        except: 
            for i in xrange(len(self._individuos)):
                indiv = self._individuos[i]
                ng = self._next_generation[i]
                try:
                    assert indiv == ng
                except:
                    raise AssertionError, ("index %i" % i)
                #print "Individuo", i
                #print "Genes:", str(indiv._genes)
                #print "Genes:", str(ng._genes)
                #print "Ext:", str(indiv._extended_cromosom)
                #print "Ext:", str(ng._extended_cromosom)
        self._next_generation = []

####

    def _cruza(self, a, b):
        assert(isinstance(a, int) and isinstance(b, int))
        parent_a = self._individuos[a]
        parent_b = self._individuos[b]

        if random() < self._prob_cruza:
            genes_child1, genes_child2 =  parent_a.crossover(parent_b)
            child1 = self._create_crom(genes_child1)
            child2 = self._create_crom(genes_child2)
        else:
            child1 = copy(self._individuos[a])
            child2 = copy(self._individuos[b])
        
        if random() < self._prob_mutacion:
            child1.mutate()
        if random() < self._prob_mutacion:
            child2.mutate()

        self._next_generation.append(choice([child1, child2]))
        
####

    def _create_crom(self, genes_crom):
        return Crom(0, self._max_length, self._problem, self._grammar, 
                    cross_meth = self._crossover_method, genes=genes_crom)

####

    def show_all_indivs(self):
        pprint(self._individuos)
    def _get_median(self):
        return self._fitness_list[self._n/2].fitness
    def plot_medians(self):
        if len(self._medians) > 0:
            pylab.plot([math.log(me) for me in self._medians])
            pylab.show()
    def _average_length(self):
        return sum(i.length() for i in self._individuos) / float(len(self._individuos))

#### XXX borrable de aca pa'delante

#bnf = ''.join(open("bnf_lineal.txt", "r").readlines())
#bnf = parse_bnf(bnf)
#grammar = Grammar(bnf)
#assert(isinstance(grammar, Grammar))
#eq1 = "4*x**2*_y''_ + 17*_y_&_y(1)_+1&_y'(1)_+0.5"
#eq2 = "_y''_ - _y_&_y(0)_-1&_y'(0)_-1"
#eq3 = "_y''_ + _y_&_y(0)_-1&_y'(0)_-2"
#eq4 = "_y'_ - ((2*x-_y_)/x)&_y(0.1)-20.1"
#eq5 = "_y_ - 2*math.sin(x) - math.cos(x)"
#problem = Problem(eq3, ff=100.0,li=0.0,ls=1.0,step=0.05)
#poblacion = Poblacion(250, 50, problem, grammar, ml=50, cm='one-point')
#poblacion.ev_and_print(1000, 0.001)

#bnf = ''.join(open("bnf_paper.txt", "r").readlines())
#grammar = Grammar(parse_bnf(bnf))
#probl = Problem(grammar,ec=eq3,li=0,ps=10.0)
#pobl = Poblacion(100, 20, probl, pm=0.1, tm='s')
