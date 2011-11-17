import sys
import random as random_module
import math
import pylab
from random import randint, random, choice
from copy import copy
from pprint import pprint
from problem import Problem
from grammar import Grammar
from crom import Crom

class Pair:
    def __init__(self, i, f):
        self.indice, self.fitness = i, f
    def __cmp__(self, p2):
        return -1 if self.fitness < p2.fitness else 1
    def __str__(self):
        return "%i-%f" % (self.indice, self.fitness)

class Poblacion:
    def __init__(self, n, l, problem, grammar, dict_meta, 
                 ml=None, pc=0.8, pm=0.05, tm='simple', bg=0.1, elit=True, cm="homologous"):
        random_module.seed()
        
        self._n = n
        self._l = l
        if ml == None:
            self._max_length = 2*l
        else:
            self._max_length = ml
        self._problem = problem
        self._grammar = grammar
        self._dict_meta = dict_meta
        self._prob_cruza = pc
        self._prob_mutacion = pm
        self._tipo_mutacion = tm
        assert(isinstance(self._tipo_mutacion, str))
        self._brecha_gen = bg
        self._elitismo = elit
        self._crossover_method = cm
        
        self._individuos = []
        self._next_generation = []
        for i in xrange(n):
            self._individuos.append(Crom(self._l, self._max_length, self._problem, 
                                         self._grammar, self._dict_meta,
                                         cross_meth = self._crossover_method))


####

    def _compute_fitness_list(self, individuos):
        fitness_list = []
        for i, indiv in enumerate(individuos):
            fitness_list.append(Pair(i, indiv.eval_fitness()))
        return fitness_list

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
        
        self._medians = []
        self._bests = []
        i = 0
        while i < maxit:
            self._fitness_list = sorted(self._compute_fitness_list(self._individuos))
            mejor_fitness = self.get_best_fitness()
            self._medians.append(self._get_median())
            self._bests.append(mejor_fitness)
            print "Generacion", i
            print "mejor fitness:", mejor_fitness, \
                  "fitness mediana:", self._fitness_list[self._n/2].fitness,\
                  "invalidos:", sum(1 for indiv in self._individuos if not indiv._valid),\
                  "promedio longitud", self._average_length()
            print "Mejor individuo:"
            print self.get_best_member()
            if mejor_fitness <= tol:
                break
            self.evolucionar()
            i += 1

        if i == maxit:
            print "No hubo convergencia"
        else:
            print "Convergencia en %i generaciones" % i

####

    def evolucionar(self):
        n = self._n
        cant_padres = int(round(n * self._brecha_gen))
        aptitudes = sorted(self._fitness_list)


        if (self._elitismo):
            self._next_generation.append(self._individuos[aptitudes[0].indice])

        i = 0
        while i < cant_padres:
            fin_ventana = n-i-1
            a = aptitudes[randint(0, fin_ventana)].indice
            b = aptitudes[randint(0, fin_ventana)].indice
            self._cruza(a, b)
            self._next_generation.append(self._individuos[a if randint(0,1) else b])
            i += 1
        
        while len(self._next_generation) < n:
            fin_ventana = n-i-1
            a = aptitudes[randint(0, fin_ventana)].indice
            b = aptitudes[randint(0, fin_ventana)].indice
            self._cruza(a, b)
            i += 1
        
        if len(self._next_generation) > n:
            ind_eliminar = randint(1 if self._elitismo else 0, n-1)
            self._next_generation = self._next_generation[:ind_eliminar] + self._next_generation[ind_eliminar+1:]
            
        
        
        self._individuos = [copy(i) for i in self._next_generation]
        self._next_generation = []

####

    def _cruza(self, a, b):
        parent_a = self._individuos[a]
        parent_b = self._individuos[b]

        if random() < self._prob_cruza:
            genes_child1, genes_child2 =  parent_a.crossover(parent_b)
        else:
            genes_child1 = copy(self._individuos[a]._genes)
            genes_child2 = copy(self._individuos[b]._genes)
        
        #if random() < self._prob_mutacion:
        self._mutate(genes_child1)
            #index = randint(0, len(genes_child1)-1) # XXX
            #genes_child1[index] = randint(0, 255)
        #if random() < self._prob_mutacion:
        self._mutate(genes_child2)
            #index = randint(0, len(genes_child2)-1) # XXX
            #genes_child2[index] = randint(0, 255)

        child1 = self._create_crom(genes_child1)
        child2 = self._create_crom(genes_child2)
        self._next_generation.append(child1)
        self._next_generation.append(child2)
        
####
    
    def _mutate(self, genes):
        if self._tipo_mutacion == 'simple':
            if random() < self._prob_mutacion:
                index = randint(0, len(genes)-1)
                genes[index] = randint(0, 255)
        elif self._tipo_mutacion == 'multiple':
            for index in xrange(len(genes)):
                if random() < self._prob_mutacion:
                    genes[index] = randint(0, 255)
        else:
            raise Exception, "Tipo de mutacion invalido: %s" % self._tipo_mutacion



####

    def _create_crom(self, genes_crom):
        return Crom(0, self._max_length, self._problem, 
                    self._grammar, self._dict_meta,
                    cross_meth = self._crossover_method, genes=genes_crom)

####

    def plot_medians(self):
        if len(self._medians) > 0:
            pylab.plot([math.log(me) for me in self._medians])
            pylab.show()
    
    def plot_stats(self):
        if len(self._medians) > 0:
            pylab.plot([math.log(me+1) for me in self._medians])
        if len(self._bests) > 0:
            pylab.plot([math.log(be+1) for be in self._bests])
        pylab.show()

####

    def show_all_indivs(self):
        pprint(self._individuos)
    def _get_median(self):
        median_index = self._n/2
        median = self._fitness_list[median_index].fitness
        return median
    def _average_length(self):
        return sum(i.length() for i in self._individuos) / float(len(self._individuos))
    def __getitem__(self, index):
        return self._individuos[index]
