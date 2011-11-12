import random as random_module
import math
from random import randint, random
from copy import deepcopy
from pprint import pprint
from problema import Problema

class Pair:
    def __init__(self, i, f):
        self.indice, self.fitness = i, f
    def __cmp__(self, p2):
        return -1 if self.fitness < p2.fitness else 1
    def __str__(self):
        return "%i-%f" % (self.indice, self.fitness)

class Poblacion:
    def __init__(self, n, l, problema, ml=None, pc=0.8, pm=0.05, tm='s', bg=0.1, elit=True, ff=1000.0):
        random_module.seed()

        self._n = n
        self._l = l
        if ml == None:
            self._max_length = 2*l
        else:
            self._max_length = ml
        self._problema = problema
        self._prob_cruza = pc
        self._prob_mutacion = pm
        self._tipo_mutacion = tm
        self._brecha_gen = bg
        self._elitismo = elit
        self._fitness_fail = ff
        self._next_generation = []

        self._individuos = []
        for i in xrange(n):
            indiv = [randint(0, 255) for j in xrange(l)]
            self._individuos.append(indiv)


    def _compute_fitness_list(self, individuos):
        fitness_list = []
        for i in xrange(len(individuos)):
            fitness_list.append(Pair(i, self._compute_fitness(i)))
        return fitness_list

    def _compute_fitness(self, ind_indiv):
        fin, fitness = self._problema.eval_fitness(self._individuos[ind_indiv])
        if fin != -1: # individuo valido
            self._contract_indiv(ind_indiv, fin)
        return fitness

    def _contract_indiv(self, ind_indiv, fin):
        self._individuos[ind_indiv] = self._individuos[ind_indiv][:fin]

    def get_best_fitness(self):
        return min(self._fitness_list).fitness

    def get_best_index(self):
        return min(self._fitness_list).indice

    def get_best_member(self):
        best = self._individuos[self.get_best_index()]
        return self._problema.compute_program(best)

    def ev_and_print(self, maxit, tol):
        assert(isinstance(maxit, int))
        
        i = 0
        mejor_fitness_anterior = 1e5
        while i < maxit:
            self._fitness_list = sorted(self._compute_fitness_list(self._individuos))
            mejor_fitness = self.get_best_fitness()
            print "mejor fitness: ", mejor_fitness, "fitness mediana: ", self._fitness_list[self._n/2].fitness
            #print "fitness medio: ", sum(j.fitness for j in self._fitness_list)/float(self._n)
            #print [str(j) for j in sorted(self._fitness_list)]
            assert(mejor_fitness <= mejor_fitness_anterior)
            if mejor_fitness <= tol:
                break
            self.evolucionar()
            i += 1
            mejor_fitness_anterior = mejor_fitness

        if i == maxit:
            print "No hubo convergencia"
        else:
            print "Convergencia en %i generaciones" % i


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
        
        self._individuos = [deepcopy(i) for i in self._next_generation]
        assert(self._individuos == self._next_generation)
        self._next_generation = []

    def _cruza(self, a, b):
        assert(isinstance(a, int) and isinstance(b, int))

        if random() < self._prob_cruza:
            hijo1, hijo2 = self._reproducir(a, b)
        else:
            hijo1 = deepcopy(self._individuos[a])
            hijo2 = deepcopy(self._individuos[b])
        
        if self._tipo_mutacion == 's':
            if random() < self._prob_mutacion:
                    self._mutar(hijo1)
            if random() < self._prob_mutacion:
                    self._mutar(hijo2)
        elif self._tipo_mutacion == 'm':
            self._mutar(hijo1)
            self._mutar(hijo2)
        
        self._next_generation.append(hijo1)
        self._next_generation.append(hijo2)
        
    def _reproducir(self, a, b):
        punto_de_cruza = randint(1, self._l-1)
        hijo1 = self._individuos[a][:punto_de_cruza] + self._individuos[b][punto_de_cruza:]
        hijo2 = self._individuos[b][:punto_de_cruza] + self._individuos[a][punto_de_cruza:]

        return hijo1, hijo2

    def _mutar(self, indiv):
        assert(isinstance(indiv, list))
        if self._tipo_mutacion == 's':
            punto_de_mutacion = randint(0, self._l-1)
            indiv[punto_de_mutacion] = randint(0, 255)
        elif self._tipo_mutacion == 'm':
            for i in xrange(len(indiv)):
                if random() <= self._prob_mutacion:
                    indiv[i] = randint(0, 255)

eq1 = "4*x**2*_y''_ + 17*_y_&_y(1)_+1&_y'(1)_+0.5"
eq2 = "_y''_ - _y_&_y(0)_-1&_y'(0)_-1"
eq3 = "_y''_ + _y_&_y(0)_-1&_y'(0)_-2"
probl = Problema(bnf_filename="bnf.txt",ec=eq3,li=0,ps=10.0)
pobl = Poblacion(100, 20, probl, pm=0.1, tm='s')
