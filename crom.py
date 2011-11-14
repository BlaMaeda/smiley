import re
from random import randint, choice
from problem import Problem
from grammar import Grammar

START_SYMBOL = "<S>"
VARIABLE_FORMAT = '(\<[^\>|^\s]+\>)'
LEFT_DEL = '<'
RIGHT_DEL = '>'

class Crom:
    def __init__ (self, length, max_length, problem, grammar, **kwargs):
        self._length = length
        self._max_length = max_length
        self._problem = problem
        self._grammar = grammar

        if kwargs.has_key('genes'):
            self._genes = kwargs['genes']
            self._length = len(self._genes)
        else:
            self._genes = []
            for i in xrange(length):
                self._genes.append(randint(0, 255))

        self._crossover_method = kwargs.get('cross_meth', 'homologous')

        self._generate_program()

####

    def _generate_program(self):
        """ Genera el programa con los genes. Si termina
            antes, corta los genes que sobran; si no alcanza
            a terminar con lo que tiene, aumenta el cromosoma.
        """
        complete = False
        extended_cromosom = []

        program = self._grammar[START_SYMBOL, 0]
        prg_list = re.split(VARIABLE_FORMAT, program)
        i = 0 # prg_list
        j = 0 # genes
        while True:
            if i == len(prg_list):
                complete = True
                break

            item = prg_list[i]
            if item != '' and item[0] == LEFT_DEL and item[-1] == RIGHT_DEL:
                if j == len(self._genes):
                    if j == self._max_length:
                        break
                    else:
                        self._genes.append(randint(0, 255))
                choice = self._genes[j]
                replacement = self._grammar[item, choice]
                replacement = re.split(VARIABLE_FORMAT, replacement)
                prg_list = prg_list[0:i] + replacement + prg_list[i+1:]
                extended_cromosom.append("T" if len(replacement) == 1 else "S")
                j += 1
            else:
                i += 1
        
        self._genes = self._genes[:j]
        self._program = ''.join(prg_list)
        self._extended_cromosom = extended_cromosom
        self._valid = complete

####

    def eval_fitness(self):
        if self._valid:
            return self._problem.eval_fitness(self._program)
        else:
            return self._problem.get_fitness_fail()

####

    def mutate(self):
        mutable_indexes = [i for i, v in enumerate(self._extended_cromosom) 
                              if v == 'T']
        index = choice(mutable_indexes)
        index = randint(0, len(self._genes)-1)
        self._genes[index] = randint(0, 255)
        self._generate_program()

####

    def crossover(self, partner):
        if self._crossover_method == 'homologous':
            min_length = min(len(self._genes), len(partner._genes))
            if min_length > 1: # atadura con alambres para cuando los dos padrestienen longitud = 1
                crossover_point = randint(1, min_length-1)
            else:
                crossover_point = 1 # los hijos son iguales

            child1, child2 = self._genes[:crossover_point] +\
                             partner._genes[crossover_point:],\
                             partner._genes[:crossover_point] +\
                             self._genes[crossover_point:]
        elif self._crossover_method == 'one-point':
            cross_point_a = randint(1, max(len(self._genes)-1, 1))
            cross_point_b = randint(1, max(len(partner._genes)-1, 1))
            child1, child2 = self._genes[:cross_point_a] +\
                             partner._genes[cross_point_b:],\
                             partner._genes[:cross_point_b] +\
                             self._genes[cross_point_a:]

        return child1, child2

####

    def __eq__(self, crom):
        return self._genes             == crom._genes and\
               self._extended_cromosom == crom._extended_cromosom and\
               self._program           == crom._program and\
               self._grammar           == crom._grammar and\
               self._max_length        == crom._max_length and\
               self._valid             == crom._valid

####
    
    def __repr__(self):
        return self._program
    def __str__(self):
        return self._program
    def length(self):
        return len(self._genes)
