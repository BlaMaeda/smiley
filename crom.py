import re
from random import randint, choice
from problem import Problem
from grammar import Grammar
from copy import deepcopy

START_SYMBOL = "<S>"
VARIABLE_FORMAT = '(\<[^\>|^\s]+\>)'
LEFT_DEL = '<'
RIGHT_DEL = '>'

class Crom:
    def __init__ (self, length, max_length, problem, grammar, dict_meta,
                  genes=None, cross_meth=None):
        self._length = length
        self._max_length = max_length
        self._problem = problem
        self._grammar = grammar
        self._dict_meta = dict_meta

        if genes is None:
            self._genes = []
            for i in xrange(length):
                self._genes.append(randint(0, 255))
        else:
            self._genes = deepcopy(genes)
            self._length = len(self._genes)

        if cross_meth is None:
            self._crossover_method = 'homologous'
        else:
            self._crossover_method = cross_meth

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
                if j == self._max_length:
                    break
                if j == len(self._genes):
                    self._genes.append(randint(0, 255))
                choice = self._genes[j]
                replacement = self._grammar[item, choice]
                replacement = re.split(VARIABLE_FORMAT, replacement)
                prg_list = prg_list[0:i] + replacement + prg_list[i+1:]
                extended_cromosom.append(self._dict_meta[item])
                j += 1
            else:
                i += 1
        
        self._genes = self._genes[:j]
        self._program = ''.join(prg_list)
        self._extended_cromosom = extended_cromosom
        self._valid = complete
        self._length = len(self._genes)

####

    def eval_fitness(self):
        if self._valid:
            fitness = self._problem.eval_fitness(self._program)
        else:
            fitness = self._problem.get_fitness_fail()

        return fitness

####

    def crossover(self, partner):
        if self._crossover_method == 'homologous':
            min_length = min(len(self._genes), len(partner._genes))
            crossover_point = randint(0, min_length-1)
            cross_point_a = cross_point_b = crossover_point
        elif self._crossover_method == 'one-point':
            cross_point_a = randint(0, len(self._genes)-1)
            cross_point_b = randint(0, len(partner._genes)-1)
        elif self._crossover_method == 'analogous':
            cross_point_a = randint(0, len(self._genes)-1)
            meta_a = self._extended_cromosom[cross_point_a]
            crossable_indexes = [i for i, v in enumerate(partner._extended_cromosom) 
                                  if v == meta_a]
            if crossable_indexes:
                cross_point_b = choice(crossable_indexes)
            else:
                # TODO Hay otras formas de hacer esto: no hacer la cruza,
                # elegir otro cross_point_a (dudoso), etc. Por ahora
                # usamos esto.
                cross_point_b = randint(0, len(partner._genes)-1)
        else:
            raise Exception, "unknown crossover method: %s" % self._crossover_method

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
               self._dict_meta         == crom._dict_meta and\
               self._length            == crom._length and\
               self._max_length        == crom._max_length and\
               self._valid             == crom._valid

####
    
    def __repr__(self):
        return self._program
    def __str__(self):
        return self._program
    def length(self):
        return len(self._genes)
    def __getitem__(self, index):
        return self._genes[index]
