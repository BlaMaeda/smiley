#!/usr/bin/env python
#
#   Copyright (C) 2008  Don Smiley  ds@sidorof.com

#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.

#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.

#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.

#   See the LICENSE file included in this archive
#

"""
This module implements the components for grammatical evolution.

"""
import datetime

from copy import deepcopy
from threading import Thread
from Queue import Queue
from random import shuffle
import gc

from genotypes import Genotype
from fitness import FitnessList, Fitness, Replacement
from utilities import rand_int

STATEMENT_FORMAT = '<S'

class GrammaticalEvolution(object):
    """
    This class comprises the overall process of generating genotypes,
    expressing the genes as programs using grammer and evalating the fitness of
    those members.
    """

    def __init__(self):
        self._population_size = 0
        self._crossover_rate = 0.2
        self._children_per_crossover = 2
        self._mutation_type = 's'
        self._mutation_rate = 0.02
        self._wrap = True
        self._extend_genotype = True

        self._completion_type = None
        self._start_gene_length = None
        self._max_gene_length = None
        self._max_program_length = None
        self._generation = 0
        self._max_generations = 0
        self.fitness_list = FitnessList('center')
        self._fitness_fail = -1000.0
        self._max_fitness_rate = .5
        self._completion_criteria = "g"
        self._maintain_history = True

        self._queue_size = 5
        self._garbage_collection = 30
        self._timeouts = [20, 3600]
        self.current_g = None
        self._fitness_selections = []
        self._replacement_selections = []

        self.bnf = {}
        self.population = []
        self.que = None

        self._history = []

    def set_population_size(self, size):
        """
        This function sets the total number of genotypes in the population.
        This program uses a fixed size population.

        """

        if isinstance(size, int) and size > 0:
            self._population_size = size
            i = len(self.fitness_list)
            while i < size:
                self.fitness_list.append([0.0, i])
                i += 1
        else:
            raise ValueError("""
                population size, %, must be an int above 0""" % (size))

    def get_population_size(self):
        """
        This function returns total number of genotypes in the
        population.

        """

        return self._population_size

    def set_genotype_length(self, start_gene_length,
                                    max_gene_length=None):
        """
        This function sets the initial size of the binary genotype.  An
        optional max_gene_length can be entered as well.  This permits the
        genotype to grow during the mapping process of the genotype to a
        program.  The lengths are the length of the decimal genotypes, which
        are therefor 8 times longer the binary genotypes created.

        """

        if max_gene_length is None:
            start_gene_length = max_gene_length
        if not isinstance(start_gene_length, int):
            raise ValueError("start_gene_length, %s must be an int" % (
                start_gene_length))
        if not isinstance(max_gene_length, int):
            raise ValueError("max_gene_length, %s must be an int" % (
                max_gene_length))
        if start_gene_length > max_gene_length:
            raise ValueError("""max_gene_length, %s cannot be smaller
                than start_gene_length%s""" % (
                    max_gene_length, start_gene_length))

        #   Silently adjusts lengths to a multiple of 8
        self._start_gene_length = start_gene_length
        self._max_gene_length = max_gene_length

    def get_genotype_length(self):
        """
        This function returns a tuple with the length the initial genotype and
        the maximum genotype length permitted.

        """

        return (self._start_gene_length, self._max_gene_length)

    def set_extend_genotype(self, true_false):
        """
        This function sets whether the genotype is extended during the gene
        mapping process.

        """

        if isinstance(true_false, bool):
            self._extend_genotype = true_false
        else:
            raise ValueError("Extend genotype must be True or False")

    def get_extend_genotype(self):
        """
        This function returns whether the genotype is extended during the gene
        mapping process.

        """

        return self._extend_genotype

    def set_wrap(self, true_false):
        """
        This function sets whether the genotype is wrapped during the gene
        mapping process.  Wrapping would occur in the iterative process of
        getting the next codon is the basis for the variable selection process.
        If wrapped, when all the codons in the genotype are exhausted, the
        position marker is wrapped around to the first codon in the sequence
        and goes on.

        """

        if isinstance(true_false, bool):
            self._wrap = true_false
        else:
            raise ValueError("Wrap must be True or False")

    def get_wrap(self):
        """
        This function returns whether the genotype is wrapped during the gene
        mapping process.  Wrapping would occur in the iterative process of
        getting the next codon is the basis for the variable selection process.
        If wrapped, when all the codons in the genotype are exhausted, the
        position marker is wrapped around to the first codon in the sequence
        and goes on.

        """

        return self._wrap

    def set_bnf(self, bnf):
        """
        This function parses up a bnf and builds a dictionary. The incoming
        format is designed to follow a format of:  <key> ::= value1 | value2
        \n. The following lines can also hold additional values to accommodate
        longer choices.

        In addition, a set of statements are marked with a key
        starting with "<S".  These are treated differently in that spaces are
        not automatically stripped from the front.  This enables python
        oriented white space to be honored.

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
        self.bnf = bnf_dict

    def get_bnf(self):
        """
        This function returns the Backus Naur form of variables that are used
        to map the genotypes to the generated programs.

        """

        return self.bnf

    def set_completion_criteria(self, completion_type, *params):
        """
        This function sets the completion criteria.  The choices are either g
        for a maximum generation count and f for a fitness threshhold.

        If the completion type is 'g', then the parameter following should be
        an int that is the maximum number of generations that will be run.

        If the completion type is 'f', then there must be a parameter
        following that 'f'.  The parameter is the threshhold value that
        must be reached in order to complete evolutionary process.

        Finally, if a combination of factors can be helpful, it is possible
        to set a completion type of 'F' and the appropriate criteria, and then
        manually set the max generations that could be reached before halting
        the process as a failure, using set_max_generations(maxgen).

        """

        if completion_type not in ['g', 'f']:
            raise ValueError("""
                completion type must be either g(enerations) or
                    f(itness value)""")
        if completion_type == "g":
            if params is not None:
                self.set_max_generations(params[0])
            self._completion_type = "g"
        else:
            self._completion_type = "f"
            if params:
                if len(params) == 1:
                    self.set_fitness_type(params[0])
                elif len(params) == 2:
                    self.set_fitness_type(params[0], params[1])

    def set_maintain_history(self, true_false):
        """
        This function sets a flag to maintain a history of fitness_lists.

        """
        if isinstance(true_false, bool):
            self._maintain_history = true_false
        else:
            raise ValueError("Maintain history must be True or False")

    def get_maintain_history(self):
        """
        This function returns a flag indicating whether the fitness list is
        retained for each generation.

        """

        return self._maintain_history

    def set_max_program_length(self, max_program_length):
        """
        This function sets the maximum length that a program can attain before
        the genotype is declared a failure.

        """

        errmsg = """The maximum program length, %s must be an int value
                    """ % (max_program_length)
        if not isinstance(max_program_length, int):
            raise ValueError(errmsg)

        self._max_program_length = max_program_length

    def get_max_program_length(self):
        """
        This function gets the maximum length that a program can attain before
        the genotype is declared a failure.

        """

        return self._max_program_length

    def set_fitness_fail(self, fitness_fail):
        """
        This function sets the fitness fail value that will be applied to
        fitness functions that are deemed failure.  Failure would be programs
        that fail due to overflows, or programs that grow to greater than
        maximum program length, syntax failures, or other reasons.

        """

        errmsg = """The fitness_fail, %s must be a float value
                    """ % (fitness_fail)
        if not isinstance(fitness_fail, float):
            raise ValueError(errmsg)

        self._fitness_fail = fitness_fail

    def get_fitness_fail(self):
        """
        This function returns the value of fitness if the program is a failure.

        """

        return self._fitness_fail

    def set_mutation_type(self, mutation_type):
        """
        This function sets the mutation type.  The choices are s(ingle),
        m(ultiple).  If the choice is "s", then the mutation rate is applied
        as a choice of whether to alter 1 bit on a gene or not.  If the choice
        is "m", then the process applies the rate as the probability that a bit
        will be changed as it walks the gene.  In short, "s", means that if the
        gene is mutated, it will take place once.  Otherwise, the gene could be
        mutated multiple times.

        """

        errmsg = "The mutation type must be either 's' or 'm'."
        if mutation_type not in ["s", "m"]:
            raise ValueError(errmsg)

        self._mutation_type = mutation_type

    def get_mutation_type(self):
        """
        This function returns the mutation type.  See set_mutation_type for a
        more complete explanation.

        """

        return self._mutation_type

    def set_mutation_rate(self, mutation_rate):
        """
        This function sets the mutation rate that will be applied to members
        selected into the fitness pool and to newly generated children.  Note
        that the mutation rate should be vastly different depending upon the
        mutation type that you have selected.  If the mutation type is 's',
        then the rate is the probability that the genotype will be mutated.  If
        the mutation type is 'm', then the rate is the probability that the any
        given bit in the genotype will be altered.  Because of that, the
        mutation rate should be significantly lower than the rate used with a
        mutation type of 's'.

        """

        errmsg = """The mutation rate, %s must be a float value
                    from 0.0 to 1.0""" % (mutation_rate)
        if not isinstance(mutation_rate, float):
            raise ValueError(errmsg)
        if not (0.0 <= mutation_rate <= 1.0):
            raise ValueError(errmsg)

        self._mutation_rate = mutation_rate

    def get_mutation_rate(self):
        """
        This function gets the mutation rate that will be applied to members
        selected into the fitness pool and to newly generated children.  Note
        that the mutation rate should be vastly different depending upon the
        mutation type that you have selected.  If the mutation type is 's',
        then the rate is the probability that the genotype will be mutated.  If
        the mutation type is 'm', then the rate is the probability that the any
        given bit in the genotype will be altered.  Because of that, the
        mutation rate should be significantly lower than the rate used with a
        mutation type of 's'.

        """

        return self._mutation_rate

    def set_crossover_rate(self, crossover_rate):
        """
        This function sets the probablity that will be
        applied to members selected into the fitness pool.

        """

        errmsg = """The crossover rate, %s must be a float value
                    from 0.0 to 1.0""" % (crossover_rate)
        if not isinstance(crossover_rate, float):
            raise ValueError(errmsg)
        if not (0.0 <= crossover_rate <= 1.0):
            raise ValueError(errmsg)

        self._crossover_rate = crossover_rate

    def get_crossover_rate(self):
        """
        This function gets the probablity that will be applied to members
        selected into the fitness pool.

        """

        return self._crossover_rate

    def set_children_per_crossover(self, children_per_crossover):
        """
        This function sets the number of children that will generated from two
        parents.  The choice is one or two.

        """

        if children_per_crossover not in [1, 2]:
            raise ValueError(
                "The children per crossovermust be either 1 or 2.")
        self._children_per_crossover = children_per_crossover

    def get_children_per_crossover(self):
        """
        This function gets the number of children that will generated from two
        parents.

        """

        return self._children_per_crossover

    def set_max_generations(self, generations):
        """
        This function sets the maximum number of generations that will be run.

        """

        if isinstance(generations, int) and generations >= 0:
            self._max_generations = generations
        else:
            raise ValueError("""
                generations, %, must be an int 0 or greater""" % (generations))

        self._max_generations = generations

    def get_max_generations(self):
        """
        This function gets the maximum number of generations that will be run.

        """

        return self._max_generations

    def set_fitness_type(self, fitness_type, target_value=0.0):
        """
        This function sets whether the objective is to achieve as large a
        fitness value possible, small, or hit a target_value.  Therefor the
        choices are 'max', 'min', or 'center'.  If center is used, then a
        target value should be entered as well.  For example, suppose that you
        wanted to hit a target somewhere near zero.  Setting the target_value
        at .001 would cause the process to complete if a fitness value achieved
        .001 or less.

        """

        self.fitness_list.set_fitness_type(fitness_type)
        self.fitness_list.set_target_value(target_value)

    def get_fitness_type(self):
        """
        This function gets whether the objective is to achieve as large a
        fitness value possible, small, or hit a target_value.  Therefor the
        choices are 'max', 'min', or 'center'.  If center is used, then a
        target value should be entered as well.  For example, suppose that you
        wanted to hit a target somewhere near zero.  Setting the target_value
        at .001 would cause the process to complete if a fitness value achieved
        .001 or less.

        """

        return self.fitness_list.get_fitness_type()

    def set_max_fitness_rate(self, max_fitness_rate):
        """
        This function sets a maximum for the number of genotypes that can be
        in the fitness pool.  Since some fitness selection approaches can have
        a varying number selected, and since multiple selection approaches can
        be applied consequentially, there needs to be an ultimate limit on the
        total number.  The max fitness rate must be a value greater than zero
        and less than 1.0.

        """

        errmsg = """The max fitness rate, %s must be a float value
                    from 0.0 to 1.0""" % (max_fitness_rate)
        if not isinstance(max_fitness_rate, float):
            raise ValueError(errmsg)
        if not (0.0 <= max_fitness_rate <= 1.0):
            raise ValueError(errmsg)

        self._max_fitness_rate = max_fitness_rate

    def get_max_fitness_rate(self):
        """
        This function gets a maximum for the number of genotypes that can be
        in the fitness pool.  Since some fitness selection approaches can have
        a varying number selected, and since multiple selection approaches can
        be applied consequentially, there needs to be an ultimate limit on the
        total number.  The max fitness rate must be a value greater than zero
        and less than 1.0.

        """

        return self._max_fitness_rate

    def set_fitness_selections(self, *params):
        """
        This function loads the fitness selections that are to be used to
        determine genotypes worthy of continuing to the next generation.  There
        can be multiple selections, such as elites and tournaments.  See the
        section Fitness Selection for further information.

        """

        for fitness_selection in params:
            if isinstance(fitness_selection, Fitness):
                self._fitness_selections.append(fitness_selection)
            else:
                raise ValueError("Invalid fitness selection")

    def set_replacement_selections(self, *params):
        """
        This function loads the replacement selections that are used to
        determine genotypes are to be replaced.  Basically, it is the grim
        reaper. Multiple replacement types can be loaded to meet the criteria.
        The number replaced is governed by the fitness selection functions to
        ensure that the population number stays constant.

        """

        for replacement_selection in params:
            if isinstance(replacement_selection, Replacement):
                self._replacement_selections.append(replacement_selection)
            else:
                raise ValueError("Invalid replacement selection")

    def get_fitness_history(self, statistic='best_value'):
        """
        This funcion returns a list of values that represent historical values
        from the fitness history.  While there is a default value of
        'best_value', other values are 'mean', 'min_value', 'max_value',
        'worst_value', 'min_member', 'max_member', 'best_member', and
        'worst_member'. The order is from oldest to newest.

        """

        hist_list = []
        for fitness_list in self._history:
            hist_list.append(fitness_list.__getattribute__(statistic)())
        return hist_list

    def get_best_member(self):
        """
        This function returns the member that it is most fit according to the
        fitness list.  Accordingly, it is only functional after at least one
        generation has been completed.

        """

        return self.population[self.fitness_list.best_member()]

    def get_worst_member(self):
        """
        This function returns the member that it is least fit according to the
        fitness list.  Accordingly, it is only functional after at least one
        generation has been completed.

        """

        return self.population[self.fitness_list.worst_member()]

    def set_timeouts(self, preprogram, program):
        """
        This function sets the number of seconds that the program waits until
        declaring that the process is a runaway and cuts it off.  During the
        mapping process against the preprogram, due to random chance a function
        can be calling another function, which calls another, until the process
        becomes so convoluted that the resulting program will be completely
        useless. While the total length of a program can be guide to its
        uselessnes as well, this is another way to handle it. Since variables
        can also be generated during the running of the program there is a
        second variable for the running program. Clearly, the second value must
        be in harmony with the nature of the program that you are actually
        running. Otherwise, you will be cutting of your program prematurely.
        Note that the second timeout is checked only if the running program
        requests an additional variable.  Otherwise, it will not be triggered.

        """

        if isinstance(preprogram, int) and preprogram >= 0:
            self._timeouts[0] = preprogram
        else:
            raise ValueError("""
                timeout, %, must be an int 0 or above""" % (preprogram))

        if isinstance(preprogram, int) and preprogram >= 0:
            self._timeouts[1] = program
        else:
            raise ValueError("""
                timeout, %, must be an int 0 or above""" % (program))

    def get_timeouts(self):
        """
        This function returns the number of seconds that must elapse before
        the mapping process cuts off the process and declares that the genotype
        is a failure.  It returns a tuple for the number of seconds for the
        preprogram and the program itself.

        """

        return self._timeouts

    def set_queue_size(self, size):
        """
        This function sets the queue size of genotype fitness functions that
        would be run simultaneously.  The default is zero, which indicates that
        it is not used.  It should be regarded as experimental, and used only
        on that basis.

        """

        if isinstance(size, int) and size >= 0:
            self._queue_size = size
        else:
            raise ValueError("""
                queue size, %, must be an int 0 or above""" % (size))

    def get_queue_size(self):
        """
        This function gets the queue size of genotype fitness functions that
        would be run simultaneously.  The default is zero, which indicates that
        it is not used.  It should be regarded as experimental, and used only
        on that basis.

        """

        return self._queue_size

    def set_garbage_collection(self, interval):
        """
        This function sets how often forced garbage collection takes place.
        Although python normally handles this automatically, there is a lot of
        string creation, deletion, replacements, this makes enforces this
        process to hopefully improve performance.  For example, by setting this
        interval to 10, for every 10 times that a genotype is mapped and a
        fitness function is calculated, garbage collection is fored.

        """

        if isinstance(interval, int) and interval >= 0:
            self._garbage_collection = interval
        else:
            raise ValueError("""
                interval, %, must be an int 0 or above""" % (interval))

    def get_garbage_collection(self, interval):
        """
        This function gets how often forced garbage collection takes place.
        Although python normally handles this automatically, there is a lot of
        string creation, deletion, replacements, this makes enforces this
        process to hopefully improve performance.  For example, by setting this
        interval to 10, for every 10 times that a genotype is mapped and a
        fitness function is calculated, garbage collection is fored.

        """

        if isinstance(interval, int) and interval >= 0:
            self._garbage_collection = interval
        else:
            raise ValueError("""
                interval, %, must be an int 0 or above""" % (interval))

    def _start_thread_pool(self):
        """
        This function, an experimental one, starts the queue that holds a pool
        of threads for processing.  If used, it starts each fitness function in
        a separate thread, adding an addition genotype as one finishes until
        all have been processes.

        """

        def fitness_processor():
            """
            This gets a thread from the queue and computes fitness.

            """
            while True:
                self.current_g = self.que.get()

                self.current_g.compute_fitness()
                self.fitness_list[self.current_g.member_no][0] = \
                    self.current_g.get_fitness()
                self.que.task_done()

        self.que = Queue()
        count = 0
        while count < self._queue_size:
            thrd = Thread(target=fitness_processor)
            thrd.start()
            count += 1

    def _compute_fitness(self):
        """
        This function runs the process of computing fitness functions for each
        genotype and calculating the fitness function.

        """

        if self._queue_size > 0:
            if self.que is None:
                self._start_thread_pool()

            for gene in self.population:
                starttime = datetime.datetime.now()
                #print "Starting member G %s: %s at %s" % (
                #    self._generation, gene.member_no,
                #    starttime.strftime('%m/%d/%y %H:%M'))
                self.que.put(gene)
                self._check_garbage(gene.member_no, self._garbage_collection)

            self.que.join()
        else:
            for gene in self.population:
                starttime = datetime.datetime.now()
                #print "Starting member G %s: %s at %s" % (
                #    self._generation, gene.member_no,
                #    starttime.strftime('%m/%d/%y %H:%M'))
                gene.starttime = starttime
                self.current_g = gene
                gene.compute_fitness()
                #print "fitness=%s" % (gene.get_fitness())
                self.fitness_list[gene.member_no][0] = gene.get_fitness()
                self._check_garbage(gene.member_no, self._garbage_collection)

        gc.collect()
        #   end of generation

    @staticmethod

    def _check_garbage(member_no, interval):
        """
        This function checks the interval for garbage collection and
        initiates it if it is time.

        """

        if interval > 0:
            if member_no % interval == 0:
                gc.collect()

    def run(self):
        """
        Once the parameters have all been set governing the course of the
        evolutionary processing, this function starts the process running.  It
        will continue until it the completion criteria have been set.

        """
        self._generation = 0
        while True:
            self._compute_fitness()
            if self._maintain_history:
                self._history.append(deepcopy(self.fitness_list))

            if self._continue_processing():
                self._perform_endcycle()
                #print "mean=", self.fitness_list.mean()
                #print "median=", self.fitness_list.median()
                print "best_value=", self.fitness_list.best_value()
                #print "stddev=", self.fitness_list.stddev()
                self._generation += 1
            else:
                print "Finished at generation ", self._generation
                break

        return self.fitness_list.best_member()

    def create_genotypes(self):
        """
        This function creates a genotype using the input parameters for each
        member of the population, and transfers operating parameters to the
        genotype for running the fitness functions.

        """

        member_no = 0
        while member_no < self._population_size:
            gene = Genotype(self._start_gene_length,
                        self._max_gene_length,
                        member_no)
            #   a local copy is made because variables
            #   can be saved within the local_bnf
            gene.local_bnf = deepcopy(self.bnf)
            gene.local_bnf['<member_no>'] = [gene.member_no]
            gene._max_program_length = self._max_program_length
            gene._fitness = self._fitness_fail
            gene._fitness_fail = self._fitness_fail
            gene._extend_genotype = self._extend_genotype
            gene._max_program_length = self._max_program_length
            gene._timeouts = self._timeouts
            gene._wrap = self._wrap
            self.population.append(gene)
            member_no += 1

    def _perform_endcycle(self):
        """
        This function runs after each member of the population has computed
        their fitness function.  Then, the fitness selection objects will
        evaluate those members according to their respective criteria and
        develop a pool of members that will potentially survive to the next
        generation. Crossovers will take place from that pool and each member
        will be subject to the possibility of mutatuting.  Finally, a
        replacement process will find which members should be replaced. The
        fitness pool will then replace those members.

        """

        fitness_pool = self._evaluate_fitness()
        child_list = self._perform_crossovers(fitness_pool)

        fitness_pool.extend(child_list)
        #shuffle(fitness_pool)
        self._perform_mutations(fitness_pool)
        self._perform_replacements(fitness_pool)

    def _evaluate_fitness(self):
        """
        This function evaluates the fitness of the members in the light of the
        fitness criteria functions.  It returns a list of members that will be
        used for crossovers and mutations.

        """

        flist = []
        total = self._max_fitness_rate * float(self._population_size)
        count = 0
        for fsel in self._fitness_selections:
            fsel.set_fitness_list(self.fitness_list)
            for i in fsel.select():
                flist.append(i)
                count += 1
                if count == total:
                    #   Done
                    break

        flist1 = []
        for member_no in flist:
            flist1.append(deepcopy(self.population[member_no]))

        return flist1

    def _perform_crossovers(self, flist):
        """
        This function accepts a list genotypes that are to be crossed.  The
        list is processed two at a time, and a child list holding the offspring
        is returned.  The _children_per_crossover indicator governs whether two
        children are produced or one.

        """

        child_list = []
        length = len(flist)
        if length % 2 == 1:
            length -= 1

        if length >= 2:
            for i in xrange(0, length, 2):
                parent1 = flist[i]
                parent2 = flist[i + 1]

                child1, child2 = self._crossover(parent1, parent2)
                if self._children_per_crossover == 2:
                    child_list.append(child1)
                    child_list.append(child2)
                else:
                    child_list.append(child1)

        return child_list

    def _crossover(self, parent1, parent2):
        """
        This function accepts two parents, randomly selects which is parent1
        and which is parent2.  Then, executes the crossover, and returns two
        children.

        """

        if not isinstance(parent1, Genotype):
            raise ValueError("Parent1 is not a genotype")
        if not isinstance(parent2, Genotype):
            raise ValueError("Parent2 is not a genotype")

        if rand_int(0, 1):
            child1 = deepcopy(parent1)
            child2 = deepcopy(parent2)
        else:
            child1 = deepcopy(parent2)
            child2 = deepcopy(parent1)

        child1_binary, child2_binary = self._crossover_function(
            child1.binary_gene, child2.binary_gene)

        child1.set_binary_gene(child1_binary)
        child1.generate_decimal_gene()
        child2.set_binary_gene(child2_binary)
        child2.generate_decimal_gene()

        return (child1, child2)

    @staticmethod

    def _crossover_function(child1_binary, child2_binary):
        """
        This function performs the actual crossover of material at a random
        point.

        """

        minlength = min(len(child1_binary), len(child2_binary))
        crosspoint = rand_int(1, minlength - 1)

        child1_binary, child2_binary = child1_binary[0:crosspoint] + \
                                       child2_binary[crosspoint:], \
                                       child2_binary[0:crosspoint] + \
                                       child1_binary[crosspoint:]
        return (child1_binary, child2_binary)

    def _perform_mutations(self, mlist):
        """
        This functions accepts a list of genotypes that are subject to
        mutation.  Each genotype is then put at risk for mutation and may or
        may not be mutated.

        """
        for gene in mlist:
            gene.mutate(self._mutation_rate, self._mutation_type)

    def _perform_replacements(self, fitness_pool):
        """
        This function accepts a list of members that will replace lesser
        performers.  The replacement process then applies the fitness pool to
        the population.

        """

        position = 0
        for rsel in self._replacement_selections:
            rsel.set_fitness_list(self.fitness_list)

            for replaced_no in rsel.select():
                replaced_g = self.population[replaced_no]
                if position < len(fitness_pool):
                    new_g = fitness_pool[position]
                    new_g.member_no = replaced_g.member_no
                    new_g._generation = self._generation + 1

                    #   update local bnf
                    new_g.local_bnf['<member_no>'] = [new_g.member_no]

                    self.population[new_g.member_no] = new_g
                    position += 1
                else:
                    break

    def _continue_processing(self):
        """
        This function, using the criteria for ending the evolutionary process
        after each generation, returns a flag of whether to continue or not.

        """

        status = True
        fitl = self.fitness_list
        if self._max_generations != 0:
            if self._generation > self._max_generations:
                status = False

        if fitl.get_target_value() is not None:
            if fitl.get_fitness_type() == 'max':
                if fitl.best_value() >= fitl.get_target_value():
                    status = False
            elif fitl.get_fitness_type() == 'min':
                if fitl.best_value() <= fitl.get_target_value():
                    status = False
            elif fitl.get_fitness_type() == 'center':
                if fitl.best_value() <= fitl.get_target_value():
                    status = False
        return status
