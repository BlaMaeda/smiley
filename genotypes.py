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
This module implements genotypes for grammatical evolution.

"""
import datetime
import random
import re

from utilities import base10tobase2, base2tobase10, rand_int

STOPLIST = ['runtime_resolve', 'set_bnf_variable']
VARIABLE_FORMAT = '(\<([^\>|^\s]+)\>)'

class Genotype(object):
    """
    The Genotype class holds the genetic material.  It has the ability to run
    fitness functions and mutate.  It is an internal object, and so few aspects
    of it would be regarded as public.

    """

    def __init__(self, start_gene_length,
                        max_gene_length,
                        member_no):
        """
        This function initiates the genotype.  It must open with the starting
        gene length and the maximum gene length.  These lengths are the decimal
        lengths not the binary lengths.  In addition, the member number is
        needed, since the genotype creation process is controlled by the
        grammatic evolution class.

        """

        self.member_no = member_no
        self.local_bnf = {}
        self._max_program_length = None
        self._fitness = None
        self._fitness_fail = None
        self._wrap = True
        self._extend_genotype = True
        self.starttime = None
        self._timeouts = (0, 0)

        self._gene_length = start_gene_length
        self._max_gene_length = max_gene_length

        self.binary_gene = None
        self.decimal_gene = None
        self._generate_binary_gene(self._gene_length)
        self.generate_decimal_gene()

        self._position = (0, 0)

    def _generate_binary_gene(self, length):
        """
        This function creates a random set of bits.

        """

        geno = []
        count = 0
        while count < length * 8:
            geno.append(str(rand_int(0, 1)))
            count += 1
        self.binary_gene = ''.join(geno)

    def set_binary_gene(self, binary_gene):
        """
        This function sets the value of the binary gene directly.  This is
        used in the crossover and mutation functions.  There is an automatic
        adjustment to trim the length to a multiple of 8.

        """

        length = len(self.binary_gene)
        self.binary_gene = binary_gene[:length - length % 8]
        self._gene_length = len(self.binary_gene) / 8

    def generate_decimal_gene(self):
        """
        This function converts the binary gene to a usable decimal gene.

        """

        if self._gene_length == 0:
            raise ValueError("Invalid gene length")
        dec_geno = []
        for i in xrange(0, self._gene_length * 8, 8):
            item = self.binary_gene[i:i + 8]
            str_trans = base2tobase10(item)
            dec_geno.append(int(str_trans))

        self.decimal_gene = dec_geno
        self._position = (0, 0)

    @staticmethod

    def _dec2bin_gene(dec_gene):
        """
        This is a utility function that converts a decimal list to binary
        string.

        """

        bin_gene = []
        for item in dec_gene:
            bin_gene.append(base10tobase2(item, zfill=8))
        return ''.join(bin_gene)

    @staticmethod

    def _place_material(program, item, start_pos, end_pos):
        """
        This is a utility function that replaces a part of a string in a
        specific location.

        """

        if end_pos > len(program) - 1:
            raise ValueError("end_pos greater than len(program)")
        if start_pos > end_pos:
            raise ValueError("starting position > end postion")
        if start_pos == 0:
            if end_pos == len(program) - 1:
                program = item
            else:
                program = item + program[end_pos + 1:]
        else:
            if end_pos == len(program) - 1:
                program = program[:start_pos] + item
            else:
                program = program[:start_pos] + item + \
                            program[end_pos:]
        return program

    def runtime_resolve(self, item, return_type):
        """
        This function is callable by the generated program to enable
        additional values be pulled from genotype and BNF as the need arises
        during execution of the program.

        Usage is self.runtime_resolve('<myvariable>', return_type);

        The return type casts the result back to the format needed.  Supported
        return types are: 'int', 'float', 'str', and 'bool'.

        """
        value = self._map_variables(item, False)
        value = self._fmt_resolved_vars(value, return_type)
        return value

    @staticmethod

    def _fmt_resolved_vars(value, return_type):
        """
        This method formats the result for a resolved variable for use
        during runtime so that the information can fit into the context of what
        is running.

        Note that if the execute code was to be subclassed to a parser to avoid
        the use of exec, then this funtion should also be done as well, since
        it uses eval.

        """

        return_types = ['int', 'float', 'str', 'bool']

        def _conv_int(str_value):
            """
            This method attempts to convert string value to an int

            """
            try:
                value = int(str_value)
            except:
                value = None
            return value

        if return_type == 'str':
            return value
        elif return_type == 'int':
            int_value = _conv_int(value)
            if int_value is None:
                try:
                    value = eval(value)
                except:
                    msg = "Invalid evaluation of %s" % (value)
                    raise ValueError(msg)
            else:
                value = int_value
        elif return_type == 'float':
            try:
                value = float(value)
            except:
                #   allow normal error message to bubble up
                value = eval(value)
        elif return_type == 'bool':
            if value == 'True':
                value = True
            elif value == 'False':
                value = False
        else:
            msg = "return_type, %s must be in %s" % (value, return_types)
            raise ValueError(msg)

        return value

    def set_bnf_variable(self, variable_name, value):
        """
        This function adds a variable to the bnf.  The format is the name,
        typically bounded by <>, such as "<variable_name>", and the parameters
        are in the form of a list. The items in the list will be converted to
        strings, if not already.

        """

        if isinstance(value, list):
            self.local_bnf[variable_name] = value
        else:
            self.local_bnf[variable_name] = [str(value)]

    def resolve_variable(self, variable):
        """
        This function receives a variable and using the variable as a key
        looks it up in the local_bnf.  The list of possible values available
        are then used by the genotype via a codon to select a final value that
        would be used.

        """

        values = self.local_bnf[variable]
        try:
            value = self._select_choice(self._get_codon(), values)
        except:
            raise ValueError("""
                Failure to resolve variable: variable: %s values: %s
                """ % (variable, values))

        return str(value)

    def _map_variables(self, program, check_stoplist):
        """
        This function looks for a variable in the form of <variable>.  If
        check_stoplist is True, then there will be a check to determine if it
        is a run-time variable, and therefor will be resolved later.

        This process runs until all of the variables have been satisfied, or a
        time limit on the process has been reached.

        """

        def on_stoplist(item):
            """
            Checks the stop list for runtime variables

            """

            status = False
            for stopitem in STOPLIST:
                if item.find(stopitem) > -1:
                    status = True

            return status

        incomplete = True
        prg_list = re.split(VARIABLE_FORMAT, program)
        while incomplete:
            position = 0
            continue_map = False
            while position < len(prg_list):
                item = prg_list[position]
                if item.strip() == '':
                    del(prg_list[position])
                else:
                    if item[0] == "<" and item[-1] == ">":
                        #   check stop list
                        status = True
                        if check_stoplist and position > 0:
                            if on_stoplist(prg_list[position - 1]):
                                status = False
                        if status:
                            prg_list[position] = self.resolve_variable(item)
                            continue_map = True

                        del(prg_list[position + 1])
                    position += 1

            program = ''.join(prg_list)
            prg_list = re.split(VARIABLE_FORMAT, program)
            elapsed = datetime.datetime.now() - self.starttime

            #   Reasons to fail the process
            if check_stoplist:
                #   Program already running
                if elapsed.seconds > self._timeouts[1]:
                    continue_map = False
            else:
                #   Preprogram
                if elapsed.seconds > self._timeouts[0]:
                    continue_map = False

            if len(program) > self._max_program_length:
                #   Runaway process, let it fail on syntax
                continue_map = False

            if continue_map is False:
                return program

    def _get_codon(self):
        """
        This function gets the next decimal codon from the genotype.  If the
        end of the genotype is reached, and the wrap flag is True, then the
        position for the next codon is taken from the front again.
        Additionally, if wrapping has taken place and the extend_genotype flag
        is set, then the genotype will continue to grow in length until the
        max_gene_length is reached.

        If the wrap flag is not set, when the end of the genotype is
        reached, an error is raised.

        """

        i, position = self._position

        length = len(self.decimal_gene)
        if position < self._max_gene_length:
            codon = self.decimal_gene[i]

            if self._extend_genotype:
                if i != position:
                    self.decimal_gene.append(codon)
                    self._gene_length += 1

            i += 1
            position += 1
            if i == length:
                if self._wrap:
                    i = 0
                else:
                    raise StandardError("End of genotype without wrapping")
            self._position = (i, position)
            return codon
        else:
            #   Roll around position
            position = 0
            i = 0
            codon = self.decimal_gene[i]
            return codon

    def _update_genotype(self):
        """
        This function updates the binary genotype from the decimal gene if the
        genotype is extended.

        """

        self.set_binary_gene(self._dec2bin_gene(self.decimal_gene))

    def compute_fitness(self):
        """
        This function computes the fitness function.  The process consists
        mapping the codon to the program variables and running the resulting
        program and computing the fitness value.  In addition, the binary gene
        is updated if the decimal gene has been extended.

        """

        self._map_gene()
        if self._extend_genotype:
            self._update_genotype()

        return self._fitness

    def _map_gene(self):
        """
        This function applies the genotype information to build a program.
        Mapping the variables into the search space is an initial load, and can
        also iteratively accept values as the program that has been created
        executes via the runtime_resolve function.

        If for any reason the mapping process fails to create a viable
        program, or it takes too long, then the process is aborted and the
        fitness_fail value is applied instead.

        This function uses the print command to show the program that has been
        generated as well as print the fitness value.  It is expected that this
        will be converted to a more useful logging system at some point.

        """

        try:
            #print "==================================================="
            program = self._map_variables(self.local_bnf['<S>'][0], True)
            #print program
            self._execute_code(program)
            #print "==================================================="
        except:
            self.local_bnf['<fitness>'] = [str(self._fitness_fail)]

        self._fitness = float(self.local_bnf['<fitness>'][0])

    def _execute_code(self, program):
        """
        This function executes code that has been generated. This function
        would be subclassed if executing the code on a remote server, or
        swapping in a custom parser.

        """

        self.local_bnf['program'] = program
        program_comp = compile(program, '<program>', 'exec')
        exec program_comp

    def mutate(self, mutation_rate, mutation_type):
        """
        This is function randomly mutates a binary genotype by changing 1 to 0
        and vice versa.  It is not context-perserving at this level.

        """

        if mutation_type == 's':
            if random.random() < mutation_rate:
                self._single_mutate()
        else:
            self._multiple_mutate(mutation_rate)

    def _multiple_mutate(self, mutation_rate):
        """
        This function walks the gene and based upon the mutation rate will
        alter a bit.

        """

        length =  len(self.binary_gene)
        gene = [i for i in self.binary_gene]
        for i in xrange(length):
            if random.random() < mutation_rate:
                #   Mutate
                if gene[i] == '0':
                    gene[i] = '1'
                else:
                    gene[i] = '0'

        self.set_binary_gene(''.join(gene))
        self.generate_decimal_gene()

    def _single_mutate(self):
        """
        This function with a randomly selects a mutation point within the gene
        and changes a 1 to 0, or vice versa.

        """

        position = rand_int(0, self._gene_length * 8 - 1)
        if self.binary_gene[position] == '0':
            bit = '1'
        else:
            bit = '0'
        self.set_binary_gene(self._place_material(self.binary_gene, bit,
                                position, position))
        self.generate_decimal_gene()

    @staticmethod

    def _select_choice(codon, selection):
        """
        This function, based upon the codon, makes a choice from the list.
        The determination is based upon the module of the codon to the length
        of the list of choices.  For example, if the codon is 10 and the list
        is 7 choices long, then the selection would be from selection[3].  This
        ensures that for every codon for every selection there is some choice.

        """

        if isinstance(selection, list):
            return selection[codon % len(selection)]
        else:
            msg = "selection. %s, must be a list" % (selection)
            raise ValueError(msg)

    def get_program(self):
        """
        This function returns the program that has been generated.  It is only
        valid after the gene has been mapped.

        """

        return self.local_bnf['program']

    def get_preprogram(self):
        """
        This function returns the prototype program to which the variables
        will be applied.

        """

        return self.local_bnf['<S>']

    def get_fitness(self):
        """
        This function returns the fitness value that has been created as a
        result of running the fitness function.

        """

        return self._fitness
