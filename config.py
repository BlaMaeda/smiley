import os
from poblacion import Poblacion
from problem import Problem
from grammar import Grammar
from parser import parse_bnf

SECTION = 'Parametros'

SIZE_PARAMETER = 'size'
LENGTH_PARAMETER = 'length'
MAX_LENGTH_PARAMETER = 'max_length'
CROSSOVER_RATE_PARAMETER = 'probabilidad_de_cruza'
MUTATION_RATE_PARAMETER = 'probabilidad_de_mutacion'
GENERATION_GAP_PARAMETER = 'brecha_generacional'
ELITISM_PARAMETER = 'elitismo'
CROSSOVER_METHOD_PARAMETER = 'crossover_method'
EQUATION_PARAMETER = 'ecuacion'
FITNESS_FAIL_PARAMETER = 'fitness_fail'
LIMINF_PARAMETER = 'lim_inf'
LIMSUP_PARAMETER = 'lim_sup'
STEP_PARAMETER = 'step'
ADJUSTMENT_WEIGHT_PARAMETER = 'peso_ajuste'
SATISFACTION_WEIGHT_PARAMETER = 'peso_satisfaccion'
BNF_FILENAME_PARAMETER = 'bnf_filename'
BNF_META_FILENAME_PARAMETER = 'bnf_meta_filename'
META_SUFFIX = '_meta'

def config_to_population(config):
    # Default values
    n = 50
    l = 20
    ml = 50
    pc = 0.8
    pm = 0.05
    bg = 0.1
    elit = True
    cm = 'homologous'
    equation = '_y_ - x'
    li = 0.0
    ls = 5.0
    step = 0.1
    pa = 1.0
    ps = 1.0

    # Poblacion parameters
    if config.has_option(SECTION, SIZE_PARAMETER):
        n = config.getint(SECTION, SIZE_PARAMETER)
    if config.has_option(SECTION, LENGTH_PARAMETER):
        l = config.getint(SECTION, LENGTH_PARAMETER)
    if config.has_option(SECTION, MAX_LENGTH_PARAMETER):
        ml = config.getint(SECTION, MAX_LENGTH_PARAMETER)
    if config.has_option(SECTION, CROSSOVER_RATE_PARAMETER):
        pc = config.getfloat(SECTION, CROSSOVER_RATE_PARAMETER)
    if config.has_option(SECTION, MUTATION_RATE_PARAMETER):
        pm = config.getfloat(SECTION, MUTATION_RATE_PARAMETER)
    if config.has_option(SECTION, GENERATION_GAP_PARAMETER):
        bg = config.getfloat(SECTION, GENERATION_GAP_PARAMETER)
    if config.has_option(SECTION, ELITISM_PARAMETER):
        elit = config.getboolean(SECTION, ELITISM_PARAMETER)
    if config.has_option(SECTION, CROSSOVER_METHOD_PARAMETER):
        cm = config.get(SECTION, CROSSOVER_METHOD_PARAMETER)

    # Problem parameters
    if config.has_option(SECTION, EQUATION_PARAMETER):
        equation = config.get(SECTION, EQUATION_PARAMETER)
    if config.has_option(SECTION, FITNESS_FAIL_PARAMETER):
        ff = config.getfloat(SECTION, FITNESS_FAIL_PARAMETER)
    if config.has_option(SECTION, LIMINF_PARAMETER):
        li = config.getfloat(SECTION, LIMINF_PARAMETER)
    if config.has_option(SECTION, LIMSUP_PARAMETER):
        ls = config.getfloat(SECTION, LIMSUP_PARAMETER)
    if config.has_option(SECTION, STEP_PARAMETER):
        step = config.getfloat(SECTION, STEP_PARAMETER)
    if config.has_option(SECTION, ADJUSTMENT_WEIGHT_PARAMETER):
        pa = config.getfloat(SECTION, ADJUSTMENT_WEIGHT_PARAMETER)
    if config.has_option(SECTION, SATISFACTION_WEIGHT_PARAMETER):
        ps = config.getfloat(SECTION, SATISFACTION_WEIGHT_PARAMETER)
    problem = Problem(equation, ff, li, ls, step, pa, ps)
    
    # Grammar parameters
    if config.has_option(SECTION, BNF_FILENAME_PARAMETER):
        bnf_filename = config.get(SECTION, BNF_FILENAME_PARAMETER)
    else:
        raise Exception, "%s parameter is required" % BNF_FILENAME_PARAMETER
    if not os.path.isfile(bnf_filename):
        raise Exception, "%s doesn't exist" % bnf_filename

    bnf = ''.join(open(bnf_filename, 'r').readlines())
    grammar = Grammar(parse_bnf(bnf))

    # Meta grammar parameters
    if config.has_option(SECTION, BNF_META_FILENAME_PARAMETER):
        bnf_meta_filename = config.get(SECTION, BNF_META_FILENAME_PARAMETER)
        if not os.path.isfile(bnf_meta_filename):
            raise Exception, "%s doesn't exist" % bnf_meta_filename
    else:
        bnffn = os.path.splitext(bnf_filename)
        bnf_meta_filename = bnffn[0] + META_SUFFIX + bnffn[1]
        if not os.path.isfile(bnf_meta_filename):
            raise Exception, "%s parameter is required" % BNF_META_FILENAME_PARAMETER

    bnf_meta = ''.join(open(bnf_meta_filename, 'r').readlines())
    dict_meta = parse_bnf(bnf_meta)
    dict_meta = dict((k,v[0]) for (k,v) in dict_meta.items())

    popul = Poblacion(n, l, problem, grammar, dict_meta, ml, pc, pm, bg, elit, cm)
    return popul
