#!/usr/bin/python

import ConfigParser 
import getopt, os, sys
from config import config_to_population
from poblacion import Poblacion
from problem import Problem
from grammar import Grammar
from parser import parse_bnf

DEFAULT_PARAMS_FILENAME = 'parameters.cfg'

if __name__=='__main__':
    config = ConfigParser.RawConfigParser()
    opts, args = getopt.getopt(sys.argv[1:], "f:p:")

    params_filename = None
    stats_to_plot = []
    for o, a in opts:
        if o == "-f":
            params_filename = a
        if o == "-p": # TODO or --plot
            stats_to_plot.append(a)

    if params_filename is None:
        params_filename = DEFAULT_PARAMS_FILENAME
        if not os.path.isfile(params_filename):
            print("No filename specified")
            sys.exit(1)
    else:
        if not os.path.isfile(params_filename):
            print("Invalid filename: %s" % params_filename)
            sys.exit(1)

    read_result = config.read(params_filename)
    if read_result == []:
        print("Invalid filename")
        sys.exit(1)
    
    poblacion = config_to_population(config)

    #------------------
    poblacion.ev_and_print(10000, 0.01)
    print poblacion.get_best_member()
    if "medians" in stats_to_plot:
        poblacion.plot_medians()
