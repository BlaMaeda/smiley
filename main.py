#!/usr/bin/python

import ConfigParser 
import getopt, sys
from config import config_to_population
from poblacion import Poblacion
from problem import Problem
from grammar import Grammar
from parser import parse_bnf


if __name__=='__main__':
    config = ConfigParser.RawConfigParser()
    opts, args = getopt.getopt(sys.argv[1:], "f:")

    params_filename = None
    for o, a in opts:
        if o == "-f":
            params_filename = a

    if params_filename is None:
        print("No filename specified")
        sys.exit(1)

    read_result = config.read(params_filename)
    if read_result == []:
        print("Invalid filename")
        sys.exit(1)
    
    poblacion = config_to_population(config)
    poblacion.ev_and_print(10000, 0.01)
    print poblacion.get_best_member()
