from grammatical_evolution import GrammaticalEvolution
from fitness import *
from fitness import ReplacementTournament
import random

bnf =   """
<expr>              ::= <expr> <biop> <expr> | <uop> <expr> | <real> | <int-const>
                        math.log(abs(x)) | <pow> | math.sin(<exprf> ) | math.cos(<exprf> )
                        (<exprf> ) | math.exp(<exprf> ) 
<exprf>             ::= <real> | math.log(abs(x)) | <powf> | math.sin(x) | x 
                        <int-const> | math.exp(x) | math.cos(x)
<biop>              ::= + | - | * | /
<uop>               ::= + | -
<pow>               ::= pow(abs(<exprf> ), <real>) | pow(abs(<exprf> ), <int-const>)
<powf>              ::= pow(abs(x), <real>)
<plus>              ::= +
<minus>             ::= -
<real>              ::= <int-const>.<int-const>
<int-const>         ::= <int-const> | 1 | 2 | 3 | 4 | 5 | 6 |
                        7 | 8 | 9 | 0
<S>                 ::=
import math
import scipy
def f(x):
    import math
    return <expr>
total = 0.0
for i in xrange(100):
    x = value = float(i) / float(100)
    self.set_bnf_variable('<value>', value)
    total += abs(f(x) - (math.sin(x) + math.cos(x)))
fitness = total
self.set_bnf_variable('<fitness>', fitness)
        """


random.seed()
ges = GrammaticalEvolution()

ges.set_bnf(bnf)
ges.set_genotype_length(start_gene_length=20,
                        max_gene_length=50)
ges.set_population_size(100)
ges.set_wrap(True)

ges.set_completion_criteria("g", 200)
ges.set_completion_criteria("f", "center", 0.01)

ges.set_max_program_length(500)
ges.set_timeouts(1, 1)
ges.set_fitness_fail(100.0)

ges.set_max_fitness_rate(.5)
ges.set_fitness_selections(
    FitnessElites(ges.fitness_list, .05),
    FitnessProportionate(ges.fitness_list, 'linear'))
    #FitnessTournament(ges.fitness_list, tournament_size=2))

ges.set_children_per_crossover(2)
ges.set_mutation_type('s')
ges.set_mutation_rate(.075)

ges.set_replacement_selections(
        ReplacementTournament(ges.fitness_list, tournament_size=3))

ges.set_queue_size(0)
ges.set_garbage_collection(5)
ges.set_maintain_history(False)
ges.create_genotypes()
ges.run()
print ges.fitness_list.sorted()
print
print
gene = ges.population[ges.fitness_list.best_member()]
print gene.get_program()
