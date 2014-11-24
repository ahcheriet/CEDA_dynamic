from deap import benchmarks, algorithms, base, creator, tools
from deap.benchmarks.tools import diversity, convergence
import random
'''
Created on 06-06-2012

@author: pita
'''

N=100
GEN=200
U=0
V=1


def mainnsga2(GEN,pb):
    def my_rand():
        return random.random()*(V-U) - (V+U)/2
    
    file = "pareto_front/"+pbs+"_front.json"
    optimal_front = json.load(open(file))
    optimal_front = [optimal_front[i] for i in range(0, len(optimal_front), 2)]
    first_point_opt = optimal_front[0]
    last_point_opt = optimal_front[-1]

    creator.create("FitnessMax", base.Fitness, weights=(-1.0,-1.0,))
    creator.create("Individual", list, fitness=creator.FitnessMax) #@UndefinedVariable
    toolbox = base.Toolbox()
    toolbox.register("attr_float", my_rand)
    toolbox.register("individual", tools.initRepeat, creator.Individual,toolbox.attr_float, n=30) #@UndefinedVariable
    toolbox.register("evaluate", pb) #
    toolbox.register("population", tools.initRepeat, list, toolbox.individual)
    toolbox.register("mate", tools.cxSimulatedBinaryBounded, eta=0.5, low=U, up=V)
    toolbox.register("mutate", tools.mutPolynomialBounded, eta=0.5, low=U, up=V, indpb=1)
    toolbox.register("select", tools.selNSGA2)
    toolbox.register("selectTournament", tools.selTournamentDCD)
    obj_count = 0
    dvrst_list = []

    
    # init population
    pop = toolbox.population(n=N)
    for ind in pop:
        ind.fitness.values = toolbox.evaluate(ind)
    
    # sort using non domination sort (k is the same as n of population - only sort is applied)
    pop = toolbox.select(pop, k=N)
    
    for _ in xrange(GEN):
        #select parent pool with tournament dominated selection
        parent_pool = toolbox.selectTournament(pop, k=N)
        offspring_pool = map(toolbox.clone, parent_pool)
    
        # Apply crossover and mutation on the offspring
        for child1, child2 in zip(offspring_pool[::2], offspring_pool[1::2]):
            if random.random() < 0.9:
                toolbox.mate(child1, child2)
        for mutant in offspring_pool:
            if random.random() < 0.1:
                toolbox.mutate(mutant)
    
        # evaluate offsprings
        for ind in offspring_pool:
            ind.fitness.values = toolbox.evaluate(ind)
        
        # extend base population with offsprings, pop is now 2N size
        pop.extend(offspring_pool)
        
        # sort and select new population
        pop = toolbox.select(pop, k=N)
        
        obj_count = tools.spea2_count + obj_count
        dvrst =  diversity(pop, first_point_opt, last_point_opt)
        dvrst_list.append([obj_count , dvrst])

    
    first_front = tools.sortFastND(pop, k=N)[0]
    return first_front,