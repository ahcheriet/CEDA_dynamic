from deap import benchmarks, algorithms, base, creator, tools
from deap.benchmarks.tools import diversity, convergence
import random
import json
from hv import  *

'''
Created on 13-11-2012

@author: wysek
'''

N = 80
Nbar = 40
GEN = 5
U = 0
V = 1

referencePoint = [11,11]
hv = HyperVolume(referencePoint)      



def mainspea2(GEN,pb,pbs):
    def my_rand():
        return random.random() * (V - U) - (V + U) / 2

#-------------------------------------------------
    file = "pareto_front/"+pbs+"_front.json"
    optimal_front = json.load(open(file))
    optimal_front = [optimal_front[i] for i in range(0, len(optimal_front), 2)]
    first_point_opt = optimal_front[0]
    last_point_opt = optimal_front[-1]
# -----------------------------------------------
    creator.create("FitnessMax", base.Fitness, weights=(-1.0, -1.0,))
    creator.create("Individual", list, fitness=creator.FitnessMax) #@UndefinedVariable
    toolbox = base.Toolbox()
    toolbox.register("attr_float", my_rand)

    toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.attr_float, n=30) #@UndefinedVariable
    toolbox.register("evaluate", pb ) #benchmarks.zdt1
    toolbox.register("population", tools.initRepeat, list, toolbox.individual)
    toolbox.register("mate", tools.cxSimulatedBinaryBounded, eta=0.5, low=U, up=V)
    toolbox.register("mutate", tools.mutPolynomialBounded, eta=0.5, low=U, up=V, indpb=1)
    toolbox.register("select", tools.selSPEA2)
    
    #binary tournament selection
    toolbox.register("selectTournament", tools.selTournament, tournsize=2)


    # Step 1 Initialization
    pop = toolbox.population(n=N)
    archive = []
    curr_gen = 1
    obj_count = 0
    dvrst_list = []
    
    while True:
        # Step 2 Fitness assignement
        for ind in pop:
            ind.fitness.values = toolbox.evaluate(ind)

        for ind in archive:
            ind.fitness.values = toolbox.evaluate(ind)
        obj_count = tools.spea2_count + obj_count
        # Step 3 Environmental selection
        archive = toolbox.select(pop + archive, k=Nbar)
        aa = []
        for ee in archive:
            aa.append([ee.fitness.values[0],ee.fitness.values[1]])
        dvrst =  diversity(archive, first_point_opt, last_point_opt)
        hyprv = hv.compute(aa) 
        dvrst_list.append([obj_count , float(hyprv)])#  /float(dvrst) ])
        # Step 4 Termination
        if curr_gen >= GEN:
            final_set = archive
            break
        
        # Step 5 Mating Selection
        mating_pool = toolbox.selectTournament(archive, k=N)
        offspring_pool = map(toolbox.clone, mating_pool)

        # Step 6 Variation
        # crossover 100% and mutation 6%
        for child1, child2 in zip(offspring_pool[::2], offspring_pool[1::2]):
            toolbox.mate(child1, child2)

        for mutant in offspring_pool:
            if random.random() < 0.06:
                toolbox.mutate(mutant)

        pop = offspring_pool

        curr_gen += 1
        


   
 
    final_set.sort(key=lambda x: x.fitness.values)
    dvrst =  diversity(final_set, first_point_opt, last_point_opt)


    return final_set,dvrst,dvrst_list,obj_count
        
