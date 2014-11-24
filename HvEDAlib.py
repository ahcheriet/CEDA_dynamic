#
#   utilisation de Fast Sorting avec demi choix
#
from math import *
from copy import *
import random as rd
from decimal import *
from scipy.misc.common import *
import os
import numpy as np
from deap import benchmarks
import Gnuplot
import pickle

from itertools import chain
from operator import attrgetter, eq
from collections import Sequence, defaultdict

import json
#from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from openturns import *
from hv import HyperVolume
from copulalib.copulalib import Copula
#import multiprocessing as mp


rho = 0


Contrainte01 = []
FDA01 = []
Contrainte01.append('lambda x: True')
FDA01.append('lambda x:x[0]')
g = lambda x, nt, rhot: 1 + \
    sum(map(lambda a: pow(a - (sin(0.5 * pi * (1 / nt) * (floor(rho / rhot)))), 2), x[1:]))
h = lambda x, nt, rhot: (1 - sqrt(x[0] / g(x, nt, rhot)))
# rho generation counter
FDA01.append('lambda x:g(x,nt = 10.0,rhot =5.0)*h(x,nt = 10.0,rhot =5.0)')

Contrainte01 = []
FDA02 = []
Contrainte01.append('lambda x: True')
FDA02.append('lambda x:x[0]')
g = lambda x, nt, rhot: 1 + 9 * sum(map(lambda a: pow(a, 2), x[1:]))
h = lambda x, nt, rhot: (
    1 - pow(x[0] / g(x, nt, rhot), (1.25 + 0.75 * sin(0.5 * pi * (1 / nt) * (floor(rho / rhot))))))
# rho generation counter
FDA02.append('lambda x:g(x,nt = 10.0,rhot =5.0)*h(x,nt = 10.0,rhot =5.0)')



class HvEDA:
    def __init__(self,nbits,fonction,xmin,xmax,ndimension,method,estimate_methode):
        " testing HV With EDA "
        self.ndimension = ndimension
        self.xmin = xmin
        self.xmax = xmax
        self.ind = []
        self.C = range(self.ndimension)
        self.ListDominated = []
        self.fonctionObject = fonction
        self.X = []
        self.method = method
        self.estimate_methode = estimate_methode
        self.Tbit = nbits
        self.nsga2_count = 0
        self.spea2_count = 0
        self.time = 0
      #  self.selected_number = 50
    def IsDominate( self,Lfonction, x, y , MinOrMax):
        fonction = eval(Lfonction)
        if MinOrMax == 1:
            if min(fonction(x),fonction(y)) == fonction(x):
                return True
            else:
                return False
        else:
            if max(fonction(x),fonction(y)) == fonction(x):
                return True
            else:
                return False           
    def Choice(self):
        List_Sorted = []
        list_one = map(eval(self.fonctionObject[0]),self.ListNonDominated)#fonction1)
        list_two = sorted(list_one)
        for i in range(len(self.ListNonDominated)):
            index_one = list_one.index(list_two[i])
            element = deepcopy(self.ListNonDominated[index_one])
            List_Sorted.append(element)
        return List_Sorted
    def isContrainte(self,X):
        Zvar = deepcopy(X)
        Indexs = []
        for index in range(len(Zvar)): # test du contraintes
#             if ( ( Zvar[index][0]<0 or Zvar[index][0]>1)  or  any(map(lambda x:True if (x < -5 or x > 5) else False ,Zvar[index][1:])) == True ):
             if ( any(map(lambda x:True if (x < self.xmin or x > self.xmax) else False ,Zvar[index])) == True ):
                Indexs.append(Zvar[index])
        for index in range(len(Indexs)):
            Zvar.remove(Indexs[index])
        return Zvar
    
    def Dominating(self,x,y):
        if (self.fitnesss[x][0] < self.fitnesss[y][0] and self.fitnesss[x][1] < self.fitnesss[y][1]):
            return 1
        if (self.fitnesss[y][0] < self.fitnesss[x][0] and self.fitnesss[y][1] < self.fitnesss[x][1]):
            return 2
        
 #   @profile
    def selSPEA2(self, individuals , k):
        N = len(individuals)
        L = len(self.fonctionObject)
        K = sqrt(N)
        strength_fits = [0] * N
        fits = [0] * N
        dominating_inds = [list() for i in xrange(N)]
        for i, ind_i in enumerate(individuals):
            for j, ind_j in enumerate(individuals[i+1:], i+1):
                kobject = 0
                CountI = 0
                CountJ = 0
                self.spea2_count = self.spea2_count + 4
#                while kobject < len ( self.fonctionObject ):                    
#                    if  self.IsDominate( self.fonctionObject[kobject], ind_i , ind_j ,1): # minimisation problem
#                        CountI = CountI + 1
#                    if  self.IsDominate( self.fonctionObject[kobject], ind_j , ind_i ,1):
#                        CountJ = CountJ + 1
#                    kobject = kobject +1
#                if CountI == len( self.fonctionObject ):
#                    strength_fits[i] += 1
#                    dominating_inds[j].append(i)
#                if CountJ == len( self.fonctionObject ):
#                    strength_fits[j] += 1
#                    dominating_inds[i].append(j)
                if  self.Dominating(  i , j ) == 1: # minimisation problem
                    strength_fits[i] += 1
                    dominating_inds[j].append(i)
                if  self.Dominating(  j , i ) == 1:
                    strength_fits[j] += 1
                    dominating_inds[i].append(j)
        
        for i in xrange(N):
            for j in dominating_inds[i]:
                fits[i] += strength_fits[j]
        
        #print "Choose all non-dominated individuals"
        chosen_indices = [i for i in xrange(N) if fits[i] < 1]
        self.front_0 = deepcopy([individuals[i] for i in chosen_indices])
        if len(chosen_indices) < k:     # The archive is too small
            for i in xrange(N):
                distances = [0.0] * N
                for j in xrange(i + 1, N):
                    dist = 0.0
                    for l in xrange(L):
#                        val = eval(self.fonctionObject[l])(individuals[i]) - eval(self.fonctionObject[l])(individuals[j]) # euclidienne distence**2 = a - b 
                        val = self.fitnesss[i][l] - self.fitnesss[j][l]
                        dist += val * val
                    distances[j] = dist
                kth_dist = self.randomizedSelect(distances, 0, N - 1, K)
                density = 1.0 / (kth_dist + 2.0)
                fits[i] += density
                
            next_indices = [(fits[i], i) for i in xrange(N)
                            if not i in chosen_indices]
            next_indices.sort()
            #print next_indices
            chosen_indices += [i for _, i in next_indices[:k - len(chosen_indices)]]
                    
        elif len(chosen_indices) > k:   # The archive is too large
            N = len(chosen_indices)
            distances = [[0.0] * N for i in xrange(N)]
            sorted_indices = [[0] * N for i in xrange(N)]
            for i in xrange(N):
                for j in xrange(i + 1, N):
                    dist = 0.0
                    for l in xrange(L):
#                        val = eval(self.fonctionObject[l])(individuals[chosen_indices[i]]) - eval(self.fonctionObject[l])(individuals[chosen_indices[j]])
                        val = self.fitnesss[chosen_indices[i]][l] - self.fitnesss[chosen_indices[j]][l]
                        dist += val * val
                    distances[i][j] = dist
                    distances[j][i] = dist
                distances[i][i] = -1
            
           # print "Insert sort is faster than quick sort for short arrays"
            for i in xrange(N):
                for j in xrange(1, N):
                    l = j
                    while l > 0 and distances[i][j] < distances[i][sorted_indices[i][l - 1]]:
                        sorted_indices[i][l] = sorted_indices[i][l - 1]
                        l -= 1
                    sorted_indices[i][l] = j
            
            size = N
            to_remove = []
            while size > k:
                #print "Search for minimal distance"
                min_pos = 0
                for i in xrange(1, N):
                    for j in xrange(1, size):
                        dist_i_sorted_j = distances[i][sorted_indices[i][j]]
                        dist_min_sorted_j = distances[min_pos][sorted_indices[min_pos][j]]
                        
                        if dist_i_sorted_j < dist_min_sorted_j:
                            min_pos = i
                            break
                        elif dist_i_sorted_j > dist_min_sorted_j:
                            break
                
                #print "Remove minimal distance from sorted_indices"
                for i in xrange(N):
                    distances[i][min_pos] = float("inf")
                    distances[min_pos][i] = float("inf")
                    
                    for j in xrange(1, size - 1):
                        if sorted_indices[i][j] == min_pos:
                            sorted_indices[i][j] = sorted_indices[i][j + 1]
                            sorted_indices[i][j + 1] = min_pos
                
                #print "Remove corresponding individual from chosen_indices"
                to_remove.append(min_pos)
                size -= 1
            
            for index in reversed(sorted(to_remove)):
                del chosen_indices[index]
        
        return [individuals[i] for i in chosen_indices]
        
    def randomizedSelect(self,array, begin, end, i):
        """Allows to select the ith smallest element from array without sorting it.
        Runtime is expected to be O(n).
        """
        if begin == end:
            return array[begin]
        q = self.randomizedPartition(array, begin, end)
        k = q - begin + 1
        if i < k:
            return self.randomizedSelect(array, begin, q, i)
        else:
            return self.randomizedSelect(array, q + 1, end, i - k)
    
    def randomizedPartition(self,array, begin, end):
        i = rd.randint(begin, end)
        array[begin], array[i] = array[i], array[begin]
        return self.partition(array, begin, end)
        
    def partition(self,array, begin, end):
        x = array[begin]
        i = begin - 1
        j = end + 1
        while True:
            j -= 1
            while array[j] > x:
                j -= 1
            i += 1
            while array[i] < x:
                i += 1
            if i < j:
                array[i], array[j] = array[j], array[i]
            else:
                return j

    def diversity(self,first_front, first, last):
        """Given a Pareto front `first_front` and the two extreme points of the
        optimal Pareto front, this function returns a metric of the diversity
        of the front as explained in the original NSGA-II article by K. Deb.
        The smaller the value is, the better the front is.
        """
        
        df = hypot(first_front[0][0] - first[0],
                   first_front[0][1] - first[1])
        dl = hypot(first_front[-1][0] - last[0],
                   first_front[-1][1] - last[1])
        dt = [hypot( first[0] - second[0],
                     first[1] - second[1])
              for first, second in zip(first_front[:-1], first_front[1:])]

#        dt = [hypot(eval(self.fonctionObject[0])(first) - eval(self.fonctionObject[0])(second),
#                    eval(self.fonctionObject[1])(first) - eval(self.fonctionObject[1])(second))
#              for first, second in zip(first_front[:-1], first_front[1:])]
    
        if len(first_front) == 1:
            return df + dl
    
        dm = sum(dt)/len(dt)
        di = sum(abs(d_i - dm) for d_i in dt)
        delta = (df + dl + di)/(df + dl + len(dt) * dm )
        return delta
    
    def convergence(self,first_front, optimal_front):
        """Given a Pareto front `first_front` and the optimal Pareto front,
        this function returns a metric of convergence
        of the front as explained in the original NSGA-II article by K. Deb.
        The smaller the value is, the closer the front is to the optimal one.
        """
        distances = []
       
        for ind in first_front:
            distances.append(float("inf"))
            for opt_ind in optimal_front:
                dist = 0.
                for i in xrange(len(opt_ind)):
                    dist += ( ind[i] - opt_ind[i])**2
                if dist < distances[-1]:
                    distances[-1] = dist
            distances[-1] = sqrt(distances[-1])
           
        return sum(distances) / len(distances)
###################    
    nsga2_count = 0
    def selNSGA2(self,individuals, k):
        """Apply NSGA-II selection operator on the *individuals*. Usually, the
        size of *individuals* will be larger than *k* because any individual
        present in *individuals* will appear in the returned list at most once.
        Having the size of *individuals* equals to *k* will have no effect other
        than sorting the population according to a non-domination scheme. The list
        returned contains references to the input *individuals*. For more details
        on the NSGA-II operator see [Deb2002]_.
        
        :param individuals: A list of individuals to select from.
        :param k: The number of individuals to select.
        :returns: A list of selected individuals.
        
        .. [Deb2002] Deb, Pratab, Agarwal, and Meyarivan, "A fast elitist
           non-dominated sorting genetic algorithm for multi-objective
           optimization: NSGA-II", 2002.
        """
        pareto_fronts = self.sortFastND(individuals, k)
        for front in pareto_fronts:
            self.assignCrowdingDist(front)
        
        chosen = list(chain(*pareto_fronts[:-1]))
        k = k - len(chosen)
        if k > 0:
            sorted_front = sorted(pareto_fronts[-1], key=attrgetter("self.crowding_dist"), reverse=True)
            chosen.extend(sorted_front[:k])
            
        return chosen
        
    def sortFastND(self,individuals, k, first_front_only=False):
        """Sort the first *k* *individuals* according the the fast non-dominated
        sorting algorithm.
        
        :param individuals: A list of individuals to select from.
        :param k: The number of individuals to select.
        :param first_front_only: If :obj:`True` sort only the first front and
                                 exit.
        :returns: A list of Pareto fronts (lists), with the first list being the
                  true Pareto front.
        """
        nsga2_count = 0
        if k == 0:
            return []
    
        unique_fits = defaultdict(list)
        for i,ind in enumerate(individuals):
            unique_fits[self.fitnesss[i]].append(i)
        fits = unique_fits.keys()
    
        N = len(fits)
        pareto_fronts = []
        
        pareto_fronts.append([])
        pareto_sorted = 0
        dominating_fits = [0] * N
        dominated_fits = [list() for i in xrange(N)]
        
        # Rank first Pareto front
        for i, fit_i in enumerate(fits):
            for j, fit_j in enumerate(fits[i+1:], i+1):
             #   nsga2_count = nsga2_count + 4 
                if self.Dominating(fit_i, fit_j):
                    dominating_fits[i] += 1
                    dominated_fits[j].append(i)
                elif self.Dominating(fit_j, fit_i):
                    dominating_fits[j] += 1
                    dominated_fits[i].append(j)
            if dominating_fits[i] == 0:
                pareto_fronts[-1].append(i)
                pareto_sorted += 1
        
        # Rank the next front until all individuals are sorted or the given
        # number of individual are sorted
        if not first_front_only:
            N = min(N, k)
            while pareto_sorted < N:
                pareto_fronts.append([])
                for indice_p in pareto_fronts[-2]:
                    for indice_d in dominated_fits[indice_p]:
                        dominating_fits[indice_d] -= 1
                        if dominating_fits[indice_d] == 0:
                            pareto_fronts[-1].append(indice_d)
                            pareto_sorted += 1
        
        total = min(len(individuals), k)
        fronts = list()
        for front in pareto_fronts:
            fronts.append([])       
            for index in front:
                fronts[-1].extend(unique_fits[fits[index]])
            total -= len(fronts[-1])
            if total <= 0:
                break
    
        return fronts
    
    def assignCrowdingDist(self,individuals):
        """Assign a crowding distance to each individual of the list. The 
        crowding distance is set to the :attr:`crowding_dist` attribute of
        each individual.
        """
        if len(individuals) == 0:
            return
        
        distances = [0.0] * len(individuals)
        crowd = [(self.fitnesss[i], i) for i, ind in enumerate(individuals)]
        
        nobj = len(self.fitnesss[0])
        
        for i in xrange(nobj):
            crowd.sort(key=lambda element: element[0][i])
            distances[crowd[0][1]] = float("inf")
            distances[crowd[-1][1]] = float("inf")
            if crowd[-1][0][i] == crowd[0][0][i]:
                continue
            norm = nobj * float(crowd[-1][0][i] - crowd[0][0][i])
            for prev, cur, next in zip(crowd[:-2], crowd[1:-1], crowd[2:]):
                distances[cur[1]] += (next[0][i] - prev[0][i]) / norm
        self.crowding_dist = deepcopy(self.fitnesss)
        for i, dist in enumerate(distances):
            self.crowding_dist[i] = dist
    
    
    
    
    
    

####################
    def Fast_Sorting(self,X):
        tmpX = deepcopy(X)
        n = [-1]* len(tmpX)
        S = [0]* len(tmpX)
        rank = [-1]* len(tmpX)
        F = [[]]
        for x_i in range(len(tmpX)):            
            S[x_i] = []
            n[x_i] = 0
            for x_j in  range(len(tmpX)):               
                if (x_j != x_i):
                    self.nsga2_count = self.nsga2_count + 2
#                    k = 0
#                    CountI = 0
#                    CountJ = 0
#                    self.nsga2_count = self.nsga2_count + 1      
#                    while k < len ( self.fonctionObject ):                    
#                        if  self.IsDominate( self.fonctionObject[k], tmpX[x_i], tmpX[x_j],1): # minimisation problem
#                            CountI = CountI + 1
#                        if  self.IsDominate( self.fonctionObject[k], tmpX[x_j], tmpX[x_i],1):
#                            CountJ = CountJ + 1
#                        k = k +1
#                    if CountI == len( self.fonctionObject ):
#                        S[x_i].append(x_j)
#                    if CountJ == len( self.fonctionObject ):
#                        n[x_i] = n[x_i] + 1
                    if self.Dominating(x_i,x_j) == 1:
                        S[x_i].append(x_j)
                    if self.Dominating(x_j,x_i) == 1:
                        n[x_i] = n[x_i] + 1
            if n[x_i] == 0:
                rank[x_i] = 0
                F[0]=F[0]+[x_i]
        i = 0
        ri = 1
        while F[i] != []:
            Q = []
            for p in F[i]:
                for q in S[p]:
                    n[q] = n[q]-1
                    if n[q] == 0:
                        rank[q] = ri + 1
                        Q.append(q)
            i = i + 1
            ri = ri + 1
            F.append(deepcopy(Q))
        return F
    
    def Load_ListND(self,X,n):
        tmpND = []
        f = self.Fast_Sorting(X)
        for i in range(len(f)):
            for j in range(len(f[i])):
                tmpND.append(X[f[i][j]])
                if len(tmpND) == n:
                    return tmpND ,tmpND[:len(f[0])]
        return tmpND , tmpND[:len(f[0])]
    def Evaluate(self):        
        i = 0
        self.X = []
        self.ListNonDominated = self.Make_First_Pop()# self.isContrainte(self.Make_First_Pop())
    #    self.ListNonDominated = deepcopy(self.X)
        self.fitnesss =  self.ApplyFunctions(self.ListNonDominated)
      #  self.filtrate()
        self.front_0 = self.ListNonDominated
    def filtrate(self):
        i = 0
        self.filtrate_count = 0
        while i < len(self.ListNonDominated)-1:
                j = i+1
                while j < len(self.ListNonDominated):
                        CountI = 0
                        CountJ = 0
                        k = 0
                        self.filtrate_count = self.filtrate_count + 4
#                        while k < len ( self.fonctionObject ):                    
#                            if  self.IsDominate( self.fonctionObject[k], self.ListNonDominated[i], self.ListNonDominated[j],1):
#                                CountI = CountI + 1
#                            if  self.IsDominate( self.fonctionObject[k], self.ListNonDominated[j], self.ListNonDominated[i],1):
#                                CountJ = CountJ + 1                               
#                            k = k +1
                        if  self.Dominating(  i , j ) == 1: # minimisation problem
                                x = self.ListNonDominated[j]
                                j = j - 1
                                self.ListNonDominated.remove(x)
                        if  self.Dominating(  j , i ) == 1:
                                x = self.ListNonDominated[i]
                                j = j - 1
                                i = i - 1
                                self.ListNonDominated.remove(x)
                                break
#
#                        if CountI == len( self.fonctionObject ):
#                                x = self.ListNonDominated[j]
#                                j = j - 1
#                                self.ListNonDominated.remove(x)
#                        if CountJ == len( self.fonctionObject ):
#                                x = self.ListNonDominated[i]
#                                j = j - 1
#                                i = i - 1
#                                self.ListNonDominated.remove(x)
#                                break
                        j = j + 1
                i = i +1
    def Make_First_Pop(self):
        ligne = []
        i = 1
        ligne.append((np.random.uniform(self.xmin,self.xmax,self.Tbit)).tolist())
        while i < self.ndimension:
            ligne.append((np.random.uniform(self.xmin,self.xmax,self.Tbit)).tolist()) # 100 individu
            i = i + 1       # for zdt1,zdt2,zdt3
        return zip(*ligne)

#        while i < self.ndimension:
#            ligne.append((np.random.uniform(-5,5,self.Tbit)).tolist()) # zdt4
#            i = i + 1
#        return zip(*ligne)


    def Make_X_UsingCopula_half(self):
        Xvar = range(self.ndimension)        # pour que la liste puisse contenir 30 elements
#        Varrtmp = self.Choice()
        varr = list(zip(*self.ListNonDominated))
        indexs = rd.sample(range(self.ndimension),self.ndimension)
        xx = rd.sample( self.ListNonDominated , len(self.ListNonDominated)/2) # demi de meilleur sol
        yy = [item for item in self.ListNonDominated if item not in xx]
        if len(xx)>len(yy):
            xx.pop()            # test egality
        if len(yy)>len(xx):
            yy.pop()
        k = 0        
        for j in range(self.ndimension/2):
            xm = deepcopy(varr[k])
            ym = deepcopy(varr[k+1])
            x = np.array(xm)
            y = np.array(ym)
            foo = Copula(x, y, family='frank') # gumbel, clayton, frank
            self.C[j] = deepcopy(foo)            
            XX, YY = foo.generate_xy(self.Tbit)
            X1 = XX.tolist()
            Y1 = YY.tolist()
            Xvar[k] = deepcopy(X1)
            Xvar[k+1] = deepcopy(Y1)
            k = k + 2
        Zvar = zip(*Xvar)
        Indexs = []
        for index in range(len(Zvar)): # test du contraintes
#            if ( ( Zvar[index][0]<0 or Zvar[index][0]>1)  or  any(map(lambda x:True if (x < self.xmin or x > self.xmax) else False ,Zvar[index][1:])) == True ):
            if ( any(map(lambda x:True if (x < self.xmin or x > self.xmax ) else False ,Zvar[index])) == True ):
                Indexs.append(Zvar[index])
        for index in range(len(Indexs)):
            Zvar.remove(Indexs[index])
        return Zvar
        
  #  @profile    
    def Make_X_UsingCopula_xi(self):
        Xvar = []
        varr = list(zip(*self.ListNonDominated))
        xx = rd.sample( self.ListNonDominated , len(self.ListNonDominated)/2) # demi de meilleur sol and random
        yy = [item for item in self.ListNonDominated if item not in xx]
        if len(xx)>len(yy):
            xx.pop()            # test egality
        if len(yy)>len(xx):
            yy.pop()        
        for j in range(self.ndimension):
            xm = []
            ym = []
            for i in range(len(yy)):   # youo are here 
                xm.append(xx[i][j])
                ym.append(yy[i][j])
            x = np.array(xm)
            y = np.array(ym)
            foo = Copula(x, y, family='frank') # gumbel, clayton, frank
            self.C[j] = deepcopy(foo)            
            XX, YY = foo.generate_xy(self.Tbit)
            X1 = XX.tolist()
            Y1 = YY.tolist()
            Xvar.append(X1+Y1)
        Zvar = zip(*Xvar)
        Indexs = []
        for index in range(len(Zvar)): # test du contraintes
#            if ( ( Zvar[index][0]<0 or Zvar[index][0]>1)  or  any(map(lambda x:True if (x < -5 or x > 5) else False ,Zvar[index][1:])) == True ):
            if ( any(map(lambda x:True if (x < self.xmin or x > self.xmax) else False ,Zvar[index])) == True ):
                Indexs.append(Zvar[index])
        for index in range(len(Indexs)):
            Zvar.remove(Indexs[index])
        return Zvar
    
    def Make_X_UsingCopula_openturns(self):  # look at the filter
        Zvar  = []
        Xvar = []
#        varr = list(zip(*self.ListNonDominated))
        print len(self.ListNonDominated)
        xx = rd.sample( self.ListNonDominated , len(self.ListNonDominated)/2) # demi de meilleur sol and random
        yy = [item for item in self.ListNonDominated if item not in xx]
        if len(xx)>len(yy):
            xx.pop()            # test egality
        if len(yy)>len(xx):
            yy.pop()  
        copulaColl = CopulaCollection(30) 
        estimatedDistribution = []
        for j in range(self.ndimension):
            xm = []
            ym = []
            x = NumericalSample(0,2)
            for i in range(len(yy)):    
                xm.append(deepcopy(xx[i][j]))
                ym.append(deepcopy(yy[i][j]))
                x.add([xx[i][j],yy[i][j]])  
            copulaColl[j] = ClaytonCopulaFactory().build(x)
        self.C = ComposedCopula(copulaColl)
        list_new = self.C.getNumericalSample(30) # sampling
        print list(list_new[1])[0::2]
 
        for i in range(self.ndimension):
            X = list(list_new[i])[0::2]
            print X
            Y = list(list_new[i])[1::2]
            Zvar.append(X+Y)
        Zvar = zip(*Zvar)
        Indexs = []

        for index in range(len(Zvar)): # test du contraintes
            if ( any(map(lambda x:True if (x < self.xmin or x >self.xmax) else False ,Zvar[index])) == True ):
                Indexs.append(Zvar[index])
        for index in range(len(Indexs)):
            Zvar.remove(Indexs[index])
        return Zvar

    def Evaluate_After(self):  
        " evaluter avec hypervolume et Copula "
        referencePoint = [6.0,6.0]
        ListNonDominatedtmp = []      
        self.indexbest = 0
        if (self.estimate_methode == "xi"):
            self.ListNonDominated = self.Make_X_UsingCopula_xi() + self.ListNonDominated
        if (self.estimate_methode =="half"):
            self.ListNonDominated = self.Make_X_UsingCopula_half() + self.ListNonDominated
        self.fitnesss = self.ApplyFunctions(self.ListNonDominated)
        if (self.method == "spea2"):
            self.ListNonDominated  = self.selSPEA2(self.ListNonDominated, self.selected_number)
        if (self.method == "nsga2"):
#            self.ListNonDominated = self.selNSGA2(self.ListNonDominated, self.selected_number)
            self.ListNonDominated , self.front_0 = self.Load_ListND(self.ListNonDominated, self.selected_number)
#            self.ListNonDominated = self.front_0 
            #self.filtrate()
            #self.front_0 = self.ListNonDominated
            
    def ApplyFunctions(self,ListX):
        tmpResult = []
        for i in range(len( ListX )):               
            tmpResult.append([eval(fonction)(ListX[i]) for fonction in self.fonctionObject])
        return tmpResult
    
    def Tofile( self , Str_File , X  , AorW):
        if AorW == 1:
            fileData = open(Str_File+".txt","a")
        else:
            fileData = open(Str_File+".txt","w")            
        i = 0
        while i < len( X ):
            fonction1 = eval(self.fonctionObject[0])
            fonction2 = eval(self.fonctionObject[1])
            fileData.write(str(fonction1(X[i]))+'\t'+str(fonction2(X[i]))+'\n')
            i = i +1
        fileData.close()
    
    def Calculed_Optimaml_front(self,problem):
        file = "pareto_front/"+problem+"_front.json"
        optimal_front = json.load(open(file))
        # Use 500 of the 1000 points in the json file
        optimal_front = sorted(optimal_front[i] for i in range(0, len(optimal_front), 2))
        return optimal_front
    
    def helper3(self,lst):
        lst1, lst2 ,lst3= [], [],[]
        for el in lst:
            lst1.append(el[0])
            lst2.append(el[1])
            lst3.append(el[2])
        return lst1, lst2 ,lst3
    def helper(self,lst):
        lst1, lst2= [], []
        for el in lst:
            lst1.append(el[0])
            lst2.append(el[1])
        return lst1, lst2 
    def plotting3d(self,list,xlabel,ylabel,title,format,figure):
        l1,l2,l3 = self.helper(list)
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.scatter(l1, l2, l3, c='b', marker='+')
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title(title)
        plt.savefig(title+".pdf")
        

    def plotting(self,list,xlabel,ylabel,title,formats,labelf,figure,sub):
        l1,l2 = self.helper(list)
        plt.figure(figure)
        if sub!=0 :
            plt.subplot(sub)
            plt.plot(l1,l2,formats,label=labelf) 
            plt.ticklabel_format(axis='x',style='sci',scilimits=(1,2))
            plt.ylabel(ylabel)
            plt.xlabel(xlabel) 
            plt.legend()
            plt.annotate('best = '+ format(min(l2),'.3f') ,xy = ( l1[l2.index( min(l2) )] , min(l2) ),
                           xytext=(-50, +20), textcoords='offset points', fontsize=8,
                           arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2"))

  #          plt.title(title)
        else:
            plt.plot(l1,l2,formats,label=labelf)
            plt.ticklabel_format(axis='x',style='sci',scilimits=(1,2))
            plt.xlabel(xlabel)
            plt.ylabel(ylabel)
            plt.legend()
            plt.title(title)
   #     plt.xscale('log')
        plt.savefig(title+".pdf")
        
    def GetFromFile_zdt(self,file,algorithm):
        data = np.genfromtxt("./pareto_front/alg_results/"+file+"/"+algorithm+".30", delimiter=' ')#, names = ['x','y'])
        return data
    def Regenerate(self,mode,nbit):      
        if mode == "half":
            Xvar = range(30)
            k = 0        
            for j in range(self.ndimension/2):           
                XX, YY = self.C[j].generate_xy(nbit)
                X1 = XX.tolist()
                Y1 = YY.tolist()
                Xvar[k] = deepcopy(X1)             #half
                Xvar[k+1] = deepcopy(Y1)
                k = k + 2
            Zvar = zip(*Xvar)
            Indexs = []
            for index in range(len(Zvar)): # test du contraintes
                if ( any(map(lambda x:True if (x < self.xmin or x > self.xmax ) else False ,Zvar[index])) == True ):
                    Indexs.append(Zvar[index])
            for index in range(len(Indexs)):
                Zvar.remove(Indexs[index])
            return Zvar
        else:
            Xvar = []
            for i in range(self.ndimension):
                XX, YY = self.C[i].generate_xy(nbit)
                X1 = XX.tolist()
                Y1 = YY.tolist()
                Xvar.append(X1+Y1)
            Zvar = zip(*Xvar)
            Indexs = []
            for index in range(len(Zvar)):
                if ( any(map(lambda x:True if (x < self.xmin or x >self.xmax) else False ,Zvar[index])) == True ):
                    Indexs.append(Zvar[index])
            for index in range(len(Indexs)):
                Zvar.remove(Indexs[index])
            return Zvar
 #   @profile            
    def Run(self,ngen):    
        volume = []
        spea2_count_hypervolume = []
        nsga2_count_hypervolume = []
        referencePoint = [11,11]
        hv = HyperVolume(referencePoint)      
        self.Evaluate()
        gplot = Gnuplot.Gnuplot()
        print "Starting"
        for i in xrange(ngen):
            if i % 1 ==0 :
                print "Iteration  ",i
            self.Evaluate_After()
            d=Gnuplot.Data(self.ApplyFunctions(self.ListNonDominated))
            gplot.plot(d)
            front = self.ApplyFunctions(self.front_0)
            volume.append (hv.compute(front))
            nsga2_count_hypervolume.append([self.nsga2_count,hv.compute(front)])
            spea2_count_hypervolume.append([self.spea2_count,hv.compute(front)])
        return nsga2_count_hypervolume,spea2_count_hypervolume
##-----------------------------------------------------------------------------
# Run Dynamic

    def d_Run(self, ngen):
        volume = []
        spea2_count_hypervolume = []
        nsga2_count_hypervolume = []
        referencePoint = [11, 11]
        hv = HyperVolume(referencePoint)
        self.Evaluate()
        gplot = Gnuplot.Gnuplot()
        print "Starting"
        for i in xrange(ngen):
            print "Iteration  ", i
            self.selected_set = 50
            if i % 10 == 0:
                print "Change maded"
                self.selected_set = 100
    #            gg = self.Make_First_Pop()
   #             self.ListNonDominated = self.ListNonDominated + self.Make_First_Pop()
            global rho
            rho = rho + 1
            self.Evaluate_After()
            d = Gnuplot.Data(self.ApplyFunctions(self.ListNonDominated))
            gplot.plot(d)
            front = self.ApplyFunctions(self.front_0)
            volume.append(hv.compute(front))
            nsga2_count_hypervolume.append(
                [self.nsga2_count, hv.compute(front)])
            spea2_count_hypervolume.append(
                [self.spea2_count, hv.compute(front)])
        return nsga2_count_hypervolume, spea2_count_hypervolume

