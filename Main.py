from HvEDAlib import *
from SPEA2 import mainspea2
from NSGA2 import mainnsga2


xmin = 0
xmax = 1
dimension = 20

rho = 0



Contrainte01 = []
dtlz1 = []
Contrainte01.append('lambda x: True')
dtlz1.append('lambda x:benchmarks.dtlz1(x,3)[0]')
dtlz1.append('lambda x:benchmarks.dtlz1(x,3)[1]')
dtlz1.append('lambda x:benchmarks.dtlz1(x,3)[2]')

Contrainte01 = []
kursawe = []
Contrainte01.append('lambda x: True')
kursawe.append('lambda x:benchmarks.kursawe(x)[0]')
kursawe.append('lambda x:benchmarks.kursawe(x)[1]')

Contrainte01 = []
poloni = []
Contrainte01.append('lambda x: True')
poloni.append('lambda x:benchmarks.poloni(x)[0]')
poloni.append('lambda x:benchmarks.poloni(x)[1]')

Contrainte01 = []
fonseca = []
Contrainte01.append('lambda x: True')
fonseca.append('lambda x:benchmarks.fonseca(x)[0]')
fonseca.append('lambda x:benchmarks.fonseca(x)[1]')


Contrainte01 = []
schaffer_mo = []
Contrainte01.append('lambda x: True')
schaffer_mo.append('lambda x:benchmarks.schaffer_mo(x)[0]')
schaffer_mo.append('lambda x:benchmarks.schaffer_mo(x)[1]')


Contrainte01 = []
dtlz6 = []
Contrainte01.append('lambda x: True')
dtlz6.append('lambda x:benchmarks.dtlz1(x,2)[0]')
dtlz6.append('lambda x:benchmarks.dtlz1(x,2)[1]')


Contrainte01 = []
zdt1 = []
Contrainte01.append('lambda x: True')
zdt1.append('lambda x:benchmarks.zdt1(x)[0]')
zdt1.append('lambda x:benchmarks.zdt1(x)[1]')

Contrainte01 = []
zdt2 = []
Contrainte01.append('lambda x: True')
zdt2.append('lambda x:benchmarks.zdt2(x)[0]')
zdt2.append('lambda x:benchmarks.zdt2(x)[1]')


Contrainte01 = []
zdt3 = []
Contrainte01.append('lambda x: True')
zdt3.append('lambda x:benchmarks.zdt3(x)[0]') # I Add a comment
zdt3.append('lambda x:benchmarks.zdt3(x)[1]')

Contrainte01 = []
zdt4 = []
Contrainte01.append('lambda x: True')
zdt4.append('lambda x:benchmarks.zdt4(x)[0]')
zdt4.append('lambda x:benchmarks.zdt4(x)[1]')


Contrainte01 = []
zdt6 = []
Contrainte01.append('lambda x: True')
zdt6.append('lambda x:benchmarks.zdt6(x)[0]')
zdt6.append('lambda x:benchmarks.zdt6(x)[1]')


Contrainte01 = []
kursawe = []
Contrainte01.append('lambda x: True')
kursawe.append('lambda x:benchmarks.kursawe(x)[0]')
kursawe.append('lambda x:benchmarks.kursawe(x)[1]')



xmin = 0   #1
xmax = 1   #2
dimension = 20   #3 30 for zdt1, zdt2,zdt3 and 10 for zdt4 and zdt6
n_indv = 100
volume = []
spea2_count_hypervolume = []
nsga2_count_hypervolume = []

test = HvEDA(n_indv,FDA02,xmin,xmax,dimension,"spea2","half")  #4 pb alg method
test.selected_number = 50
#optimal_front = test.Calculed_Optimaml_front("kursawe") # don't forget to change the problem
optimal_front = test.GetFromFile_zdt("zdt1","real")  #5 file alg
Bool = 0
referencePoint = [11,11]
hv = HyperVolume(referencePoint)

def difference(list1,list2):
    liste = list1 + list2
    d = len(liste) - sum([ 1 for i,ele in enumerate(liste) if ele in liste[:i]+liste[(i+1):]  ])
    return d
      
if ( Bool == 0):
#    first_front_spea2,dvrst,hyp_spea2,spea2_cc = mainspea2(100 ,benchmarks.zdt1,"zdt1",xmin,xmax,dimension) #6 comapre
 #   test.plotting(hyp_spea2,  "Objective Functions Evaluation","Hypervolume", \
  #                "zdt1 problem, SPEA2 with the Copula-based EDA",'g-',"SPEA2",1, 0)


    nsga2_count_hypervolume,spea2_count_hypervolume = test.d_Run(40)
    PFront = test.ApplyFunctions(test.front_0)
#    filehandler = open('spea2_xi_kursawe_100', 'w')  # Save Copula to 
#    pickle.dump(test.C, filehandler) 
#    filehandler.close()

    print "diversity ",test.diversity( PFront , optimal_front[0], optimal_front[-1])
    print "convergance ",test.convergence( PFront , optimal_front )
    print "hypervolume ",spea2_count_hypervolume[-1:][0][1]
    
    test.plotting(spea2_count_hypervolume,  "Objective Functions Evaluation","Hypervolume", \
                  "zdt1 problem, SPEA2 with the Copula-based EDA",'r-',"CEDA",1, 0) #7 list methods
  
    test.plotting(PFront , "Objective function  1", "Objective function 2", "zdt1 Benchmark Pareto Front", 'b^',"Copula-based EDA",2, 0)
#    test.plotting(test.GetFromFile_zdt("zdt3","SPEA"),  "Objective function 1","Objective function 2","ZDT3 Benchmark Pareto Front",'g+',"SPEA2",2,0)
#    test.plotting(test.GetFromFile_zdt("zdt3","NSGA"),  "Objective function 1","Objective function 2","ZDT3 Benchmark Pareto Front",'rx',"NSGA2",2,0)

    
else :
    print "Starting"
    filehandler = open('./Objects/spea2_xi_zdt1_600', 'r')   # load Copula from file
    test.C  = pickle.load(filehandler) 
    filehandler.close()
    test.ListNonDominated = []
    front_nsga2 = []
    front_spea2 = []   
    nsga2_count_hypervolume = []
    spea2_count_hypervolume = []
    filter_count_hypervolume = []
    nsga2_count_diversity = []
    spea2_count_diversity = []
    filtrate_count_diversity = []
    nsga_front0 = []
    spea_front0 = []
    x = 0
    from SPEA2 import *
    dvv = []
    count_ss = 0
    last_spea = []
    count_nb = 0
    aa = []
    list_nb = []
    for lm in range(5):
        print lm
        last_spea = aa
        first_front_spea2,dvrst,dvrst_spea2,spea2_cc = mainspea2(150 ,benchmarks.zdt1,"zdt1",xmin,xmax,dimension)
        for ee in first_front_spea2:
            aa.append([ee.fitness.values[0],ee.fitness.values[1]])
        x = len(aa) - difference( last_spea , aa)
        count_nb = count_nb + spea2_cc
        list_nb.append([count_nb,x])
    print count_nb,x
    test.plotting(list_nb,  "Objective Functions Evaluation","New Solutions", \
                      "New best solutions Over Time",'b-',"SPEA2",1,0) #211)
        
    
    test.ListNonDominated = []
    print "Starting My Algorithm"
    count_f = 0
    gplot = Gnuplot.Gnuplot()
    All = []
    for mm in xrange(10):
        
        filehandler = open('./Objects/spea2_xi_zdt1_600', 'r')   # load Copula from file
        test.C  = pickle.load(filehandler) 
        filehandler.close()
        

        Zvar = test.Regenerate("xi",300)
        test.ListNonDominated = deepcopy(Zvar )# + test.ListNonDominated  #+ test.ListNonDominated
        test.fitnesss = test.ApplyFunctions(test.ListNonDominated)
#        test.filtrate()


#        count_f = count_f + test.filtrate_count
#        a,nsga_front0 = test.Load_ListND((test.ListNonDominated ), 100)
        spea_front0= test.selSPEA2((test.ListNonDominated ), 100)
        spea_front0 = test.front_0
        
        
        if mm % 5==0:
            print mm 

#        front_filtrate = test.ApplyFunctions(test.ListNonDominated)



        Last = front_spea2
#        front_nsga2 = test.ApplyFunctions(nsga_front0)# + front_nsga2
        front_spea2 = test.ApplyFunctions(spea_front0) + front_spea2
        
        x = len(front_spea2) - difference(Last , front_spea2) 
        hyprv = hv.compute(front_spea2)
#        nsga2_count_hypervolume.append([ test.nsga2_count, hyprv ])
        spea2_count_hypervolume.append([ test.spea2_count, hyprv ])
#        front_nsga2.sort()
        front_spea2.sort()

#        hyprv = hv.compute(front_filtrate)
#        filter_count_hypervolume.append([ count_f , hyprv ])
#        front_filtrate.sort()
 #       dvnfltr = test.diversity(front_filtrate , optimal_front[0], optimal_front[-1]) 
 
#        dvnsga2 = test.diversity(front_nsga2, optimal_front[0], optimal_front[-1]) 
        dvspea2 = test.diversity(front_spea2, optimal_front[0], optimal_front[-1])

#        filtrate_count_diversity.append([ count_f , dvnfltr ])
#        nsga2_count_diversity.append([ test.nsga2_count , dvnsga2 ])
        spea2_count_diversity.append([ test.spea2_count , x]) #float(hyprv) ])#  /float(dvspea2) ])

# -------------------------- After Loop ---------------------------




    test.plotting(spea2_count_diversity,  "Objective Functions Evaluation","New Solutions", \
                      "New best solutions Over Time",'r-',"Copula-Based EDA",1,0) #212)
    test.plotting(spea2_count_hypervolume,  "Objective Functions Evaluation","Hypervolume", \
                      "Hypervolume Over Time",'r-',"Copula-Based EDA",3,0)
    test.plotting(test.ApplyFunctions(spea_front0),  "Objective function 1","Objective function 2", \
                     "zdt1 problem",'g+',"Copula-Based EDA",2,0)

   # test.plotting(test.ApplyFunctions(first_front_spea2),  "Objective function 1","Objective function 2", \
    #                  "ZDT3 problem,",'r+',2)
   # test.plotting(test.GetFromFile_zdt("zdt1","real"),  "Objective function 1","Objective function 2", \
    #                  "zdt1 problem",'b+',"Optimal front",2,0)

    
   # test.plotting(nsga2_count_hypervolume,  "Number of evaluation functions","Diversity", \
     #                 "ZDT1 problem, Approximation of the Pareto front NSGA2",'r-',3)


    print "stop"

##-----------------------------------------------------------------------------
