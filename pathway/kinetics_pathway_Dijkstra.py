# There is still a bug of taking the same reverse reaction back!
# e.g. CO => H2O 

""""
Similar to fc_pathway_Dijkstra.py but taking the mixing ratios from a certain level in Vulcan's output instead.
"""
import sys, os
os.chdir('../')
sys.path.insert(0, '.') # including the upper level of directory for the path of modules
#sys.path.insert(0, '/Users/tsai/Dropbox/Bern/python/Coding_test/vulcan_photo_benchmarking')

import numpy as np 
from collections import defaultdict
from scipy import interpolate
import itertools
import operator
import pydot
from PIL import Image 
import pickle

# import VULCAN modules
import store, op
import vulcan_cfg
from phy_const import kb, Navo
import chem_funs

np.set_printoptions(threshold=np.inf)

# read Vulcan pickle files
vul_data = 'output/v2718-HD189.vul'

with open(vul_data, 'rb') as handle:
  data = pickle.load(handle)
  
# the pressure to analize
p_range = [1e2,1e1,5e1] # in cgs
#p_range = [5e-2,1e-1,1e0]

#p_range = [1.e6,1e7,1e8,5e8]

remove_fw_re = [593,15,17   ,327]
remove_fw_re = []

# the pathway (conversion of the 2 species) to analize
conv_sp = ('C2H2','CH4')

species = data['variable']['species']

#species = chem_funs.spec_list
re_dict = chem_funs.re_dict
re_wM_dict = chem_funs.re_wM_dict

# options for adding more paths
use_add_more = False

    
for p_ana in p_range:
    
    # Find the index of pco closest to p_ana
    p_indx = min( range(len(data['atm']['pco'])), key=lambda i: abs(data['atm']['pco'][i]-p_ana))

    # # pressure and T from the vulcan data
    p_ana = data['atm']['pco'][p_indx]
    T_ana = data['atm']['Tco'][p_indx]
    print ("Using the kinetics results at P = " + "{:>10.2e}".format(p_ana/1e6) + " bar, T = " + str(T_ana) )
    
    
    # the third body
    mm = p_ana*1.e6/kb/T_ana

    # species to analize
    ana_sp = list(species)
    ana_sp.remove('He')
    ana_sp.remove('H')
    ana_sp.remove('H2')
#
    ana_sp.remove('OH')
    ana_sp.remove('H2O')
    
    #ana_sp.remove('CH4')
    #ana_sp.remove('CH3')
    # ana_sp.remove('C2H2')
    ana_sp.remove('O')
    # ana_sp.remove('CO')
    # ana_sp.remove('CH2')
    # ana_sp.remove('CH')
    # ana_sp.remove('C')
    # ana_sp.remove('C2H4')
    # ana_sp.remove('C2H6')
    
    # all C species
    #ana_sp = [ 'CH4', 'CH3', 'CH2', 'CH', 'C', 'CO', 'CO2','C2',  'C2H3', 'C2H5', 'C2H6', 'C2H4', 'H2CO', 'HCO', 'CH2OH', 'CH3OH', 'CH3O',  'CH3CO', 'C2H', 'C2H2']

    # all N species
    #ana_sp = [ 'CH4', 'CH3',  'N', 'NH', 'CN', 'HCN', 'NO', 'NH2', 'N2', 'NH3', 'N2H2', 'N2H', 'N2H3', 'N2H4', 'CH3NH', 'CH2NH', 'CH2NH2', 'CH3NH2', 'HNO', 'H2CN', 'HNCO', 'CH3NO', 'NO2', 'HNO2']

    #remove_fw_re = [523,197,383,207,193,127,153, 309,  251] # 251: CH4 + CN
    #remove_fw_re = [251,  207,213, 523, 193,127,153,197,199,157, 333, 321.   ,369,329] #333: NH3 + CN
    #remove_fw_re = [251, 207, 193,  369] #,321
    # R473 is new from Moses, I added it after the relaxation paper.
    #remove_fw_re = [329,333,369, 345,309,373] #, 247, 329, 333, 309,  371
    
    #remove_fw_re = [17,173]
   
    
    remove_re = list(remove_fw_re)
    [remove_re.append(re+1) for re in remove_fw_re]

    
    # all the combination from ana_sp
    ana_comb = list(itertools.combinations(ana_sp, 2))
    remove_pair = []
    for pair in remove_pair:  
        if pair in ana_comb: ana_comb.remove(pair)
        if pair[::-1] in ana_comb: ana_comb.remove(pair[::-1])

    # store the combination of prod and reac
    re_comb = {}

    # store all the connections
    edges = defaultdict(list)

    # store the net reaction rates of a pair of species, and the fastest reaction channel
    netRate_pair, maxRate = {}, {}
    # store the contribution (percentage) of the fastest reaction between the two species
    max_contri = {}
    # initializing zero
    for pair in ana_comb: 
        netRate_pair[pair] = 0.
        # reverse
        netRate_pair[pair[::-1] ] = 0.
        # first element stores the rate and the 2md element stores the reaction number 
        maxRate[pair] = [0. , np.nan]
        maxRate[pair[::-1]] = [0., np.nan]

        edges[pair[0]].append(pair[1]) 
        edges[pair[1]].append(pair[0])
        edges[pair[0]] = list(set(edges[pair[0]]))
        edges[pair[1]] = list(set(edges[pair[1]]))

    # Setting the current working directory to the script location
    abspath = os.path.abspath(sys.argv[0])
    dname = os.path.dirname(abspath)
    os.chdir(dname)

    from chem_funs import ni, nr  # number of species and reactions in the network
    np.set_printoptions(threshold='nan')  # print all for debuging


    # Only to print the reaction names
     
    # ### reading rates using vulcan modules
    # ## creat the instances for storing the variables and parameters
    data_var = store.Variables()
    data_atm = store.AtmData()
    # ##
    data_atm.Tco = T_ana
    data_atm.pco = p_ana
    # data_atm.M = mm

    # for reading rates
    rate = op.ReadRate()
    # read-in the network
    data_var = rate.read_rate(data_var, data_atm)
    # reversing rates
 #    data_var = rate.rev_rate(data_var, data_atm)
 #    ### reading rates using vulcan modules




    # Reading rates
    for re in range(1,nr+1): 
        if not re in remove_re:
            re_prod = [prod for prod in re_dict[re][0] if prod in ana_sp]
            re_reat = [reat for reat in re_dict[re][1] if reat in ana_sp]

            comb = []
            for prod in re_prod:
                for reat in re_reat:
                    comb.append((prod, reat)) #tuple
            re_comb[re] = comb

            # collecting rates (not rate constants) for each combination in ana_comb
            for pair in re_comb[re]:
                
                if pair in ana_comb or pair[::-1] in ana_comb:
                    
                    
                    rate_const = data['variable']['k'][re][p_indx]
                    #rate_const = data_var.k[re]
                    
                    
                    this_rate = rate_const
                    for sp in re_wM_dict[re][0]: # Looping all the "reactants" including M
                        # sp needs to be in ana_sp?
                        # if sp in ana_sp:
                        if not sp == 'M':
                            
                            #if sp == 'CH3CO': sp = 'C2H3O'
                             
                            #this_rate *= float(fEQ(sp, p_ana, T_ana)) * mm # # density
                            this_rate *= data['variable']['y'][p_indx,species.index(sp)] # # density
                            
                            
                        elif sp == 'M': this_rate *= mm
                        else: raise IOError ('\nUnknow species name in the network.')
                    
                    #if re == 307: print (str(re) + "{:>10.2e}".format(this_rate)) 
                    #if re == 291: print (str(re) + "{:>10.2e}".format(this_rate)) 
                    if re == 341: print (str(re) + "{:>10.2e}".format(this_rate)) 
                    if re == 605: print (str(re) + "{:>10.2e}".format(this_rate)) 
                    
                    if np.abs(this_rate) > maxRate[pair][0]: 
                        maxRate[pair][0] = np.abs(this_rate)
                        maxRate[pair][1] = re
                    netRate_pair[pair] +=  this_rate

    # Calculating the contribution of the main reactions in maxRate
    for pair in ana_comb:
        if netRate_pair[pair] == 0:
            max_contri[pair] = np.nan
        elif netRate_pair[pair[::-1] ] == 0:
            max_contri[pair[::-1] ] = np.nan
        else:
            max_contri[pair] = maxRate[pair][0]/netRate_pair[pair]
            max_contri[pair[::-1] ] = maxRate[pair[::-1] ][0]/netRate_pair[pair[::-1] ]



    # Graph Data
    '''
    ana_sp = nodes  
    edges_dict = edges
    maxRate = distance
    '''

    class Graph:
      def __init__(self):
        self.nodes = ana_sp
        self.edges = edges
        self.distance = maxRate

    g = Graph()

    def shortest_path(graph, ini, end):
        edges = graph.edges.copy()
        distance = graph.distance.copy()

        # Initialization
        distance_indx = defaultdict(lambda:np.inf)
        distance_indx[ini] = 0
        unvisited_distance = {node: np.inf for node in list(graph.nodes)}
        visited, path = [], []
        prev = {}
        nodes = set(graph.nodes)
        # the current node
        now = ini

        while unvisited_distance: # while unvisited_distance is not empty
            for neighbor in edges[now]:
                if neighbor in unvisited_distance.keys(): # never check the visited nodes again!
                    rate = distance[(now, neighbor)][0]
                    if not rate == 0:
                        temp_distance = distance_indx[now] + 1./rate 

                        if temp_distance < distance_indx[neighbor]: 
                            # distance_indx is only for recording footprints; unvisited_distance is for the algorithm 
                            prev[neighbor] = now # recording the pathway
                            distance_indx[neighbor] = temp_distance  
                            unvisited_distance[neighbor] = temp_distance 
    
            # marking the node "JUST VISITED"
            visited.append(now) 
            del unvisited_distance[now]

            # select the node that is marked with the smallest distance_indx from ALL unvisited nodes
            next_now = min(unvisited_distance,key=unvisited_distance.get) # key = xx.get gives the access to the values
            now = next_now
    
            # check if stops 
            if now == end: 
                visited.append(now)
                path.append(now)
                goback = now
                while not goback==ini:
                    path.append(prev[goback])
                    goback = prev[goback]
                return (path[::-1], distance_indx)
            elif unvisited_distance[next_now] == np.inf:
                print ('Not connected!')
                break
        

    # Plotting     
    # Creating path graph
    gdot = pydot.Dot(graph_type='graph')

    # to weed out taking the same reaction twice
    # e.g. R1: A + B -> C + D 
    # A (R1)-> C (R1)-> B is not allowed

    # also exclude the reaction that produces the straing species again
    # e.g. Path A -> B -> C 
    # the dominant reaction for B -> C to be B + X -> C + "A" should not be allowed!!!



    path, time_lables = shortest_path(g, conv_sp[0], conv_sp[1])
    pathway = {}
    min_rate = np.inf
    max_rate = 0.
    for sp in path:
        if not sp==path[-1]:
            #print ( (sp, (path[path.index(sp)+1]) ) )
            pair = (sp, (path[path.index(sp)+1]))
            pathway[pair] = maxRate[pair]
            if maxRate[pair][0] < min_rate:
                min_rate = maxRate[pair][0]
                rls = pair
                rls_re = maxRate[pair][1]
            
            # for scaling the plot
            if maxRate[pair][0] > max_rate:
                max_rate = maxRate[pair][0] 


    # maximum and minumum rates used to scale the edge width
    w_min = np.log10(min_rate)
    w_max = np.log10(max_rate)


    for sp in path:
        node = pydot.Node(sp, color="k",fontcolor='k',style='bold')
        gdot.add_node(node)

    # RLS
    if rls_re % 2 == 0: rls_re -= 1
    re_exp = ' ' + data_var.Rf[rls_re] + ' '
    re_exp = re_exp.replace("->", "=")
    edge = pydot.Edge( rls[0], rls[1], label="{:0.2E}".format(pathway[rls][0]) + '\nR' + str(pathway[rls][1])\
        + ' ('+"{:0.2f}".format(max_contri[rls]*100.) +'%)\n' + re_exp, penwidth=1, color='red', fontsize=13, fontcolor='red')
    gdot.add_edge(edge)

    for pair in pathway.keys():
        if pair is not rls:
            lwidth = 1. + 4.5*(np.log10(pathway[pair][0]) - w_min)/(w_max-w_min)
            main_re = pathway[pair][1]
            if main_re % 2 == 0: main_re -= 1
            re_exp = ' ' + data_var.Rf[main_re] + " "
            re_exp = re_exp.replace("->", "=")
            edge = pydot.Edge( pair[0], pair[1], label="{:0.2E}".format(pathway[pair][0]) + '\nR' + str(pathway[pair][1])\
            + ' ('+"{:0.2f}".format(max_contri[pair]*100.) +'%)\n' + re_exp, penwidth=lwidth, color='k', fontsize=12, fontcolor='k', nodecolor='red')

            gdot.add_edge(edge)


    outfile = 'plot/pathways/kinetics_' + conv_sp[0]+"-"+conv_sp[1]+'_' + "{:0.1E}".format(p_ana/1e6) + '_bar_T' + str(int(T_ana))+'_path.png'
    gdot.write_png(outfile)
    plot = Image.open(outfile)
    #plot.show()
    #gdot.write_pdf('plot/pathways/' + conv_sp[0]+"-"+conv_sp[1]+'_' + "{:0.1E}".format(p_ana) + '_bar_T' + str(int(T_ana))+'_path.pdf')

