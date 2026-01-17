#!/usr/bin/env python3
# -*- coding: utf-8 -*-

##################################################################################
## Python script to simulate 3D tumour growth and multi-region sequencing data ##
## using a deme-based agent-based model with flexible spatial constrain        ##
##                                                                             ##
## Deme model structure is adapted from:                                       ##
##   - R. Sun, Z. Hu, A. Sottoriva, T. A. Graham, A. Harpak, Z. Ma, et al.,    ##
##     "Polyclonal-to-monoclonal transition in colorectal precancerous         ##
##      evolution."                                                            ##
##   - Z. Lu, S. Mo, D. Xie, X. Zhai, S. Deng, K. Zhou, et al.,                ##
##     "Between-region genetic divergence reflects the mode and tempo of       ##
##      tumor evolution."                                                      ##
##                                                                             ##
## The pushing algorithm is adapted from:                                      ##
##   - J. Househam, T. Heide, G. D. Cresswell, I. Spiteri, C. Kimberley,       ##
##     L. Zapata, et al., "Phenotypic plasticity and genetic control in        ##
##     colorectal cancer evolution."                                           ##
##   - K. Chkhaidze, T. Heide, B. Werner, M. J. Williams, W. Huang,            ##
##     G. Caravagna, et al., "Spatially constrained tumour growth affects      ##
##     the patterns of clonal selection and neutral drift in cancer            ##
##     genomic data."                                                          ##
##                                                                             ##
###################################################################################


import sys, math, random
import numpy as np
from collections import Counter
import heapq
import subprocess


class deme():
    def __init__(self):
        self.present = 0         ## whether the deme is empty or occupied: 0-empty;1-occupied
        self.background = []     ## the background founder lineage after tumor tranformation
        self.advant = []         ## the advantageous cells

        
def schedule_event(coord, t):
    """
    Schedule the next event time for a deme at coord in the global event heap and map.
    coord - lattice coordinate of the deme (2D: (x,y) or 3D: (x,y,z))
    t - scheduled global time
    """
    event_time_map[coord] = t
    heapq.heappush(event_heap, (t, coord))

def pop_next_event():
    """ 
    Pop the earliest valid event from the global event heap and keep event_time_map consistent.
    """
    while event_heap:
        t, coord = heapq.heappop(event_heap)
        if event_time_map.get(coord) != t:
            continue
        del event_time_map[coord]
        return coord, t
    raise KeyError("No more events!")


def createLattice(rd, dim):
    """
    Create a 2D or 3d cubic lattice with side length of 2rd+1 where each site contains a empty deme.
    rd - radius of the lattice (half-side length of the grid)
    dim - simulation dimensionality (2 for 2D, 3 for 3D)
    """
    lattice = {}
    if dim == 3:
        for x in range(0, 2*rd + 1):
            for y in range(0, 2*rd + 1):
                for z in range(0, 2*rd + 1):
                    lattice[(x, y, z)] = deme()
    elif dim == 2:
        for x in range(0, 2*rd + 1):
            for y in range(0, 2*rd + 1):
                lattice[(x, y)] = deme()
    else:
        raise ValueError("Only support 2D or 3D simualtion")
    return lattice


def neighbor8(x, y):
    """
    2D Moore neighbourhood: 8 neighbour sites of (x,y).
    """
    neighbors = []
    for dx in [-1, 0, 1]:
        for dy in [-1, 0, 1]:
            if dx == 0 and dy == 0:
                continue
            neighbors.append((x + dx, y + dy))
    return neighbors

def neighbor26(a, b, c):
    """
    3D Moore neighbourhood: 26 neighbour sites of (a,b,c).
    """
    neighbor = [(a-1, b-1, c-1),(a-1, b-1, c),(a-1, b-1, c+1),
                (a-1, b, c-1),(a-1, b, c),(a-1, b, c+1),
                (a-1, b+1, c-1),(a-1, b+1, c),(a-1, b+1, c+1),
                (a, b-1, c-1),(a, b-1, c),(a, b-1, c+1),
                (a, b, c-1),(a, b, c+1),
                (a, b+1, c-1),(a, b+1, c),(a, b+1, c+1),
                (a+1, b-1, c-1),(a+1, b-1, c),(a+1, b-1, c+1),
                (a+1, b, c-1),(a+1, b, c),(a+1, b, c+1),
                (a+1, b+1, c-1),(a+1, b+1, c),(a+1, b+1, c+1)]
    return neighbor


def sign(x):
    """
    Return the sign of a scalar x as -1, 0, or 1.
    """
    return 1 if x>0 else (-1 if x<0 else 0)


def choose_push_direction_deme(space, pos, max_dist, random_push=True):
    """
    Choose a random push direction for a deme at pos and search along it for the nearest empty site within max_dist.
    space - lattice dictionary mapping coordinates to deme objects
    pos - coordinate of the deme to be pushed (2D: (x,y); 3D: (x,y,z))
    max_dist - maximum distance (in lattice steps) to search for an empty site
    random_push - randomly choose a discrete direction before searching
    """
    dim = len(pos)
    # Build the discrete direction list (26 directions)
    if dim == 2:
        dirs = [(dx,dy,0) for dx in (-1,0,1) for dy in (-1,0,1) if not (dx==0 and dy==0)]
    else:
        dirs = [(dx,dy,dz)
                for dx in (-1,0,1)
                for dy in (-1,0,1)
                for dz in (-1,0,1)
                if not (dx==dy==dz==0)]
    # randomly choose a direction and search empty site along this direction           
    if random_push:
        dx_cell, dy_cell, dz_cell = random.choice(dirs)
        x0, y0 = pos[0], pos[1]
        z0 = pos[2] if dim==3 else 0
        for step in range(1, int(max_dist)+1):
            x = x0 + dx_cell*step
            y = y0 + dy_cell*step
            z = z0 + dz_cell*step if dim==3 else 0
            cand = (round(x), round(y)) if dim==2 else (round(x),round(y),round(z))
            if cand not in space:
                break
            if space[cand].present == 0:
                return dx_cell, dy_cell, dz_cell, cand
        return 0, 0, 0 if dim==3 else 0, None
    



def can_divide_near_edge_deme(space, pos, max_division_distance):
    """
    Test whether a deme at pos lies in a boundary region where division with pushing is allowed.
    space - lattice dictionary mapping coordinates to deme objects
    pos - coordinate of the deme (2D: (x,y); 3D: (x,y,z))
    max_division_distance - maximum allowed distance to an empty site for permitting fission
    """
    random_push = True
    dx, dy, dz, target_coord = choose_push_direction_deme(space, pos, max_dist=max_division_distance, random_push=random_push)
    if target_coord is not None:
        if SIMULATION_DIM == 2:
            return True,((dx,dy), target_coord)
        else:
            return True, ((dx, dy, dz), target_coord)
    return False, None



def traceLineage(mlineage, mutid):
    """
    A function to obtain the mutational lineage of a cell from the mutation id of the most recently occurred mutation in the cell. 
    For example, the input ID (most recently occurred mutation) of target cell is "100" and the output is "1-12-35-56-100", which is the mutation lineage of the cell
    
    mlineage - the list that could be used to recover the mutational lineage given the most recent mutation id of a lineage
    mutid - the mutation ID of the most recently occurred mutation in the cell
    """
    recent_muts = mutid.split(',') # it is possible that multiple mutations occur during in a cell division. For instance, the mutation id of most recently occurred mutations is "100,101"
    recent_muts = [int(t) for t in recent_muts]
    first_mut = recent_muts[0] # the first mutation id in a multi-mutation event
    trace = []
    while first_mut > 0:
        trace += recent_muts
        recent_muts = mlineage[first_mut].split(',')
        recent_muts = [int(t) for t in recent_muts]
        first_mut = recent_muts[0]
    return trace

def lowerORupper(value):
    """
    Stochastically round a real value to either its lower or upper integer based on the fractional part.
    value - real-valued number to be rounded
    """
    lower_int = int(value)
    upper_int = lower_int + 1
    if random.random() < value - lower_int:
        return upper_int
    else:
        return lower_int


def initiateFirstDeme(maxsize, lineage, current_id):
    """
    The growth of the initial deme from a single transformed tumor cell via a random discrete-time birth-death process
    maxsize - size limit of a deme
    lineage - a list that stores the lineage information of mutations
    current_id - the starting mutation ID
    """
    global adv_clone, adv_id

    sfit = birth_rate * s_coef

    neu_list = [str(current_id) for i in range(10)] # ensure enough cells 
    adv_list = []
    current_deme_size = 1

    while current_deme_size < maxsize:
        n1, n2 = len(neu_list), len(adv_list)    #n1 and n2 are the current number of neutral founder cells and advantageous cells, respectively
        if n1 > 0:
            neu_qui_number = lowerORupper(n1 * quies_rate)  #number of quiescent cells (remain viable but stop dividing) of neutral lineage in this generation. 
            neu_div_number = lowerORupper(n1 * birth_rate)  #number of dividing cells of neutral lineage in this generation.
            neu_pass_number = neu_qui_number + neu_div_number
            if neu_pass_number > n1:
                neu_qui_number = max(0, n1 - neu_div_number)
                neu_pass_number = neu_qui_number + neu_div_number

            random.shuffle(neu_list)
            neu_qui_list = neu_list[0:neu_qui_number]
            neu_div_list = neu_list[neu_qui_number:neu_pass_number]
            neu_div_list_double = neu_div_list * 2
        else:
            neu_qui_list, neu_div_list_double = [], []

        if n2 > 0:
            adv_qui_number = lowerORupper(n2 * (quies_rate - sfit))     #number of quiescent cells (remain viable but stop dividing) of selected lineage in this generation.
            adv_div_number = lowerORupper(n2 * (birth_rate + sfit))     #number of dividing cells of selected lineage in this generation.
            adv_pass_number = adv_qui_number + adv_div_number
            if adv_pass_number > n2:
                adv_qui_number = max(0, n2 - adv_div_number)
                adv_pass_number = adv_qui_number + adv_div_number

            random.shuffle(adv_list)
            adv_qui_list = adv_list[0:adv_qui_number]
            adv_div_list = adv_list[adv_qui_number:adv_pass_number]
            adv_div_list_double = adv_div_list * 2
        else:
            adv_qui_list, adv_div_list_double = [], []

        n1_double = len(neu_div_list_double)
        n2_double = len(adv_div_list_double)

        if n1_double > 0:
            new_mut1 = np.random.poisson(mut_rate * n1_double)  # the total number of mutations occurring in a generation follows Poission distribution with lambda=u*n
            if new_mut1 > 0:
                mut_assig1 = Counter(np.random.choice(n1_double, new_mut1))
                for x1 in mut_assig1.keys():
                    nmut = mut_assig1[x1]
                    new_ids = range(current_id + 1, current_id + 1 + nmut)
                    mut_str = ",".join(map(str, new_ids))
                    for _ in range(nmut):
                        current_id += 1
                        lineage += [neu_div_list_double[x1]]
                    neu_div_list_double[x1] = mut_str

        if n2_double > 0:
            new_mut2 = np.random.poisson(mut_rate * n2_double)
            if new_mut2 > 0:
                mut_assig2 = Counter(np.random.choice(n2_double, new_mut2))
                for x2 in mut_assig2.keys():
                    nmut = mut_assig2[x2]
                    new_ids = range(current_id + 1, current_id + 1 + nmut)
                    mut_str = ",".join(map(str, new_ids))
                    for _ in range(nmut):
                        current_id += 1
                        lineage += [adv_div_list_double[x2]]
                    adv_div_list_double[x2] = mut_str

        if n1_double > 0 and adv_clone < 1 and (random.random() < adv_rate * n1_double):  # occurence of advantageous mutation on the neutral lineage (only allowed one advantageous mutation)
            adv_clone += 1
            current_id += 1
            adv_id = current_id
            lineage += [str(neu_div_list_double[n1_double - 1])]
            adv_div_list_double += [str(current_id)]
            neu_div_list_double = neu_div_list_double[0:n1_double - 1]

        neu_list = neu_qui_list + neu_div_list_double
        adv_list = adv_qui_list + adv_div_list_double
        print(adv_list)
        current_deme_size = len(neu_list) + len(adv_list)

    return neu_list, adv_list, current_id, lineage



def demeGrowthFission(neu_list, adv_list, lineage, current_id, current_deme_number):
    """
    A function to simulate deme expansion and fission and keep track of the mutational lineages.
    neu_list - list of mutation IDs for neutral cells in the deme
    adv_list - list of mutation IDs for advantaged cells in the deme
    lineage - list storing the lineage information of all mutations
    current_id - current maximum mutation ID used in the simulation
    current_deme_number - current total number of demes in the tumor
    """
    global adv_clone, adv_id
    sfit = birth_rate * s_coef 

    current_deme_size = len(neu_list) + len(adv_list)
    while current_deme_size < 2 * deme_size:        # when the deme size doubles, it will split into two offspring demes
        n1, n2 = len(neu_list), len(adv_list)
        neu_qui_list, neu_div_list_double = [], []
        if n1 > 0:
            neu_qui_number = lowerORupper(n1 * quies_rate)      # number of quiescent cells (remain viable but stop dividing) of neutral lineage in this generation.
            neu_div_number = lowerORupper(n1 * birth_rate)      #number of dividing cells of neutral lineage in this generation.
            neu_pass_number = neu_qui_number + neu_div_number
            if neu_pass_number > n1:
                neu_qui_number = max(0, n1 - neu_div_number)
                neu_pass_number = neu_qui_number + neu_div_number

            random.shuffle(neu_list)
            neu_qui_list = neu_list[0:neu_qui_number]
            neu_div_list = neu_list[neu_qui_number:neu_pass_number]
            neu_div_list_double = neu_div_list * 2  

        adv_qui_list, adv_div_list_double = [], []
        if n2 > 0:
            adv_qui_number = lowerORupper(n2 * (quies_rate - sfit))      #number of quiescent cells (remain viable but stop dividing) of selected lineage in this generation.
            adv_div_number = lowerORupper(n2 * (birth_rate + sfit))      #number of dividing cells of neutral lineage in this generation.
            adv_pass_number = adv_qui_number + adv_div_number
            if adv_pass_number > n2:
                adv_qui_number = max(0, n2 - adv_div_number)
                adv_pass_number = adv_qui_number + adv_div_number
            random.shuffle(adv_list)
            adv_qui_list = adv_list[0:adv_qui_number]
            adv_div_list = adv_list[adv_qui_number:adv_pass_number]
            adv_div_list_double = adv_div_list * 2

        n1_double = len(neu_div_list_double)
        n2_double = len(adv_div_list_double)
        if current_deme_number < 5*pow(10,7)/deme_size:     #stop mutation occurring when the tumor size is larger than 5*10^7 cells. The reason is that late occuring mutations have very small chance to present at detectable frequency even under selection.
            if n1_double > 0:
                new_mut_count1 = np.random.poisson(mut_rate * n1_double)
                if new_mut_count1 > 0:
                    mut_assig1 = Counter(np.random.choice(n1_double, new_mut_count1))
                    for x1 in mut_assig1.keys():
                        nmut = mut_assig1[x1]
                        new_ids = range(current_id + 1, current_id + 1 + nmut)
                        mut_str = ",".join(map(str, new_ids))
                        for _ in range(nmut):
                            current_id += 1
                            lineage += [neu_div_list_double[x1]]
                        neu_div_list_double[x1] = mut_str

            if n2_double > 0:
                new_mut_count2 = np.random.poisson(mut_rate * n2_double)
                if new_mut_count2 > 0:
                    mut_assig2 = Counter(np.random.choice(n2_double, new_mut_count2))
                    for x2 in mut_assig2.keys():
                        nmut = mut_assig2[x2]
                        new_ids = range(current_id + 1, current_id + 1 + nmut)
                        mut_str = ",".join(map(str, new_ids))
                        for _ in range(nmut):
                            current_id += 1
                            lineage += [adv_div_list_double[x2]]
                        adv_div_list_double[x2] = mut_str

            if adv_clone < 1:
                if random.random() < adv_rate * n1_double:
                    adv_clone += 1
                    current_id += 1
                    adv_id = current_id
                    lineage += [str(neu_div_list_double[n1_double - 1])]
                    adv_div_list_double += [str(current_id)]
                    neu_div_list_double = neu_div_list_double[0:n1_double - 1]

        neu_list = neu_qui_list + neu_div_list_double
        adv_list = adv_qui_list + adv_div_list_double
        current_deme_size = len(neu_list) + len(adv_list)

    random.shuffle(neu_list)
    if len(neu_list) > 0:
        offspring_neu = np.random.binomial(len(neu_list), 0.5) # the offpring deme size is determined by a Binomial distribution B(n,0.5)
    else:
        offspring_neu = 0
    neu_list1 = neu_list[0:offspring_neu]
    neu_list2 = neu_list[offspring_neu:len(neu_list)]

    random.shuffle(adv_list)
    if len(adv_list) > 0:
        offspring_adv = np.random.binomial(len(adv_list), 0.5)
    else:
        offspring_adv = 0
    adv_list1 = adv_list[0:offspring_adv]
    adv_list2 = adv_list[offspring_adv:len(adv_list)]

    return neu_list1, neu_list2, adv_list1, adv_list2, current_id, lineage


def seqProcessing(sp,sample_keys,mlineage,size_par,mean_depth,purity):
    """
    Model the random sampling process in NGS and report the sequencing allele frequencies in a sample of cells
    
    sp- the lattice space
    sample_keys- the locations for the demes in a bulk sample
    size_par- variance parameter for negative-binomial distribution
    mean_depth- the mean depth of the sequencing
    purity- tumor purity
    """
    all_cur_id = []                                     # all most recently occurred mutations
    all_mut_id = []                                     # all mutations in the sampled cells
    for key in sample_keys:
        smuts = list(sp[key].background + sp[key].advant)
        all_cur_id += smuts
    sample_size = len(all_cur_id)
    sample_id = all_cur_id
    id_count = Counter(sample_id)
    for x in id_count.keys():
        xlineage = traceLineage(mlineage,x)
        all_mut_id += xlineage*id_count[x]
    mut_count = Counter(all_mut_id)
    prob_par=size_par*1.0/(size_par+mean_depth)
    sampleAF = {}                                       # a dictionary storing the mutation IDs and corresponding depth and allele frequency the seq data
    for x in mut_count.keys():
        true_af = mut_count[x]*0.5*purity/sample_size   # the true allele frequency in the sample
        if true_af > 0.001:                             # filter mutations with very low frequency that is not detectable by ~100X sequencing depth
            site_depth = np.random.negative_binomial(size_par,prob_par)
            if site_depth >= 2:                        # seq depth cutoff for "calling" a mutation
                var_reads = np.random.binomial(site_depth,true_af)
                seq_af = var_reads*1.0/site_depth
                if var_reads >= 2:                      # variant reads cutof for "calling" a mutation
                    sampleAF[str(x)] = (site_depth,seq_af)
    return sampleAF


def pubMutGenerator(n,size_par,mean_depth,purity):
    """
    A function to generate the public clonal mutations occurred during the multi-step tumorigenesis before transformation.
    
    n- number of clonal mutations
    size_par- variation parameter in the negative binomial distribution
    mean_death- mean seq depth
    """
    prob_par=size_par*1.0/(size_par+mean_depth)
    mean_af = 0.5*purity
    depth_pub = []
    maf_pub = []
    for k in range(0,n):
        correct = 0
        while correct == 0:
            site_depth = np.random.negative_binomial(size_par,prob_par)
            if site_depth >= 15:
                correct =1
        var_reads = np.random.binomial(site_depth,mean_af)
        site_maf = var_reads*1.0/site_depth
        depth_pub += [site_depth]
        maf_pub += [site_maf]
    return depth_pub,maf_pub


def bulkTissueSampling(sp,location,radius):
    """
    A function to sampling a bulk sample in a local region.
    """
    if SIMULATION_DIM == 2:
        location = (location[0], location[1], 0)
    local_region = localNeighbor_square(location[0],location[1],location[2],radius)
    bulk_tissue = []
    for x in local_region:
        try: 
            sp[x]
        except KeyError:
            continue
        if sp[x].present == 1:
            bulk_tissue += [x]
    return bulk_tissue
        
def missingDepth(vafdata,absent_muts,mean_depth):
    """
    Randomly generate the sequencing depth for the mutation-absent sites across samples
    """
    for x in absent_muts:
        done = 0
        while done == 0:
            missing_depth = np.random.negative_binomial(2,2.0/(2+mean_depth))
            if missing_depth >= 15:
                done = 1
        vafdata[str(x)] = (missing_depth,0)
    return vafdata

def record_snapshot(label, space, SIMULATION_DIM, rd):
    """
    Record a spatial snapshot of all present demes to a text file (2D slice for 3D, full grid for 2D).
    """
    filename = f"{path}/{snapname}_prop{label}.txt"
    with open(filename, "w") as f:
        if SIMULATION_DIM == 3:
            f.write("x\ty\tcategory\tmutation_count\n")
            for coord, deme in space.items():
                # 只记录 z==rd 的平面
                if len(coord) >= 3 and coord[2] == rd and deme.present:
                    x, y = coord[0], coord[1]
                    category = 2 if len(deme.advant) > 0 else 1
                    mut_count = get_mutation_count(deme)
                    f.write(f"{x}\t{y}\t{category}\t{mut_count}\n")
        else:
            f.write("x\ty\tcategory\tmutation_count\n")
            for coord, deme in space.items():
                if deme.present:
                    x, y = coord[0], coord[1]
                    category = 2 if len(deme.advant) > 0 else 1
                    mut_count = get_mutation_count(deme)
                    f.write(f"{x}\t{y}\t{category}\t{mut_count}\n")
    print(f"Snapshot {label} saved to {filename}.")
    

def localNeighbor_square(a, b, c, r):
    """
    Enumerate lattice sites within a square (2D) or cube (3D) neighbourhood of radius r around (a, b, c).
    """
    neighbor = []
    if SIMULATION_DIM == 2:
        # 仅返回二维邻域，忽略 c 的影响
        for dx in range(-r, r):
            for dy in range(-r, r):
                neighbor.append((a+dx, b+dy))
    elif SIMULATION_DIM == 3:
        for dx in range(-r, r):
            for dy in range(-r, r):
                for dz in range(-r, r):
                    neighbor.append((a+dx, b+dy, c+dz))
    else:
        raise ValueError("SIMULATION_DIM必须为2或3")
    return neighbor

def generate_samples(punch_diameter, spacing):
    """
    Generate sampling locations across the tumor, forming a square grid in 2D or a cubic grid in 3D.
    punch_diameter - effective diameter of each biopsy region, used to set sampling radius
    spacing - gap between adjacent sampling centers along x and y axis
    """
    sample_list = []

    if SIMULATION_DIM == 2:
        x_coords = range(min_x + punch_diameter, max_x, punch_diameter + spacing)
        y_coords = range(min_y + punch_diameter, max_y, punch_diameter + spacing)

        for x in x_coords:
            for y in y_coords:
                tissue_tmp = bulkTissueSampling(space, (round(x), round(y)), int(punch_diameter / 2))
                if len(tissue_tmp) >= puch_density * pow(int(punch_diameter / 2) * 2, 2):       # Add this sample if it covers enough tumor demes (puch_density of sampling cube)
                    sample_list.append((round(x), round(y), 0))

    elif SIMULATION_DIM == 3:
        x_coords = range(min_x + punch_diameter, max_x, punch_diameter + spacing)
        y_coords = range(min_y + punch_diameter, max_y, punch_diameter + spacing)
        z_coords = [rd]

        for x in x_coords:
            for y in y_coords:
                for z in z_coords:
                    tissue_tmp = bulkTissueSampling(space, (round(x), round(y), round(z)), int(punch_diameter / 2))
                    if len(tissue_tmp) >= puch_density * pow(int(punch_diameter / 2) * 2, 3):
                        sample_list.append((round(x), round(y), round(z)))
    else:
        raise ValueError("SIMULATION_DIM must be 2 or 3")

    return sample_list



def get_mutation_count(deme):
    """
    Count the total number of mutations carried by all cells in a given deme.
    """
    all_cur_id = []                                     # all most recently occurred mutations
    all_mut_id = []                                     # all mutations in the sampled cells
    all_cur_id = list(deme.background + deme.advant)
    for x in all_cur_id:
        xlineage = traceLineage(mutlineage,x)
        all_mut_id += xlineage
    mut_count = len(all_mut_id)

    return mut_count


def _min_gens_to_fission_quies(n_neu, n_adv, birth_rate, s_coef, quies_rate, deme_size):
    """
    Compute the minimal integer number of generations needed for a deme to reach size 2*deme_size under quiescent/dividing dynamics.
    n_neu - initial number of neutral cells
    n_adv - initial number of advantaged cells
    birth_rate - baseline birth rate for neutral cells
    s_coef - selective advantage coefficient scaling the birth rate of advantaged cells
    quies_rate - probability that a cell becomes quiescent instead of dividing per generation
    deme_size - target deme size (before fission) used to define the 2*deme_size threshold
    """
    target = 2 * deme_size
    n0 = n_neu + n_adv
    if n0 >= target:
        return 1  

    if n0 == 0:
        return None

    b = float(birth_rate)
    q = float(quies_rate)
    s = float(birth_rate * s_coef)  # sfit

    a_n = q + 2.0 * b          # neutral per-generation growth factor
    a_a = q + 2.0 * b + s      # advantageous per-generation growth factor (q-s)+2(b+s) = q+2b+2s

    # Check reachability: if all existing types have growth factor ≤ 1, the target cannot be reached
    if (n_neu > 0 and a_n <= 1.0) and (n_adv > 0 and a_a <= 1.0):
        return None
    if n_neu == 0 and a_a <= 1.0:
        return None
    if n_adv == 0 and a_n <= 1.0:
        return None

    # Get a reasonable lower bound g0 from the fastest-growing type
    fastest = max(a_n if n_neu > 0 else 0.0, a_a if n_adv > 0 else 0.0)
    g0 = 1
    if fastest > 1.0:
        need = target / float(n0)
        if need > 1.0:
            g0 = max(1, int(math.floor(math.log(need, fastest))))

    # Linearly scan from g0 upwards until reaching the target size, or give up at MAX_G
    MAX_G = 100000
    for g in range(g0, MAX_G + 1):
        val = (n_neu * (a_n ** g) if n_neu > 0 else 0.0) + (n_adv * (a_a ** g) if n_adv > 0 else 0.0)
        if val >= target:
            return g
    return None  # 认为不可达




def _schedule_next_with_expected_divisions(coord, now):
    """
    Estimate the expected generations to fission for a deme and schedule its next event as an exponential waiting time.
    coord - lattice coordinate of the deme whose next event is being scheduled
    now - current simulation time 
    """
    n_neu = len(space[coord].background)
    n_adv = len(space[coord].advant)
    n_tot = n_neu + n_adv


    if n_tot == 0 or n_tot >= 2 * deme_size:
        return

    g = _min_gens_to_fission_quies(
        n_neu, n_adv,
        birth_rate=birth_rate,
        s_coef=s_coef,
        quies_rate=quies_rate,
        deme_size=deme_size
    )
    if g is None or g <= 0:
        return  # if fission is unreachable, do not schedule

    rate = 1.0 / g
    dt = np.random.exponential(1.0 / rate)  # scale = g
    schedule_event(coord, now + dt)

def move_event_timer_simple(src, dst):
    """
    Move the scheduled next-event time from coord src to coord dst and reinsert it into the event heap.
    src - original lattice coordinate whose event time should be moved
    dst - new lattice coordinate that should inherit the event time
    """
    old_time = event_time_map.pop(src, None)
    if old_time is not None:
        event_time_map[dst] = old_time
        heapq.heappush(event_heap, (old_time, dst))

def migrate_one_random_push(space, coord, max_push):
    """
    Attempt a single random-direction migration for the deme at coord, allowing a push chain along that direction.
    space - lattice dictionary mapping coordinates to deme objects
    coord - coordinate of the deme to be migrated (2D: (x,y); 3D: (x,y,z))
    max_push - maximum lattice distance allowed when searching for an empty site to push into
    """
    # Skip if this coord is out of bounds or no longer occupied
    if coord not in space or space[coord].present != 1:
        return None

    #  pick a direction and find the nearest empty site within max_push
    dx, dy, dz, target = choose_push_direction_deme(
        space, coord, max_dist=max_push, random_push=True
    )
    if target is None:
        return None  

    # Calculate empty neighbor coord created by push
    dim = len(coord)
    if dim == 2:
        step = (dx, dy)
        nextkey = (coord[0] + dx, coord[1] + dy)
    else:
        step = (dx, dy, dz)
        nextkey = (coord[0] + dx, coord[1] + dy, coord[2] + dz)

    # psuh algorithm
    push_target = target
    while True:
        push_before = tuple(c - d for c, d in zip(push_target, step))
        if push_before == coord:
            break
        space[push_target] = space[push_before]
        space[push_before] = deme()
        move_event_timer_simple(push_before, push_target)  
        push_target = push_before

    # Move the original deme from coord into nextkey
    space[nextkey] = space[coord]
    space[coord]   = deme()
    space[nextkey].present = 1
    move_event_timer_simple(coord, nextkey)

    return nextkey

def perform_migration_batch(space, p, max_push):
    """
    Perform one batch migration step where a fraction p of occupied demes attempt random-push migration.
    space - lattice dictionary mapping coordinates to deme objects
    p - fraction of currently occupied demes to select for migration attempts (0–1)
    max_push - maximum lattice distance allowed when searching for an empty site during pushes
    """
    present_coords = [k for k, v in space.items() if v.present == 1]
    N = len(present_coords)
    if N == 0 or p <= 0:
        return 0, 0

    # select floor(N * p) demes
    k = int(N * p)
    if k <= 0:
        return 0, 0

    # Randomly sample k distinct demes for migration attempts
    selected = random.sample(present_coords, min(k, N))
    success = 0
    for coord in selected:
        # If this coord has already been pushed earlier in this batch, it may now be empty
        if coord not in space or space[coord].present != 1:
            continue
        moved_to = migrate_one_random_push(space, coord, max_push)
        if moved_to is not None:
            success += 1
    return k, success


# ========================= Parameter initialisation =========================
# Global simulation settings (can be adjusted)
adv_clone = 0           # count for whether an advantageous clone
npub = 175              # number of public (clonal) mutations to generate
seq_depth = 100         # mean sequencing depth

# Command-line input parameters
sim_dim = int(sys.argv[1])          # 2 for 2D simulations, 3 for 3D simulations
deme_size = int(sys.argv[2])        # the deme size
mut_rate = float(sys.argv[3])       # the neutral mutation rate at whole exonic region
adv_rate = float(sys.argv[4])       # the advantageous mutation rate 
s_coef = float(sys.argv[5])         # the selection coefficient
repl = int(sys.argv[6])             # replication of simulation
path = str(sys.argv[7])             # output directory base path
rd = int(sys.argv[8])               # the side length of the 2/3D space
birth_rate = float(sys.argv[9])     # birth probability
death_rate = float(sys.argv[10])     # death probability (birth_rate + death_rate <= 1)
push_prop = float(sys.argv[11])     # proportion of radius allowed as max push distance
mig_rate = float(sys.argv[12])      # migration probability for a deme 
punch_diameter = int(sys.argv[13])  # diameter of the biopsy / sampling region
puch_density = float(sys.argv[14])  # minimum fraction of tumor demes that must be present within the sampling cube for the punch biopsy to be considered valid
spacing = int(sys.argv[15])         # spacing between bulk sampling centres
title = str(sys.argv[16])           # prefix for output file names                
snapname = str(sys.argv[17])        # prefix for snapshot file names

# deme_size=100
# mut_rate=0.6
# adv_rate=0.1
# s_coef=0.4
# repl=4
# path="./csv"
# rd=10
# birth_rate=0.4
# death_rate=0.3
# push_prop=0.1
# mig_rate=0
# punch_diameter=3
# puch_density=0.1
# spacing=3
# title="example"
# snapname="snap_example"


# Derived parameters for pushing and quiescence
max_push = 2*rd*push_prop                   # maximum push distance in lattice units
quies_rate = 1 - birth_rate - death_rate    # probability of becoming quiescent
MIGRATION_PROB = mig_rate
SIMULATION_DIM = sim_dim
# Random seeds
random_seed = repl
random.seed(random_seed)
np.random.seed(random_seed)

# Output directory
subprocess.call(f"mkdir -p {path}", shell=True)

######################################################################################
# Initial mutation lineage and first deme
mut_id = 0
mutlineage = ['0']
first_neu, first_adv, mut_id, mutlineage = initiateFirstDeme(deme_size, mutlineage, mut_id)

# Lattice / space initialisation
space = createLattice(rd, SIMULATION_DIM)
if SIMULATION_DIM == 3:
    initial_coord = (rd, rd, rd)
    final_tumor_size = pow(rd,3) * (4/3) * 3.1415926 * deme_size # the number of cells in the final tumor
    final_deme_number = final_tumor_size/deme_size    # the final number of demes in the tumor
else:  # 2D
    initial_coord = (rd, rd)
    final_tumor_size = pow(rd,2)  * 3.1415926 * deme_size  # the number of cells in the final tumor
    final_deme_number = final_tumor_size/deme_size    # the final number of demes in the tumor
space[initial_coord].present = 1
space[initial_coord].background = list(first_neu)
space[initial_coord].advant = list(first_adv)
space[initial_coord].activate = True
current_keys = [initial_coord]
current_deme_number = 1
current_time = 0
boundary = False

# Predefined thresholds at which to record spatial snapshots (fractions of final deme number)
threshold_fractions = [1/32, 1/16, 1/8, 2/8, 4/8, 6/8]
thresholds = [final_deme_number * 1/2 * frac for frac in threshold_fractions] 
threshold_labels = ['1_32', '1_16', '1_8', '2_8', '4_8', "6_8"]
snapshots_recorded = {label: False for label in threshold_labels}

# Initialize the event heap for deme-fission timing events priority queue `event_heap` to store `(time, coord)` pairs
event_heap = []            # list-based min-heap storing tuples (time, coord)
event_time_map = {}        # dict mapping coord -> latest scheduled event time
schedule_event(initial_coord, current_time + 1.0)  # schedule first event at the initial deme

# ========================= Main tumour growth loop =========================
while current_deme_number < final_deme_number:
    # Check and record snapshots when tumour size crosses predefined thresholds  
    for idx, thresh in enumerate(thresholds):
        label = threshold_labels[idx]
        if (current_deme_number >= thresh) and (not snapshots_recorded[label]):
            record_snapshot(label, space, SIMULATION_DIM, rd)
            snapshots_recorded[label] = True
    if MIGRATION_PROB > 0:
        # Batch migration step (random deme migration with pushing) --- 
        _ksel, _kdone = perform_migration_batch(space, MIGRATION_PROB, max_push)

    # Pop the least scheduled deme event 
    ckey, event_time = pop_next_event()

    current_time = event_time
    
    # Collect empty neighbour sites
    if SIMULATION_DIM == 3:
        rx, ry, rz = ckey
        nei_sites = neighbor26(rx, ry, rz)
        center_coord = (rd, rd, rd)
    else:
        rx, ry = ckey
        nei_sites = neighbor8(rx, ry)
        center_coord = (rd, rd)
    empty_sites = []
    for key in nei_sites:
        if key not in space:
            boundary = True
            break
        if space[key].present == 0:
            empty_sites.append(key)

    if len(empty_sites) > 0:
        # Standard deme fission into an empty neighbour
        (post_neu_l1, post_neu_l2, 
         post_adv_l1, post_adv_l2, 
         mut_id, mutlineage) = demeGrowthFission(space[ckey].background, space[ckey].advant, 
                                                mutlineage, mut_id, current_deme_number)
        space[ckey].background = list(post_neu_l1)
        space[ckey].advant = list(post_adv_l1)
        nextkey = random.choice(empty_sites)
        space[nextkey].background = list(post_neu_l2)
        space[nextkey].advant = list(post_adv_l2)
        space[nextkey].present = 1
        current_deme_number += 1

        # Reschedule events time for two daughter deme
        _schedule_next_with_expected_divisions(ckey, current_time)
        _schedule_next_with_expected_divisions(nextkey, current_time)


    else:
        # deme fission without empty neighbour 
        can_divide, push_params = can_divide_near_edge_deme(space, ckey, max_division_distance=max_push) # check whether pushing division is allowed; if allowed, return push direction and target empty coord
        if not can_divide:
            _schedule_next_with_expected_divisions(ckey, current_time) # Deme cannot divide via pushing; just reschedule its next event
        else:
            psuh_direct, target_coord = push_params
            
            (post_neu1, post_neu2, 
             post_adv1, post_adv2, 
             mut_id, mutlineage) = demeGrowthFission(space[ckey].background, space[ckey].advant, 
                                      mutlineage, mut_id, current_deme_number)
            space[ckey].background = list(post_neu1)
            space[ckey].advant = list(post_adv1)
            current_deme_number += 1
            nextkey = tuple(c + d for c, d in zip(ckey, psuh_direct)) # calculate empty neighbor coord created by push 
            current_keys.append(nextkey) 

            # Push chain: push demes along the direction until the target_coord is reached (target empty coord)
            push_target = target_coord 
            while True:
                push_before = tuple(c - d for c, d in zip(push_target, psuh_direct)) # One step back along the push direction
                if push_before == ckey:
                    break
                space[push_target] = space[push_before] # move deme from push_before -> push_target
                space[push_before] = deme() # reset source to a new empty deme
                old_time = event_time_map.get(push_before) # move scheduled event time along with the deme
                schedule_event(push_target, old_time)
                push_target = push_before
            
            
            space[nextkey].background = list(post_neu2)
            space[nextkey].advant = list(post_adv2)
            space[nextkey].present = 1
            
            
            # deme fission without empty neighbour 
            _schedule_next_with_expected_divisions(ckey, current_time)
            _schedule_next_with_expected_divisions(nextkey, current_time)
                     
    if boundary:
        print("Reach boundary", ckey)
        break

# ========================= Final snapshot at full size (prop8_8) =========================
x_coords_list = []
y_coords_list = []
filename =f"{path}/{snapname}_prop8_8.txt"
with open(filename, "w") as f:
    # For 3D, only record the central z-slice (z == rd); for 2D, record all present demes
    if SIMULATION_DIM == 3:
        f.write("x\ty\tcategory\tmutation_count\n")
        for coord, deme in space.items():
            if len(coord) >= 3 and coord[2] == rd and deme.present:
                x_coords_list.append(coord[0])
                y_coords_list.append(coord[1])
                x, y = coord[0], coord[1]
                category = 2 if len(deme.advant) > 0 else 1
                mut_count = get_mutation_count(deme)
                f.write(f"{x}\t{y}\t{category}\t{mut_count}\n")
    else:
        f.write("x\ty\tcategory\tmutation_count\n")
        for coord, deme in space.items():
            if deme.present:
                x_coords_list.append(coord[0])
                y_coords_list.append(coord[1])
                x, y = coord[0], coord[1]
                category = 2 if len(deme.advant) > 0 else 1
                mut_count = get_mutation_count(deme)
                f.write(f"{x}\t{y}\t{category}\t{mut_count}\n")


# ========================= Bulk sampling region generation =========================
# Compute tumour bounding box in the recorded slice for defining sampling limits
max_x, max_y, min_x, min_y = max(x_coords_list), max(y_coords_list), min(x_coords_list), min(y_coords_list)
sample49 = generate_samples(punch_diameter, spacing)
sample49.sort(key = lambda x:(x[0], x[1]))
location_file = open(path+"/"+title+"_deme_location.txt", "w")
location_file.write("x\ty\tz\tsample\n")
for k in range(len(sample49)):
    if SIMULATION_DIM == 3:
        location_file.write(f"{sample49[k][0]}\t{sample49[k][1]}\t{sample49[k][2]}\tS{str(k+1)}\n")
    else:
        location_file.write(f"{sample49[k][0]}\t{sample49[k][1]}\t{0}\tS{str(k+1)}\n")
location_file.close()

# Build tissue lists for each bulk sample
tissue_list = []
deme_num = 0
for k in range(len(sample49)):
    if SIMULATION_DIM == 3:
        tissue_tmp = bulkTissueSampling(space, sample49[k], int(punch_diameter / 2))
    else: 
        tissue_tmp = bulkTissueSampling(space, (sample49[k][0], sample49[k][1]), int(punch_diameter / 2))
    if len(tissue_tmp) >= puch_density * pow(int(punch_diameter / 2) * 2, 3 if SIMULATION_DIM == 3 else 2):
        tissue_list.append(tissue_tmp)
        deme_num += len(tissue_list[-1])
print("Average # of demes in the samples: " + str(deme_num / len(sample49)))

# VAF matrix construction
maf_list = []
maf_name = []
muts_all = []
MAF_file = open(path+"/"+title+"_deme_vaf.txt", "w")  # write MAF file header
MAF_file.write("mut_id" + "\t")
for k in range(len(sample49)):
    maf_tmp = seqProcessing(space, tissue_list[k], mutlineage, 2, seq_depth, 1)

    if len(maf_tmp) == 0:
        continue
    maf_list.append(maf_tmp)
    maf_name.append("S" + str(k + 1))
    MAF_file.write("depth" + str(len(maf_list)) + "\t" + "S" + str(len(maf_list)) + "\t")  # write MAF file header
    muts_all = set(muts_all) | set(maf_list[-1].keys())  # collect mutations from all the samples
MAF_file.write("\n")

# generate public mutations
for k in range(0, npub):
    pdepth, pmaf = pubMutGenerator(len(sample49), 2, seq_depth, 1)
    MAF_file.write(f"-{k}\t")
    for t in range(len(maf_list)):
        MAF_file.write(str(pdepth[t]) + "\t" + str(pmaf[t]) + "\t")  # write depth and VAF of public mutations into MAF file
    MAF_file.write("\n")

# generate missing depth for private mutations
for k in range(len(maf_list)):
    absent_tmp = muts_all - set(maf_list[k].keys())
    maf_list[k] = missingDepth(maf_list[k], absent_tmp, seq_depth)

# write depth and VAF of private mutations into MAF file
for mt in list(muts_all):
    MAF_file.write(str(mt) + "\t")
    for k in range(len(maf_list)):
        maf_tmp = maf_list[k]
        MAF_file.write(str(maf_tmp[mt][0]) + "\t" + str(maf_tmp[mt][1]) + "\t")
    MAF_file.write("\n")
MAF_file.close()

# Record driver mutation ID
adv_file = open(path+"/"+title+"_adv_mutation.txt", "w")
if adv_clone >=1:
    adv_file.write(f"{adv_id}")
else:
    adv_file.write("none")
adv_file.close()

