#!/usr/bin/env python2.7

"""
    ###spec.py###

    Counts premature stop codon-causing mutations given a multiple alignment file.
    The codon frame must begin at the first base position of the multiple alignment FASTA input.

    This new version of spec.py has a debugged version of the normalize_sbm_counts() function and prints out the list of consensus codon frequency counts. 

"""

__version__ = "16.0"
__author__ = "Jimmy Chan and Jay W. Kim"

import sys
import argparse
import itertools
import random
from random import randrange
from collections import Counter
from collections import defaultdict
import numpy as np

def main(args):
    pargs = parse_args(args)
    fastas = parse_fasta(pargs.input)
    column_seqs = get_column_seqs(fastas)
    if pargs.print_depth:
        mean_depth,std_depth = calc_depth(column_seqs)
        print "{}\t{}".format(mean_depth,std_depth)
        exit(0)
    if pargs.all_muts:
        muts = count_all_muts(column_seqs,pargs.mut_call_threshold)
        spec = get_spectrum(muts)
    else:
        # check the arguments are valid
        if (pargs.apobec_switch != "on" and pargs.apobec_switch != "off"): print "error in argument '-as' usage: '-as on' or '-as off'"; exit(0)
        if (pargs.norm_switch != "on" and pargs.norm_switch != "off"): print "error in argument '-ns' usage: '-ns on' or '-ns off'"; exit(0)
        if (pargs.mut_type != "single" and pargs.mut_type != "double" and pargs.mut_type != "both"): print "error in argument '-m' usage: '-m single' or '-m double' or -m both'"; exit(0)
        if (pargs.base_mut_option != "TAA" and pargs.base_mut_option != "TAG" and pargs.base_mut_option != "all"): print "error in argument '-bm' usage: '-bm TAA' or '-bm TAG' or '-bm all'"; exit(0)
        sbm_count = count_single_muts(pargs.base_mut_option)
        psc_pos = get_pstop_pos(fastas, pargs.base_mut_option)
        psc_ids = get_psc_ids_dict(psc_pos, fastas, pargs.base_mut_option)
        single_muts, breakList = get_muts(psc_pos,psc_ids,column_seqs,fastas,pargs.mut_call_threshold)
        apobec = get_apobec_muts(fastas, breakList, pargs.base_mut_option, pargs.mut_call_threshold)
        cons_codon = count_cons_codons(fastas, pargs.mut_call_threshold) # compiles frequency of consensus codons in a fasta file
        sbm_spec = get_sbm_spectrum(apobec, single_muts, pargs.apobec_switch, pargs.norm_switch, pargs.base_mut_option, cons_codon)
        # fasta_simulation_PSCMs = fasta_simulation(pargs.seq_length, pargs.num_of_seq_entries, pargs.num_of_base_muts) # creates test fasta file
        # for SC_type in fasta_simulation_PSCMs: print SC_type, fasta_simulation_PSCMs[SC_type] # prints out PSCM dictionary containing the PSCMs obtained from the simulated fasta file
        spec = sbm_spec # for the print_spectrum() function to work properly regardless of whether you enter the if pargs.all_muts: or the following else loop 
    if pargs.print_spectrum:
        print_spectrum(spec)
    
def calc_depth(column_seqs):
    depths = []
    for col in column_seqs.values():
        depths.append(col.count('A') + col.count('C') + col.count('G') + col.count('T'))
    depths = np.array(depths)
    return np.mean(depths),np.std(depths)

def capitalizeBases(seq):
    ''' capitalizes the bases of a sequence '''
    capSeq = ""
    switch={'A':'A','C':'C','G':'G','T':'T','a':'A','c':'C','g':'G','t':'T','N':'N','n':'N', '.':'-', '*':'-', '-':'-'}
    return ''.join([switch[b] for b in seq])     

def count_all_muts(column_seqs,threshold): # I've checked that this function works fine (time stamp: 2/7/18)
    muts = defaultdict(dict)
    id = 1
    f = {}
    total = 0
    mutBase = ''
    maxBase = ''
    for column in column_seqs.values(): # for each column
        for c in column: # for each char in a column; going through each char to record each base mutation, i.e., base that differs from the consensus base (base that meets a given variable threshold value)
            if (c == "-" or c == "N"): continue 
            mutBase = c # base that is the mutated/deviated base from the consensus base in a column
            check = -1
            total = 0
            for b in ['A', 'C', 'T', 'G']:
                f[b] = column.count(b)
                total += f[b] 
            for b in ['A', 'C', 'T', 'G']:
                if (f[b] >= total * threshold):
                    f['M'] = f[b]/float(total) # 'M' is the key for the max base 
                    maxBase = b # consensus base in a column
                    check = 1 
            if (mutBase != maxBase and check != -1): 
                muts[id]['freq'] = f['M']
                muts[id]['maxbase'] = maxBase
                muts[id]['mutbase'] = mutBase 
                id += 1
    return muts

def count_cons_codons(fastas, threshold): 
    '''
    keeps track of the frequency of every consensus codon at every AA position in a dictionary (key: cons codon + value: frequency)
    The resulting dict will be used to normalize the SBM spectrum according to the new normalization method described at the start of the script 
    '''
    seq_length = 0 # sequence length
    AAlist = list()
    codons_freq = defaultdict(int) # records the frequency for each consensus codon at every AA position
    consCodon = "" # consensus codon at a certain AA position
    for id, seq in fastas.iteritems(): seq_length = len(seq); break
    for bpos in xrange(1, seq_length-3, 3): #iterating correctly
        skip = 'no' # boolean variable for later use
        for pos in xrange (bpos, bpos+3, 1): # iterate through each codon
            if (get_cons_base(pos, fastas, threshold) == ""):
                skip = 'yes' # skip this codon because it contains one or more null consensus bases
                break # break out of for loop or else skip string may get overwritten from 'yes' back to 'no'
        if (skip == 'no'): # only continue if all three base positions have distinct consensus bases, i.e., no non-nucleotide chars
            consCodon = get_cons_base(bpos, fastas, threshold) + get_cons_base(bpos+1, fastas, threshold) + get_cons_base(bpos+2, fastas, threshold)          
            codons_freq[consCodon] += 1
    print 'mutating codons counts:'
    for codon_type in sorted(codons_freq): print codon_type + ':', codons_freq[codon_type]
    print ''
    return codons_freq

def count_single_muts(base_mut_option): 
    '''
    Counts all possible single-base mutations (SBMs) that cause PSCs
    '''
    table = defaultdict(Counter) # 2D dictionary - dict: base that mutates, key: mutBase, value: number of chances that SBM can cause a PSC
    # must have this block of code before the latter because we're zero-ing out G-A SBMs that cause TAG and TGA
    if (base_mut_option == "TAA" or base_mut_option == "all"):
        for sb in ['T', 'A']: # the following counts SBMs that cause the PSC: TAA
            all = ['A', 'C', 'T', 'G'] 
            all.remove(sb)
            for b in all:
                if (sb == "A"): table['TAA'][b + sb] += 2 # give value of 2 since there's two base positions where the SBM is to A
                else: table['TAA'][b + sb] += 1 # give value of 1 for the SBM to T, the first base position
        table['TAA']['GA'] = 0

    if (base_mut_option == "TAG" or base_mut_option == "all"):
        for sb in ['T', 'A', 'G']: # the following counts SBMs that cause the PSCs: TAG and TGA
            all = ['A', 'C', 'T', 'G'] 
            all.remove(sb)
            for b in all:
                table['TAG'][b + sb] += 1
                table['TGA'][b + sb] += 1 
        table['TAG']['AG'] = 0 # we must zero out the A-G SBMs since they do not occur 
        table['TGA']['AG'] = 0 # we must zero out the A-G SBMs since they do not occur
    
    SBM_norm = defaultdict(Counter)  # key: psc, value: product of DBM chances of occuring
    # the following represents the lowest common product of each SBM count for each psc
    if (base_mut_option == "TAA" or base_mut_option == "all"):
        SBM_norm['TAA']['TAA'] = 2 
    if (base_mut_option == "TAG" or base_mut_option == "all"):
        SBM_norm['TAG']['TAG'] = 1
        SBM_norm['TGA']['TGA'] = 1
    # the following calculates the normalization factor
    # (number needed to multiply the base mut chance to achieve the product for that psc)
    # for each base mutation 
    for psc in table:
        for base_mut in table[psc]:
            if (table[psc][base_mut] != 0):
                SBM_norm[psc][base_mut] = SBM_norm[psc][psc]/table[psc][base_mut]
                SBM_norm[psc]['tot'] += table[psc][base_mut]
    return SBM_norm

def fasta_simulation(seq_length, num_of_seq_entries, num_of_base_muts):
    '''
    Creates a fasta file containing artificial HIV id/seq pairs with a variable number of id-seq pair entries,
    each of a (uniform) variable sequence length containing the same number of variable random base mutations
    '''
    fasta = open('simulation_fasta.txt', 'w')
    # spec_file = open('spectrum_simulation', 'a') # outputs all mutations throughout the simulation_fasta file
    base = ['A', 'C', 'G', 'T']
    SC = ['TAA', 'TAG', 'TGA']
    id_num = 1
    WT_seq = ""
    codon = list() # list representation of a given codon; allows for random base change
    codon_string = "" # string representation of a given codon
    old_codon = "" # original codon that mutates
    new_codon = "" # new codon that contains a SBM (which is possible a PSCM)
    PSCM = defaultdict(Counter) # 2D dictionary containing the PSCM counts separated by the PSC type
    all_muts_spectrum = defaultdict(int) # contains the complete base mutation counts

    # creates the WT sequence
    fasta.write('>id_WT \n')
    for AA in range(seq_length/3):
        base = ['A', 'C', 'G', 'T']
        codon = list(SC[randrange(0, 3)])
        
        # mutate randomly selected PSC; this will increase the likelihood of a PSCM occurring in the later steps
        random_bpos = randrange(0, 3)
        old_base = codon[random_bpos]
        base.remove(old_base)
        new_base = random.choice(base)
        codon[random_bpos] = new_base
        
        # keep mutating the random SC until it's not another SC type
        codon_string = codon[0] + codon[1] + codon[2]
        while codon_string in SC:
            base.remove(new_base)
            new_base = random.choice(base)
            codon[random_bpos] = new_base
            codon_string = codon[0] + codon[1] + codon[2]
        fasta.write(codon_string)
        WT_seq += codon_string
    fasta.write('\n\n')

    # the following block of code causes a base mutation (num_of_seq_entries number of times per sequence)
    for id_seq in range(num_of_seq_entries): # for each new mutated sequence
        codon = list(WT_seq)
        for i in range(num_of_base_muts): # for a variable number of induced base mutations per sequence entry
            base = ['A', 'C', 'G', 'T']
            random_AA = randrange(1, (seq_length-1)/3) # any random AA excluding the last codon
            bpos = (random_AA -1)*3 + 1 # convert AA pos to leading bpos
            random_bpos = bpos + randrange(0, 3) # random bpos of the selected AA for there to be a base mutation
            old_base = codon[random_bpos - 1] # original base; random_bpos - 1 because list indices start at 0 so have to adjust for that
            old_codon = codon[bpos - 1] + codon[bpos] + codon[bpos + 1] # old original codon
            # implement the base change
            base.remove(old_base)
            new_base = random.choice(base) # mutated base
            codon[random_bpos -1] = new_base # save mutated base at the randomly selected bpos; overwriting codon base (old_base to new_base)
            new_codon = codon[bpos -1] + codon[bpos] + codon[bpos + 1] # new mutated codon
            if new_codon in ['TAG', 'TGA']: # if the mutated codon is a PSC, save the base mutation as a PSCM in the PSCM dictionary
                PSCM[new_codon][old_base + new_base] += 1 # record the PSCM in the PSCM dictionary
                PSCM[new_codon]['tot'] += 1
            all_muts_spectrum[old_base + new_base] += 1
            all_muts_spectrum['tot'] += 1
            old_base, new_base = "", ""
        fasta.write('>id_' + str(id_num) + '\n' + ''.join(codon) + '\n\n')
        id_num += 1
        codon_string = ''
    fasta.close()

    # print_spectrum(all_muts_spectrum) # prints out all_muts_spectrum 
    return PSCM

def freq_dict(fastas): # haven't called this function yet
    '''
    Loop through each key(id) and value(seq) in fastas.
    For each seq, iterate through each char c and keep count of the char in a 2D dictionary.
    Then divide each count by the total number of valid base chars {A, T, C, G} in each respective column.
    Put those decimal frequencies in a tuple as the value for each key (base position) in the dictionary.
    '''
    freq = defaultdict(Counter) 
    bpos = 1 
    L = float(len(fastas))
    for id, seq in fastas.iteritems():
        bpos = 1 
        for c in seq:
            for char in ['A', 'C', 'T', 'G', 'N', '-']:
                if (c == char): freq[bpos][c] += 1 
            bpos += 1 # iterate to next base position
    # divides the count by the column length for each base position, where {N, -} are excluded  
    for bpos, d in freq.iteritems(): 
        for base, count in sorted(d.iteritems()):
            length = L - freq[bpos]['N'] - freq[bpos]['-']
            if (freq[bpos][base] <= length):
                for b in ['A', 'C', 'T', 'G', 'N']:
                    if (base == b): freq[bpos][b] = count/length
    return freq

def get_apobec_muts(fastas, breakList, base_mut_option, threshold): 
    '''
    Counts the APOBEC muts causing SBMs and DBMs and places those counts into a 2D dictionary,
    with the following setup: dict[id = SBM or DBM][key = A3G or A3DFH] = (value = number of occurrences).
    This is counting each instance of apobec in every input sequence so there may be more than 1 apobec mut at
    a codon column and get_muts() would thus get a negative value since get_muts() would only get one mut at that column.
    breaklist is the list of codon positions that should be skipped when considering apobec muts
    because get_muts() skipped those positions due to none of the bases meeting the designated threshold.
    '''
    apobec = defaultdict(lambda : defaultdict(Counter))
    psc_pos = list()
    for id, seq in fastas.iteritems():
        pos = 0
        bpos = 0
        consCodon = ""
        for i in xrange(0, len(seq), 3):
            codon = seq[i:i+3]
            pos += 1
            if (pos in breakList or "-" in codon): continue
            if (codon=="TAA" or codon=="TGA" or codon=="TAG"):
                if (pos not in psc_pos and i+3!=len(seq)): # appends if 1) not a duplicate AA pos and 2) not the last stop codon
                    psc_pos.append(pos)
                    bpos = (pos-1)*3 + 1 # base position of starting base in the pSC
                    consCodon = get_cons_base(bpos, fastas, threshold) + get_cons_base(bpos+1, fastas, threshold) + get_cons_base(bpos+2, fastas, threshold) # consensus codon for the pSC AA position
                    if (consCodon== "TGG"): # check for typtophan
                        base4 = seq[(pos-1)*3 + 4 - 1] # post-codon base
                        '''
                        Date Stamp 2-7-18
                        What if base4 is a non-nucleotide char? Then the GA mutation would automatically be attributed to RT even though there's a chance it could be APOBEC-induced
                        Should we consider these GA mutations separately in a class of RT_or_APOBEC?
                        '''
                        if (codon=="TAG"):
                            if (base_mut_option == "TAG" or base_mut_option == "all"):
                                apobec["SBM"]["TAG"]["A3G"] += 1
                        elif (codon=="TGA"):
                            if (base_mut_option == "TAG" or base_mut_option == "all"):
                                if (base4=="A"): apobec["SBM"]["TGA"]["A3DFH"] += 1
                                if (base4=="G"): apobec["SBM"]["TGA"]["A3G"] += 1
                        else: # codon == "TAA" so it's a DBM
                            if (base_mut_option == "TAA" or base_mut_option == "all"): # base_mut_option must be checked to include TAA muts
                                apobec["DBM"]["TAA"]["A3G"] += 1 # automatic A3G mutation 
                                if (base4=="A"): apobec["DBM"]["TAA"]["A3DFH"] += 1
                                if (base4=="G"): apobec["DBM"]["TAA"]["A3G"] += 1
    return apobec
    
def get_codons(pstop_pos, fastas): # function has not been used so far in main()
    '''
    Get all codons from all isolates, corresponding to the columns containing premature STOPs.
    '''
    codons = {}
    for pos in pstop_pos:
        for id, seq in fastas.iteritems():
            codon = seq[(pos-1)*3:(pos*3)]
            if ("-" in codon): continue
            codons[pos] = codon
    return codons

def get_column_seqs(fastas):
    '''
    returns dictionary of column seqs, where key: base position and value: corresponding column seq
    '''
    column_seqs = {}
    matrix = np.array([list(seq) for seq in fastas.values()])
    for j in xrange(matrix.shape[1]):
        bpos = j+1
        column_seqs[bpos] = ''.join(matrix[:,j])
    return column_seqs

def get_cons_base(bpos, fastas, threshold): 
    '''
    Returns the consensus base for a specific bpos 
    '''
    columnSeq = ""
    consBase = ""
    freq = defaultdict(int)
    for id, seq in sorted(fastas.iteritems()): columnSeq += seq[bpos - 1] # bpos - 1 to give corrected char in seq
    for b in ['A', 'C', 'G', 'T']: freq[b] = columnSeq.count(b) # find frequency of each nucleotide base
    total = freq['A'] + freq['C'] + freq['T'] + freq['G'] # total frequency of nucleotide base chars in that column for the given bpos
    # find the nucleotide that appears at least 90% (or whatever the threshold parameter is) of the time out of the four nucleotide bases (i.e., excluding dashes and N chars)
    for b in ['A', 'C', 'T', 'G']: 
        x = ['A', 'C', 'T', 'G']
        x.remove(b)
        if (freq[b] != 0 and freq[b] >= threshold * total): # freq[b] != 0 check is for when the consensus char is non-nucleotide, i.e., N, -, ., *
            consBase = b
            break # found the consensus base so break out of for loop; no need to continue searching for it
    return consBase

def get_muts(psc_pos, psc_ids, column_seqs, fastas, threshold): # DEBUG; I drastically changed this function so make sure it still works by testing it on a simple test FASTA file
    '''
    Returns the SBMs in a 2D dict (depending on given parameter), in the following format: dict[PSC][baseMut] = count
    '''
    muts = defaultdict(Counter) # 2D dict containing the DBM counts; this is returned at the end of the function
    codon = defaultdict(Counter) # 2D dict for storing the maxBase and mutBase at a given base position
    f = defaultdict(int) # key/vaue pairs are the column count of each base and all bases; rename this var to be more descriptive
    column = "" # seq of a column of a given base positon
    maxBase = mutBase = "" 
    mut = "" # two char string for reprensenting the mutation
    check = 1 # boolean variable used for checking threshold
    x = 1 # keeps track of base position
    plist = defaultdict(dict) # checks to make sure the DBMs aren't counted more than once
    breakList = list() # list used in apobec() so that apobec muts at this codon pos are skipped like in get_apobec_muts()
    for pos in psc_pos:
        for id, seq in sorted(fastas.iteritems()):
            if id in psc_ids[pos]:
                pSCbpos = (pos-1)*3 + 1 # pSC AA position converted to base position
                column = ""
                base = "" 
                x = 0
                mutCount = 0 # boolean var to keep track of number of base muts in a codon
                for c in seq: # i could optimize this search by making seq a list of nucleotide base elements and going through the indexes for the amino acid of interest
                    x += 1
                    if (x < pSCbpos): continue # immediately skip the rest of the for loop to continue searching for the bpos of interest
                    if (x >= pSCbpos and x < pSCbpos + 3): 
                        column = column_seqs[x] 
                        if (c == "-" or c == "N"): continue 
                        mutBase = c
                        f['tot']= 0 # reset 
                        for b in ['A', 'C', 'T', 'G']: f[b] = column.count(b); f['tot'] += f[b]
                        check = -1
                        for b in ['A', 'C', 'T', 'G']: # determine which base freq is >= threshold; only consider it a mut if mutbase != maxbase
                            if (f[b] >= f['tot'] * threshold): f['M'] = f[b]/float(f['tot']); maxBase = b; check = 1
                        if (check == -1): breakList.append(pos); break 
                        codon[x]['maxBase'] = maxBase
                        codon[x]['mutBase'] = mutBase
                        codon[x]['check'] = check
                    if (x == pSCbpos + 3): # save the SBM in muts dictionary if it's a valid mutation
                        mutType = ''
                        mut = ''
                        # make sure all three base positions met the threshold (i.e. check == 1 for all); if not, skip the codon and don't save this SBM in the muts dictionary
                        if (codon[pSCbpos]['check'] == 1 and codon[pSCbpos + 1]['check'] == 1 and codon[pSCbpos + 2]['check'] == 1):
                            if ((codon[pSCbpos]['maxBase'] + codon[pSCbpos +1]['maxBase'] + codon[pSCbpos]['maxBase']) != 'TAA' ): 
                                for basepos in codon:
                                    # there's a better way of implementing the following without iterating through the entire codon; change later 
                                    if (basepos == pSCbpos or basepos == pSCbpos + 1 or basepos == pSCbpos +2): # check for specific codon
                                        # bug detected: get_muts is not differentiating between SBMs and DBMs so it's saving both types when we want it to only save SBMs
                                        if (codon[basepos]['maxBase'] != codon[basepos]['mutBase']):
                                            if (mutType == ''): # first base mutation detected; may be a SBM or the first SBM of a DBM
                                                mutType = 'SBM' # if another SBM is detected, then this is a DBM, so do not save it into the muts dictionary
                                                mut = codon[basepos]['maxBase'] + codon[basepos]['mutBase']
                                                continue
                                            if (mutType == 'SBM'):
                                                mutType = 'DBM' # this change will mean get_muts() will not save the base mut (which is necessarily the second SBM of the DBM detected here)
                                                # print 'DBM' # for debugging
                                if (mutType == 'SBM'):
                                    # save the mutation into the mut dictionary
                                    psc = seq[pSCbpos -1: pSCbpos +2]
                                    # psc = codon[pSCbpos]['mutBase'] + codon[pSCbpos +1]['mutBase'] + codon[pSCbpos]['mutBase']
                                    # mut = codon[basepos]['maxBase'] + codon[basepos]['mutBase'] # incorrect to put this here because the basepos may not correspond to the mut basepos
                                    muts[psc][mut] += 1
                    if (x > pSCbpos + 3): break # done with this amino acid position so go on to the next one
    '''
    for bm in muts:
        print 'muts', muts[bm]
    '''
    return muts, breakList

def get_psc_ids_dict(psc_pos, fastas, base_mut_option):
    '''
    Keeps track of which fastas sequences have premature stop codons and
    puts their id into a dictionary as the corresponding value to the key (pSC index).
    '''
    psc_ids = defaultdict(list) 
    for pos in psc_pos:
        for id, seq in fastas.iteritems():
            codon = seq[(pos-1)*3:(pos*3)]
            if (codon == "TAA"):
                if (base_mut_option == "TAA" or base_mut_option == "all"):
                    psc_ids[pos].append(id)
            if (codon == "TAG" or codon == "TGA"):
                if (base_mut_option == "TAG" or base_mut_option == "all"):
                    psc_ids[pos].append(id)
    return psc_ids


def get_pstop_pos(fastas, base_mut_option): 
    """
    Creates a list containing the AA positions (no repeat positions) of the premature stop codons.
    """
    psc_pos = list()
    for id, seq in fastas.iteritems():
        pos = 0
        bpos = 0
        consCodon= ""
        for i in xrange(0, len(seq), 3):
            codon = seq[i:i+3]
            pos += 1
            if ("-" in codon): continue
            if (codon=="TGA" or codon=="TAG"): 
                if (base_mut_option == 'TAG' or base_mut_option == "all") :
                    if (pos not in psc_pos and i+3!=len(seq)): # appends if 1) not a duplicate AA pos and 2) not the last stop codon
                        psc_pos.append(pos)
            if (codon=="TAA"):
                if (base_mut_option == "TAA" or base_mut_option == "all"):
                    if (pos not in psc_pos and i+3!=len(seq)): # appends if 1) not a duplicate AA pos and 2) not the last stop codon
                        psc_pos.append(pos)
    psc_pos.sort(key=int)  # sort list numerically
    return psc_pos

def get_sbm_spectrum(apobec, single_muts, apobec_switch, norm_switch, base_mut_option, cons_codon): 
    '''
    returns the sbm spectrum
    '''
    if (apobec_switch == "on"): # implement apobec filter
        single_muts['TAG']['GA'] -= apobec["SBM"]['TAG']["A3G"]
        single_muts['TGA']['GA'] -= (apobec["SBM"]['TGA']["A3G"] + apobec["SBM"]['TGA']["A3DFH"])
    sbm_spec = defaultdict(float)
    
    # print 'sbm_counts_pre_normalization', single_muts
    if (norm_switch == "on"): # new normalization of SBMs is implemented below
        sbm_spec = normalize_sbm_counts(cons_codon, single_muts, base_mut_option) # this function is outlined in pseudo-code in my notebook
        # print 'sbm_count_post_normalization', sbm_spec #debugging purposes
    elif (norm_switch == "off"):
        if (base_mut_option == "TAA" or base_mut_option == "all"):
            for base_mut in single_muts['TAA']: sbm_spec[base_mut] += single_muts['TAA'][base_mut]
        if (base_mut_option == "TAG" or base_mut_option == "all"):
           for psc in ['TAG', 'TGA']:
               for base_mut in single_muts[psc]: sbm_spec[base_mut] += single_muts[psc][base_mut]
        # print 'sbm_counts_normalization_off', sbm_spec
    total = 0
    for base_mut in sbm_spec: total += sbm_spec[base_mut]
    sbm_spec['tot'] = total
    return sbm_spec

def get_spectrum(muts): 
    spec = Counter()
    fail = 0
    maxBase, mutBase =  ['A', 'C', 'T', 'G'],  ['A', 'C', 'T', 'G']
    for m,props in muts.iteritems():
        # is m not needed in this for loop because it's not used; the m used in the for loop pertains to mutBase, not muts
        for b in maxBase:
            mutBase.remove(b)
            if (props['maxbase'] == b):
                for m in mutBase:
                    if (props['mutbase'] == m): spec[b + m] += 1; fail = 1
            mutBase = ['A', 'C', 'T', 'G']
        if (fail == 0): print "'maxbase' is not A,C,G or T. Exiting."; exit(1)
    spec['tot'] = sum(spec.values())
    return spec

def identify_pstop(file, threshold): 
    """
    This function has not been used thus far (12/04/17)
    Returns a list whose elements are the pSCs at different AA positions
    this function's purpose is meant for when there are more than one type of pSC at a AA position
    when that happens, we need to select the pSC that occurs most at that AA position to be added to the list
    """
    pSC = list()
    TAA = TAG = TGA = 0
    boolean = 0
    for line in file:
        line = line.rstrip()
        char = line.split()
        if (char[0][0]=="p" and boolean==0):
            boolean = 1
            continue
        if (char[0][0]=="p" and boolean==1):
            # threshold check of >= 90% in choosing the pSC for a AA position
            for codon in [TAA, TAG, TGA]:
                other = [TAA, TAG, TGA]
                other.remove(codon)
                if (codon >= (other[0]/threshold) and codon >= (other[1]/threshold)):
                    pSC.append(codon)
            TAA = TAG = TGA = 0 # reseting the vars
        for codon in [TAA, TAG, TGA]:
            if (line== codon):
                codon += 1
    # threshold check for the last/only AA position pSC
    for codon in [TAA, TAG, TGA]:
        other = [TAA, TAG, TGA]
        other.remove(codon)
        if (codon >= (other[0]/threshold) and codon >= (other[1]/threshold)): pSC.append(codon)
    return pSC

def normalize_sbm_counts(cons_codon, single_muts, base_mut_option): # new function added to version 12; go through and debug just in case 
    '''
    New normalization method that normalizes SBM counts via the following formula:
    normalized BM count = # BM counts / total possibilities of that type of BM occurring.
    The denominator is the frequency of corresponding consensus mutating codons (MCs) that the selected SBM type 
    could cause to mutate into a PSC (i.e. the frequency of consensus codons (one count per AA position) that could be mutated to a PSC). 
    '''
    SBM_TAA = defaultdict(list) # list dictionary for TAA mutating codons
    SBM_TAG = defaultdict(list) # list dictionary for TGA mutating codons
    norm_SBM = defaultdict(float) # normalized values of inputed SBM counts 
    # Manually fill in the dict: key = SBM type, value = list of corresponding MCs
    # TAA stop codon
    SBM_TAA['CA'] = ['TAC', 'TCA']
    SBM_TAA['TA'] = ['TAT', 'TTA']
    SBM_TAA['AT'] = ['AAA']
    SBM_TAA['CT'] = ['CAA'] 
    SBM_TAA['GT'] = ['GAA']
    # TAG and TGA stop codons
    SBM_TAG['CA'] = ['TCG', 'TGC']
    SBM_TAG['GA'] = ['TGG']
    SBM_TAG['TA'] = ['TTA', 'TTG']
    SBM_TAG['CG'] = ['TAC', 'TCA']
    SBM_TAG['TG'] = ['TAT', 'TTA']
    SBM_TAG['AT'] = ['AAG', 'AGA']
    SBM_TAG['CT'] = ['CAG', 'CGA'] 
    SBM_TAG['GT'] = ['GAG', 'GGA']

    # 1. Iterate through each SBM (key)
    # 2. Sum up the counts of different MCs a SBM corresponds to
    # 3. Divide the corresponding SBM count (from single_muts) by the total count of MCs
    total = 0 # tallies up the counts of each MC for all SBMs that induce TAA, TAG & TGA, or all three stop codons

    if (base_mut_option == 'TAA'):
        for base_mut in sorted(SBM_TAA): # iterates through each dict key
            total = 0
            for mut_codon in SBM_TAA[base_mut]: # iterates through each dict value (list) 
                total += cons_codon[mut_codon]
            if (total > 0): norm_SBM[base_mut] = float (single_muts['TAA'][base_mut]) / (total)
            # the following print statement prints for each SBM type the corresponding MC, count, and norm values
            # print 'SBM:', base_mut, '\tMCs:', total, '\t\tcount:', single_muts['TAA'][base_mut], '\tnorm', norm_SBM[base_mut]
            # MCs is the frequency of the corresponding mutating codons for a given SBM, count is the orignal count for a given SBM, norm is the normalized value of the count
        
    elif (base_mut_option == 'TAG'):
        for base_mut in sorted(SBM_TAG): 
            total = 0
            for mut_codon in SBM_TAG[base_mut]:
                total += cons_codon[mut_codon]
                if (base_mut == 'GA'):
                    total += cons_codon[mut_codon] # double the count because this mutating codon applies to TAG and TGA
            if (total > 0):
                norm_SBM[base_mut] = float (single_muts['TAG'][base_mut] + single_muts['TGA'][base_mut]) / (total)
            # print 'SBM:', base_mut, '\tMCs:', total, '\t\tcount:', single_muts['TAG'][base_mut] + single_muts['TGA'][base_mut], '\tnorm', norm_SBM[base_mut]
               
    elif (base_mut_option == 'all'):
        for base_mut in sorted(SBM_TAA): 
            total = 0
            for mut_codon in SBM_TAA[base_mut]:
                total += cons_codon[mut_codon]
            if (total > 0):
                norm_SBM[base_mut] = float (single_muts['TAA'][base_mut]) / (total)
        for base_mut in sorted(SBM_TAG): 
            total = 0
            for mut_codon in SBM_TAG[base_mut]: 
                total += cons_codon[mut_codon]
                if (base_mut == 'GA'):
                    total += cons_codon[mut_codon] 
            if (total > 0):
                norm_SBM[base_mut] += float (single_muts['TAG'][base_mut] + single_muts['TGA'][base_mut]) / (total)
            # print 'SBM:', base_mut, '\tMCs:', total, '\t\tcount:', single_muts['TAA'][base_mut] + single_muts['TAG'][base_mut] + single_muts['TGA'][base_mut], '\tnorm', norm_SBM[base_mut]      
    return norm_SBM

def parse_args(args):
    """
    Command-line arguments, some of which act as parameters for different functions 
    """
    parser = argparse.ArgumentParser(description=__doc__,
                                    epilog='count_pscm.py v{}'.format(__version__))
    
    parser.add_argument('input',
        nargs="?",
        type=argparse.FileType('r'),
        help="Input FASTA file containing multiple alignment; alignment must be in frame (codon).")
    parser.add_argument('-t','--mut_call_threshold',
        type=float,
        default=0.9,
        help="Consensus base threshold for calling mutation vs polymorphism.")
    parser.add_argument('-a','--all_muts',
        action='store_true',
        default=False,
        help="Count all mutations and not just PSCM.")
    parser.add_argument('-s','--print_spectrum',
        action='store_true',
        default=True,
        help="Output: print the calculated mutation spectrum.")
    parser.add_argument('-d','--print_depth',
        action='store_true',
        default=False,
        help="Output: print the average depth of the input.")
    parser.add_argument('-as', '--apobec_switch',
        type=str,
        default = "off",
        help="Turns off the apobec filter by typing '-as off'.")
    parser.add_argument('-ns', '--norm_switch',
        type=str,
        default = "off",
        help="Turns off the normalization setting by typing '-ns off'.")
    parser.add_argument('-m', '--mut_type',
        type=str,
        default = "single",
        help="Counts SBMs, DBMs or both by typing '-m single/double/both'.")
    parser.add_argument('-bm', '--base_mut_option',
        type=str,
        default = "TAG",
        help="Counts (1) TAA muts, (2) TAG/TGA muts or all muts by typing '-bm TAA, TAG, all'.")
    parser.add_argument('-sl', '--seq_length',
        type=int,
        default = "10000",
        help="Runs fasta_simulation() with n sequence length per entry by typing '-sl n', where n is a positive integer.")
    parser.add_argument('-se', '--num_of_seq_entries',
        type=int,
        default = "1000",
        help="Runs fasta_simulation() with n number of sequence entries by typing '-se n', where n is a positive integer.")
    parser.add_argument('-nm', '--num_of_base_muts',
        type=int,
        default = "1",
        help="Runs fasta_simulation() with n number of base muts per sequence entry by typing '-nm n', where n is a positive integer.")
    pargs = parser.parse_args()
    return pargs

def parse_fasta(fasta):
    """
    parses FASTA file into a dictionary (keys: sequence IDs and value: corresponding sequence)
    """
    fastas = {}
    id = "" # fasta id var
    seq = "" # fasta sequence var
    chars = ['A', 'C', 'G', 'T', 'N', 'a', 'c', 'g', 't', 'n', '-', '*', '.']
    for line in fasta:
        line = line.rstrip()
        if (len(line) == 0): continue 
        if line[0] == ">": # if first character of the line is a > sign
            if seq!="": # value is now the sequence string
                seq = capitalizeBases(seq) # optional: capitalize the bases in the sequence
                fastas[id] = seq 
            id = line[1:] # id is every character after > sign of that line or you can do id = line[0:] so you can just do print id in main()
            seq = ""
        else:
            for c in line:
                if (c not in chars): seq += '-' # converts any non-readable chars to '-'
                else: seq+= c
    seq = capitalizeBases(seq) # optional: capitalize the bases in the sequence
    fastas[id] = seq
    return fastas

def print_spectrum(spec):
    '''
    different method of printing out spectrum
    print ('AC      AG      AT      CA      CG      CT      GA      GC      GT      TA      TC      TG      tot') # prints the order for reference 
    print "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(spec['AC'],spec['AG'],spec['AT'],spec['CA'],spec['CG'],spec['CT'],spec['GA'],spec['GC'],spec['GT'],spec['TA'],spec['TC'],spec['TG'],spec['tot'])
    '''
    print ('AC : {} \nAG : {} \nAT : {} \nCA : {} \nCG : {} \nCT : {} \nGA : {} \nGC : {} \nGT : {} \nTA : {} \nTC : {} \nTG : {} \ntot: {} \n'.format(spec['AC'],spec['AG'],spec['AT'],spec['CA'],spec['CG'],spec['CT'],spec['GA'],spec['GC'],spec['GT'],spec['TA'],spec['TC'],spec['TG'],spec['tot']))
 
if __name__ == "__main__" :
    try:
        sys.exit(main(sys.argv))
    except EnvironmentError as (errno,strerr):
        sys.stderr.write("ERROR: " + strerr + "\n")
        sys.exit(errno)

