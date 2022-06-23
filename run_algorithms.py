import math
import subprocess
import pandas as pd
import numpy as np
import os
import argparse
from copy import deepcopy
import time
import shutil
from copy import copy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import colorConverter as cc
import matplotlib.ticker as mticker
from pathlib import Path
from sklearn.metrics import r2_score
from multiprocessing import Pool
from matplotlib import rcParams
rcParams['figure.dpi'] = 300

folder_name = 'plots'
if not os.path.exists(folder_name):
    os.makedirs(folder_name)

parser = argparse.ArgumentParser(description='parse parameters')
parser.add_argument('ALGORITHM_NAME', metavar='ALGORITHM_NAME', type=str,
                    help='ALGORITHM_NAME')
parser.add_argument('PROBLEM_NUMBER', metavar='PROBLEM_NUMBER', type=int,
                    help='PROBLEM_NUMBER')
parser.add_argument('FIRST_RUN', metavar='FIRST_RUN', type=int,
                    help='number of first optimization run')
parser.add_argument('N_RUNS', metavar='N_RUNS', type=int,
                    help='number of optimization runs')
args = parser.parse_args()

ALGORITHM_NAME = args.ALGORITHM_NAME
PROBLEM_NUMBER = int(args.PROBLEM_NUMBER)
FIRST_RUN = int(args.FIRST_RUN)
N_RUNS = int(args.N_RUNS)

PROBLEM_NAMES = {0: 'MaxCut_Sparse',
                 1: 'Concatenated_Deceptive_Trap_4_4',
                 2: 'Concatenated_Deceptive_Trap_6_6',
                 3: 'Concatenated_Deceptive_Trap_8_8',               
                 4: 'NK_3_1',
                 5: 'NK_4_1',
                 6: 'NK_5_1',
                 7: 'Concatenated_Deceptive_Trap_4_1',
                 8: 'Concatenated_Deceptive_Trap_6_1',
                 9: 'MAXSAT'
                 }

PROBLEMS_DIMS = {
    0: [25, 49, 100, 196, 400, 784, 1600],
    1: [60, 120, 240, 480, 960, 1920],
    2: [60, 120, 240, 480, 960, 1920], 
    3: [40, 120, 240, 480, 960, 1920], 
    4: [20, 40, 80, 160, 320, 640, 1280],
    5: [20, 40, 80, 160, 320, 640, 1280],
    6: [20, 40, 80, 160, 320, 640, 1280],
    7: [40, 80, 160, 320, 640, 1280], 
    8: [40, 80, 160, 320, 640, 1280],
    9: [20, 50, 100]}
    
def problem_name_to_index(name):
    if 'Trap' in name:
        return 1
    if 'NK' in name:
        return 2
    if 'MaxCut' in name:
        return 3
    if 'MAXSAT' in name:
        return 4
    


TIME_LIMIT = 3600 #1 hour

def get_trap_vtr(L, k, s):
    num_subfunctions = L // s
    vtr = num_subfunctions * k
    return vtr

def get_vtr(problem_name, L, instance_number):

    if 'NK' in problem_name:     
        k = int(problem_name.split('_')[-2])
        s = int(problem_name.split('_')[-1])    
        vtr = np.loadtxt('problem_data/%s/L%d/%d_vtr.txt' % ('n%d_s%d' %(k,s), L, instance_number)) - 1e-6  

    elif problem_name == 'Spinglass':     
        vtr = open('problem_data/SPIN/%d/%d_%d' % (L, L, instance_number), 'r').readlines()[1]
        vtr = -float(vtr) - 1e-6  

    elif 'Concatenated_Deceptive_Trap' in problem_name:
        k = int(problem_name.split('_')[-2])
        s = int(problem_name.split('_')[-1])
        vtr = get_trap_vtr(L, k, s)

    elif problem_name == 'MaxCut_Sparse':       
        vtr = np.loadtxt('problem_data/maxcut/set0b/n%.7di%.2d.bkv' % (L, instance_number))

    elif problem_name == 'MaxCut_Dense':       
        vtr = np.loadtxt('problem_data/maxcut/set0a/n%.7di%.2d.bkv' % (L, instance_number))

    elif problem_name =='MAXSAT':
        vtr = 0.0
    
    else:
        print ('problem not implemented!')
        exit(0)

    return vtr

def get_instance_name(problem_name, L, instance_number):

    if 'NK' in problem_name:     
        k = int(problem_name.split('_')[-2])
        s = int(problem_name.split('_')[-1])    
        return 'problem_data/%s/L%d/%d.txt' % ('n%d_s%d' %(k,s), L, instance_number)

    elif problem_name == 'Spinglass':     
        return 'problem_data/SPIN/%d/%d_%d' % (L, L, instance_number)

    elif problem_name == 'MaxCut_Sparse':       
        return 'problem_data/maxcut/set0b/n%.7di%.2d.txt' % (L, instance_number)

    elif problem_name == 'MaxCut_Dense':       
        return 'problem_data/maxcut/set0a/n%.7di%.2d.txt' % (L, instance_number)

    elif problem_name =='MAXSAT':
        return 'problem_data/SAT/uf%d/uf%d-0%d.cnf' % (L, L, instance_number)
    
    else:
        print ('problem not implemented!')
        exit(0)

    return vtr

def write_vtr(vtr, folder):
    f = open('%s/vtr.txt' % folder, 'w')
    f.write(str(vtr))
    print(vtr, '%s/vtr.txt' % folder)
    f.close()

def get_data(folder_name, first_run, last_run, check_linkage_discovery):
    files = []
    all_evals, all_times, all_elitists, all_times_opt = [], [], [], []
    
    for run in range(first_run, last_run):
        try:
            for file in os.listdir(folder_name+'/'+str(run)):
                files.append(folder_name+'/'+str(run)+'/'+str(file))
        except Exception as e:
            return [], [], [], [-1]
            
    files_elitists = [file for file in files if 'elitist.dat' in file]
    
    if check_linkage_discovery:
        files_times = [file for file in files if 'linkage_discovery.dat' in file]
    
    for file in files_elitists:
        #print file
        lines = open(file, 'r').readlines()
        line = lines[-1].replace('\n', '').split()
        n_evals = int(line[0])
        time = float(line[1])
        elitist = float(line[2])
        
        if not check_linkage_discovery:
            all_evals.append(n_evals)
            all_times.append(time)

        all_elitists.append(elitist)
    
    if check_linkage_discovery:
        for file in files_times:

            lines = open(file, 'r').readlines()
            line = lines[-1].replace('\n', '').split()
            
            n_evals = int(line[0])
            all_evals.append(n_evals)

            print (file, line)
            time_full = float(line[1])
            time_opt = float(line[2])            
            all_times.append(time)
            all_times_opt.append(time_opt)

    print(all_evals, all_times, all_times_opt, all_elitists)

    stat_times = np.median(all_times), np.sort(all_times)[min(len(all_times)-1, 1)], np.sort(all_times)[max(len(all_times)-2,0)]
    stat_evals = np.median(all_evals), np.sort(all_evals)[min(len(all_times)-1, 1)], np.sort(all_evals)[max(len(all_times)-2,0)]
    all_elitists = np.array(all_elitists)

    if check_linkage_discovery:
        stat_times_opt = np.median(all_times_opt), np.sort(all_times_opt)[min(len(all_times_opt)-1, 1)], np.sort(all_times_opt)[max(len(all_times_opt)-2,0)]    
        return stat_evals, stat_times, stat_times_opt, all_elitists
    else:
        return stat_evals, stat_times, stat_times, all_elitists

def check_data(folder_name, max_run):
    #print folder_name
    files, files_main, files_times = [], [], []
    for run in range(max_run+1):
        try:
            for file in os.listdir(folder_name+'/'+str(run)):
                files.append(folder_name+'/'+str(run)+'/'+str(file))
        except Exception as e:
            return 0

    files_main = [file for file in files if 'elitist' in file]
    files_times = [file for file in files if 'linkage_discovery.dat' in file]    
        
    if len(files_main) == 0 or len(files_times) == 0:
        return 0
    return 1
    

gomea_time = []

def run_algorithm_for_problem(args):
    problem_index, algorithm = args['problem_index'], args['algorithm']
    DIMS = PROBLEMS_DIMS[problem_index]
    PROBLEM_NAME = PROBLEM_NAMES[problem_index]

    folder_name = 'results/%s/%s' % (algorithm, PROBLEM_NAME)    
    os.makedirs(folder_name, exist_ok=True)

    if 'Trap' in PROBLEM_NAME or 'NK' in PROBLEM_NAME:
        k,s = PROBLEM_NAME.split('_')[-2], PROBLEM_NAME.split('_')[-1]

    failed = False

    print(folder_name)

    for i in range(len(DIMS)):

        if algorithm == 'GOMEA':
            folder_name = 'results/GOMEA/%s/%d' % (PROBLEM_NAME, DIMS[i])    

            if os.path.exists(folder_name):
                shutil.rmtree(folder_name)     
            os.makedirs(folder_name, exist_ok=True)
            
            solved = True
            for run in range(FIRST_RUN, N_RUNS):
                
                run_folder = folder_name + '/%d' % run
                if os.path.exists(run_folder):
                    shutil.rmtree(run_folder)    
                os.mkdir(run_folder)
                
                vtr = get_vtr(PROBLEM_NAMES[problem_index], DIMS[i], run+1)
                print(vtr, run_folder)
                write_vtr(vtr, run_folder)
                seed = (DIMS[i] * N_RUNS + run) % 10**9
                probing_call = ['./BinaryGOMEA', '-i', '-r', str(problem_name_to_index(PROBLEM_NAME)), str(DIMS[i]), str(100), run_folder, str(TIME_LIMIT * 1000), str(seed)]
                if 'Trap' in PROBLEM_NAME or 'NK' in PROBLEM_NAME:
                    probing_call += [k,s]
                if 'NK' in PROBLEM_NAME or 'MaxCut' in PROBLEM_NAME or 'MAXSAT' in PROBLEM_NAME:
                    probing_call += [get_instance_name(PROBLEM_NAMES[problem_index], DIMS[i], run+1)]    
                print (probing_call)
                subprocess.call(probing_call)
                
                stat_evals, _, _, elitists = get_data(folder_name, run, run+1, check_linkage_discovery=False)
                print (elitists, stat_evals)
                if elitists[-1] < vtr:
                    solved = False
                    gomea_failed = 1
                    break
                if elitists[-1] > vtr+1e-6:
                    print('WRONG VTR!', elitists[-1], vtr)
                    exit(0)

            print ('GOMEA failed', failed)
            if failed:
                break
        #####################################################################################################
        folder_name = 'results/GOMEA/%s/%d' % (PROBLEM_NAME, DIMS[i])    
        stat_evals, times, _, elitists = get_data(folder_name, FIRST_RUN, N_RUNS, check_linkage_discovery=False)
        print('GOMEA stats', stat_evals, times, elitists)
        GOMEA_time = times[0]

        if algorithm == 'ELL':
            folder_name = 'results/ELL/%s/%d' % (PROBLEM_NAME, DIMS[i])
            os.makedirs(folder_name, exist_ok=True)
            
            for run in range(FIRST_RUN, N_RUNS):

                run_folder = folder_name + '/%d' % run
                if os.path.exists(run_folder):
                    shutil.rmtree(run_folder)    
                os.mkdir(run_folder)
                
                vtr = get_vtr(PROBLEM_NAMES[problem_index], DIMS[i], run+1)
                print(vtr, run_folder)
                write_vtr(vtr, run_folder)
                seed = (DIMS[i] * N_RUNS + run) % 10**9


                probing_call = ['./LinkageDiscovery', str(problem_name_to_index(PROBLEM_NAME)), str(DIMS[i]), '-1', '1.0', run_folder, '0', str(TIME_LIMIT * 1000), str(GOMEA_time), str(seed)]
                if 'Trap' in PROBLEM_NAME or 'NK' in PROBLEM_NAME:
                    probing_call += [k,s]
                if 'NK' in PROBLEM_NAME or 'MaxCut' in PROBLEM_NAME or 'MAXSAT' in PROBLEM_NAME:
                    probing_call += [get_instance_name(PROBLEM_NAMES[problem_index], DIMS[i], run+1)]    
                print (probing_call)
                subprocess.call(probing_call)

                stat_evals, _, _, elitists = get_data(folder_name, run, run+1, check_linkage_discovery=True)
                print (elitists, vtr, stat_evals)
                if elitists[-1] < vtr:
                    failed = 1
                    break
                if elitists[-1] > vtr+1e-6:
                    print('WRONG VTR!', elitists[-1], vtr)
                    exit(0)

            if failed:
                break

            
        # #####################################################################################################
        if algorithm == 'PLL':    
            folder_name = 'results/PLL/%s/%d' % (PROBLEM_NAME, DIMS[i])
            os.makedirs(folder_name, exist_ok=True)
            
            for run in range(FIRST_RUN, N_RUNS):

                run_folder = folder_name + '/%d' % run
                if os.path.exists(run_folder):
                    shutil.rmtree(run_folder)    
                os.mkdir(run_folder)
                
                vtr = get_vtr(PROBLEM_NAMES[problem_index], DIMS[i], run+1)
                print(vtr, run_folder)
                write_vtr(vtr, run_folder)
                seed = (DIMS[i] * N_RUNS + run) % 10**9


                probing_call = ['./LinkageDiscovery', '-H', str(problem_name_to_index(PROBLEM_NAME)), str(DIMS[i]), '-1', '1.0', run_folder, '0', str(TIME_LIMIT * 1000), str(GOMEA_time), str(seed)]
                if 'Trap' in PROBLEM_NAME or 'NK' in PROBLEM_NAME:
                    probing_call += [k,s]
                if 'NK' in PROBLEM_NAME or 'MaxCut' in PROBLEM_NAME or 'MAXSAT' in PROBLEM_NAME:
                    probing_call += [get_instance_name(PROBLEM_NAMES[problem_index], DIMS[i], run+1)]   
                print(probing_call)
                subprocess.call(probing_call)

                if check_data(folder_name, run) == 0:
                    failed = 1
                    break

                stat_evals, _, _, elitists = get_data(folder_name, run, run+1, check_linkage_discovery=True)
                print (elitists, vtr, stat_evals)
                if elitists[-1] < vtr:
                    failed = 1
                    break
                if elitists[-1] > vtr+1e-6:
                    print('WRONG VTR!', elitists[-1], vtr)
                    exit(0)

            if failed:
                break

        # #########################################################################################
        if algorithm == 'LARSLL':
            folder_name = 'results/LARSLL/%s/%d' % (PROBLEM_NAME, DIMS[i])
            os.makedirs(folder_name, exist_ok=True)
            os.environ['OMP_NUM_THREADS'] = '8'

            for run in range(FIRST_RUN, N_RUNS):

                run_folder = folder_name + '/%d' % run
                if os.path.exists(run_folder):
                    shutil.rmtree(run_folder)    
                os.mkdir(run_folder)
                
                vtr = get_vtr(PROBLEM_NAMES[problem_index], DIMS[i], run+1)
                print(vtr, run_folder)
                write_vtr(vtr, run_folder)
                seed = (DIMS[i] * N_RUNS + run) % 10**9
                
                lasso_call = ['./Lasso', str(problem_name_to_index(PROBLEM_NAME)), str(DIMS[i]),  '-1', '0.99999', run_folder, str(TIME_LIMIT * 1000), str(GOMEA_time), str(seed)]
                if 'Trap' in PROBLEM_NAME or 'NK' in PROBLEM_NAME:
                    lasso_call += [k,s]
                if 'NK' in PROBLEM_NAME or 'MaxCut' in PROBLEM_NAME or 'MAXSAT' in PROBLEM_NAME:
                    lasso_call += [get_instance_name(PROBLEM_NAMES[problem_index], DIMS[i], run+1)]   
                print(lasso_call)
                subprocess.call(lasso_call)

                if check_data(folder_name, run) == 0:
                    failed = 1
                    break

                stat_evals, _, _, elitists = get_data(folder_name, run, run+1, check_linkage_discovery=True)
                print (elitists, vtr, stat_evals)
                if elitists[-1] < vtr:
                    failed = 1
                    break
                if elitists[-1] > vtr+1e-6:
                    print('WRONG VTR!', elitists[-1], vtr)
                    exit(0)

            if failed:
                break


run_algorithm_for_problem({'algorithm':ALGORITHM_NAME, 'problem_index':PROBLEM_NUMBER})
