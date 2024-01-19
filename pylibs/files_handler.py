#!/usr/bin/env python3

import os
from pathlib import Path

import numpy as np

def build_directories(Par):

    gamma, delta, u_EE, u_PP, L, RotMax, TranMax, N_max = Par

    fold = '.'

    res_fold = '/res'
    fold = build_fold(fold + res_fold)

    geometry_fold = '/gamma%g_delta%g' %(gamma, delta)
    fold = build_fold(fold + geometry_fold)

    energy_fold = '/eEE%g_ePP%g' %(u_EE, u_PP)
    fold = build_fold(fold + energy_fold)

    L_fold = '/L%d' %L
    fold = build_fold(fold + L_fold)

    move_fold = '/RotMax%g_TranMax%g' %(RotMax, TranMax)
    fold = build_fold(fold + move_fold)

    Nmax_fold = '/Nmax%d' %N_max
    fold = build_fold(fold + Nmax_fold)
    
    return fold

def build_fold(fold):

    if not os.path.isdir(fold):
        Path(fold).mkdir(parents=False, exist_ok=False)

    return fold


def build_inner_directories(fold, T, mu, Nparal):

    infold = fold + '/T%.4f_mu%.3f' %(T, mu)
    _ = build_fold(infold)

    for j in range(Nparal):
        confs_fold = '/confs_num%d' %(j+1)
        _ = build_fold(infold + confs_fold)

    return infold


def single_input_file(Par, MC_settings, t, n, mu, fold, run_num):

    gamma, delta, sigma_p, a, n_patches, e_EE, e_EP, e_PP, L, RotMax, TranMax, N_max = Par

    MCsteps, print_every, save_every, Ndays, Nhours, Nminutes, restart, seed = MC_settings
    seed = seed[run_num-1]
    
    numeric_ens = 0
    numeric_move = 0
    
    log_file = fold + '/log_num%d.txt' %run_num
    energy_file = fold + '/energy_num%d.txt' %run_num
    density_file = fold + '/density_num%d.txt' %run_num
    acceptance_file = fold + '/acceptance_num%d.txt' %run_num

    confs_fold = fold + '/confs_num%d' %run_num
    init_cond_file = confs_fold + '/init_cond_num%d.txt' %run_num
    last_conf_file = confs_fold + '/last_num%d.txt' %run_num

    settings_file = fold + '/settings_num%d.txt' %run_num

    Ftmp = open('settings_template.txt', 'r')
    F = open(settings_file, 'w')

    for _ in range(4): update_line(Ftmp, F, '')

    update_line(Ftmp, F, MCsteps)
    update_line(Ftmp, F, print_every)
    update_line(Ftmp, F, save_every)

    update_line(Ftmp, F, '')

    update_line(Ftmp, F, numeric_move)
    update_line(Ftmp, F, numeric_ens)
    update_line(Ftmp, F, t)

    update_line(Ftmp, F, '')

    update_line(Ftmp, F, RotMax)
    update_line(Ftmp, F, TranMax)

    update_line(Ftmp, F, '')

    update_line(Ftmp, F, L)
    update_line(Ftmp, F, n)

    update_line(Ftmp, F, '')

    update_line(Ftmp, F, Ndays)
    update_line(Ftmp, F, Nhours)
    update_line(Ftmp, F, Nminutes)

    update_line(Ftmp, F, '')

    update_line(Ftmp, F, init_cond_file)

    for _ in range(5): update_line(Ftmp, F, '')

    update_line(Ftmp, F, mu)
    update_line(Ftmp, F, N_max)
    
    for _ in range(5): update_line(Ftmp, F, '')

    update_line(Ftmp, F, restart)
    update_line(Ftmp, F, seed)

    update_line(Ftmp, F, '')

    ### FILES

    update_line(Ftmp, F, log_file)
    update_line(Ftmp, F, energy_file)
    update_line(Ftmp, F, density_file)
    update_line(Ftmp, F, acceptance_file)
    update_line(Ftmp, F, confs_fold)
    update_line(Ftmp, F, last_conf_file)

    for _ in range(5): update_line(Ftmp, F, '')

    update_line(Ftmp, F, 2)
    update_line(Ftmp, F, delta)
    update_line(Ftmp, F, 0.5)
    update_line(Ftmp, F, delta)
    update_line(Ftmp, F, sigma_p)
    update_line(Ftmp, F, a)

    update_line(Ftmp, F, '')

    update_line(Ftmp, F, e_EE)
    update_line(Ftmp, F, e_EP)
    update_line(Ftmp, F, e_PP)

    Ftmp.close()
    F.close()

def update_line(Ftmp, F, var):

    l = next(Ftmp)[:-2]
    l += ' ' + str(var) + '\n'
    F.write(l)
    
    
def executers_generator(T, mu, Nparal, Ncores, fold, program):
    
    Ftmp = open('executers_template.sh', 'r')
    F = open('executer.sh', 'w')
    
    p = build_params(T, mu, Nparal)
        
    line_A = build_params_line(p)
    
    for _ in range(2): copy_line(Ftmp, F)
    
    F.write('max_tasks=%d\n\n' %Ncores)

    F.write(line_A)
    
    F.write('Nruns=%d\n\n' %len(p))
    
    F.write('path=%s \n' %fold)
    
    for _ in range(4): copy_line(Ftmp, F)
    
    line_B = '  running_tasks=`ps -C %s --no-headers | wc -l`' %program
    F.write(line_B+'\n')
    
    for _ in range(15): copy_line(Ftmp, F)
    
    line_C = '    ./IPPcode/%s $f &' %program
    F.write(line_C + '\n\n')
    F.write('  ' + line_B+'\n')

    for _ in range(5): copy_line(Ftmp, F)
            
    Ftmp.close()
    F.close()

        
def build_params(p1, p2, num):
    
    p = []
    for j in range(1, num+1):
        p.append((p1, p2, j))
        
    return p


def build_params_line(p):
    
    line = 'declare -a _params=('

    i=1
    for c in p:
        p1, p2, n = c
        s='\"%.4f %.3f %d\" ' %(p1, p2, n)
            
        line += s
        i += 1
        if i == 8:
            line = line[:-1]
            line += '\n\t\t    '
            i = 1


    line += ')'

    line = line[:-2] + ')\n\n'
    
    return line


def copy_line(Ftmp, F):
    l = next(Ftmp)
    F.write(l)
    
    
def load_results(fold, num_of_runs, equilibration):

    TS = {}
    for m in range(1, num_of_runs+1):
        
        t, e, n, acc = [], [], [], []

        # Energy
        f = fold + '/energy_num%d.txt' %m
        with open(f) as F:
            for line in F:
                l = line.split()
                t.append(float(l[0]))
                e.append(float(l[1]))


        # density
        f = fold + '/density_num%d.txt' %m
        with open(f) as F:
            for line in F:
                l = line.split()
                n.append(float(l[1]))
                
        # acceptance rate
        f = fold + '/acceptance_num%d.txt' %m
        with open(f) as F:
            for line in F:
                l = line.split()
                acc.append(float(l[1]))

        TS[m] = {'t' : np.array(t).astype(int), 'e' : e, 'n' : n, 'acc' : acc}
  
    t = np.array(t).astype(int)
    ind = np.where(t == equilibration)[0][0]

    data = {'e' : [], 'n' : [], 'acc' : []}

    for m in TS:
        data['e'].extend(TS[m]['e'][ind:])
        TS[m]['e'] = np.array(TS[m]['e'])
        data['n'].extend(TS[m]['n'][ind:])
        TS[m]['n'] = np.array(TS[m]['n'])

        data['acc'].extend(TS[m]['acc'][ind:])
        
        print('The average acceptance rate after equilibration of run %d is %.3f' 
              %(m, sum(TS[m]['acc'][ind:]) / len(TS[m]['acc'][ind:])))

    print('\nThe total acceptance rate after equilibration is %.3f' %(sum(data['acc']) / len(data['acc'])))
    
    data.pop('acc', None)
    for k in data:
        data[k] = np.array(data[k])

    I = np.where(data['n'] == 0)
    data['n'] = np.delete(data['n'], I)
    data['e'] = np.delete(data['e'], I)

    return {'TS' : TS, 'data' : data}, ind


