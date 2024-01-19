#!/usr/bin/env python3

import numpy as np

## Basic particles functions: build the datastructure representing the particle and 
## check that the parameters settings is a legitimate one

def build_particles(sigma_c, delta, gamma):

    cg = np.cos(gamma)

    a = delta*(2*sigma_c + delta) / (2*(sigma_c+delta-sigma_c*cg))
    sigma_p = sigma_c + delta - delta*(2*sigma_c + delta) / (2*(sigma_c+delta-sigma_c*cg))

    IPP1 = {
        'n_patches' : 2,
        'a' : a,
        'delta' : delta,            
        'sigma_c' : sigma_c,
        'sigma_p' : sigma_p,
        'r' : np.array([0, 0, 0]),
        'rp1' : np.array([0, 0, a]),
        'rp2' : np.array([0, 0, -a]) }

    IPP2 = {
        'n_patches' : 2,
        'a' : a,
        'delta' : delta,
        'sigma_c' : sigma_c,
        'sigma_p' : sigma_p,
        'r' : np.array([2.0*sigma_c, 0, 0]),
        'rp1' : np.array([0, 0, a]),
        'rp2' : np.array([0, 0, -a]) }


    info_about_IPP_setting(IPP1)
    
    return IPP1, IPP2, a, sigma_p
    
def info_about_IPP_setting(IPP):

    delta = IPP['delta']
                
    a = IPP['a']
    sigma_c = IPP['sigma_c']
    sigma_p = IPP['sigma_p']

    cos_g = (sigma_c**2 + a**2 - sigma_p**2) / (2*a*sigma_c)
    gamma = np.arccos(cos_g)
    print('Patch extension = %.4f (degrees):' %(gamma*180/np.pi))
    print('Patch eccentricity = %.4f' %a)
    print('Patch radius = %.4f' %sigma_p)
    
    
    print('\n --- Parameters consistency ---')
    
    if a <= sigma_c: print('Eccentricity is smaller than the hard-core radius.')
    else: print('ERROR! The eccentricity is larger than the hard-core radius.')
        
    if a+sigma_p > sigma_c: print('Patches extend outside the hard-core.')
    else: print('ERROR! Patches do not extend outside the hard-core.')

    
    print('\n --- Conditions for isolation of interactions ---')
    
    cond1 = sigma_p < sigma_c
    if cond1: print('The EE configuration isolates the center-center interaction.')
    else: print('PROBLEM! Patches are too large for the EE interaction to be isolated.')
    
    cond2 = (sigma_p + sigma_c + delta)**2 <= a**2 + 4*sigma_c**2
    if cond2: print('The EP configuration isolates the center-off center interaction')
    else: print('PROBLEM! Patches are too large for the EP interaction to be isolated')
    
    cond3 = 4*(sigma_p**2) <= a**2 + (2*sigma_c - a)**2
    if cond3: print('The PP configuration isolates the off center-off center interaction')
    else: print('PROBLEM! Patches are too large for the PP interaction to be isolated')


def reference_config(IPP1, IPP2, a, sigma_c, conf):
    
    IPP1['r'] = np.array([0, 0, 0])
    IPP1['rp1'] = np.array([0, 0, a])
    IPP1['rp2'] = np.array([0, 0, -a])
    
    IPP2['r'] = np.array([0, 2*sigma_c, 0])
    IPP2['rp1'] = IPP2['r'] + np.array([0, 0, a])
    IPP2['rp2'] = IPP2['r'] + np.array([0, 0, -a])
    
    T = np.pi/2
    v = np.array([1, 0, 0])

    if conf == 'EE':
        pass
    elif conf == 'EP':
        IPP2 = rotate_IPP_around_its_center(v, T, IPP2)
    elif conf == 'PP':
        IPP1 = rotate_IPP_around_its_center(v, -T, IPP1)
        IPP2 = rotate_IPP_around_its_center(v, T, IPP2)

    return IPP1, IPP2

# Energy calculation

def energy_given_config(IPP1, IPP2, e):
    O = overlap_vector(IPP1, IPP2)
    return np.dot(e, O)

def overlap_vector(IPP1, IPP2):
        
    R = IPP1['sigma_c'] + IPP1['delta']
    r = IPP1['sigma_p']
    
    rab = IPP1['r'] - IPP2['r']
    rab = np.sqrt(np.dot(rab, rab))
    if rab < 2*IPP1['sigma_c']:
        return np.array([1e12, 1e12, 1e12])
    
    # EE
    O_EE = bonding_volume(R, R, IPP1['r'], IPP2['r'])
    
    # EP
    O_EP = 0
    O_EP += bonding_volume(R, r, IPP1['r'], IPP2['rp1'])       # E_1, P_1
    O_EP += bonding_volume(R, r, IPP1['r'], IPP2['rp2'])       # E_1, P_2
    O_EP += bonding_volume(r, R, IPP1['rp1'], IPP2['r'])       # E_2, P_1
    O_EP += bonding_volume(r, R, IPP1['rp2'], IPP2['r'])       # E_2, P_2

    # PP
    O_PP = 0
    O_PP += bonding_volume(r, r, IPP1['rp1'], IPP2['rp1'])       # P_1, P_1
    O_PP += bonding_volume(r, r, IPP1['rp1'], IPP2['rp2'])       # P_1, P_2
    O_PP += bonding_volume(r, r, IPP1['rp2'], IPP2['rp1'])       # P_2, P_1
    O_PP += bonding_volume(r, r, IPP1['rp2'], IPP2['rp2'])       # P_2, P_2
    
    O = (4/3) * np.pi * (IPP1['sigma_c']**3)
    
    return np.array([O_EE, O_EP, O_PP]) / O


def bonding_volume(Ra, Rb, ra, rb):
    
    rab = ra - rb
    rab = np.sqrt(np.dot(rab, rab))
    
    Rm = min(Ra, Rb)
    
    rmax = Ra + Rb
    rmin = np.fabs(Ra - Rb)
    
    if rab >= rmax:
        return 0
    elif rab <= rmin:
        return (4/3)*np.pi*(Rm**3)
    elif rmin <= rab and rab <= rmax:
        
        f1 = 2*Ra + (Ra**2 - Rb**2 + rab**2) / (2*rab) 
        f2 = Ra - (Ra**2 - Rb**2 + rab**2) / (2*rab) 
        
        f3 = 2*Rb - (Ra**2 - Rb**2 - rab**2) / (2*rab)
        f4 = Rb + (Ra**2 - Rb**2 - rab**2) / (2*rab)
        
        f = f1*f2*f2 + f3*f4*f4
        
        return f * np.pi/3
        
    else:
        print('Something wierd is going on')
        return np.nan


# These functions switch representation between energy vectors 

def e_from_u(W, u):
    WI = np.linalg.inv(W)
    e = WI.dot(u)
    return e

def u_from_e(W, e):
    u = W.dot(e)
    return u


# This function computes the potential along radial and angular paths
        
def model_potential(IPP1, IPP2, e, num):

    Potential = {}

    a = IPP1['a']
    sigma_c = IPP1['sigma_c']

    dt = 2.0 * IPP1['delta'] / num

    ### Radial potential from EE configuration
    IPP1, IPP2 = reference_config(IPP1, IPP2, a, sigma_c, 'EE')
    x, y = radial_potential(IPP1, IPP2, e, dt, num)
    Potential['EE'] = [x, y]

    ### Radial potential from EP configuration
    IPP1, IPP2 = reference_config(IPP1, IPP2, a, sigma_c, 'EP')
    x, y = radial_potential(IPP1, IPP2, e, dt, num)
    Potential['EP'] = [x, y]

    ### Radial potential from PP configuration
    IPP1, IPP2 = reference_config(IPP1, IPP2, a, sigma_c, 'PP')
    x, y = radial_potential(IPP1, IPP2, e, dt, num)
    Potential['PP'] = [x, y]

    dt = np.pi / num
    v = np.array([1, 0, 0])

    ### Angular potential from EE configuration
    IPP1, IPP2 = reference_config(IPP1, IPP2, a, sigma_c, 'EE')
    x, y = angular_potential(IPP1, IPP2, e, v, dt, num)
    Potential['radial_EE'] = [x, y]

    ### Angular potential from PP configuration
    IPP1, IPP2 = reference_config(IPP1, IPP2, a, sigma_c, 'PP')
    x, y = angular_potential(IPP1, IPP2, e, v, dt, num)
    Potential['radial_PP'] = [x, y]
    
        
    return Potential

# This function computes the potential along the radial path

def radial_potential(IPP1, IPP2, e, dt, num):

    x, y = np.zeros(num), np.zeros(num)
    for i in range(num):
        E = energy_given_config(IPP1, IPP2, e)
        IPP2['r'][1] += dt
        IPP2['rp1'][1] += dt
        IPP2['rp2'][1] += dt
        
        x[i] = i*dt
        y[i] = E
    
    return x, y

# This function computes the potential along the angular path

def angular_potential(IPP1, IPP2, e, v, dt, num):

    x, y = np.zeros(num), np.zeros(num)    
    for i in range(num):
        E = energy_given_config(IPP1, IPP2, e)
        IPP2 = rotate_IPP_around_its_center(v, dt, IPP2)

        x[i] = i*dt
        y[i] = E

    return x, y


# Rotations

def rotate_IPP_around_its_center(v, dtheta, IPP):
    R = rotation_matrix(v, dtheta)

    dr1 = IPP['rp1'] - IPP['r']
    dr2 = IPP['rp2'] - IPP['r']

    IPP['rp1'] = IPP['r'] + R.dot(dr1)
    IPP['rp2'] = IPP['r'] + R.dot(dr2)

    return IPP


def rotation_matrix(v, theta):
    
    ct = np.cos(theta)
    omct = 1 - ct
    st = np.sin(theta)
    omst = 1 - st
    
    R = np.array((
        [ct + v[0]*v[0]*omct, v[0]*v[1]*omct - v[2]*st, v[0]*v[2]*omct + v[1]*st], 
        [v[0]*v[1]*omct + v[2]*st, ct + v[1]*v[1]*omct, v[1]*v[2]*omct - v[0]*st],
        [v[0]*v[2]*omct - v[1]*st, v[1]*v[2]*omct + v[0]*st, ct + v[2]*v[2]*omct] ))
    
    return R
