#!/usr/bin/env python3

import random
import numpy as np
from datetime import datetime
from scipy import optimize

import sys
sys.path.insert(0, '.')
from plot_settings import *

def running_time(MC, fold, Nparal):

    for m in range(1, Nparal+1):
        f = fold+'/log_num%d.txt'%m
    
        dt = []
        flag = 0
        with open(f) as F:
            for line in F:
                if line.startswith('START:'):
                    if flag == 1:
                        de = di + timedelta(days=0, hours=72)
                        dt_partial = de - di
                        dt.append(dt_partial)
                    l = line.split()
                    day = l[1]
                    time = l[2]
                    init = day + ' ' + time
                    di = datetime.strptime(init, '%Y-%m-%d %H:%M:%S')
                    flag = 1
                if line.startswith('STOP:'):
                    l = line.split()
                    day = l[1]
                    time = l[2]
                    end = day + ' ' + time
                    de = datetime.strptime(end, '%Y-%m-%d %H:%M:%S')

                    dt_partial = de - di
                    dt.append(dt_partial)

                    flag = 0

        s = sum([dt_val.total_seconds() for dt_val in dt])

        nd = s / (24*60*60)
        nh = (nd - int(nd)) * 24
        nm = (nh - int(nh)) * 60
        ns = (nm - int(nm)) * 60
        
        t = max(MC['TS'][m]['t'])
        rate = t/s
        
        print('RUN %d' %m)
        print('Runtime: %d - %d:%d:%d' %(nd, nh, nm, ns))
        print('Average number of MC steps per second: %.2f\n' %(rate))

        
        
def all_histos_from_one_conf(res):

    pdf = {}

    x, y, dx = lin_hist(res['e'], 20)
    pdf['e'] = [x, y, dx]

    x, y, dx = lin_hist(res['n'], 20)
    pdf['n'] = [x, y, dx]

    return pdf


def lin_hist(data, bins):

    xM, xm = max(data) + 1e-6, min(data) - 1e-6

    dx = (xM - xm) / bins

    x = np.zeros(bins)
    I = np.arange(bins)
    y = np.zeros(bins)

    ind = ((data - xm) / dx).astype(int)

    np.add.at(x, ind, data)
    np.add.at(y, ind, 1)

    ind_on = np.where(y>0)
    ind_off = np.delete(I, ind_on)
    x[ind_on] = x[ind_on] / y[ind_on]
    x[ind_off] = xm + (ind_off + 0.5) * dx

    y = y / (len(data)*dx)

    return x, y, dx



def single_run_plot(TS, pdfs, ind, V):

    fig, ax = plt.subplots(nrows=1, ncols=2)
    fig.set_size_inches(14, 5)

    m = ['o', 's', 'p', 'v']

    L = list(TS.keys())
    I = random.sample(L, min([len(L), 4]))

    for j, i in enumerate(I):
        t, N = TS[i]['t'], TS[i]['n']
        rho = N / V

        ax[0].plot(t, rho)

        x, y, dx = lin_hist(rho[ind:], 20)
        ax[1].plot(x, y, '-%s' %m[j])

    ax[0].axvline(t[ind], color='grey', linestyle='--')
    ax[0].set_xlabel('MC step')
    ax[1].set_xlabel('Particle density')

    ax[0].set_ylabel('Particle density')
    ax[1].set_ylabel('Probability')

    plt.tight_layout()
    plt.show()
    
    fig, ax = plt.subplots(nrows=1, ncols=2)
    fig.set_size_inches(14, 5)

    for j, i in enumerate(I):
        t, E, N = TS[i]['t'], TS[i]['e'], TS[i]['n']
        e = E / N

        ax[0].plot(t, e)

        x, y, dx = lin_hist(e[ind:], 20)
        ax[1].plot(x, y, '-%s' %m[j])

    ax[0].axvline(t[ind], color='grey', linestyle='--')
    ax[0].set_xlabel('MC step')
    ax[1].set_xlabel('System energy')

    ax[0].set_ylabel('System energy')
    ax[1].set_ylabel('Probability')

    plt.tight_layout()
    plt.show()
    
    
def target_tmp(to_fit, not_to_fit):

    mu_prime, beta_prime, s = to_fit
    N, E, beta, mu, b = not_to_fit

    c, m1, renorm, Mc, PM, Ising = chi_tmp(mu_prime, beta_prime, s, N, E, beta, mu, b)

    return 1e3*c


def chi_tmp(mu_prime, beta_prime, s, N, E, beta, mu, b):

    M = N + s*E
    XM = np.linspace(min(M), max(M), num=b, endpoint=1)

    Mm = min(M)
    MM = max(M)
    dM = (MM - Mm) / (b-1)
    ind = ((M - Mm) / dM).astype(int)
    PM = build_histo(N, E, b, ind, beta_prime, mu_prime, beta, mu)
    PM /= sum(PM)

    m1 = sum(XM*PM)
    m2 = sum(XM*XM*PM)
    renorm = np.sqrt(m2 - m1**2)

    Mc = (XM -m1)/renorm
    dMc = Mc[1] - Mc[0]
    PM /= dMc

    Ising =  Ising_scaled_tmp(Mc)

    c = np.linalg.norm(PM - Ising)

    return c, m1, renorm, Mc, PM, Ising

    
    
def Ising_scaled_tmp(xx):

    nnav = 0
    alpha = 1
    x=(xx-nnav)/alpha;
    bimodal= 0.19155+0.26293*(x*x)-0.17003*(x*x)**2+0.47086*(x*x)**3-0.52411*(x*x)**4+0.22116*(x*x)**5\
            -0.033064*(x*x)**6-0.0031301*(x*x)**7+0.0016974*(x*x)**8-0.00020888*(x*x)**9+0.0000088492*(x*x)**10
    bimodal=bimodal/alpha;

    return bimodal

    
    
def target(to_fit, not_to_fit):

    mu_prime, beta_prime, s = to_fit
    N, E, beta, mu, b = not_to_fit

    c, m1, renorm, Mc, PM, Ising = chi(mu_prime, beta_prime, s, N, E, beta, mu, b)

    return 1e3*c


def chi(mu_prime, beta_prime, s, N, E, beta, mu, b):

    M = N + s*E
    XM = np.linspace(min(M), max(M), num=b, endpoint=1)

    Mm = min(M)
    MM = max(M)
    dM = (MM - Mm) / (b-1)
    ind = ((M - Mm) / dM).astype(int)
    PM = build_histo(N, E, b, ind, beta_prime, mu_prime, beta, mu)
    PM /= sum(PM)

    m1 = sum(XM*PM)
    m2 = sum(XM*XM*PM)
    renorm = np.sqrt(m2 - m1**2)

    Mc = (XM -m1)/renorm
    dMc = Mc[1] - Mc[0]
    PM /= dMc

    Ising =  Ising_scaled(Mc)

    c = np.linalg.norm(PM - Ising)

    return c, m1, renorm, Mc, PM, Ising



def build_histo(N, E, b, ind, beta_prime, mu_prime, beta, mu):

    y = np.zeros(b)

    f1 = (beta_prime - beta) * E
    f2 = (beta_prime*mu_prime - beta*mu) * N
    f = f2 - f1
    fact = np.exp(f)

    np.add.at(y, ind, fact)

    return y

def single_reweighted_hist(mu_prime, beta_prime, X, N, E, beta, mu, b):

    XM, Xm = max(X), min(X)
    dX = (XM - Xm) / b
    ind = ((X - Xm) / dX).astype(int)

    y = build_histo(N, E, b+1, ind, beta_prime, mu_prime, beta, mu)

    Xvals = np.linspace(Xm, XM, endpoint=1, num=b+1)
    edges = np.linspace(Xm, XM+dX, endpoint=1, num=b+2)
    widths = np.diff(edges)

    PX = y / (sum(y)*widths)

    return Xvals, PX, widths


def Ising_scaled(x):

    c1 = 0.1582
    c2 = 0.7762
    x0 = 1.13203953

    f1 = ((x**2)/(x0**2) - 1)**2
    f2 = (c1*(x**2)/(x0**2) + c2)

    I = np.exp(-f1*f2)

    i = int(len(I)/2)

    dx = x[i] - x[i-1]
    I /= (sum(I)*dx)

    return I
