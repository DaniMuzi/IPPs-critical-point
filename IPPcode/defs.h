#ifndef DEFS_H_
#define DEFS_H_

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

// SYSTEM

#define N_MOVES 3
#define ROTO_TRASL 0
#define ADD 1
#define REMOVE 2

#define GC 0

#define RTMC 0

#define OVERLAP -1
#define NO_BOND 0
#define BOND 1

// INPUT

#define UNPARSED 0
#define PARSED 1
#define ERROR 2
#define OPT_MAX_LENGTH 255
#define KEY_NOT_FOUND -1
#define KEY_FOUND 0

// ALL LIBS

#include "cells.h"

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>

#include <string.h>
#include <signal.h>
#include <assert.h>
#include <complex.h>

#include <ctype.h>


// SYSTEM


typedef double vector[3];
typedef double quaternion[4];
typedef double matrix[3][3];
typedef long long int llint;
typedef struct Output Output;

typedef struct PatchyParticle {
	vector r, r_old;
	int index;

	int n_patches;
	vector *patches, *patches_old;

	int cell, cell_old;
} PatchyParticle;

typedef struct System {
	int N, N_min, N_max;
	vector box;
	double T;
	double z;
	double V;
	double energy;
	int n_patches;

	int dynamics;
	int ensemble;
	void (*do_dynamics)(struct System *, Output *);
	void (*do_ensemble)(struct System *, Output *);

	int overlap;

	int tries[N_MOVES];
	int accepted[N_MOVES];

	double disp_max;
	double theta_max;

	double sigma_c, sigma_p, delta_c, delta_p, a;
	double r_cut;
	double sqr_rcut;
	vector e;

	Cells *cells;

	int seed;
	PatchyParticle *particles;
} System;


// OUTPUT


typedef struct input_file input_file;
typedef struct System System;
typedef long long int llint;

typedef struct Output {
	llint start_from;
	llint save_every;
	llint print_every;

	double ndays;
	double nhours;
	double nminutes;

	int restart_step_counter;

	char configuration_folder[256];
	char configuration_last[256];


	FILE *log;
	FILE *energy;
	FILE *density;
	FILE *acc;
} Output;


// INPUT



typedef char input_string[OPT_MAX_LENGTH];

typedef struct input_file {
	int state;
	int N_opts;
	int N_alloc;
	int N_unread_keys;
	input_string *keys;
	input_string *values;
	int *reads;
	input_string *unread_keys;
} input_file;


#endif /* DEFS_H_ */
