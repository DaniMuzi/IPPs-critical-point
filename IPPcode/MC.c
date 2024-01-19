#include "MC.h"
#include "utils.h"
#include "basic_functions.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>


////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////


void MC_init(input_file *input, System *syst, Output *IO) {

	// Here we set the pointer to the function that will be used to make a Monte Carlo step
	// according to the ensemble specified in the input file
	switch(syst->ensemble) {
		case GC:
			syst->do_ensemble = &do_GC;
			break;
		default:
			output_exit(IO, "Ensemble %d not supported\n", syst->ensemble);
			break;

	}


	// Here we set the pointer to the function that will be used to perform the type of dynamics
	// specified in the input file
	switch(syst->dynamics) {
		case RTMC:
			syst->do_dynamics = &MC_move_rototranslate;
			break;
		default:
			output_exit(IO, "Dynamics %d not supported\n", syst->dynamics);
			break;
	}


}

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
// Ensembles

void do_GC(System *syst, Output *output_files) {
	int i;
	for(i = 0; i < syst->N_max; i++) {
		double R = drand48();
		if(R < 0.01) {
			MC_add_remove(syst, output_files);
		}
		else if(syst->N > 0) {
			syst->do_dynamics(syst, output_files);
		}
	}
}


////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
// Energy calculation

double MC_energy(System *syst, PatchyParticle *p) {

	vector O;
	double E = 0.;
	int j, k, l, val;
	int ind[3], loop_ind[3], loop_index;

	syst->overlap = 0;

	cells_fill_and_get_idx_from_particle(syst, p, ind);

	for(j = -1; j < 2; j++) {
		loop_ind[0] = (ind[0] + j + syst->cells->N_side[0]) % syst->cells->N_side[0];
		for(k = -1; k < 2; k++) {
			loop_ind[1] = (ind[1] + k + syst->cells->N_side[1]) % syst->cells->N_side[1];
			for(l = -1; l < 2; l++) {
				loop_ind[2] = (ind[2] + l + syst->cells->N_side[2]) % syst->cells->N_side[2];
				loop_index = (loop_ind[0] * syst->cells->N_side[1] + loop_ind[1]) * syst->cells->N_side[2] + loop_ind[2];

				PatchyParticle *q = syst->cells->heads[loop_index];
				while(q != NULL) {
					if(q->index != p->index) {
						val = MC_interact(syst, p, q, O);

						if(val == BOND) {
							E += SCALAR(syst->e, O);
						}
						else if(val == OVERLAP) {
							syst->overlap = 1;
							return 1e8;
						}
					}
					q = syst->cells->next[q->index];
				}
			}
		}
	}

	return E;
}



inline int MC_interact(System *syst, PatchyParticle *p, PatchyParticle *q, vector O) {
	return overlap_volume(syst, p, q, O);
}


int overlap_volume(System *syst, PatchyParticle *p, PatchyParticle *q, vector O){

	int i, j;
	double R, r, d;

	R = syst->sigma_c + syst->delta_c;
	r = syst->sigma_p;

	vector rab = {q->r[0] - p->r[0], q->r[1] - p->r[1], q->r[2] - p->r[2]};

	rab[0] -= syst->box[0] * rint(rab[0] / syst->box[0]);
	rab[1] -= syst->box[1] * rint(rab[1] / syst->box[1]);
	rab[2] -= syst->box[2] * rint(rab[2] / syst->box[2]);
	d = sqrt(SCALAR(rab, rab));

	if (d - 2.0*syst->sigma_c < -1e-12) return OVERLAP;
	if (2.0*R - d < -1e-12) return NO_BOND;                                              // This MUST be changed if the IPC constraint is OFF!!!


	double O_EE, O_EP, O_PP, norm;

	// EE
	O_EE = 0;
	O_EE += bonding_volume(syst, R, R, p->r, q->r);

	// EP
	O_EP = 0;
	for(i=0; i<syst->n_patches; i++) {
		O_EP += bonding_volume(syst, R, r, p->r, q->patches[i]);                           // E1, P_i
		O_EP += bonding_volume(syst, r, R, p->patches[i], q->r);                           // E2, P_i
	}


	// PP
	O_PP = 0;
	for (i=0; i<syst->n_patches; i++) {
		for (j=0; j<syst->n_patches; j++) {
			O_PP += bonding_volume(syst, r, r, p->patches[i], q->patches[j]);                 // P_i, P_j
		}
	}

	norm = (4.0/3.0) * M_PI * (syst->sigma_c*syst->sigma_c*syst->sigma_c);

	O[0] = O_EE / norm;
	O[1] = O_EP / norm;
	O[2] = O_PP / norm;

	return BOND;

}

double bonding_volume(System *syst, double Ra, double Rb, vector ra, vector rb) {

	double Rm, rmax, rmin, d, f;

	vector rab = {ra[0] - rb[0], ra[1] - rb[1], ra[2] - rb[2]};

	rab[0] -= syst->box[0] * rint(rab[0] / syst->box[0]);
	rab[1] -= syst->box[1] * rint(rab[1] / syst->box[1]);
	rab[2] -= syst->box[2] * rint(rab[2] / syst->box[2]);

	d = sqrt(SCALAR(rab, rab));

	Rm = (Ra <= Rb) ? Ra : Rb;
	rmax = Ra + Rb;
	rmin = fabs(Ra - Rb);

	if (d >= rmax) { f = 0; }
	else if (d <= rmin) { f = (4.0/3.0)*M_PI*Rm*Rm*Rm; }
	else {                                                                        	// THAT IS: else if (rmin <= d && d <= rmax){

		double f1, f2, f3 ,f4;

		f1 = 2*Ra + (SQR(Ra) - SQR(Rb) + SQR(d)) / (2*d);
		f2 = Ra - (SQR(Ra) - SQR(Rb) + SQR(d)) / (2*d);

		f3 = 2*Rb - (SQR(Ra) - SQR(Rb) - SQR(d)) / (2*d);
		f4 = Rb + (SQR(Ra) - SQR(Rb) - SQR(d)) / (2*d);

		f = (f1*f2*f2 + f3*f4*f4) * M_PI / 3.0;

	}

	return f;


}


////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
// MC moves

// RTMC
void MC_move_rototranslate(System *syst, Output *IO) {

	matrix R;
	double theta;
	vector disp, axis;


	PatchyParticle *p = syst->particles + (int) (drand48() * syst->N);
	int type = ROTO_TRASL;
	syst->tries[type]++;

	disp[0] = (drand48() - 0.5) * syst->disp_max;
	disp[1] = (drand48() - 0.5) * syst->disp_max;
	disp[2] = (drand48() - 0.5) * syst->disp_max;

	// new orientation
	random_vector_on_sphere(axis);
	// theta = drand48() * syst->theta_max;   THIS IS WRONG! The angle must come from a distribution ~sin^2 theta, wich we approximate with theta^2 for small theta.
	theta = biased_angle(syst->theta_max);
	get_rotation_matrix(axis, theta, R);

	double E0, E1;
	E0 = MC_energy(syst, p);
	MC_rototraslate_particle(syst, p, disp, R);
	E1 = MC_energy(syst, p);
	double deltaE = E1 - E0;

	if(!syst->overlap && (deltaE < 0. || drand48() < exp(-deltaE / syst->T))) {
		syst->energy += deltaE;
		syst->accepted[type]++;
	}
	else {
		MC_rollback_particle(syst, p);
		syst->overlap = 0;
	}

}



void MC_rototraslate_particle(System *syst, PatchyParticle *p, vector disp, matrix R) {

	int i;

	_store_dof(p);

	p->r[0] += disp[0];
	p->r[1] += disp[1];
	p->r[2] += disp[2];

	vector dp, dp_tmp;
	for(i=0; i<p->n_patches; i++) {
		normalized_diff_vector(syst, p->r_old, p->patches[i], dp_tmp);              // Patch to center relative vector
		MATRIX_VECTOR_MULTIPLICATION(R, dp_tmp, dp);                                // Rotation of relative vector
		p->patches[i][0] = p->r[0] + dp[0];                               					// New absolute position = (p->r_old + disp) + dp = new_center + rotated_relative_vector =
		p->patches[i][1] = p->r[1] + dp[1];                               					//                       = (p->r_old + dp) + disp = rotated_absolute_vector + rigid_translation
		p->patches[i][2] = p->r[2] + dp[2];
	}

	MC_change_cell(syst, p);
}



void MC_rollback_particle(System *syst, PatchyParticle *p) {

	_restore_dof(p);

	Cells *cells = syst->cells;

	// bring the particle back in the old cell
	if(p->cell != p->cell_old) {
		if(cells->heads[p->cell]->index == p->index) cells->heads[p->cell] = cells->next[p->index];
		else {
			PatchyParticle *q = cells->heads[p->cell];
			while(cells->next[q->index] != p) {
				q = cells->next[q->index];
			}
			cells->next[q->index] = cells->next[p->index];
		}

		PatchyParticle *old = cells->heads[p->cell_old];
		cells->heads[p->cell_old] = p;
		cells->next[p->index] = old;
		int c_old = p->cell;
		p->cell = p->cell_old;
		p->cell_old = c_old;
	}
}

void MC_change_cell(System *syst, PatchyParticle *p) {

	int ind[3];
	int cell_index = cells_fill_and_get_idx_from_particle(syst, p, ind);
	if(cell_index == p->cell) {
		p->cell_old = p->cell;
		return;
	}

	// Remove the particle from the old cell
	Cells *cells = syst->cells;
	PatchyParticle *previous = NULL;
	PatchyParticle *current = cells->heads[p->cell];
	assert(current != NULL);
	while(current->index != p->index) {
		previous = current;
		current = cells->next[current->index];
		assert(cells->next[previous->index]->index == current->index);
	}
	if(previous == NULL) cells->heads[p->cell] = cells->next[p->index];
	else cells->next[previous->index] = cells->next[p->index];

	// Add the particle to the new cell
	cells->next[p->index] = cells->heads[cell_index];
	cells->heads[cell_index] = p;
	p->cell_old = p->cell;
	p->cell = cell_index;
}


void MC_add_remove(System *syst, Output *IO) {

	int i, cell_index;
	int ind[3];
	double delta_E, acc;

	// Try to add a particle
	if(drand48() < 0.5) {

		if(syst->N == syst->N_max) {
			if(syst->ensemble == GC) output_exit(IO, "The system contains the maximum number of particles set in the input file (%d), increase GC_N_max\nNMAX\n", syst->N_max);
			return;
		}
		syst->tries[ADD]++;

		PatchyParticle *p = syst->particles + syst->N;
		p->index = syst->N;

		random_orientation(syst, p);

		p->r[0] = drand48() * syst->box[0];
		p->r[1] = drand48() * syst->box[1];
		p->r[2] = drand48() * syst->box[2];

		for(i=0; i<syst->n_patches; i++) {
			p->patches[i][0] += p->r[0];
			p->patches[i][1] += p->r[1];
			p->patches[i][2] += p->r[2];
		}


		delta_E = MC_energy(syst, p);
		acc = exp(-delta_E / syst->T) * syst->z * syst->V / (syst->N + 1.);

		if(!syst->overlap && drand48() < acc) {

			syst->energy += delta_E;

			// Add the particle to the new cell
			cell_index = cells_fill_and_get_idx_from_particle(syst, p, ind);
			syst->cells->next[p->index] = syst->cells->heads[cell_index];
			syst->cells->heads[cell_index] = p;
			p->cell = p->cell_old = cell_index;

			syst->N++;
			syst->accepted[ADD]++;
		}
		else { syst->overlap = 0; }
	}
	// Try to remove a particle
	else {

		if(syst->N == syst->N_min) return;
		syst->tries[REMOVE]++;

		PatchyParticle *p = syst->particles + (int) (drand48() * syst->N);
		delta_E = -MC_energy(syst, p);
		acc = exp(-delta_E / syst->T) * syst->N / (syst->V * syst->z);

		if(drand48() < acc) {
			syst->energy += delta_E;
			syst->N--;

			// Remove the particle from the old cell
			PatchyParticle *previous = NULL;
			PatchyParticle *current = syst->cells->heads[p->cell];
			assert(current != NULL);
			while(current->index != p->index) {
				previous = current;
				current = syst->cells->next[current->index];
				assert(syst->cells->next[previous->index]->index == current->index);
			}
			if(previous == NULL) syst->cells->heads[p->cell] = syst->cells->next[p->index];
			else syst->cells->next[previous->index] = syst->cells->next[p->index];


			if(p->index != syst->N) {
				PatchyParticle *q = syst->particles + syst->N;
				assert(q->index != p->index);

				// We have to change the last particle's identity
				// First we remove it from its cell
				previous = NULL;
				current = syst->cells->heads[q->cell];
				assert(current != NULL);
				while(current->index != q->index) {
					previous = current;
					current = syst->cells->next[current->index];
					assert(syst->cells->next[previous->index]->index == current->index);
				}
				if(previous == NULL) { syst->cells->heads[q->cell] = syst->cells->next[q->index]; }
				else { syst->cells->next[previous->index] = syst->cells->next[q->index]; }

				// Copy its type, position and patches onto p's memory position
				p->r[0] = q->r[0];
				p->r[1] = q->r[1];
				p->r[2] = q->r[2];

				for(i=0; i<syst->n_patches; i++) {
					p->patches[i][0] = q->patches[i][0];
					p->patches[i][1] = q->patches[i][1];
					p->patches[i][2] = q->patches[i][2];
				}

				// And finally add it back to its former cell
				p->cell = q->cell;
				syst->cells->next[p->index] = syst->cells->heads[p->cell];
				syst->cells->heads[p->cell] = p;
			}

			syst->accepted[REMOVE]++;
		}
	}

}



////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////


void _store_dof(PatchyParticle * p) {
	int i, j;
	for(i = 0; i < 3; i++) {
		p->r_old[i] = p->r[i];
		for (j=0; j<p->n_patches; j++) {
			p->patches_old[j][i] = p->patches[j][i];
		}
	}
}


void _restore_dof(PatchyParticle * p) {
	int i, j;
	for(i = 0; i < 3; i++) {
		p->r[i] = p->r_old[i];
		for (j=0; j<p->n_patches; j++) {
			p->patches[j][i] = p->patches_old[j][i];
		}
	}
}


void normalized_diff_vector(System *syst, vector a, vector b, vector c) {

	c[0] = b[0] - a[0];
	c[1] = b[1] - a[1];
	c[2] = b[2] - a[2];
	c[0] -= syst->box[0] * rint(c[0] / syst->box[0]);
	c[1] -= syst->box[1] * rint(c[1] / syst->box[1]);
	c[2] -= syst->box[2] * rint(c[2] / syst->box[2]);

}
