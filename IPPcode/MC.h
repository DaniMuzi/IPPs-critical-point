#ifndef MC_H_
#define MC_H_

#include "defs.h"

typedef struct input_file input_file;


void MC_init(input_file *input, System *syst, Output *IO);                 // This exists in the MC.c file
void do_GC(System *syst, Output *output_files);                            // Originally this line was not present at all in this file. How?

double MC_energy(System *syst, PatchyParticle *p);
int MC_interact(System *syst, PatchyParticle *p, PatchyParticle *q, vector O);

int overlap_volume(System *syst, PatchyParticle *p, PatchyParticle *q, vector O);
double bonding_volume(System *syst, double Ra, double Rb, vector ra, vector rb);

void MC_move_rototranslate(System *syst, Output *output_files);
void MC_rototraslate_particle(System *syst, PatchyParticle *p, vector disp, matrix R);
void MC_rollback_particle(System *syst, PatchyParticle *p);
void MC_change_cell(System *syst, PatchyParticle *p);

void MC_add_remove(System *syst, Output *IO);


void _store_dof(PatchyParticle * p);
void _restore_dof(PatchyParticle * p);
void normalized_diff_vector(System *syst, vector a, vector b, vector c);

#endif /* MC_H_ */
