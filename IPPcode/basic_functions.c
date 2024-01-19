#include <time.h>
#include <math.h>
#include <float.h>
#include <assert.h>
#include <stdlib.h>
#include <unistd.h>

#include "MC.h"
#include "defs.h"
#include "utils.h"
#include "basic_functions.h"



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//
//                           SYSTEM INITIALIZATION
//
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


void system_init(input_file *input, System *syst, Output *output_files) {

	double mu;
	int res, i, j;

	// MC parameters
	syst->ensemble = GC;
	syst->dynamics = RTMC;

	getInputDouble(input, "Disp_max", &syst->disp_max, 1);
	getInputDouble(input, "Theta_max", &syst->theta_max, 1);
	getInputDouble(input, "Temperature", &syst->T, 1);


	// IPP parameters
	syst->n_patches = 2;

	getInputDouble(input, "delta_c", &syst->delta_c, 1);
	getInputDouble(input, "delta_p", &syst->delta_p, 1);
	getInputDouble(input, "sigma_c", &syst->sigma_c, 1);
	getInputDouble(input, "sigma_p", &syst->sigma_p, 1);
	getInputDouble(input, "a", &syst->a, 1);

	getInputDouble(input, "e_EE", &syst->e[0], 1);
	getInputDouble(input, "e_EP", &syst->e[1], 1);
	getInputDouble(input, "e_PP", &syst->e[2], 1);


	// Initial condition
	char name[256];
	getInputString(input, "Initial_conditions_file", name, 1);

	if(getInputInt(input, "Seed", &syst->seed, 0) == KEY_NOT_FOUND) {
		syst->seed = time(NULL);
		output_log_msg(output_files, "Using seed %d\n", syst->seed);
	}
	srand48(syst->seed);

	FILE *conf = fopen(name, "r");
	if(conf == NULL) output_exit(output_files, "Initial_conditions_file '%s' is not readable\n", name);
	res = fscanf(conf, "%*d %d %lf %lf %lf\n", &syst->N, syst->box, syst->box + 1, syst->box + 2);
	if(res != 4) output_exit(output_files, "The initial configuration file '%s' is empty or its headers are malformed\n", name);
	syst->V = syst->box[0] * syst->box[1] * syst->box[2];

	// Ensemble related inputs pt. 1
	getInputInt(input, "GC_N_max", &syst->N_max, 1);
	getInputDouble(input, "Mu", &mu, 1);
	syst->N_min = 0;
	syst->z = exp(mu/syst->T);


	// IPP initialization
	syst->particles = malloc(syst->N_max * sizeof(PatchyParticle));
	syst->energy = 0;
	syst->overlap = 0;

	for(i = 0; i < syst->N_max; i++) {
		PatchyParticle *p = syst->particles + i;
		p->index = i;

		p->n_patches = syst->n_patches;
		p->patches = malloc(sizeof(vector) * p->n_patches);
		p->patches_old = malloc(sizeof(vector) * p->n_patches);
	}

	i = 0;
	vector p0, p1;
	char myline[256];
	char *s_res = fgets(myline, 256, conf);
	while(s_res != NULL) {

		PatchyParticle *p = syst->particles + i;

		sscanf(myline, "%lf %lf %lf\n", p0, p0 + 1, p0 + 2);
		set_vector(p->r, p0[0], p0[1], p0[2]);
		for (j=0; j<syst->n_patches; j++){
			res = fscanf(conf, "%lf %lf %lf\n", p1, p1 + 1, p1 + 2);
			set_vector(p->patches[j], p1[0], p1[1], p1[2]);
		}

		i++;

		s_res = fgets(myline, 256, conf);
	}
	fclose(conf);
	if(i != syst->N) output_exit(output_files, "Number of particles found in configuration (%d) is different from the value found in the header (%d)\n", i, syst->N);

	utils_reset_acceptance_counters(syst);


	// Cells
	syst->r_cut = 1. + 2*syst->delta_p;
	syst->sqr_rcut = SQR(syst->r_cut);
	cells_init(syst, output_files, syst->r_cut);
	cells_fill(syst);
}



void system_free(System *syst) {
	cells_free(syst->cells);

	int i;
	for(i = 0; i < syst->N_max; i++) {
		PatchyParticle *p = syst->particles + i;
		if(p->patches != NULL) {
			free(p->patches);
		}
	}

	free(syst->particles);

}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//
//                           OUTPUT MANAGEMENT
//
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////



void output_init(input_file *input, Output *output_files) {


	output_files->log = stderr;
	getInputInt(input, "Restart_step_counter", &output_files->restart_step_counter, 1);
	const char *mode = (output_files->restart_step_counter) ? "a" : "w";


	char name[256];
	if(getInputString(input, "Log_file", name, 0) == KEY_FOUND) {
		if(strcmp("none", name) != 0) {
			FILE *mylog = fopen(name, mode);
			if(mylog == NULL) output_exit(output_files, "Log file '%s' is not writable\n", name);
			output_files->log = mylog;
		}
	}


	// Load the initial step from the configuration file, if the user requested to not restart the step counter
	if(!output_files->restart_step_counter) { output_files->start_from = 0; }
	else {
		getInputString(input, "Initial_conditions_file", name, 1);
		FILE *conf = fopen(name, "r");
		int res = fscanf(conf, "%lld %*d %*f %*f %*f\n", &output_files->start_from);
		if(res != 1) output_exit(output_files, "Invalid initial configuration: the first value in the first row should be the time step of the configuration");
		fclose(conf);
	}


	getInputLLInt(input, "Print_every", &output_files->print_every, 1);
	getInputLLInt(input, "Save_every", &output_files->save_every, 1);


	getInputString(input, "Acceptance_file", name, 1);
	FILE *myacc = fopen(name, mode);
	if(myacc == NULL) output_exit(output_files, "Acceptance file '%s' is not writable\n", name);
	output_files->acc = myacc;

	getInputString(input, "Energy_file", name, 1);
	FILE *myener = fopen(name, mode);
	if(myener == NULL) output_exit(output_files, "Energy file '%s' is not writable\n", name);
	output_files->energy = myener;


	getInputString(input, "Density_file", name, 1);
	FILE *myden = fopen(name, mode);
	if(myden == NULL) output_exit(output_files, "Density file '%s' is not writable\n", name);
	output_files->density = myden;


	getInputString(input, "Configuration_folder", output_files->configuration_folder, 1);
	if(access(output_files->configuration_folder, W_OK) != 0) {
		output_exit(output_files, "Cannot create files in directory '%s': please make sure that the directory exists and it is writable\n", output_files->configuration_folder);
	}
	char name2[1024];
	sprintf(name2, "%s/confs_all.txt", output_files->configuration_folder);
	FILE *out = fopen(name, "w");
	if(out == NULL) output_exit(output_files, "File '%s' is not writable\n", name);
	fclose(out);

	sprintf(output_files->configuration_last, "last.rrr");
	getInputString(input, "Configuration_last", output_files->configuration_last, 0);

	getInputDouble(input, "ndays", &output_files->ndays, 1);
	getInputDouble(input, "nhours", &output_files->nhours, 1);
	getInputDouble(input, "nminutes", &output_files->nminutes, 1);

}




void output_print(Output *output_files, System *syst, llint step) {

	// Print the energy
	double E = (syst->N > 0) ? syst->energy : 0;
	fprintf(output_files->energy, "%lld %lf\n", step, E);
	fflush(output_files->energy);

	// Print the density, if we are simulating in a non-canonical ensemble
	fprintf(output_files->density, "%lld %d\n", step, syst->N );
	fflush(output_files->density);

	// Print acceptances for the different moves, according to the chosen dynamics
	fprintf(output_files->acc, "%lld", step);
	fprintf(output_files->acc, " %e", syst->accepted[ROTO_TRASL] / (double) syst->tries[ROTO_TRASL]);

	fprintf(output_files->acc, " %e", syst->accepted[ADD]/ (double) syst->tries[ADD]);
	fprintf(output_files->acc, " %e", syst->accepted[REMOVE]/ (double) syst->tries[REMOVE]);
	fprintf(output_files->acc, "\n");
	fflush(output_files->acc);

	utils_reset_acceptance_counters(syst);

	// Print the current configuration
	output_save(output_files, syst, step, output_files->configuration_last, 1);


}



void output_save(Output *output_files, System *syst, llint step, char *name, int new_file_flag) {

	FILE *out;
	if (new_file_flag == 0) {
		out = fopen(name, "a");
		if(out == NULL) output_exit(output_files, "File '%s' is not writable\n", name);
		fprintf(out, "##################################################################################\n");
	} else {
		out = fopen(name, "w");
		if(out == NULL) output_exit(output_files, "File '%s' is not writable\n", name);
	}
	fprintf(out, "%lld %d %.12lf %.12lf %.12lf\n", step, syst->N, syst->box[0], syst->box[1], syst->box[2]);

	int i, j;
	PatchyParticle *p = syst->particles;
	for(i = 0; i < syst->N; i++) {
		fprintf(out, "%.12lf %.12lf %.12lf\n", p->r[0], p->r[1], p->r[2]);
		for (j=0; j<syst->n_patches; j++) {
			fprintf(out, "%lf %lf %lf\n", p->patches[j][0], p->patches[j][1], p->patches[j][2]);
		}
		p++;
	}
	fclose(out);
}



void output_free(Output *output_files) {
	fclose(output_files->acc);
	if(output_files->log != stderr) fclose(output_files->log);
	if(output_files->energy != NULL) fclose(output_files->energy);
	if(output_files->density != NULL) fclose(output_files->density);
}



void output_exit(Output *output_files, char *format, ...) {

	va_list args;
	if(output_files->log != stderr) {
		va_start(args, format);
		vfprintf(output_files->log, format, args);
		va_end(args);
	}

	va_start(args, format);
	vfprintf(stderr, format, args);
	va_end(args);

	fflush(output_files->log);
	exit(1);
}


void output_log_msg(Output *output_files, char *format, ...) {

	va_list args;
	if(output_files->log != stderr) {
		va_start(args, format);
		vfprintf(output_files->log, format, args);
		va_end(args);
	}
	else{
		va_start(args, format);
		vfprintf(stderr, format, args);
		va_end(args);
	}

	fflush(output_files->log);
}


void output_exit_stderr(char *format, ...) {
	va_list args;
	va_start(args, format);
	vfprintf(stderr, format, args);
	va_end(args);

	exit(1);
}

void output_current_time(Output *output_files, int flag) {

	time_t t0 = time(NULL);
  struct tm tm0 = *localtime(&t0);
	if (flag == 0) output_log_msg(output_files, "START: %d-%02d-%02d %02d:%02d:%02d\n",tm0.tm_year+1900, tm0.tm_mon+1, tm0.tm_mday, tm0.tm_hour, tm0.tm_min, tm0.tm_sec);
	if (flag == 1) output_log_msg(output_files, "STOP: %d-%02d-%02d %02d:%02d:%02d\nSTOPPED\n",tm0.tm_year+1900, tm0.tm_mon+1, tm0.tm_mday, tm0.tm_hour, tm0.tm_min, tm0.tm_sec);
	fflush(output_files->log);
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//
//                           INPUT MANAGEMENT
//
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////



void loadInputFile(input_file *inp, const char *filename) {

	FILE *inp_file = fopen(filename, "r");
	if(inp_file == NULL) {
	fprintf(stderr, "Input file not found\n");
		inp->state = ERROR;
		return;
	}

	loadInput (inp, inp_file);
	fclose (inp_file);

	return;
}



void loadInput(input_file *inp, FILE * inp_file) {

	int eof, i;
	size_t alloc_size;
	char *delim, *option = NULL;
	char t_option[256];

	inp->keys = NULL;
	inp->values = NULL;
	inp->reads = NULL;
	inp->unread_keys = NULL;
	inp->N_unread_keys = 0;
	inp->state = UNPARSED;

	inp->N_alloc = inp->N_opts = eof = 0;
	while(!eof) {
		if(getline(&option, &alloc_size, inp_file) == -1) {
			eof = 1;
			free(option);
			continue;
		};
		if(strlen(option) == 0) continue;

		if(strlen(option) > 0 && option[strlen(option)-1] == '\n') option[strlen(option)-1] = '\0';
		getTrimmedString(option, t_option);

		if(strlen(t_option) > 0 && t_option[0] != '#') {
			delim = strchr(t_option, '=');
			if(delim == NULL) {
				fprintf(stderr, "WARNING: malformed line '%s'\n", option);
				free(option);
				option = NULL;
				continue;
			}
			*delim = '\0';

			// We could need a realloc
			if(inp->N_opts == inp->N_alloc) {
				inp->N_alloc += 10;
				inp->keys = (input_string *) realloc(inp->keys, inp->N_alloc * sizeof(input_string));
				inp->unread_keys = (input_string *) realloc(inp->unread_keys, inp->N_alloc * sizeof(input_string));
				inp->values = (input_string *) realloc(inp->values, inp->N_alloc * sizeof(input_string));
				inp->reads = (int *) realloc(inp->reads, inp->N_alloc * sizeof(int));
				for(i = inp->N_alloc-10; i < inp->N_alloc; i++) inp->reads[i] = 0;
			}

			// Split the option in key and value and trim them both
			getTrimmedString(t_option, inp->keys[inp->N_opts]);
			getTrimmedString(delim+1, inp->values[inp->N_opts]);

			inp->N_opts++;
		}

		free(option);
		option = NULL;
	}

	inp->state = PARSED;
}


////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////


void getTrimmedString(const char *src, char *dest) {
	int start = 0;
	int end = strlen(src)-1;
	if(end < 0) {
		dest[0] = '\0';
		return;
	}

	while(isspace(src[start])) start++;
	while(isspace(src[end])) end--;

	strncpy(dest, src + start, end+1);
	dest[end+1] = '\0';
}



int getInputKeyIndex(input_file *inp, const char *skey, int mandatory)  {
	int i;
	for(i = 0; i < inp->N_opts; i++)
		if(!strncmp(skey, inp->keys[i], OPT_MAX_LENGTH) && strlen(inp->values[i]) > 0) {
			inp->reads[i]++;
			return i;
		}

	if(mandatory) {
		fprintf(stderr, "Mandatory key '%s' not found, exiting\n", skey);
		exit(-1);
	}
	// else fprintf(stderr, "INFO: optional key '%s' not found\n", skey);

	return KEY_NOT_FOUND;
}

int getInputString(input_file *inp, const char *skey, char *dest, int mandatory) {
	int key = getInputKeyIndex(inp, skey, mandatory);
	if(key != KEY_NOT_FOUND) strncpy(dest, inp->values[key], OPT_MAX_LENGTH);
	else return KEY_NOT_FOUND;

	return KEY_FOUND;
}

int getInputInt(input_file *inp, const char *skey, int *dest, int mandatory) {
	int key = getInputKeyIndex(inp, skey, mandatory);
	if(key != KEY_NOT_FOUND) *dest = atoi(inp->values[key]);
	else return KEY_NOT_FOUND;

	return KEY_FOUND;
}

int getInputLLInt(input_file *inp, const char *skey, long long int *dest, int mandatory) {
	int key = getInputKeyIndex(inp, skey, mandatory);
	if(key != KEY_NOT_FOUND) *dest = atol(inp->values[key]);
	else return KEY_NOT_FOUND;

	return KEY_FOUND;
}

int getInputUInt(input_file *inp, const char *skey, unsigned int *dest, int mandatory) {
	int key = getInputKeyIndex(inp, skey, mandatory);
	if(key != KEY_NOT_FOUND) *dest = atoi(inp->values[key]);
	else return KEY_NOT_FOUND;

	return KEY_FOUND;
}

int getInputDouble(input_file *inp, const char *skey, double *dest, int mandatory) {
	int key = getInputKeyIndex(inp, skey, mandatory);
	if(key != KEY_NOT_FOUND) *dest = atof(inp->values[key]);
	else return KEY_NOT_FOUND;

	return KEY_FOUND;
}

int getInputFloat(input_file *inp, const char *skey, float *dest, int mandatory) {
	int key = getInputKeyIndex(inp, skey, mandatory);
	if(key != KEY_NOT_FOUND) *dest = atof(inp->values[key]);
	else return KEY_NOT_FOUND;

	return KEY_FOUND;
}

int getInputChar(input_file *inp, const char *skey, char *dest, int mandatory) {
	int key = getInputKeyIndex(inp, skey, mandatory);
	if(key != KEY_NOT_FOUND) *dest = inp->values[key][0];
	else return KEY_NOT_FOUND;

	return KEY_FOUND;
}

void setUnreadKeys(input_file *inp) {
	int i;

	inp->N_unread_keys = 0;
	for(i = 0; i < inp->N_opts; i++) {
		if(inp->reads[i] == 0) {
			sprintf(inp->unread_keys[inp->N_unread_keys], "%s", inp->keys[i]);
			inp->N_unread_keys++;
		}
	}
}

void cleanInputFile(input_file *inp) {
	free(inp->keys);
	free(inp->values);
	free(inp->reads);
	free(inp->unread_keys);
	inp->state = UNPARSED;
}
