#include "MC.h"
#include "time.h"
#include "defs.h"
#include "utils.h"
#include "basic_functions.h"

int stop = 0;
void gbl_terminate(int arg) {

	// The next time the signal is intercepted, the default handler will be invoked
	// in order to avoid making the program unkillable.

	signal(arg, SIG_DFL);
	fprintf(stderr, "# Caught SIGNAL %d; setting stop = 1\n", arg);
	stop = 1;
}


////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////


int main(int argc, char *argv[]) {


	// Here we handle a few SIG* signals: whenever we intercept one of these signals
	// the program will quit the main loop and die as gracefully as possible.
	signal(SIGTERM, gbl_terminate);
	signal(SIGABRT, gbl_terminate);
	signal(SIGINT, gbl_terminate);
	signal(SIGUSR2, gbl_terminate);

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

	// Initialise the data structures
	Output output_files;
	input_file input;
	System syst;

	output_files.log = stderr;
	if(argc == 1) output_exit(&output_files, "Usage is %s input\n", argv[0]);

	loadInputFile(&input, argv[1]);
	if(input.state == ERROR) exit(1);
	output_init(&input, &output_files);
	system_init(&input, &syst, &output_files);
	MC_init(&input, &syst, &output_files);

////////////////////////////////////////////////////////////////////////////////////

	// Handling duration of the simulation

	time_t t0,t1;
	int max_time;
	double dt, ndays, nhours, nminutes;

	time(&t0);
	output_current_time(&output_files, 0);

	ndays = output_files.ndays;
	nhours = output_files.nhours;
	nminutes = output_files.nminutes;
	max_time = ndays*24.0*60.0*60.0 + nhours*60.0*60.0 + nminutes*60.0;

////////////////////////////////////////////////////////////////////////////////////

// Compute the initial energy and check whether there are overlaps in the initial configuration

	int i;
	syst.energy = 0;
	for(i = 0; i < syst.N; i++)	{
		syst.energy += MC_energy(&syst, syst.particles + i);
		if(syst.overlap == 1) output_exit(&output_files, "Initial configuration contains an overlap, aborting\n");
	}
	syst.energy *= 0.5;

	// Get the number of steps to be run from the input file
	llint steps, curr_step;
	getInputLLInt(&input, "Steps", &steps, 1);

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

// Main loop


	char name[1024];
	for(curr_step = output_files.start_from; curr_step < steps && !stop; curr_step++) {

		// Print the output (energy, density, acceptance, etc.) every "print_every" steps
		if((curr_step % output_files.print_every) == 0) { output_print(&output_files, &syst, curr_step); }

		// Print the configuration every "save_every" steps
		if(curr_step > 0 && (curr_step % output_files.save_every) == 0) {
			sprintf(name, "%s/confs_all.txt", output_files.configuration_folder);
			output_save(&output_files, &syst, curr_step, name, 0);
		output_save(&output_files, &syst, curr_step, output_files.configuration_last, 1);
		}

		// // Perform a Monte Carlo sweep
		syst.do_ensemble(&syst, &output_files);

		// Stop if needed
		time(&t1);
		dt = difftime(t1, t0);
		if (dt > max_time) break;
	}

	output_current_time(&output_files, 1);

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

// Print the last configuration and the last line of the output

	output_save(&output_files, &syst, curr_step, output_files.configuration_last, 1);
	if (curr_step % output_files.save_every == 0) {
		sprintf(name, "%s/confs_all.txt", output_files.configuration_folder);
		output_save(&output_files, &syst, curr_step, name, 0);
	}
	if (curr_step == steps) {
		output_print(&output_files, &syst, curr_step);
		output_log_msg(&output_files, "END OF MC SIMULATION\nENDED\n");
	}

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////


	// Cleanup
	system_free(&syst);
	output_free(&output_files);
	cleanInputFile(&input);

	return 0;
}
