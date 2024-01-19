#ifndef basic_fs_H_
#define basic_fs_H_


#include "defs.h"
#include "utils.h"
#include "basic_functions.h"

#include <time.h>
#include <math.h>
#include <float.h>
#include <assert.h>
#include <stdlib.h>
#include <unistd.h>

// #include "MC.h"



////////////////////////////////////////////////////////////////////////////////


void system_init(input_file *input, System *syst, Output *output_files);
void system_free(System *syst);

////////////////////////////////////////////////////////////////////////////////


void output_init(input_file *input, Output *output_files);

void output_print(Output *output_files, System *syst, llint step);
void output_save(Output *output_files, System *syst, llint step, char *name, int new_file_flag);
void output_free(Output *output_files);
void output_exit(Output *output_files, char *format, ...);
void output_log_msg(Output *output_files, char *format, ...);
void output_exit_stderr(char *format, ...);
void output_current_time(Output *output_files, int flag);

////////////////////////////////////////////////////////////////////////////////



void loadInputFile(input_file *inp, const char *filename);
void loadInput(input_file *inp, FILE * inp_file);

void getTrimmedString(const char *src, char *dest);
int getInputKeyIndex(input_file *inp, const char *skey, int mandatory);
int getInputString(input_file *inp, const char *skey, char *dest, int mandatory);
int getInputInt(input_file *inp, const char *skey, int *dest, int mandatory);
int getInputLLInt(input_file *inp, const char *skey, long long int *dest, int mandatory);
int getInputUInt(input_file *inp, const char *skey, unsigned int *dest, int mandatory);
int getInputDouble(input_file *inp, const char *skey, double *dest, int mandatory);
int getInputFloat(input_file *inp, const char *skey, float *dest, int mandatory);
int getInputChar(input_file *inp, const char *skey, char *dest, int mandatory);
void setUnreadKeys(input_file *inp);
void cleanInputFile(input_file *inp);





#endif /* basic_fs_H_ */
