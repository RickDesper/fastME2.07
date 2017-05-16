#ifndef OPTIONS_H_
#define OPTIONS_H_

//#include "interface_free.h"
#include "interface_utilities.h"

#ifndef BOLD
#define BOLD      "\033[00;01m"
#endif

#ifndef FLAT
#define FLAT      "\033[00;00m"
#endif

#ifndef LINE
#define LINE      "\033[00;04m"
#endif


Options *chooseSettings(int argc, char **argv);

void Usage();
void Set_Defaults_Input (Options *input);
Options *Get_Input(int argc, char **argv);
void Get_Input_Interactive (Options *input);
void Get_Input_CommandLine(Options *input, int argc, char **argv);

#endif /*OPTIONS_H_*/
