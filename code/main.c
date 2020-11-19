/***************************************************************************
 *
 *   File        : main.c
 *   Student Id  : <INSERT STUDENT ID HERE>
 *   Name        : <INSERT STUDENT NAME HERE>
 *
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <string.h>
#include "tasks.h"

int main(int argc, char *argv[]) {
	
	/* TODO: Parse Command Line Arguments
	DONOT explicitly set arguments to filenames */
	char* q1_file = NULL;
	char* q3_file = NULL;
	char* q4_file = argv[3];
	double xo=atof(argv[4]);
	char* q5_file = argv[5];

	/* TODO: Add timing for each task and output running time in ms */
    
	/* Question 1 */
	shockwave(q1_file);
	
	/* Question 3 */
	linsolve(q3_file);
	
	/* Question 4 */
	interp(q4_file,xo);
	
	/* Question 5 */
	advection(q5_file);
    
	return (EXIT_SUCCESS);
}
