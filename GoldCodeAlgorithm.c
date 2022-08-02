/*
 * ============================================================================
 * Name        : Lab03_GoldCodeGeneration.c
 * Author      : Tom West
 * Username    : twest625
 * UID         : 117659399
 * Section     : 0101
 * Description : Gold Code generator algorithm
 * ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//Constants:
const int TRUE   = 1;
const int FALSE  = 0;
const int LENGTH = 5;

int goldCodeAlgorithm(FILE *f, int stageLoad[]);

int main(void) {

	//////////
	//INPUT://
	//////////

	//.txt file heading/project info:
	FILE *fout = fopen("twest625_Lab03_GoldCodeGeneration.txt", "wt");
	fprintf(fout, "========================================================================\n");
	fprintf(fout, "Author        : Tom West\nUsername      : twest625\nUID           : 117659399\n");
	fprintf(fout, "Section       : 0101\nConfiguration : Win10, Eclipse2021\n");
	fprintf(fout, "Description : Simulates one Gold Code cycle using array representation of\n");
	fprintf(fout, "   shift register, XOR conditional statement, for loop array\n");
	fprintf(fout, "   traversal, and various programming best-practice pretty print.\n");
	fprintf(fout, "   Output array is Big Endian.\n");
	fprintf(fout, "========================================================================\n\n");
	
	//Variables:
	int sequenceLength = ((int)pow(2, LENGTH)) - 1;
	int shiftRegister[ 5 ]; //simulates the C/A-Code generator's shift register
	//since bool is not a type, I will substitute 1/0 for true/false
	int caCode[ sequenceLength ];
	int preLoad[ 5 ] = {TRUE, FALSE, FALSE, FALSE, FALSE};
	int output;

	///////////////
	//PROCESSING://
	///////////////

	//initialize shift register with preload value [1,0,0,0,0]
	fprintf(fout, "\nPRELOADING [1,0,0,0,0] TO SHIFT REGISTER...\n");
	for (int i = 0; i < 5; i++){
		shiftRegister[i] = preLoad[i];
	}//end for, shiftRegister locked and loaded!

	//Preamble:
	fprintf(fout, "\nGENERATING C/A-CODE...\n\n");

	//generate Gold Code sequence
	for (int i = 0; i < sequenceLength; i++){

		//display shift register before algorithm operation
		fprintf(fout, "   STEP %02d: ", (i + 1));
		for (int j = 0; j < LENGTH; j++){
			if (j < (LENGTH - 1) ){
				fprintf(fout, "| %d ", shiftRegister[j]);
			} else {
				fprintf(fout, "| %d |   ", shiftRegister[j]);
			}
		}

		//call algorithm, generate output
		output = goldCodeAlgorithm(fout, &shiftRegister);
	
		//update sequence array w/ new C/A-code output:
		caCode[i] = output;

		//display c/a value to user
		fprintf(fout, "OUTPUT: %d\n", output);
	}

	///////////
	//OUTPUT://
	///////////

	//Display C/A-code to user:
	fprintf(fout, "\n\nGenerated C/A-code with [1,0,0,0,0] input:\n\n   ");
	for (int i = 0; i < sequenceLength; i++){
		fprintf(fout, "%d", caCode[i]);
	}
	
	return EXIT_SUCCESS;
}

int goldCodeAlgorithm(FILE *f, int stageLoad[]){

	int output, feedback;
	output = stageLoad[LENGTH - 1];

	if ( ((stageLoad[1] == TRUE) && (stageLoad[4] == FALSE)) || ((stageLoad[1] == FALSE) && (stageLoad[4] == TRUE)) ){
		feedback = TRUE;
	} else {
		feedback = FALSE;
	}//end XOR gate

	for (int i = 4; i > 0; i--){
		stageLoad[i] = stageLoad[i - 1];
	}//end shift right

	stageLoad[0] = feedback;

	return output;
}//end GenerateGoldCode()
