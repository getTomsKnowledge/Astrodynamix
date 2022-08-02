/*

 =================================================================

 Project Name:   Lab05 -- "Two-Line Elements"

 Author      :   Tom West (UID 117659399)

 Username    :   twest625

 Date        :   11/19/2021

 Filename    :   Lab05_TwoLineElements.c

 ===================================================================

 */

//PUNCH LIST:
	//input echoes
	//avg inclination angle
	//max inclination angle
	//min inclination angle


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Satellite.h"

const int PUNCH_CARD_COLS = 80;
const char inputFileName[65] = "C://Users/14105/Desktop/UMD/FA21/ENAE380/Labs/Lab5/TLE_Input.txt";
const char outputFileName[14] = "TLE_Output.txt";

int getTotalSats(FILE *fin);
void setSatArray(FILE *fin, FILE *fout, struct Satellite satArray[]);
void parseTLELines(struct Satellite satArray[], int);
void getInclinationStats(FILE *fout, struct Satellite satArray[], double satStats[], int);


int main(void) {

	//////////////////////////
	//						//
	//		PART ONE:		//
	//				  		//
	//						//
	//////////////////////////



	////////////////////
	//open input file://
	////////////////////

//	printf("Begin main().");

	FILE *fin = fopen(inputFileName, "r");
		if(!fin){
			printf("ERROR: Could not open file.");
			return EXIT_FAILURE;
		}


	/////////////////////
	//open output file://
	/////////////////////

	FILE *fout = fopen(outputFileName, "w+");
	fprintf(fout,
			"========================================================================\n");
	fprintf(fout,
			"Author        : Tom West\nUsername      : twest625\nUID           : 117659399\n");
	fprintf(fout, "Section       : 0101\nConfiguration : Win10, Eclipse2021\n");
	fprintf(fout, "Description : Reads GPS TLE data file (.txt) and converts\n");
	fprintf(fout,
			"   various quantities to orbital elements.  Specifically calculates\n");
	fprintf(fout, "   inclinations of each satellite and organizes max/min by name.\n");
	fprintf(fout, "   Also calculates average inclination of input set.  Echoes vital input.\n");
	fprintf(fout,
			"========================================================================\n\n");

	fprintf(fout, "\nPART ONE:\n");
	fprintf(fout, "\n     1. Obtain current GPS TLE data set from celestrak.com.\n");
	fprintf(fout, "     2. Store TLE data in .txt file.\n");


	//////////////////////////////
	//read data, set satellites://
	//////////////////////////////

	fprintf(fout, "     3. Read TLE data from .txt file.\n");

	//get satellite count:
	int satTotal = 0;
	satTotal = getTotalSats(fin);
	fprintf(fout, "\n\n\t\tTelemetry data received from %d GPS satellites:\n\n", satTotal);
	rewind(fin);

	//create Satellite struct array
	struct Satellite gpsArray[satTotal];

	//store TLE lines in each Satellite struct
	setSatArray(fin, fout, gpsArray);

	//initialize Satellite fields
	parseTLELines(gpsArray, satTotal);

	//print names to user:
	for(int i = 0; i < satTotal; i++){
		fprintf(fout, "\t\t\t%s", gpsArray[i].lineZero);
	}

	///////////////////////
	//					 //
	//     PART TWO:     //
	//					 //
	//					 //
	///////////////////////

	fprintf(fout, "\n\n\nPART TWO:\n\n");
	fprintf(fout, "     1. Process TLE inclination data to obtain average, min, and max values.\n\n");
	//get inclination information:
	double satStats[5];
	getInclinationStats(fout, gpsArray, satStats, satTotal);

	fprintf(fout, "\n\n\n     2. Print results to user.\n\n");
	fprintf(fout, "\t\t\t\tMin. Inclination: %lf deg.\n\t\t\t\tSatellite: %s\n\n",
			satStats[1], gpsArray[((int) satStats[0])].lineZero);
	fprintf(fout, "\t\t\t\tMax. Inclination: %lf deg.\n\t\t\t\tSatellite: %s\n\n",
			satStats[3], gpsArray[((int) satStats[2])].lineZero);
	fprintf(fout, "\t\t\t\tAvg. Inclination: %lf deg.", satStats[4]);


	//testing data:
	printf("\n\nINCLINATION: %lf   Mean Motion: %.9lf \n\n\n", gpsArray[6].inclination, gpsArray[2].rAvg);


	/* Close files: */
	fclose(fin);
	fclose(fout);

	/* All done.*/
	printf("\n\nRun was successful.");

	return EXIT_SUCCESS;
}


int getTotalSats(FILE *fin){

	int lineCount, satelliteCount;
	char str[80];

	lineCount = 0;

	/*traverse fin */
	while (1 == 1){

		if (feof(fin)) {
			//printf("\n\nEOF reached.");
			break;
		}

		fgets(str, PUNCH_CARD_COLS, fin); //read one line of the file
		lineCount++;//increment number of lines
	}

	/* divide total lines by 3 to get number of satellites */
	satelliteCount = (int) lineCount/3;

	return satelliteCount;
}

void setSatArray(FILE *fin, FILE *fout, struct Satellite satArray[]){

//	printf("Inside setSatArray().");

	char junkLine[80];
	int index, lineNum;
	index = 0;
	lineNum = 0;

	while(1){

		if ( feof(fin) ){
			//printf("\n\n\nEnd of file reached again.");
			break;
		}

		if (index > 29) {
			fgets(junkLine, PUNCH_CARD_COLS, fin);
		} else{
			if (lineNum == 0){
				fgets(satArray[index].lineZero, PUNCH_CARD_COLS, fin);
				lineNum++;
			}
			if (lineNum == 1){
				fgets(satArray[index].lineOne, PUNCH_CARD_COLS, fin);
				lineNum++;
			}
			if (lineNum == 2){
				fgets(satArray[index].lineTwo, PUNCH_CARD_COLS, fin);
				lineNum = 0;
			}
		}

		index++;
	}
}

void parseTLELines(struct Satellite satArray[], int numSats){

	int lineOneIndex, lineTwoIndex;
	lineOneIndex = 0;
	lineTwoIndex = 0;

	for (int i = 0; i < numSats; i++){

		//parse LINE ONE
		char * token1 = strtok(satArray[i].lineOne, " ");
		while( token1 != NULL){
			switch (lineTwoIndex){
			case 0: //this is the line number identifier of TLE, do nothing
				break;
			case 1: // assign line 1 of satellite #
				strncpy(satArray[i].satNumL1, token1, 6);
				break;
			case 2: // assign International Designation Number
				strncpy(satArray[i].internationalDes, token1, 6);
				break;
			case 3: // assign Epoch in Julian Year w/ Day Fraction
				satArray[i].yearDayFrac = atof(token1);
				break;
			case 4: // assign First Derivative of the motion
				satArray[i].rAvgPrime = atof(token1);
				break;
			case 5: // assign Second Derivative of the motion
				satArray[i].rAvgDblPrime = atof(token1);
				break;
			case 6: // assign Drag/Radiation Coefficient
				satArray[i].dragRadCoeff = atof(token1);
				break;
			case 7: //assign ephemeris type
				satArray[i].ephemType = atoi(token1);
				break;
			case 8: //assign element number, checksum
				satArray[i].elemNumChecksum = atoi(token1);
				break;
			default:
				break;
			}
			token1 = strtok(NULL, " ");
			lineOneIndex++;
		}

		//parse LINE TWO
		char * token2 = strtok(satArray[i].lineTwo, " ");
		while( token2 != NULL){
			switch (lineTwoIndex){
			case 0: //this is the line number identifier of TLE, do nothing
				break;
			case 1: // assign line 2 of satellite #
				strncpy(satArray[i].satNumL2, token2, 7);
				break;
			case 2: // assign inclination angle
				satArray[i].inclination = atof(token2);
				break;
			case 3: // assign right ascension
				satArray[i].rightAscension = atof(token2);
				break;
			case 4: // assign e value for elliptical orbit
				satArray[i].eccentricity = atof(token2);
				break;
			case 5: // assign argument of perigee of orbit
				strncpy(satArray[i].argOfPerigee, token2, 5);
				break;
			case 6: // assign average anomaly of elliptical orbit
				satArray[i].avgAnomaly = atof(token2);
				break;
			case 7: //assign mean motion of satellite
				satArray[i].rAvg = atof(token2);
				break;
			case 8:
				satArray[i].revNumAtEpochChecksum = atoi(token2);
				break;
			default:
				break;
			}

			token2 = strtok(NULL, " ");
			lineTwoIndex++;
		}

		lineOneIndex = 0;
		lineTwoIndex = 0;
	}
}

void getInclinationStats(FILE *fout, struct Satellite satArray[], double statsArray[], int numSats){

	int minID, maxID;
	double min, max, sum, tempVal;

	tempVal = 0;
	min = satArray[0].inclination;
	max = min;

	for (int i = 0; i < numSats; i++) {

		tempVal = satArray[i].inclination;
		fprintf(fout, "\tInclination: %lf deg.\t\tSatellite: %s", tempVal, satArray[i].lineZero);

		sum += tempVal;//update running total

		if (min > tempVal) { //check for minimum
			min = tempVal;
			minID = i;
		}
		if (max < tempVal) { //check for maximum
			max = tempVal;
			maxID = i;
		}
	}

	statsArray[0] = (double) minID;
	statsArray[1] = min;
	statsArray[2] = (double) maxID;
	statsArray[3] = max;
	statsArray[4] = sum / numSats;

}

/*
 * https://www.celestrak.com/NORAD/documentation/tle-fmt.php
 * https://celestrak.com/publications/AIAA/2008-6770/AIAA-2008-6770.pdf
 *
 */
