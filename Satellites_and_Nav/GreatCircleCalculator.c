/*

 =================================================================

 Project Name:   Lab04 -- "Great Circle Calculation"

 Author      :   Tom West (UID 117659399)

 Username    :   twest625

 Date        :   11/05/2021

 Filename    :   Lab04_GreatCircleCalculation.c

 ===================================================================

 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Coordinate.h"

//DEFINITIONS://
const int MAX = 30; // take no more than MAX coordinate inputs
const double R_EARTH_AVG = 6371.0;
const double TO_RADIANS = M_PI/180.0;
const double FROM_RADIANS = 180.0/M_PI;
//const double

//PROTOTYPES://
int validateCoordinate(FILE*, double, int);
double getDecimalCoordinate(double, double, double);
//void calcGreatCircle(FILE*, struct Coordinate coordinates[]);

int main() {

	////////////////////
	//OPEN InPUT FILE://
	////////////////////

	//open .txt file
	FILE *fin = fopen(
			"C://Users/14105/Desktop/UMD/FA21/ENAE380/Labs/Lab4/dummyData_GreatCircle.txt",
			"r");
	if (!fin) {
		printf("Error opening file");
		return EXIT_FAILURE;
	}

	/////////////////////
	//OPEN OutPUT FILE://
	/////////////////////

	FILE *fout = fopen("twest625_Lab04_GreatCircleCalculation.txt", "w+");
	fprintf(fout,
			"========================================================================\n");
	fprintf(fout,
			"Author        : Tom West\nUsername      : twest625\nUID           : 117659399\n");
	fprintf(fout, "Section       : 0101\nConfiguration : Win10, Eclipse2021\n");
	fprintf(fout, "Description : Takes double representation of coordinate,\n");
	fprintf(fout,
			"   checks for basic error conditions, then calculates distance\n");
	fprintf(fout, "   using Great Circle method.\n");
	fprintf(fout,
			"========================================================================\n\n");

	fprintf(fout, "\nFORMAT: \n    1. COORDINATES (\"L\"atitude, lon\"G\"itude)\n");
	fprintf(fout, "    2. DISTANCE (km) & BEARING (deg.)\n\n\n");


	////////////////////
	//READ InPUT FILE://
	////////////////////

	//create coordinate array
	struct Coordinate coordinateArray[MAX];
	int counter, index, error, eof;
	counter = 0;
	index = 0;
	error = 0;
	eof = 1;

	while (1 == 1) {

		if (feof(fin)) {
			fprintf(fout, "\n\nEOF reached.");
			break;
		}

		/*
		if (counter > 20) {
			printf("oops");
			break;
		}
		*/

		//increment counter:
		counter++;

		//store base-60 information about coordinates
		double degMinSecONE[6];
		double degMinSecTWO[6];

		//GET TWO COORDINATES:
		fscanf(fin, " %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf ",
				&degMinSecONE[0], &degMinSecONE[1], &degMinSecONE[2],
				&degMinSecONE[3], &degMinSecONE[4], &degMinSecONE[5],
				&degMinSecTWO[0], &degMinSecTWO[1], &degMinSecTWO[2],
				&degMinSecTWO[3], &degMinSecTWO[4], &degMinSecTWO[5]);

		//////////////////////////////
		//CHECK COORDINATE VALIDITY://
		//////////////////////////////

		for (int i = 0; i < 6; i++) {
			if (validateCoordinate(fout, degMinSecONE[i], i) == -1) {
				fprintf(fout, "(LOCATION: Line %00d )\n",
						counter);
				if (index > 0) {
					index--; //go back to previous valid index
				}
				error = -1;
				break;
			}
			if (validateCoordinate(fout, degMinSecTWO[i], i) == -1) {
				fprintf(fout, "(LOCATION: Line %00d )\n",
						counter);
				if (index > 0) {
					index--; //go back to previous valid index
				}
				error = -1;
				break;
			}
		}

		if (error == -1){
			continue;
		}

		//////////////////////////////
		//PROCESS VALID COORDINATES://
		//////////////////////////////

		if (error == 0) {


			double latDeg1, latMin1, latSec1, lonDeg1, lonMin1, lonSec1,
					latDeg2, latMin2, latSec2, lonDeg2, lonMin2, lonSec2;

			//assign temporary values:
			latDeg1 = degMinSecONE[0];
			latMin1 = degMinSecONE[1];
			latSec1 = degMinSecONE[2];
			lonDeg1 = degMinSecONE[3];
			lonMin1 = degMinSecONE[4];
			lonSec1 = degMinSecONE[5];

			latDeg2 = degMinSecTWO[0];
			latMin2 = degMinSecTWO[1];
			latSec2 = degMinSecTWO[2];
			lonDeg2 = degMinSecTWO[3];
			lonMin2 = degMinSecTWO[4];
			lonSec2 = degMinSecTWO[5];

			//COORDINATE ONE ASSIGNMENT://

			//assign latitude ONE in base-60
			coordinateArray[index].latDeg = latDeg1;
			coordinateArray[index].latMin = latMin1;
			coordinateArray[index].latSec = latSec1;
			//assign longitude ONE in base-60
			coordinateArray[index].lonDeg = lonDeg1;
			coordinateArray[index].lonMin = lonMin1;
			coordinateArray[index].lonSec = lonSec1;
			//assign latitude ONE in decimal:
			coordinateArray[index].decimalLat = getDecimalCoordinate(
					coordinateArray[index].latDeg,
					coordinateArray[index].latMin,
					coordinateArray[index].latSec);
			//assign longitude ONE in decimal:
			coordinateArray[index].decimalLon = getDecimalCoordinate(
					coordinateArray[index].lonDeg,
					coordinateArray[index].lonMin,
					coordinateArray[index].lonSec);

			//COORDINATE TWO ASSIGNMENT://
			//shift to next Coordinate word
			index++;

			//assign latitude TWO in base-60
			coordinateArray[index].latDeg = latDeg2;
			coordinateArray[index].latMin = latMin2;
			coordinateArray[index].latSec = latSec2;
			//assign longitude TWO in base-60
			coordinateArray[index].lonDeg = lonDeg2;
			coordinateArray[index].lonMin = lonMin2;
			coordinateArray[index].lonSec = lonSec2;
			//assign latitude TWO in decimal:
			coordinateArray[index].decimalLat = getDecimalCoordinate(
					coordinateArray[index].latDeg,
					coordinateArray[index].latMin,
					coordinateArray[index].latSec);
			//assign longitude TWO in decimal:
			coordinateArray[index].decimalLon = getDecimalCoordinate(
					coordinateArray[index].lonDeg,
					coordinateArray[index].lonMin,
					coordinateArray[index].lonSec);

			/////////////////////////////////////
			//OUTPUT COORDINATE VALUES TO FILE://
			/////////////////////////////////////

			//output base-60 coordinates to file:
			fprintf(fout, "Coordinate %02d:\n", index);
			fprintf(fout,
					"     L: %.2f*%.2f\'%.2f\" , G: %.2f*%.2f\'%.2f\" \n",
					coordinateArray[index - 1].latDeg,
					coordinateArray[index - 1].latMin,
					coordinateArray[index - 1].latSec,
					coordinateArray[index - 1].lonDeg,
					coordinateArray[index - 1].lonMin,
					coordinateArray[index - 1].lonSec);
			fprintf(fout, "     (decimal):  %.6lf  ,  (decimal): %.6lf \n",
					coordinateArray[index - 1].decimalLat,
					coordinateArray[index - 1].decimalLon);
			fprintf(fout, "Coordinate %02d:\n", index + 1);
			fprintf(fout,
					"     L: %.2f*%.2f\'%.2f\" , G: %.2f*%.2f\'%.2f\" \n",
					coordinateArray[index].latDeg,
					coordinateArray[index].latMin,
					coordinateArray[index].latSec,
					coordinateArray[index].lonDeg,
					coordinateArray[index].lonMin,
					coordinateArray[index].lonSec);
			fprintf(fout, "     (decimal):  %.6lf  ,  (decimal): %.6lf \n",
					coordinateArray[index].decimalLat,
					coordinateArray[index].decimalLon);

			////////////////////////////////////
			//CALCULATE DISTANCES & BEARINGS:///
			////////////////////////////////////
			double distance, bearing, lat1, lat2, lon1, lon2;

			//convert from degrees to radians
			lat1 = coordinateArray[index - 1].decimalLat * TO_RADIANS;
			lon1 = coordinateArray[index - 1].decimalLon * TO_RADIANS;
			lat2 = coordinateArray[index].decimalLat * TO_RADIANS;
			lon2 = coordinateArray[index].decimalLon * TO_RADIANS;

			//calculate distance in nautical miles
			distance = R_EARTH_AVG * acos( sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2)*cos(lon1 - lon2) );
			bearing = atan2(sin(lon1-lon2)*cos(lat2), cos(lat1)*sin(lat2) - sin(lat1)*cos(lat2)*cos(lon1-lon2));
			bearing *= FROM_RADIANS;

			//OUTPUT RESULT://
			fprintf(fout, "\n     From point %02d to point %02d:\n          %lf km at %lf degrees\n\n\n",
					(index), (index + 1), distance, bearing);

			//MOVE TO NEXT INDEX:
			index++;
		}


	}

	/* Close files: */
	fclose(fin);
	fclose(fout);

	/* All done.*/
	printf("Run was successful.");

	return EXIT_SUCCESS;
}

int validateCoordinate(FILE* fout, double coordVal, int coordType) {

	switch (coordType) {
		case 0: //lat degrees
			if ((coordVal >= -90.0) && (coordVal <= 90.0)) {
				return 0;
			} else {
				fprintf(fout, "INVALID INPUT: Latitude degrees out of bounds. ");
				return -1;
			}
		case 1: //lat mins
			if ((coordVal >= 0.0) && (coordVal <= 59.0)) {
				return 0;
			} else {
				fprintf(fout, "INVALID INPUT: Latitude minutes out of bounds. ");
				return -1;
			}
		case 2: //lat seconds
			if ((coordVal >= 0.0) && (coordVal <= 59.99)) {
				return 0;
			} else {
				fprintf(fout, "INVALID INPUT: Latitude seconds out of bounds. ");
				return -1;
			}
		case 3: //lon degrees
			if ((coordVal >= -180.0) && (coordVal <= 180.0)) {
				return 0;
			} else {
				fprintf(fout, "INVALID INPUT: Longitude degrees out of bounds. ");
				return -1;
			}
		case 4: //lon minutes
			if ((coordVal >= 0.0) && (coordVal <= 59.0)) {
				return 0;
			} else {
				fprintf(fout, "INVALID INPUT: Longitude minutes out of bounds. ");
				return -1;
			}
		case 5: //lon seconds
			if ((coordVal >= 0.0) && (coordVal <= 59.99)) {
				return 0;
			} else {
				fprintf(fout, "INVALID INPUT: Longitude seconds out of bounds. ");
				return -1;
			}
		default:
			printf("DEFAULT REACHED\n");
			return -1;
	}

}

double getDecimalCoordinate(double deg, double min, double sec) {

	double decimalCoordinate = 0.0;

	if (sec > 0.01){
		decimalCoordinate += sec / 3600.0;
	}
	if (min > 0.01){
		decimalCoordinate += min / 60.0;
	}
	if (deg > 0.01){
		decimalCoordinate += deg;
	}

	return decimalCoordinate;
}


