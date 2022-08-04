/*
 * Satellite.h
 *
 *  Created on: Nov 19, 2021
 *      Author: 14105
 */

#ifndef SATELLITE_H_
#define SATELLITE_H_


typedef struct Satellite {

	//Store data from each line of TLE:
	//extra element added to all char[] to allow for strcopy function NULL character

	//testing:
	char lineZero[80], lineOne[80], lineTwo[80];

	//Name:
	char satName[24];

	//Line 1:
	char satNumL1[6];
	char internationalDes[6];
	double yearDayFrac;
	double rAvgPrime;
	double rAvgDblPrime;
	double dragRadCoeff;
	int ephemType;
	int elemNumChecksum;

	//Line 2:
	char satNumL2[5];
	double inclination;
	double rightAscension;
	double eccentricity;
	char argOfPerigee[7];
	double avgAnomaly;
	double rAvg;
	int revNumAtEpochChecksum;

} Satellite;


#endif /* SATELLITE_H_ */
