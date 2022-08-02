/*
 * Coordinate.h
 *
 *  Created on: Nov 5, 2021
 *      Author: 14105
 */

#ifndef COORDINATE_H_
#define COORDINATE_H_

typedef struct Coordinate{

	//base-60 representation
	double latDeg;
	double latMin;
	double latSec;
	double lonDeg;
	double lonMin;
	double lonSec;

	//decimal representation
	double decimalLat;
	double decimalLon;

} Coordinate;


#endif /* COORDINATE_H_ */