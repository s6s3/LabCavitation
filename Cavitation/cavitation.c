#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "define.h"
#include "struct.h"
#include "weight.h"
#include "extern.h"
#include "cavitation.h"
#include "Buoyancy.h" 
#include "distance.h"
#include "file.h"
#include "gravity.h"
#include "neigh.h"
#include "outflow.h"

void
CAVITATION_calculateBubble(void) {

	if (timer.iTimeStep < timer.dt_initial * 10)return;

	if (timer.simulationTime < parameter.timeToUpdateAveragePressure)return;

	CAVITATION_calculateBubbleRising();

	CAVITATION_surfaceJudgement();

	CAVITATION_calculateBuoyancy();

}

double
CAVITATION_getBodyPressure() {
	double tmp = 0;
	int iParticle, number;
	number = 0;
	for (iParticle = 0; iParticle < particle.totalNumber; iParticle++) {
		if (particle.type[iParticle] != 0)continue;

		if (OUTFLOW_judgeBounds(particle.position[XDIM][iParticle], particle.position[YDIM][iParticle], particle.position[ZDIM][iParticle])
			< -10.0 * particle.averageDistance) {
			tmp += particle.pressure[iParticle];
			number++;
		}

	}

	return tmp / number;
}

double
CAVITATION_getSaturatedVaporPressure(double T) {
	return 610.78*pow(10, (7.5*T) / (T + 237.3));
}

double
CAVITATION_calculateCosine(double distanceIJ, int iParticle, int jParticle) {

	double cosine;

	cosine = (particle.position[YDIM][jParticle] - particle.position[YDIM][iParticle]) / fabs(distanceIJ);

	return(cosine);
}

double
CAVITATION_setInfluenceRadius(void) {

	double influenceRadius;

	if (NumberOfDimensions == 2) {

		influenceRadius = 4.1*particle.averageDistance;

	}
	else {

		influenceRadius = 3.1*particle.averageDistance;

	}

	return(influenceRadius);
}

double
CAVITATION_calculateVolume(double diameter) {

	double volume;

	//volume = (4.0 / 3.0)*M_PI*pow((diameter / 2.0), 3.0);

	volume = pow(diameter, NumberOfDimensions);

	return(volume);
}

void
CAVITATION_setBetaZero(void) {

	int iParticle;
	double influenceRadius;
	double distanceIJ;
	double cosine;
	double Beta;

	Beta = 0.0;

	influenceRadius = CAVITATION_setInfluenceRadius();

	for (iParticle = 0; iParticle<particle.totalNumber; iParticle++) {

		if (particle.type[iParticle] == parameter.wallType || particle.type[iParticle] == parameter.dummyWallType)continue;

		distanceIJ = DISTANCE_calculateDistanceBetweenParticles(particle.position, iParticle, parameter.numperOfParticleForCalculatingBeta);

		if (distanceIJ>influenceRadius || iParticle == parameter.numperOfParticleForCalculatingBeta)continue;

		cosine = (particle.position[YDIM][iParticle] - particle.position[YDIM][parameter.numperOfParticleForCalculatingBeta]) / distanceIJ;

		if (cosine <= 0.0)continue;

		Beta += (particle.position[YDIM][iParticle] - particle.position[YDIM][parameter.numperOfParticleForCalculatingBeta])*cosine*WEIGHT_calculateWeightFunction(distanceIJ, influenceRadius);

	}

	parameter.betaZeroOfParticles = Beta;
}

double
CAVITATION_calculateDifferenceRadius(int iParticle, double averagePressure) {

//	double volume;
	double deltaR, deltaP;
//	double voidrate;
	/*if (particle.pressure[iParticle] + physicalProperty.headPressure > physicalProperty.saturatedVaporPressure) {
		radius = -sqrt((2 * (averagePressure - particle.pressure[iParticle])) / (3 * physicalProperty.massDensity[0]));
	}
	else {
		radius = sqrt((2 * (particle.pressure[iParticle] - averagePressure)) / (3 * physicalProperty.massDensity[0]));

	}*/

	deltaP = physicalProperty.saturatedVaporPressure - particle.bucketPressure[iParticle];
	deltaR = copysign(1, deltaP) * sqrt(2 * fabs(deltaP) / (3 * physicalProperty.massDensity[0]));

	return(deltaR * timer.dt);
}

double
CAVITATION_calculateRisingVelocity(double diameter) {

	double risingVelocity;
	double viscosityCoefficientofFluid;
	double Re;

	viscosityCoefficientofFluid = physicalProperty.massDensity[0] * physicalProperty.kinematicViscosity;

	risingVelocity = (pow(diameter, 2)*(physicalProperty.massDensityOfBubble - physicalProperty.massDensity[0]))*physicalProperty.gravity[YDIM] / (18.0*viscosityCoefficientofFluid);

	Re = (risingVelocity*diameter) / physicalProperty.kinematicViscosity;

	if (2 <= Re) {

		risingVelocity = pow(((4.0 / 225.0)*pow(((physicalProperty.massDensityOfBubble - physicalProperty.massDensity[0])*physicalProperty.gravity[YDIM]), 2.0) / (physicalProperty.massDensity[0] * viscosityCoefficientofFluid)), 1.0 / 3.0)*diameter;
		Re = (risingVelocity*diameter) / physicalProperty.kinematicViscosity;

		if (500 <= Re) {

			risingVelocity = pow((4.0 / (3 * 0.44))*(physicalProperty.massDensityOfBubble - physicalProperty.massDensity[0])*physicalProperty.gravity[YDIM] * diameter / physicalProperty.massDensity[0], 1.0 / 2.0);
			Re = (risingVelocity*diameter) / physicalProperty.kinematicViscosity;
		}
	}

	return(risingVelocity);
}


void
CAVITATION_calculateBubbleRising(void) {

	int iParticle, jParticle;
	double risingVelocity;
	double distanceIJ;
	double cosine;
	double influenceRadius;
	double beta;
	double aIJ;
//	double bubbleDiameter;
	double averagePressure;

	int    iNeigh;
	int    *numberOfNeighborParticles;
	int    **neighborTable;

	NEIGH_selectNeighborTable(&numberOfNeighborParticles
		, &neighborTable
		, parameter.radiusOfLaplacianForPressure_ratio
		);
	averagePressure = CAVITATION_getBodyPressure();

	for (iParticle = 0; iParticle<particle.totalNumber; iParticle++) {

		if (particle.type[iParticle] == parameter.wallType || particle.type[iParticle] == parameter.dummyWallType || particle.type[iParticle] == parameter.typeNumberOfRigidParticle_forForcedMotion || particle.type[iParticle] == GHOST)continue;


		particle.cavitationBubblesRadius[iParticle] += CAVITATION_calculateDifferenceRadius(iParticle, averagePressure);

		if (particle.cavitationBubblesRadius[iParticle] < 0)particle.cavitationBubblesRadius[iParticle] = 0;

		risingVelocity = CAVITATION_calculateRisingVelocity(particle.cavitationBubblesRadius[iParticle]*2);

		if (risingVelocity <= 0.0)continue;

		beta = parameter.betaZeroOfParticles / risingVelocity;
		influenceRadius = CAVITATION_setInfluenceRadius();


		for (iNeigh = 0; iNeigh < numberOfNeighborParticles[iParticle]; iNeigh++) {

			jParticle = neighborTable[iParticle][iNeigh];

			if (particle.type[jParticle] == parameter.wallType || particle.type[jParticle] == parameter.dummyWallType || particle.type[jParticle] == parameter.typeNumberOfRigidParticle_forForcedMotion || particle.type[iParticle] == GHOST)continue;

			distanceIJ = DISTANCE_calculateDistanceBetweenParticles(particle.position, iParticle, jParticle);
			cosine = CAVITATION_calculateCosine(distanceIJ, iParticle, jParticle);

			if (cosine <= 0.0 || jParticle == iParticle)continue;

			aIJ = (timer.dt / beta)*WEIGHT_calculateWeightFunction(distanceIJ, influenceRadius)*cosine;
			particle.numberOfBubbles[jParticle] += aIJ*particle.numberOfBubbles_previous[iParticle];
			particle.numberOfBubbles[iParticle] -= aIJ*particle.numberOfBubbles_previous[iParticle];
			//particle.moleOfBubbles[jParticle] += aIJ*particle.moleOfBubbles_previous[iParticle];
			//particle.moleOfBubbles[iParticle] -= aIJ*particle.moleOfBubbles_previous[iParticle];
			particle.cavitationBubblesRadius[jParticle] += pow(aIJ, 1.0 / NumberOfDimensions)*particle.cavitationBubblesRadius[iParticle];
			particle.cavitationBubblesRadius[iParticle] -= pow(aIJ, 1.0 / NumberOfDimensions)*particle.cavitationBubblesRadius[iParticle];

		}
	}
}

void
CAVITATION_surfaceJudgement(void) {
	/*
	int iParticle;
	double timeConstant;
	double risingVelocity;
	double diameter;


	for (iParticle = 0; iParticle<particle.totalNumber; iParticle++) {

		if (particle.type[iParticle] == parameter.wallType || particle.type[iParticle] == parameter.dummyWallType)continue;

		if (particle.particleNumberDensity[iParticle]< 0.97*parameter.nZeroOfParticleNumberDensity) {

			diameter = CAVITATION_calculateDifferenceRadius(iParticle);
			risingVelocity = CAVITATION_calculateRisingVelocity(diameter);
			if (risingVelocity == 0)continue;
			timeConstant = (particle.averageDistance) / risingVelocity;

			particle.numberOfBubbles[iParticle] -= (timer.dt / timeConstant)*particle.numberOfBubbles[iParticle];
			particle.moleOfBubbles[iParticle] -= (timer.dt / timeConstant)*particle.moleOfBubbles[iParticle];

			if (particle.numberOfBubbles[iParticle] <= 0.0 || particle.moleOfBubbles[iParticle] <= 0.0) {

				particle.numberOfBubbles[iParticle] = 0.0;
				particle.moleOfBubbles[iParticle] = 0.0;

			}
		}
	}*/
}


void
CAVITATION_calculateBuoyancy(void) {

	int iParticle;
	int iDim;
	double voidrate;

	for (iParticle = 0; iParticle < particle.totalNumber; iParticle++) {

		if (particle.type[iParticle] == GHOST) continue;
		if (particle.type[iParticle] == parameter.dummyWallType) continue;
		if (particle.type[iParticle] == parameter.wallType) continue;
		/*
		voidrate = (particle.moleOfBubbles[iParticle] * physicalProperty.gasConstant*physicalProperty.temperature) / (((particle.pressure[iParticle]) + physicalProperty.headPressure)*(CAVITATION_calculateVolume(particle.averageDistance)));

		if (particle.pressure[iParticle] + physicalProperty.headPressure < physicalProperty.saturatedVaporPressure) {
			//todo
			voidrate = (particle.moleOfBubbles[iParticle] * physicalProperty.gasConstant*physicalProperty.temperature) / (physicalProperty.saturatedVaporPressure*(CAVITATION_calculateVolume(particle.averageDistance)));

		}*/
		voidrate = (NumberOfDimensions == 3)? 
			4 * M_PI * pow(particle.cavitationBubblesRadius[iParticle], 3) / (3 * pow(particle.averageDistance, 3)) : 
			4 * M_PI * pow(particle.cavitationBubblesRadius[iParticle], 2) / pow(particle.averageDistance, 2);

		for (iDim = 0; iDim < NumberOfDimensions; iDim++) {

			particle.velocity[iDim][iParticle] -= ((physicalProperty.massDensity[0] - physicalProperty.massDensityOfBubble) / physicalProperty.massDensity[0]) * voidrate * physicalProperty.gravity[iDim] * timer.dt;

		}

	}

}