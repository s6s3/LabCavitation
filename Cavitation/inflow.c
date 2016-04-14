#include <stdio.h>
#include <math.h>
#include <float.h>

#include "define.h"
#include "extern.h"
#include "struct.h"
#include "neigh.h"
#include "init.h"
#include "stack.h"

double INFLOW_getSignedDistanceFromInflowBounds(int iParticle) {
	int iDim;
	double norm;
	double normal_squared, func;
	normal_squared = 0.0;
	for (iDim = 0; iDim < NumberOfDimensions; iDim++) {
		normal_squared += parameter.inflowVelocity[iDim] * parameter.inflowVelocity[iDim];
	}
	norm = sqrt(normal_squared);
	for (iDim = 0; iDim < NumberOfDimensions; iDim++) {
		func += parameter.inflowVelocity[iDim] * (particle.position[iDim][iParticle] - parameter.inflowPosition[iDim]);
	}
	return func / norm;
}

int INFLOW_getTypeOfInflowParticle(int iParticle) {
	if (parameter.flagOfInflowBoundaryCondition == OFF)return -1;
	if (particle.isInflow[iParticle] == FALSE)return -1;
	if (INFLOW_getSignedDistanceFromInflowBounds(iParticle) < -particle.averageDistance)return parameter.dummyWallType;
	if (INFLOW_getSignedDistanceFromInflowBounds(iParticle) < particle.averageDistance)return parameter.wallType;
	return 0;

}

int INFLOW_setInflowLevel() {
	int iParticle, iDim;
	double tmp, max;
	max = -1000000000;
	int no = -1;

	if (parameter.flagOfAutoSettingOfInflowLevel == OFF)return -1;

	for (iParticle = 0; iParticle < particle.totalNumber; iParticle++) {
		if (particle.isInflow[iParticle] == FALSE)continue;
		tmp = INFLOW_getSignedDistanceFromInflowBounds(iParticle);
		if (tmp > max) {
			max = tmp;
			no = iParticle;
		}
	}
	if (no == -1)return -1;
	for (iDim = 0; iDim < 3; iDim++) {
		parameter.inflowPosition[iDim] = particle.position[iDim][no];
	}
	return 0;
}

void INFLOW_setParticleType() {
	int iParticle;
	for (iParticle = 0; iParticle < particle.totalNumber; iParticle++) {
		if (particle.isInflow[iParticle] == FALSE)continue;
		particle.type[iParticle] = INFLOW_getTypeOfInflowParticle(iParticle);
	}
}

void INFLOW_setVelocity() {
	int iParticle, iDim;
	if (parameter.flagOfInflowBoundaryCondition == OFF)return;
	for (iParticle = 0; iParticle < particle.totalNumber; iParticle++) {
		if (particle.isInflow[iParticle] == FALSE)continue;
		for (iDim = 0; iDim < 3; iDim++) {
			particle.velocity[iDim][iParticle] = parameter.inflowVelocity[iDim];
		}
		if (NumberOfDimensions == 2)particle.velocity[ZDIM][iParticle] = 0.0;

	}

}

void INFLOW_changeInflowParticles() {
	int iParticle, iDim;
	int dummy = -1;
	double norm = 0;

	if (parameter.flagOfInflowBoundaryCondition == OFF)return;
	
	for (iDim = 0; iDim < NumberOfDimensions; iDim++) {
		norm += parameter.inflowVelocity[iDim] * parameter.inflowVelocity[iDim];
	}
	norm = sqrt(norm);

	for (iParticle = 0; iParticle < particle.totalNumber; iParticle++) {
		if (particle.isInflow[iParticle] == FALSE)continue;
		if (particle.type[iParticle] != INFLOW_getTypeOfInflowParticle(iParticle)) {
			if (INFLOW_getTypeOfInflowParticle(iParticle) == parameter.wallType)particle.type[iParticle] = parameter.wallType;
		}

		if (INFLOW_getTypeOfInflowParticle(iParticle) == 0) {
			STACK_pop(&ghostStack, &dummy);
			
			for (iDim = 0; iDim < 3; iDim++) {
				particle.position[iDim][dummy] = particle.position[iDim][iParticle] - 4 * particle.averageDistance * parameter.inflowVelocity[iDim] / norm;
				particle.type[dummy] = parameter.dummyWallType;
				particle.isInflow[dummy] = TRUE;
				particle.moleOfBubbles[dummy] = 0;
				particle.numberOfBubbles[dummy] = 0;
				particle.concentrationOfImpurities[dummy] = 0;
				particle.cavitationBubblesRadius[dummy] = 0;

			}
			if (NumberOfDimensions == 2)particle.velocity[ZDIM][dummy] = 0.0;

			particle.type[iParticle] = 0;
			particle.moleOfBubbles[iParticle] = parameter.inflowMoleOfBubbles;
			particle.numberOfBubbles[iParticle] = parameter.inflowNumberOfBubbles;
			particle.concentrationOfImpurities[iParticle] = 0;
			particle.isInflow[iParticle] = FALSE;

		}
	}

}
