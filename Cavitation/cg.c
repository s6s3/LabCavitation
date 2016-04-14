#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "define.h"
#include "extern.h"
#include "struct.h"

#include "other.h"
#include "init.h"
#include "memory.h"
#include "mathe.h"
#include "copy.h"

#include "solver.h"
#include "neigh.h"
#include "density.h"
#include "pressure.h"
#include "outflow.h"

int CG_solver(int totalNumber, double** matrix, double* solution, double* sourceTerm, int maxIteration, int minIteration, double error) {
	if (totalNumber <= 0)return NO;

	double *r, *p, *ax;
	double *r_n, *p_n, *solution_n;
	int iIteration,iParticle;
	double alpha, beta;
	double dot1, dot2;
	double normOfSourceTerm;

	solution_n = MEMORY_allocateMemoryFor1dimDoubleArray(totalNumber, "solution_n [in cg.c]");
	solution = MEMORY_allocateMemoryFor1dimDoubleArray(totalNumber, "solution [in cg.c]");
	r = MEMORY_allocateMemoryFor1dimDoubleArray(totalNumber, "r [in cg.c]");
	p = MEMORY_allocateMemoryFor1dimDoubleArray(totalNumber, "p [in cg.c]");
	r_n = MEMORY_allocateMemoryFor1dimDoubleArray(totalNumber, "r_n [in cg.c]");
	p_n = MEMORY_allocateMemoryFor1dimDoubleArray(totalNumber, "p_n [in cg.c]");
	ax = MEMORY_allocateMemoryFor1dimDoubleArray(totalNumber, "ax [in cg.c]");

	OTHER_fillOneDimDoubleArrayWithZero(totalNumber, solution);
	
	MATHE_multiplySquareMatrixByVector(totalNumber, ax, matrix, solution);
	MATHE_subtractVectorFromVector(totalNumber, r, sourceTerm, ax, 1);
	COPY_copy1dimDoubleArray(totalNumber, p, r);

	normOfSourceTerm = MATHE_absoluteValueOfVector(sourceTerm, totalNumber);

	for (iIteration = 0; iIteration < maxIteration; iIteration++) {
		dot1 = MATHE_innerProductOfVectors(r, r, totalNumber);
		MATHE_multiplySquareMatrixByVector(totalNumber, ax, matrix, p);
		dot2 = MATHE_innerProductOfVectors(p, ax, totalNumber);

		alpha = dot1 / dot2;

		MATHE_addVectorToVector(totalNumber, solution_n, solution, p, alpha);
		MATHE_subtractVectorFromVector(totalNumber, r_n, r, ax, alpha);
		
		if (MATHE_absoluteValueOfVector(r_n, totalNumber) <= error * normOfSourceTerm && iIteration>=minIteration) {
			COPY_copy1dimDoubleArray(totalNumber, solution, solution_n);

			SOLVER_displayStateOfConvergence(FpForLog, iIteration, minIteration);

			free(solution_n);
			free(r);
			free(p);
			free(r_n);
			free(p_n);
			free(ax);

			return YES;
		}

		dot2 = MATHE_innerProductOfVectors(r_n, r_n, totalNumber);
		beta = dot2 / dot1;
		MATHE_addVectorToVector(totalNumber, p_n, r_n, p, beta);

		COPY_copy1dimDoubleArray(totalNumber, r, r_n);
		COPY_copy1dimDoubleArray(totalNumber, p, p_n);
		COPY_copy1dimDoubleArray(totalNumber, solution, solution_n);

	}

	SOLVER_displayWarnigMessageForExcessOfIterationNumber(iIteration, maxIteration);


	free(solution_n);
	free(r);
	free(p);
	free(r_n);
	free(p_n);
	free(ax);

	return NO;
}

int CG_getDimension() {
	int count = 0;
	int iParticle;
	for (iParticle = 0; iParticle < particle.totalNumber; iParticle++) {
		if (particle.flagOfBoundaryCondition[iParticle] == INNER_PARTICLE) {
			count++;
		}
	}
	return count;
}

void CG_setMap(double* totalmap, double* calcmap, int calcNum) {
	int count = 0;
	int iParticle;
	totalmap = MEMORY_allocateMemoryFor1dimDoubleArray(particle.totalNumber, "totalmap [in cg.c]");
	calcmap = MEMORY_allocateMemoryFor1dimDoubleArray(calcNum, "calcmap [in cg.c]");
	for (iParticle = 0; iParticle < particle.totalNumber; iParticle++) {
		totalmap[iParticle] = -1;
		if (particle.flagOfBoundaryCondition[iParticle] == INNER_PARTICLE) {
			calcmap[count] = iParticle;
			totalmap[iParticle] = count;
		}
	}
}


void
CG_setSourceTerm(double* sourceTerm,double* totalmap, int calcNum) {

	int iParticle;

	double n0 = parameter.nZeroOfParticleNumberDensity;
	double dt_squared = timer.dt * timer.dt;

	sourceTerm = MEMORY_allocateMemoryFor1dimDoubleArray(calcNum, "sourceTerm [in cg.c]");

	int iCalc;

	for (iParticle = 0; iParticle < particle.totalNumber; iParticle++) {
		iCalc = totalmap[iParticle];
		if (iCalc == -1)continue;
		sourceTerm[iCalc] = 0.0;

		if (particle.flagOfBoundaryCondition[iCalc] == INNER_PARTICLE) {

			if (parameter.flagOfTanakaAndMasunagaModel == ON) {
				fprintf(FpForLog, "WARNING: Tanaka and Masunaga Model is not implemented. \n ");
				exit(1);
				/*
				sigma = 0.0;

				for (iNeigh = 0; iNeigh< numberOfNeighborParticles[iParticle]; iNeigh++) {

					jParticle = neighborTable[iParticle][iNeigh];

					if (particle.type[jParticle] == GHOST) continue;

					distanceIJ = DISTANCE_calculateDistanceBetweenParticles(particle.position, iParticle, jParticle);

					distanceIJ_squared = DISTANCE_calculateSquaredDistance(particle.position, iParticle, jParticle);

					weightIJ = WEIGHT_calculateWeightFunction(distanceIJ, parameter.radiusOfLaplacianForPressure);

					for (iDim = 0; iDim < NumberOfDimensions; iDim++) {
						sigma -= ((particle.velocity[iDim][jParticle] - particle.velocity[iDim][iParticle])*(particle.position[iDim][jParticle] - particle.position[iDim][iParticle])*weightIJ / distanceIJ_squared);
					}
				}
				particle.sourceTermOfPressure[iParticle] = (NumberOfDimensions*sigma / (dt*n0)) + parameter.valueOfGamma*((particle.particleNumberDensity_previous[iParticle] - n0) / (n0*dt_squared));
				*/
			}
			else {
				sourceTerm[iCalc] = (particle.particleNumberDensity[iParticle] - n0) / (dt_squared * n0);
			}

			if ((particle.particleNumberDensity[iParticle] / n0) < parameter.thresholdOfParticleNumberDensity_ratio && parameter.flagOfNegativePressure == ON
				&& OUTFLOW_judgeBounds(particle.position[XDIM][iParticle], particle.position[YDIM][iParticle], particle.position[ZDIM][iParticle]) > -particle.averageDistance)
				sourceTerm[iCalc] = 0.0;

		}
	}
}


/*
void
CG_setCoefficientMatrixOfPressure(double** matrix, double* totalmap, int calcNum) {

	int    iParticle, jParticle;
	int    iNeigh;

	double n0 = parameter.nZeroOfLaplacianForPressure;
	double lambda = parameter.lambdaOfLaplacianForPressure;

	double distanceIJ;
	double coefficientOfParticleJ;
	double averageMassDensity;

	double dt_squared = (timer.dt * timer.dt);
	double weightIJ;


	int    *numberOfNeighborParticles;
	int    **neighborTable;

	int iCalc;

	matrix = MEMORY_allocateMemoryFor2dimDoubleArray(iCalc
		, iCalc, "matrix [in cg.c]");


	NEIGH_selectNeighborTable(&numberOfNeighborParticles
		, &neighborTable
		, parameter.radiusOfLaplacianForPressure_ratio
		);

	OTHER_fill2dimDoubleArrayWithZero(calcNum, calcNum, matrix);

	for (iParticle = 0; iParticle < particle.totalNumber; iParticle++) {
		iCalc = totalmap[iParticle];
		if (iCalc == -1)continue;

		if (particle.flagOfBoundaryCondition[iParticle] != INNER_PARTICLE) continue;

		particle.coefficientMatrixOfPressure[iParticle][0] = 0.0;

		for (iNeigh = 0; iNeigh < numberOfNeighborParticles[iParticle]; iNeigh++) {
			jParticle = neighborTable[iParticle][iNeigh];

			if (particle.flagOfBoundaryCondition[jParticle] == GHOST_OR_DUMMY) {
				particle.coefficientMatrixOfPressure[iParticle][iNeigh + 1] = 0.0;
			}
			else {
				distanceIJ = DISTANCE_calculateDistanceBetweenParticles(particle.position, iParticle, jParticle);

				averageMassDensity = MATHE_average(physicalProperty.massDensity[particle.type[iParticle]]
					, physicalProperty.massDensity[particle.type[jParticle]]);

				weightIJ = WEIGHT_calculateWeightFunction(distanceIJ, parameter.radiusOfLaplacianForPressure);

				coefficientOfParticleJ = 2.0 * NumberOfDimensions * weightIJ / (n0 * lambda  * averageMassDensity);

				particle.coefficientMatrixOfPressure[iParticle][iNeigh + 1] = -coefficientOfParticleJ;

				if (parameter.flagOfTanakaAndMasunagaModel == ON) {
					particle.coefficientMatrixOfPressure[iParticle][0] += parameter.valueOfC*coefficientOfParticleJ;
				}
				else {
					particle.coefficientMatrixOfPressure[iParticle][0] += coefficientOfParticleJ;
				}
			}
		}

		particle.coefficientMatrixOfPressure[iParticle][0] += physicalProperty.compressibility[particle.type[iParticle]] / (dt_squared);

	}

}
*/


void CG_calculatePressure() {
	int flagOfPressureCalculation;
	int calcNumber;
	double *solution;
	double *sourceTerm;
	double *totalMap;
	double *calcMap;
	double **matrix;


	NEIGH_setNeighborTable(particle.position);

	DENSITY_calculateParticleNumberDensity(particle.position);

	PRESSURE_setDirichletBoundaryCondition();

	flagOfPressureCalculation = PRESSURE_checkNecessityOfPressureCalculation();

	if (flagOfPressureCalculation == ON) {

		calcNumber = CG_getDimension();

		CG_setMap(totalMap, calcMap, calcNumber);

		CG_setSourceTerm(sourceTerm, totalMap, calcNumber);

		/*PRESSURE_setSourceTerm();

		PRESSURE_setCoefficientMatrixOfPressure();

		SOLVER_solveSimultaniousEquations();

		if (parameter.flagOfHighViscosityCalculation == ON) {
			PRESSURE_correctPressure();
		}

		PRESSURE_setMinusPressureZero();

		PRESSURE_setMinimumPressureAmongNeighbors();*/

	}
	else {

		OTHER_fillOneDimDoubleArrayWithZero(particle.totalNumber, particle.pressure);

		fprintf(FpForLog, "WARNING: Pressure was not calculated because all particles were surface-particles or wall-particles.\n");

		PRESSURE_setMinimumPressureAmongNeighbors();

	}
}