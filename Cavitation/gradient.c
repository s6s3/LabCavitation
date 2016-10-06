#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "define.h"
#include "struct.h"
#include "extern.h"
#include "pressure.h"
#include "memory.h"
#include "weight.h"
#include "mathe.h"
#include "convection.h"
#include "neigh.h"
#include "other.h"
#include "gradient.h"



void
GRADIENT_correctParticleVelocityAndPositionUsingPressureGradient( void ){

  OTHER_fill2dimDoubleArrayWithZero( NumberOfDimensions, particle.totalNumber, particle.velocity_correction);

  GRADIENT_calculatePressureGradientAndVelocityCorrection();


  MATHE_add2dimArrayTo2dimArray(  NumberOfDimensions, particle.totalNumber
				 ,particle.velocity, particle.velocity, particle.velocity_correction );

  CONVECTION_moveParticles( particle.position, particle.velocity_correction );

}


//豊田の手法（田中益永メソッド）
//近藤越塚メソッドの人工的圧力を加える(4/19)
void
GRADIENT_calculatePressureGradientAndVelocityCorrection( void ){

	int iParticle, jParticle, iNeigh;

	double xji;
	double yji;
	double zji = 0.0;
	double d_vx;
	double d_vy;
	double d_vz = 0.0;

	double distanceIJ;
	double distanceIJ_squared;
	double absoluteValueOfVelocityCorrection;
	double averageDensity;

	int    *numberOfNeighborParticles;
	int    **neighborTable;


	NEIGH_selectNeighborTable(   &numberOfNeighborParticles
					,&neighborTable
					,parameter.radiusOfGradient_ratio
					);

  
	for(iParticle=0; iParticle < particle.totalNumber; iParticle++){

		particle.velocity_correction[XDIM][iParticle] = 0.0;
		particle.velocity_correction[YDIM][iParticle] = 0.0;

		if(NumberOfDimensions == 3){
			particle.velocity_correction[ZDIM][iParticle] = 0.0;
		}


		if( particle.flagOfBoundaryCondition[iParticle] == GHOST_OR_DUMMY )continue;
		if( particle.type[iParticle] == GHOST                   )continue;
		if( particle.type[iParticle] == parameter.dummyWallType )continue;
		if( particle.type[iParticle] == parameter.wallType      )continue;

		absoluteValueOfVelocityCorrection = 0.0;

		for(iNeigh=0; iNeigh <  numberOfNeighborParticles[iParticle]; iNeigh++){
			jParticle = neighborTable[iParticle][iNeigh];

			if(particle.type[jParticle]==GHOST || particle.type[jParticle]== parameter.dummyWallType )continue;

			xji = (particle.position[XDIM][iParticle] - particle.position[XDIM][jParticle]);
			yji = (particle.position[YDIM][iParticle] - particle.position[YDIM][jParticle]);


			distanceIJ_squared = xji*xji + yji*yji;
      
			if(NumberOfDimensions == 3){
				zji = (particle.position[ZDIM][iParticle] - particle.position[ZDIM][jParticle]);
				distanceIJ_squared += zji*zji;
			}


			distanceIJ = sqrt(distanceIJ_squared);
			averageDensity = MATHE_average( physicalProperty.massDensity[particle.type[iParticle]]
							,physicalProperty.massDensity[particle.type[jParticle]] ); 

			if (parameter.flagOfGradientTensor == ON) absoluteValueOfVelocityCorrection = timer.dt;
			else absoluteValueOfVelocityCorrection = NumberOfDimensions * timer.dt;
			//absoluteValueOfVelocityCorrection *= (particle.pressure[jParticle] - particle.minPressureAmongNeighbors[iParticle])/distanceIJ;
			/*new*/
      
			if(parameter.flagOfKondoAndKoshizukaModel == ON) absoluteValueOfVelocityCorrection *= (particle.pressure[jParticle] - particle.pressure[iParticle] + parameter.artificialPressure) / distanceIJ;
			else absoluteValueOfVelocityCorrection *= (particle.pressure[jParticle] + particle.pressure[iParticle])/distanceIJ;


			//calculating tensor.
			if (parameter.flagOfGradientTensor == ON) {
				GRADIENT_zeroGradientTensor();
				parameter.gradientTensor[0][0] = xji * xji;
				parameter.gradientTensor[0][1] = xji * yji;
				parameter.gradientTensor[1][0] = xji * yji;
				parameter.gradientTensor[1][1] = yji * yji;
			}

			absoluteValueOfVelocityCorrection *= WEIGHT_calculateWeightFunction(distanceIJ, parameter.radiusOfGradient);

			absoluteValueOfVelocityCorrection /= (averageDensity * parameter.nZeroOfGradient);

			if (parameter.flagOfGradientTensor == ON) {
				GRADIENT_getGradientTensorInverse();
				d_vx = absoluteValueOfVelocityCorrection * (xji * parameter.gradientTensorInverse[0][0] + yji * parameter.gradientTensorInverse[0][1]) / distanceIJ;
				d_vy = absoluteValueOfVelocityCorrection * (xji * parameter.gradientTensorInverse[1][0] + yji * parameter.gradientTensorInverse[1][1]) / distanceIJ;

				if (NumberOfDimensions == 3) {
					d_vz = absoluteValueOfVelocityCorrection * zji / distanceIJ;
					fprintf(FpForLog, "WARNING:Not implemented 3d tensor.");
				}
			}
			else {
				d_vx = absoluteValueOfVelocityCorrection * xji / distanceIJ;
				d_vy = absoluteValueOfVelocityCorrection * yji / distanceIJ;

				if (NumberOfDimensions == 3) {
					d_vz = absoluteValueOfVelocityCorrection * zji / distanceIJ;
				}
			}

			if (parameter.flagOfKondoAndKoshizukaModel == ON && distanceIJ < parameter.radiusOfKondoCollision * particle.averageDistance) {
				d_vx -= 0.01 * (parameter.radiusOfKondoCollision * particle.averageDistance - distanceIJ) * timer.dt * xji / distanceIJ;
				d_vy -= 0.01 * (parameter.radiusOfKondoCollision * particle.averageDistance - distanceIJ) * timer.dt * yji / distanceIJ;
				if (NumberOfDimensions == 3) {
					d_vz -= 0.01 * (parameter.radiusOfKondoCollision * particle.averageDistance - distanceIJ) * timer.dt * zji / distanceIJ;

				}
			}
				
			particle.velocity_correction[XDIM][iParticle] += d_vx;
			particle.velocity_correction[YDIM][iParticle] += d_vy;

			if(NumberOfDimensions == 3){
				particle.velocity_correction[ZDIM][iParticle] += d_vz;
			}


		}


	}
}


void GRADIENT_zeroGradientTensor(void) {
	int i, j;
	for (i = 0; i < 3; i++) {
		for (j = 0; j < 3; j++) {
			parameter.gradientTensor[i][j] = 0.0;
			parameter.gradientTensorInverse[i][j] = 0.0;

		}
	}

}

int GRADIENT_getGradientTensorInverse(void) {
	double det = 0.0;
	int i, j;

	if (NumberOfDimensions == 3) {
		return -1;

	}
	det = parameter.gradientTensor[0][0] * parameter.gradientTensor[1][1]
		- parameter.gradientTensor[0][1] * parameter.gradientTensor[1][0];

	if (det != 0) {
		for (i = 0; i < 3; i++) {
			for (j = 0; j < 3; j++) {
				parameter.gradientTensorInverse[i][j] = (i==j) ? 1.0 : 0.0;
			}
		}
		parameter.gradientTensorInverse[0][0] =  parameter.gradientTensor[1][1] / det;
		parameter.gradientTensorInverse[0][1] = -parameter.gradientTensor[0][1] / det;
		parameter.gradientTensorInverse[1][0] = -parameter.gradientTensor[1][0] / det;
		parameter.gradientTensorInverse[1][1] =  parameter.gradientTensor[0][0] / det;

		return 0;
	}
	else {
		for (i = 0; i < 3; i++) {
			for (j = 0; j < 3; j++) {
				parameter.gradientTensorInverse[i][j] = (i == j) ? NumberOfDimensions : 0.0;
			}
		}

		return 1;
	}

}