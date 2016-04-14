#include <stdio.h>
#include <math.h>

#include "define.h"
#include "extern.h"
#include "struct.h"
#include "neigh.h"
#include "other.h"
#include "copy.h"
#include "init.h"



double
OUTFLOW_judgeBounds(double x, double y, double z) {
	double norm;
	norm = sqrt(parameter.normalVector[XDIM] * parameter.normalVector[XDIM]
		+ parameter.normalVector[YDIM] * parameter.normalVector[YDIM]
		+ parameter.normalVector[ZDIM] * parameter.normalVector[ZDIM]);
	return parameter.normalVector[XDIM] * (x - parameter.positionVector[XDIM])
		+ parameter.normalVector[YDIM] * (y - parameter.positionVector[YDIM])
		+ parameter.normalVector[ZDIM] * (z - parameter.positionVector[ZDIM]) / norm;
}

int
OUTFLOW_allowedJudgingFreeSurface(int iParticle) {
	if (parameter.flagOfNegativePressure == OFF)return TRUE;
	//if (particle.type[iParticle] == parameter.wallType || particle.type[iParticle] == parameter.typeNumberOfRigidParticle_forForcedMotion)return TRUE;
	//if (OUTFLOW_judgeBounds(particle.position[XDIM][iParticle], particle.position[YDIM][iParticle], particle.position[ZDIM][iParticle]) 
	//	> -parameter.relaxationCoefficientOfBoundaryCondition * particle.averageDistance)return TRUE;

	if (OUTFLOW_judgeBounds(particle.position[XDIM][iParticle], particle.position[YDIM][iParticle], particle.position[ZDIM][iParticle])
		> -10.0 * particle.averageDistance)return TRUE;

	return FALSE;
}


int
OUTFLOW_getCurrentNumberOfFluidParticles() {
	int iParticle;
	int counter = 0;


	for (iParticle = 0; iParticle < particle.totalNumber; iParticle++) {

		if (particle.type[iParticle] == parameter.wallType ||
			particle.type[iParticle] == parameter.dummyWallType ||
			particle.type[iParticle] == parameter.typeNumberOfRigidParticle_forForcedMotion ||
			particle.type[iParticle] == GHOST)continue;

		if (OUTFLOW_judgeBounds(particle.position[XDIM][iParticle], particle.position[YDIM][iParticle], particle.position[ZDIM][iParticle]) > -parameter.lengthOfOutflow + particle.averageDistance/2) 	
		counter++;
		
	}
	return counter;
}

double
OUTFLOW_getL(int iParticle, int jParticle) {
	double l = 0;
	int iDim;
	for (iDim = 0; iDim < NumberOfDimensions; iDim++) {
		l += (particle.position[iDim][jParticle] - particle.position[iDim][iParticle]) * (-parameter.normalVector[iDim]);

	}
	return l;
}


double
OUTFLOW_getVelocityOfUpstream(int iParticle) {
	int* numberOfNeighborParticles;
	int** neighborTable;

	int jParticle;
	int totalNeigh;

	int iDim,iNeigh;
	int kParticle;

	double maxL,vel;
	double tmpL,tmpVel[3];

	NEIGH_setNeighborTable(particle.position);
	NEIGH_selectNeighborTable(&numberOfNeighborParticles, &neighborTable, parameter.radiusOfParticleNumberDensity_ratio);

	kParticle = -1;
	totalNeigh = 0;
	vel = 0;

	for (iNeigh = 0; iNeigh < numberOfNeighborParticles[iParticle]; iNeigh++) {
		jParticle = neighborTable[iParticle][iNeigh];

		if (particle.type[jParticle] == parameter.wallType ||
			particle.type[jParticle] == parameter.dummyWallType ||
			particle.type[jParticle] == parameter.typeNumberOfRigidParticle_forForcedMotion ||
			particle.type[jParticle] == GHOST)continue;


		tmpL = OUTFLOW_getL(iParticle, jParticle);
		if (iNeigh == 0)maxL = tmpL;

		if (tmpL > maxL) {
			maxL = tmpL;
			kParticle = jParticle;
		}

	}
	if (kParticle == -1) {
		for (iDim = 0; iDim < NumberOfDimensions; iDim++) {
			vel += particle.velocity[iDim][iParticle] * parameter.normalVector[iDim] 
				/ sqrt(parameter.normalVector[XDIM] * parameter.normalVector[XDIM]
				+ parameter.normalVector[YDIM] * parameter.normalVector[YDIM]
				+ parameter.normalVector[ZDIM] * parameter.normalVector[ZDIM]);
		}
		
		return vel;
	}

	for (iDim = 0; iDim < NumberOfDimensions; iDim++) {
		tmpVel[iDim] = 0;
	}
	for (iNeigh = 0; iNeigh < numberOfNeighborParticles[kParticle]; iNeigh++) {
		jParticle = neighborTable[kParticle][iNeigh];
		if (particle.type[jParticle] == parameter.wallType ||
			particle.type[jParticle] == parameter.dummyWallType ||
			particle.type[jParticle] == parameter.typeNumberOfRigidParticle_forForcedMotion ||
			particle.type[jParticle] == GHOST)continue;

		for (iDim = 0; iDim < NumberOfDimensions; iDim++) {
			tmpVel[iDim] += particle.velocity[iDim][jParticle];
		}
		totalNeigh++;

	}
	
	for (iDim = 0; iDim < NumberOfDimensions; iDim++) {
		vel += tmpVel[iDim] * parameter.normalVector[iDim] / totalNeigh;

	}
	
	return vel;
}

int
OUTFLOW_deleteOutflowParticles() {
	int number,iParticle;
	number = 0;

	if (parameter.flagOfOutflowBoundaryCondition == OFF)return -1;

	for (iParticle = 0; iParticle < particle.totalNumber; iParticle++) {
		if (particle.type[iParticle] == parameter.wallType ||
			particle.type[iParticle] == parameter.dummyWallType ||
			particle.type[iParticle] == parameter.typeNumberOfRigidParticle_forForcedMotion ||
			particle.type[iParticle] == GHOST)continue;
		if (OUTFLOW_judgeBounds(particle.position[XDIM][iParticle], particle.position[YDIM][iParticle], particle.position[ZDIM][iParticle])>0) {
			OTHER_changeParticleTypeIntoGhost(iParticle);
			number++;
		}
	}

	if (number > 0) {
		fprintf(FpForLog, "Particles overflow. -----total:               %d\n", number);

	}

	return number;

}


void
OUTFLOW_correctOutflowVelocity() {
	int iParticle, jParticle;
	int iDim;

	double velocityUpstream;
	double du;
	int numberOfParticles;
	double normalTotal;


	if (parameter.flagOfOutflowBoundaryCondition == OFF)return;

	numberOfParticles = OUTFLOW_getCurrentNumberOfFluidParticles();
	du = (parameter.initTotalOutflow - numberOfParticles) * pow(particle.averageDistance, NumberOfDimensions) /
		(parameter.widthOfOutflow * parameter.depthOfOutflow * parameter.relaxationCoefficientOfBoundaryCondition * timer.dt);
	if (du <= 0)return;
	//printf("%f\n", du);
	normalTotal = sqrt(parameter.normalVector[XDIM] * parameter.normalVector[XDIM]
		+ parameter.normalVector[YDIM] * parameter.normalVector[YDIM]
		+ parameter.normalVector[ZDIM] * parameter.normalVector[ZDIM]);

	COPY_copy2dimDoubleArray(3, particle.totalNumber, particle.velocity_upstream, particle.velocity);


	for (iParticle = 0; iParticle < particle.totalNumber; iParticle++) {
		if (particle.type[iParticle] == parameter.wallType ||
			particle.type[iParticle] == parameter.dummyWallType ||
			particle.type[iParticle] == parameter.typeNumberOfRigidParticle_forForcedMotion ||
			particle.type[iParticle] == GHOST)continue;

		if (OUTFLOW_judgeBounds(particle.position[XDIM][iParticle], particle.position[YDIM][iParticle], particle.position[ZDIM][iParticle]) > -parameter.lengthOfFixedVelocityRegion * particle.averageDistance
			&& OUTFLOW_judgeBounds(particle.position[XDIM][iParticle], particle.position[YDIM][iParticle], particle.position[ZDIM][iParticle]) <= 0) {
			velocityUpstream = OUTFLOW_getVelocityOfUpstream(iParticle);

			for (iDim = 0; iDim < NumberOfDimensions; iDim++) {
				if (parameter.normalVector[iDim] == 0)continue;
				// 法線ベクトルが座標軸と平行な場合のみ正しく動作
				//particle.velocity[iDim][iParticle] = parameter.normalVector[iDim] * (OUTFLOW_getVelocityOfUpstream(iParticle) - du) / normalTotal;
				particle.velocity_upstream[iDim][iParticle] = (velocityUpstream - du)/parameter.normalVector[iDim];
			}
			if (particle.velocity_upstream[XDIM][iParticle] * parameter.normalVector[XDIM] + particle.velocity_upstream[YDIM][iParticle] * parameter.normalVector[YDIM]
				+ particle.velocity_upstream[ZDIM][iParticle] * parameter.normalVector[ZDIM] < 0) {
				for (iDim = 0; iDim < NumberOfDimensions; iDim++) {
					if (parameter.normalVector[iDim] == 0)continue;
					particle.velocity_upstream[iDim][iParticle] = 0;
					
				}

			}
			/*if (particle.velocity_upstream[XDIM][iParticle] * parameter.normalVector[XDIM] + particle.velocity_upstream[YDIM][iParticle] * parameter.normalVector[YDIM]
				+ particle.velocity_upstream[ZDIM][iParticle] * parameter.normalVector[ZDIM]  > velocityUpstream * 2) {
				for (iDim = 0; iDim < NumberOfDimensions; iDim++) {
					if (parameter.normalVector[iDim] == 0)continue;
					particle.velocity_upstream[iDim][iParticle] = velocityUpstream * 2 * parameter.normalVector[iDim];

				}
			}*/


		}

	}

	COPY_copy2dimDoubleArray(3, particle.totalNumber, particle.velocity, particle.velocity_upstream);

}

void OUTFLOW_correctYamakawaOutflowVelocity() {
	int iParticle;// , jParticle;
	int iDim;

	double du,gamma;
	int outflowDomainNumber;
	double inflowVelocity,outflowVelocity,inflownorm;


	if (parameter.flagOfOutflowBoundaryCondition == OFF)return;
	if (parameter.flagOfInflowBoundaryCondition == OFF)return;

	inflownorm = 0;
	for (iDim = 0; iDim < NumberOfDimensions; iDim++) {
		inflownorm += parameter.inflowVelocity[iDim] * parameter.inflowVelocity[iDim];
	}
	inflownorm = sqrt(inflownorm);

	gamma = timer.dt * inflownorm / (2 * particle.averageDistance);

	if (gamma > 0.1)gamma = 0.1;

	inflowVelocity = 0; outflowVelocity = 0; outflowDomainNumber = 0;
	
	for (iParticle = 0; iParticle < particle.totalNumber; iParticle++) {
		if (particle.isInflow[iParticle] == TRUE) {
			for (iDim = 0; iDim < NumberOfDimensions; iDim++) {
				inflowVelocity += parameter.inflowVelocity[iDim] * parameter.inflowVelocity[iDim] / inflownorm;
			}
		}

		if (particle.type[iParticle] == parameter.wallType ||
			particle.type[iParticle] == parameter.dummyWallType ||
			particle.type[iParticle] == parameter.typeNumberOfRigidParticle_forForcedMotion ||
			particle.type[iParticle] == GHOST)continue;

		if (OUTFLOW_judgeBounds(particle.position[XDIM][iParticle], particle.position[YDIM][iParticle], particle.position[ZDIM][iParticle]) > -parameter.lengthOfFixedVelocityRegion * particle.averageDistance) {
			for (iDim = 0; iDim < NumberOfDimensions; iDim++) {
				outflowVelocity += particle.velocity[iDim][iParticle] * parameter.normalVector[iDim];
			}
			outflowDomainNumber++;
		}
		
	}

	inflowVelocity /= 4.0;

	fprintf(FpForLog, "OutflowDomainNumber:%d\n", outflowDomainNumber);
	fprintf(stderr, "OutflowDomainNumber:%d\n", outflowDomainNumber);
	

	if (outflowDomainNumber == 0)return;
	
	//du = gamma * ( (OUTFLOW_getCurrentNumberOfFluidParticles() - parameter.initTotalOutflow) * pow(particle.averageDistance, NumberOfDimensions) +
	//	timer.dt * particle.averageDistance * (outflowVelocity - inflowVelocity) ) / (timer.dt * particle.averageDistance * outflowDomainNumber) ;
	du = gamma * ((OUTFLOW_getCurrentNumberOfFluidParticles() - parameter.initTotalOutflow) * pow(particle.averageDistance, NumberOfDimensions) +
		timer.dt * particle.averageDistance * (outflowVelocity - inflowVelocity)) / (timer.dt * parameter.widthOfOutflow * parameter.depthOfOutflow);
	for (iParticle = 0; iParticle < particle.totalNumber; iParticle++) {
		if (particle.type[iParticle] == parameter.wallType ||
			particle.type[iParticle] == parameter.dummyWallType ||
			particle.type[iParticle] == parameter.typeNumberOfRigidParticle_forForcedMotion ||
			particle.type[iParticle] == GHOST)continue;

		if (OUTFLOW_judgeBounds(particle.position[XDIM][iParticle], particle.position[YDIM][iParticle], particle.position[ZDIM][iParticle]) > -parameter.lengthOfFixedVelocityRegion * particle.averageDistance) {
			for (iDim = 0; iDim < NumberOfDimensions; iDim++) {
				if(parameter.normalVector[iDim]!=0)particle.velocity[iDim][iParticle] += du / parameter.normalVector[iDim];

			}
		}
		//TODO
		if (particle.velocity[XDIM][iParticle] * parameter.normalVector[XDIM] + particle.velocity[YDIM][iParticle] * parameter.normalVector[YDIM] 
			+ particle.velocity[ZDIM][iParticle] * parameter.normalVector[ZDIM]
			> (particle.velocity_previous[XDIM][iParticle] * parameter.normalVector[XDIM] + particle.velocity_previous[YDIM][iParticle] * parameter.normalVector[YDIM]
				+ particle.velocity_previous[ZDIM][iParticle] * parameter.normalVector[ZDIM]) * 2.0
			&& (particle.velocity_previous[XDIM][iParticle] * parameter.normalVector[XDIM] + particle.velocity_previous[YDIM][iParticle] * parameter.normalVector[YDIM]
				+ particle.velocity_previous[ZDIM][iParticle] * parameter.normalVector[ZDIM]) > 0) {
			for (iDim = 0; iDim < NumberOfDimensions; iDim++) {
				if (parameter.normalVector[iDim] != 0)particle.velocity[iDim][iParticle] = particle.velocity_previous[iDim][iParticle] * 2;

			}
		}

		if (particle.velocity_upstream[XDIM][iParticle] * parameter.normalVector[XDIM] + particle.velocity_upstream[YDIM][iParticle] * parameter.normalVector[YDIM]
			+ particle.velocity_upstream[ZDIM][iParticle] * parameter.normalVector[ZDIM] < 0) {
			for (iDim = 0; iDim < NumberOfDimensions; iDim++) {
				if (parameter.normalVector[iDim] == 0)continue;
				particle.velocity_upstream[iDim][iParticle] = 0;

			}

		}

	}

}
	
//丼勘定やめてみる
//ノズル内で厳密に計算
void OUTFLOW_correctYamakawaOutflowVelocity2() {
	int iParticle, jParticle;
	int iDim;

	double du, gamma;
	int outflowDomainNumber;
	double inflowVelocity, outflowVelocity, inflownorm;


	if (parameter.flagOfOutflowBoundaryCondition == OFF)return;
	if (parameter.flagOfInflowBoundaryCondition == OFF)return;
	

	inflownorm = 0;
	for (iDim = 0; iDim < NumberOfDimensions; iDim++) {
		inflownorm += parameter.inflowVelocity[iDim] * parameter.inflowVelocity[iDim];
	}
	inflownorm = sqrt(inflownorm);

	//gamma = timer.dt * inflownorm / (2 * particle.averageDistance);
	gamma = 0.01;

	if (gamma > 0.1)gamma = 0.1;

	inflowVelocity = 0; outflowVelocity = 0; outflowDomainNumber = 0;

	for (iParticle = 0; iParticle < particle.totalNumber; iParticle++) {
		if (particle.type[iParticle] == parameter.wallType ||
			particle.type[iParticle] == parameter.dummyWallType ||
			particle.type[iParticle] == parameter.typeNumberOfRigidParticle_forForcedMotion ||
			particle.type[iParticle] == GHOST)continue;

		if (OUTFLOW_judgeBounds(particle.position[XDIM][iParticle], particle.position[YDIM][iParticle], particle.position[ZDIM][iParticle])
			* OUTFLOW_judgeBounds(particle.position[XDIM][iParticle] + particle.velocity[XDIM][iParticle] * timer.dt, 
				particle.position[YDIM][iParticle]  + particle.velocity[YDIM][iParticle] * timer.dt
				, particle.position[ZDIM][iParticle] + particle.velocity[ZDIM][iParticle] * timer.dt) < 0) {
			for (iDim = 0; iDim < NumberOfDimensions; iDim++) {
				outflowVelocity += particle.velocity[iDim][iParticle] * parameter.normalVector[iDim];
			}
			outflowDomainNumber++;
		}

		if (OUTFLOW_judgeBounds(particle.position[XDIM][iParticle], particle.position[YDIM][iParticle], particle.position[ZDIM][iParticle]) < -parameter.lengthOfOutflow 
			&& OUTFLOW_judgeBounds(particle.position[XDIM][iParticle] + particle.velocity[XDIM][iParticle] * timer.dt,
				particle.position[YDIM][iParticle] + particle.velocity[YDIM][iParticle] * timer.dt
				, particle.position[ZDIM][iParticle] + particle.velocity[ZDIM][iParticle] * timer.dt) > -parameter.lengthOfOutflow) {
			for (iDim = 0; iDim < NumberOfDimensions; iDim++) {
				inflowVelocity += particle.velocity[iDim][iParticle] * parameter.inflowVelocity[iDim] / inflownorm;
			}
		}

	}

	//inflowVelocity /= 4.0;

	//fprintf(FpForLog, "OutflowDomainNumber:%d\n", outflowDomainNumber);
	//fprintf(stderr, "OutflowDomainNumber:%d\n", outflowDomainNumber);


	//if (outflowDomainNumber == 0)return;

	//du = gamma * ( (OUTFLOW_getCurrentNumberOfFluidParticles() - parameter.initTotalOutflow) * pow(particle.averageDistance, NumberOfDimensions) +
	//	timer.dt * particle.averageDistance * (outflowVelocity - inflowVelocity) ) / (timer.dt * particle.averageDistance * outflowDomainNumber) ;
	du = gamma * ((OUTFLOW_getCurrentNumberOfFluidParticles() - parameter.initTotalOutflow) * pow(particle.averageDistance, NumberOfDimensions) +
		timer.dt * particle.averageDistance * (inflowVelocity - outflowVelocity)) / (timer.dt * parameter.widthOfOutflow * parameter.depthOfOutflow);

	fprintf(FpForLog, "inflowVelocity:%lf\toutflowVelocity:%lf\n", inflowVelocity, outflowVelocity);
	fprintf(stderr, "inflowVelocity:%lf\toutflowVelocity:%lf\n", inflowVelocity, outflowVelocity);
	fprintf(FpForLog, "du:%lf\n", du);
	fprintf(stderr, "du:%lf\n", du);

	if (du > 0)return;

	for (iParticle = 0; iParticle < particle.totalNumber; iParticle++) {
		if (particle.type[iParticle] == parameter.wallType ||
			particle.type[iParticle] == parameter.dummyWallType ||
			particle.type[iParticle] == parameter.typeNumberOfRigidParticle_forForcedMotion ||
			particle.type[iParticle] == GHOST)continue;

		if (OUTFLOW_judgeBounds(particle.position[XDIM][iParticle], particle.position[YDIM][iParticle], particle.position[ZDIM][iParticle]) > -parameter.lengthOfFixedVelocityRegion * particle.averageDistance) {
			for (iDim = 0; iDim < NumberOfDimensions; iDim++) {
				if (parameter.normalVector[iDim] != 0) {
					particle.velocity[iDim][iParticle] += du / parameter.normalVector[iDim];

					//逆流しないと仮定する
					if (particle.velocity[iDim][iParticle] * parameter.normalVector[iDim] < 0)
						particle.velocity[iDim][iParticle] = 0;

				}
			}

		}
		//TODO
		/*if (particle.velocity[XDIM][iParticle] * parameter.normalVector[XDIM] + particle.velocity[YDIM][iParticle] * parameter.normalVector[YDIM]
			+ particle.velocity[ZDIM][iParticle] * parameter.normalVector[ZDIM]
			>(particle.velocity_previous[XDIM][iParticle] * parameter.normalVector[XDIM] + particle.velocity_previous[YDIM][iParticle] * parameter.normalVector[YDIM]
				+ particle.velocity_previous[ZDIM][iParticle] * parameter.normalVector[ZDIM]) * 2.0
			&& (particle.velocity_previous[XDIM][iParticle] * parameter.normalVector[XDIM] + particle.velocity_previous[YDIM][iParticle] * parameter.normalVector[YDIM]
				+ particle.velocity_previous[ZDIM][iParticle] * parameter.normalVector[ZDIM]) > 0) {
			for (iDim = 0; iDim < NumberOfDimensions; iDim++) {
				if (parameter.normalVector[iDim] != 0)particle.velocity[iDim][iParticle] = particle.velocity_previous[iDim][iParticle] * 2;

			}
		}
		*/
		
	}

}

//柴田先生のを改良
void OUTFLOW_correctShibataOutflowVelocity() {
	int iParticle, jParticle;
	int iDim;

	double du;

	
	du =  (OUTFLOW_getCurrentNumberOfFluidParticles() - parameter.initTotalOutflow) * pow(particle.averageDistance, NumberOfDimensions)
		/ (parameter.relaxationCoefficientOfBoundaryCondition * timer.dt * parameter.widthOfOutflow * parameter.depthOfOutflow);


	fprintf(FpForLog, "currentN:%d\n", OUTFLOW_getCurrentNumberOfFluidParticles());
	fprintf(stderr, "currentN:%d\n", OUTFLOW_getCurrentNumberOfFluidParticles());
	fprintf(FpForLog, "du:%lf\n", du);
	fprintf(stderr, "du:%lf\n", du);

	if (du > 0)return;

	for (iParticle = 0; iParticle < particle.totalNumber; iParticle++) {
		if (particle.type[iParticle] == parameter.wallType ||
			particle.type[iParticle] == parameter.dummyWallType ||
			particle.type[iParticle] == parameter.typeNumberOfRigidParticle_forForcedMotion ||
			particle.type[iParticle] == GHOST)continue;

		if (OUTFLOW_judgeBounds(particle.position[XDIM][iParticle], particle.position[YDIM][iParticle], particle.position[ZDIM][iParticle]) > -parameter.lengthOfFixedVelocityRegion * particle.averageDistance) {
			for (iDim = 0; iDim < NumberOfDimensions; iDim++) {
				if (parameter.normalVector[iDim] != 0) {
					particle.velocity[iDim][iParticle] += du / parameter.normalVector[iDim];

					//逆流しないと仮定する
					if (particle.velocity[iDim][iParticle] * parameter.normalVector[iDim] < 0)
						particle.velocity[iDim][iParticle] = 0;

				}
			}

		}
		//TODO
		/*if (particle.velocity[XDIM][iParticle] * parameter.normalVector[XDIM] + particle.velocity[YDIM][iParticle] * parameter.normalVector[YDIM]
		+ particle.velocity[ZDIM][iParticle] * parameter.normalVector[ZDIM]
		>(particle.velocity_previous[XDIM][iParticle] * parameter.normalVector[XDIM] + particle.velocity_previous[YDIM][iParticle] * parameter.normalVector[YDIM]
		+ particle.velocity_previous[ZDIM][iParticle] * parameter.normalVector[ZDIM]) * 2.0
		&& (particle.velocity_previous[XDIM][iParticle] * parameter.normalVector[XDIM] + particle.velocity_previous[YDIM][iParticle] * parameter.normalVector[YDIM]
		+ particle.velocity_previous[ZDIM][iParticle] * parameter.normalVector[ZDIM]) > 0) {
		for (iDim = 0; iDim < NumberOfDimensions; iDim++) {
		if (parameter.normalVector[iDim] != 0)particle.velocity[iDim][iParticle] = particle.velocity_previous[iDim][iParticle] * 2;

		}
		}
		*/

	}

}
