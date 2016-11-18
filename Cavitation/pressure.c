#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "define.h"
#include "struct.h"
#include "extern.h"
#include "neigh.h"
#include "density.h"
#include "distance.h"
#include "weight.h"
#include "mathe.h"
#include "solver.h"
#include "domain.h"
#include "copy.h"
#include "neigh.h"
#include "memory.h"
#include "other.h"
#include "pressure.h"
#include "bucket.h"

#include "outflow.h"




void
PRESSURE_calculatePressure( void ){

  int flagOfPressureCalculation;


  NEIGH_setNeighborTable( particle.position );

  DENSITY_calculateParticleNumberDensity( particle.position );

  PRESSURE_setDirichletBoundaryCondition( );


  flagOfPressureCalculation = PRESSURE_checkNecessityOfPressureCalculation( );

  if( flagOfPressureCalculation == ON ){

    PRESSURE_setSourceTerm();

    PRESSURE_setCoefficientMatrixOfPressure();

    SOLVER_solveSimultaniousEquations();
    
    if(parameter.flagOfHighViscosityCalculation == ON){
       PRESSURE_correctPressure();
    }

	PRESSURE_setSurfacePressureZero();//’Ç‰Á

    PRESSURE_setMinusPressureZero();

    PRESSURE_setMinimumPressureAmongNeighbors();

  }else{

    OTHER_fillOneDimDoubleArrayWithZero( particle.totalNumber, particle.pressure );

    fprintf(FpForLog,"WARNING: Pressure was not calculated because all particles were surface-particles or wall-particles.\n");

    PRESSURE_setMinimumPressureAmongNeighbors();

  }

}




void
PRESSURE_setSourceTerm( void ){

  int iParticle;

  double n0         = parameter.nZeroOfParticleNumberDensity;
  double dt_squared = timer.dt * timer.dt;
    
  /*after*/
  int    jParticle;
  int    iNeigh;
  int    iDim;
    
  double distanceIJ;
  double distanceIJ_squared;
  double sigma;
  double s_0, s_1, s_2;
    
  double dt = timer.dt;
  double weightIJ;
    
  int    *numberOfNeighborParticles;
  int    **neighborTable;
  
  NEIGH_selectNeighborTable(  &numberOfNeighborParticles
                              ,&neighborTable
                              ,parameter.radiusOfLaplacianForPressure_ratio
                              );


  for(iParticle=0; iParticle < particle.totalNumber; iParticle++){

    if ( (particle.type[iParticle]==GHOST) || (particle.type[iParticle] == parameter.dummyWallType) ){

      particle.sourceTermOfPressure[iParticle] = 0.0;


    }else if( particle.flagOfBoundaryCondition[iParticle] == INNER_PARTICLE){
      
        if(parameter.flagOfTanakaAndMasunagaModel==ON/*&&parameter.flagOfNegativePressure==OFF*/){
           
           sigma=0.0;
           
           for(iNeigh=0; iNeigh< numberOfNeighborParticles[iParticle];iNeigh++){
            
               jParticle=neighborTable[iParticle][iNeigh];
    
               if(particle.type[jParticle]==GHOST) continue;
            
               distanceIJ = DISTANCE_calculateDistanceBetweenParticles( particle.position, iParticle, jParticle  );
            
               distanceIJ_squared = DISTANCE_calculateSquaredDistance(particle.position, iParticle, jParticle);
            
               weightIJ = WEIGHT_calculateWeightFunction( distanceIJ, parameter.radiusOfLaplacianForPressure);
            
               for(iDim= 0; iDim < NumberOfDimensions; iDim++){
                   sigma -= ((particle.velocity[iDim][jParticle]-particle.velocity[iDim][iParticle])*(particle.position[iDim][jParticle]-particle.position[iDim][iParticle])*weightIJ/distanceIJ_squared);
               }
           }
           particle.sourceTermOfPressure[iParticle]=(NumberOfDimensions*sigma/(dt*n0))+parameter.valueOfGamma*((particle.particleNumberDensity_previous[iParticle]-n0)/(n0*dt_squared));
                                                      
        }
		else if (parameter.flagOfKondoAndKoshizukaModel == ON && timer.iTimeStep>0) {
			s_2 = particle.particleNumberDensity[iParticle] - 2 * particle.particleNumberDensity_previous[iParticle] + particle.particleNumberDensity_prevstep[iParticle];
			s_1 = particle.particleNumberDensity_previous[iParticle] - particle.particleNumberDensity_prevstep[iParticle];
			s_0 = particle.particleNumberDensity_previous[iParticle] - n0;

			particle.sourceTermOfPressure[iParticle] = (s_2 * parameter.valueOfKondoAlpha + s_1 * parameter.valueOfKondoBeta + s_0 * parameter.valueOfKondoGamma) / (dt_squared * n0);

		}
		else{
           particle.sourceTermOfPressure[iParticle] = ( particle.particleNumberDensity[iParticle] - n0 )/( dt_squared * n0 ); 
        }

		/*if ((particle.particleNumberDensity[iParticle] / n0) < parameter.thresholdOfParticleNumberDensity_ratio && parameter.flagOfNegativePressure == ON 
			&& OUTFLOW_judgeBounds(particle.position[XDIM][iParticle], particle.position[YDIM][iParticle], particle.position[ZDIM][iParticle]) > -particle.averageDistance)
			particle.sourceTermOfPressure[iParticle] = 0.0;*/


    }else if(particle.flagOfBoundaryCondition[iParticle] == SURFACE_PARTICLE){

      particle.sourceTermOfPressure[iParticle] = 0.0;

    }
  }
}





void
PRESSURE_setCoefficientMatrixOfPressure( void ){

	int    iParticle, jParticle;
	int    iNeigh;

	double n0     = parameter.nZeroOfLaplacianForPressure;
	double lambda = parameter.lambdaOfLaplacianForPressure;

	double distanceIJ;
	double coefficientOfParticleJ;
	double averageMassDensity;

	double dt_squared = (timer.dt * timer.dt);
	double weightIJ;


	int    *numberOfNeighborParticles;
	int    **neighborTable;


	NEIGH_selectNeighborTable(  &numberOfNeighborParticles
					,&neighborTable
					,parameter.radiusOfLaplacianForPressure_ratio
					);

	for( iParticle=0; iParticle < particle.totalNumber; iParticle++){

		if( particle.flagOfBoundaryCondition[iParticle] != INNER_PARTICLE ) continue;

		particle.coefficientMatrixOfPressure[iParticle][0]=0.0;

		for(iNeigh=0; iNeigh <  numberOfNeighborParticles[iParticle]; iNeigh++){
			jParticle = neighborTable[iParticle][iNeigh];

			if( particle.flagOfBoundaryCondition[jParticle] == GHOST_OR_DUMMY ){
				particle.coefficientMatrixOfPressure[iParticle][iNeigh+1]=0.0;
			}else {
				distanceIJ = DISTANCE_calculateDistanceBetweenParticles( particle.position, iParticle, jParticle  );

				averageMassDensity = MATHE_average( physicalProperty.massDensity[particle.type[iParticle]]
							,physicalProperty.massDensity[particle.type[jParticle]] ); 

				weightIJ = WEIGHT_calculateWeightFunction( distanceIJ, parameter.radiusOfLaplacianForPressure);

				coefficientOfParticleJ  = 2.0 * NumberOfDimensions * weightIJ /( n0 * lambda  * averageMassDensity );

				particle.coefficientMatrixOfPressure[iParticle][iNeigh+1]  = - coefficientOfParticleJ;
    
				if(parameter.flagOfTanakaAndMasunagaModel==ON){
					particle.coefficientMatrixOfPressure[iParticle][0]        +=   parameter.valueOfC*coefficientOfParticleJ;
				}else{
					particle.coefficientMatrixOfPressure[iParticle][0]        +=   coefficientOfParticleJ; 
				}

				//Ž©—R•\–Ê—±Žq‚Ì”ñ‘ÎŠp€‚ÍŒvŽZ‚É“ü‚ê‚È‚¢
				if(particle.flagOfBoundaryCondition[jParticle]==SURFACE_PARTICLE)
					particle.coefficientMatrixOfPressure[iParticle][iNeigh + 1] = 0.0;

			}
		}

		particle.coefficientMatrixOfPressure[iParticle][0] += physicalProperty.compressibility[particle.type[iParticle]]/( dt_squared );

	}

}

void
PRESSURE_setSurfacePressureZero(void) {
	int iParticle;
	for (iParticle = 0; iParticle < particle.totalNumber; iParticle++) {
		if (particle.flagOfBoundaryCondition[iParticle] == SURFACE_PARTICLE)
			particle.pressure[iParticle] = 0.0;
	}
}

void
PRESSURE_setMinusPressureZero( void ){

  int iParticle;

  for(iParticle=0; iParticle < particle.totalNumber; iParticle++){
    
    if( particle.pressure[iParticle] < 0.0 ){
		//if (parameter.flagOfNegativePressure == OFF || particle.flagOfBoundaryCondition[iParticle] == SURFACE_PARTICLE)
		if(OUTFLOW_allowedJudgingFreeSurface(iParticle))
			particle.pressure[iParticle] = 0.0;
    }

  }

}



void
PRESSURE_setMinimumPressureAmongNeighbors( void ){

  int    iParticle,jParticle;
  int    iNeigh;

  double minimumPressure;
  double distanceIJ_squared;


  int    *numberOfNeighborParticles;
  int    **neighborTable;


  NEIGH_selectNeighborTable(   &numberOfNeighborParticles
							  ,&neighborTable
							  ,parameter.radiusOfGradient_ratio
							  );


  for(iParticle=0; iParticle < particle.totalNumber; iParticle++){

    if( particle.flagOfBoundaryCondition[iParticle] == GHOST_OR_DUMMY ) continue;
    if( particle.type[iParticle]        == parameter.wallType ) continue;

    minimumPressure = particle.pressure[iParticle];

    for(iNeigh=0; iNeigh <  numberOfNeighborParticles[iParticle]; iNeigh++){
      jParticle = neighborTable[iParticle][iNeigh];

      if(particle.type[jParticle] == GHOST              )continue;
      if(particle.type[jParticle] == parameter.dummyWallType)continue;

	  distanceIJ_squared = DISTANCE_calculateSquaredDistance(particle.position, iParticle, jParticle);
	  if(distanceIJ_squared > parameter.radiusOfGradient_squared) continue;

      if(particle.pressure[jParticle] < minimumPressure) minimumPressure = particle.pressure[jParticle];

    }

    particle.minPressureAmongNeighbors[iParticle] = minimumPressure;

  }
}



int
PRESSURE_checkNecessityOfPressureCalculation( void ){

  int iParticle;
  int flagOfPressureCalculation = OFF;

  for(iParticle=0; iParticle < particle.totalNumber; iParticle++){

    if( particle.flagOfBoundaryCondition[iParticle] == INNER_PARTICLE){
      flagOfPressureCalculation = ON;
	  return flagOfPressureCalculation;
    }

  }

  return flagOfPressureCalculation;

}





void
PRESSURE_setDirichletBoundaryCondition( void ){

  if(parameter.flagOfTanakaAndMasunagaModel==ON){
     PRESSURE_findParticlesBeingOnFreeSurfaceTheNumberOfNeighParticleBased();
  }else{
     PRESSURE_findParticlesBeingOnFreeSurface(); 
  }

  PRESSURE_checkThatDirichletBoundaryConditionIsConnected();

}





void
PRESSURE_findParticlesBeingOnFreeSurface( void ){

  int    iParticle;
  double n0 = parameter.nZeroOfParticleNumberDensity;
  
  for(iParticle=0; iParticle < particle.totalNumber; iParticle++){

    if( (particle.type[iParticle] == GHOST) || (particle.type[iParticle] == parameter.dummyWallType)){

      particle.flagOfBoundaryCondition[iParticle] = GHOST_OR_DUMMY;

    }else if( (particle.particleNumberDensity[iParticle] / n0 ) < parameter.thresholdOfParticleNumberDensity_ratio  && OUTFLOW_allowedJudgingFreeSurface(iParticle)){

      particle.flagOfBoundaryCondition[iParticle] = SURFACE_PARTICLE;

    }else{

      particle.flagOfBoundaryCondition[iParticle] = INNER_PARTICLE;

    }
  }

}

void
PRESSURE_findParticlesBeingOnFreeSurfaceTheNumberOfNeighParticleBased( void ){
    
    int    iParticle;
    double N0 = parameter.nZeroOfnumberOfNeighborParticles;
    double numberOfNeighborParticles_ratio=0.80;
    
    int    *numberOfNeighborParticles;
    int    **neighborTable;
    
    NEIGH_selectNeighborTable(   &numberOfNeighborParticles
							  ,&neighborTable
							  ,parameter.radiusOfParticleNumberDensity_ratio
							  );
    
    
    for(iParticle=0; iParticle < particle.totalNumber; iParticle++){
        
        if( (particle.type[iParticle] == GHOST) || (particle.type[iParticle] == parameter.dummyWallType)){
            
            particle.flagOfBoundaryCondition[iParticle] = GHOST_OR_DUMMY;
            
        }else if( numberOfNeighborParticles[iParticle] < N0*numberOfNeighborParticles_ratio && OUTFLOW_allowedJudgingFreeSurface(iParticle)){
            
			particle.flagOfBoundaryCondition[iParticle] = SURFACE_PARTICLE;
            
        }else{
            
            particle.flagOfBoundaryCondition[iParticle] = INNER_PARTICLE;
            
        }
    }
    
}



void
PRESSURE_checkThatDirichletBoundaryConditionIsConnected( void ){

  int iParticle, jParticle, iNeigh;
  int count;
  int *checkArray;
  char buf[256];
  int iLoop=0;

  int    *numberOfNeighborParticles = particle.numberOfNeighborParticles_large;
  int    **neighborTable            = particle.neighborTable_large;
  
  checkArray = MEMORY_allocateMemoryFor1dimIntArray( particle.totalNumber,
	                   "checkArray [in function of checkThatDirichletBoundariesAreConnected()]");

  for(iParticle=0; iParticle < particle.totalNumber; iParticle++){

    if (particle.flagOfBoundaryCondition[iParticle]== GHOST_OR_DUMMY ){
      
      checkArray[iParticle]= GHOST_OR_DUMMY;

    }else if (particle.flagOfBoundaryCondition[iParticle] == SURFACE_PARTICLE){

      checkArray[iParticle] = CONNECTED;

    }else if (particle.flagOfBoundaryCondition[iParticle] == INNER_PARTICLE){

      checkArray[iParticle] = NOT_CONNECTED;

    }else{
      
      sprintf(buf,"[in PRESSURE_checkThatDirichletBoundariesAreConnected]\n particle.surface is not adequate.\n particle.flagOfBoundaryCondition[%d] = %d", iParticle, particle.flagOfBoundaryCondition[iParticle]);

      OTHER_endProgram("buf");

    }

  }


  do {
    count=0;

    for(iParticle=0; iParticle < particle.totalNumber; iParticle++){

		if(checkArray[iParticle] == CONNECTED){
	
			for(iNeigh=0; iNeigh <  numberOfNeighborParticles[iParticle]; iNeigh++){
				jParticle = neighborTable[iParticle][iNeigh];
	  
				if(checkArray[jParticle]== NOT_CONNECTED){
					checkArray[jParticle] = CONNECTED;
				}
			}

			checkArray[iParticle] = CHECKED;
			count ++;

		}
    }

	iLoop++;

  } while (count > 0);


  PRESSURE_ifBoundaryConditionIsNotConnected_increaseDiagonalElementsOfCoefficientMatrixOfPressureToSolveSimultaneousEquationsWithoutBoundaryCondition( checkArray );

  free(checkArray);

}


void
PRESSURE_ifBoundaryConditionIsNotConnected_increaseDiagonalElementsOfCoefficientMatrixOfPressureToSolveSimultaneousEquationsWithoutBoundaryCondition( int *checkArray ){

  int iParticle;

  double increaseRate = parameter.increaseRateOfDiagonalTermInCoefficientMatrixOfPressure;


  for( iParticle=0; iParticle < particle.totalNumber; iParticle++){
    if(checkArray[iParticle] == NOT_CONNECTED ){

      particle.coefficientMatrixOfPressure[iParticle][0] = increaseRate * particle.coefficientMatrixOfPressure[iParticle][0];

      //PRESSURE_displayWarnigMessageForNoDirichletBoundaryCondition( iParticle );
    }

  }


}



void
PRESSURE_displayWarnigMessageForNoDirichletBoundaryCondition( int iParticle ){

  fprintf(FpForLog,"WARNIG: No Dirichlet boundary condition.");
  fprintf(FpForLog,"  type[%d] = %d",iParticle, particle.type[iParticle]);

  fprintf(FpForLog,"  position[XDIM][%d] = %lf",iParticle, particle.position[XDIM][iParticle]);
  fprintf(FpForLog,"  position[YDIM][%d] = %lf",iParticle, particle.position[YDIM][iParticle]);

  if(NumberOfDimensions == 3){
    fprintf(FpForLog,"  position[ZDIM][%d] = %lf",iParticle, particle.position[ZDIM][iParticle]);
  }

  fprintf(FpForLog,"\n");

}

void
PRESSURE_correctPressure( void ){
    
    int iParticle;
    double kinematicViscosity = physicalProperty.kinematicViscosity;
    
    double n0         = parameter.nZeroOfParticleNumberDensity;
    double dt_squared = timer.dt * timer.dt;
    
    for(iParticle=0; iParticle < particle.totalNumber; iParticle++){
            
            particle.pressure[iParticle] += (1.0/2.0)*kinematicViscosity*timer.dt*( particle.particleNumberDensity[iParticle] - n0 )/( dt_squared * n0 );
            
    }
}

void
PRESSURE_updateAveragePressure(void) {
	int iParticle;
	int iX, iY, iZ;
	int flag;
	int outflowsize;
	int i;
	int checkbool;

	if (parameter.flagOfAveragePressureInEachBucket == OFF || parameter.timeToUpdateAveragePressure > timer.simulationTime)return;

	for (i = 0; i < 256; i++) {
		parameter.OutflowBucketIndex[i] = -1;
	}
	parameter.OutflowAveragePressure = 0;
	parameter.OutflowBucketLength = 0;
	
	for (iParticle = 0; iParticle < particle.totalNumber; iParticle++) {
		if (particle.flagOfBoundaryCondition[iParticle] == GHOST_OR_DUMMY)continue;
		if (particle.type[iParticle] == GHOST)continue;

		BUCKET_findPressureBucketWhereParticleIsStored(&iX, &iY, &iZ, iParticle, particle.position);
		if (iX == -1) {
			particle.bucketPressure[iParticle] = 0;
			continue;
		}
		
		
		if (particle.flagOfBoundaryCondition[iParticle] == SURFACE_PARTICLE) {
			//domain.pressureBucket[iX][iY][iZ].count++;
			if(NumberOfDimensions)domain.countBucket[iX + iY * domain.pressureBucketNumber[XDIM]]++ ;
			else domain.countBucket[iX + iY * domain.pressureBucketNumber[XDIM] + iZ * domain.pressureBucketNumber[XDIM] * domain.pressureBucketNumber[YDIM]]++;
			
			if (parameter.OutflowBucketLength >= 256)break;
			checkbool = 0;
			for (i = 0; i < parameter.OutflowBucketLength; i++) {
				if (parameter.OutflowBucketIndex[i] == iX + iY * domain.pressureBucketNumber[XDIM]) {
					checkbool = 1;
					break;
				}
			}
			if (checkbool == 0) {
				parameter.OutflowBucketIndex[parameter.OutflowBucketLength] = iX + iY * domain.pressureBucketNumber[XDIM];
				parameter.OutflowBucketLength++;
			}

		}
		else if (particle.flagOfBoundaryCondition[iParticle] == INNER_PARTICLE) {
			//domain.pressureBucket[iX][iY][iZ].pressure =
				//(domain.pressureBucket[iX][iY][iZ].pressure * domain.pressureBucket[iX][iY][iZ].count + particle.pressure[iParticle]) / (domain.pressureBucket[iX][iY][iZ].count + 1);
			//domain.pressureBucket[iX][iY][iZ].count++;
			if (NumberOfDimensions)domain.pressureBucket[iX + iY * domain.pressureBucketNumber[XDIM]] 
				= (domain.pressureBucket[iX + iY * domain.pressureBucketNumber[XDIM]] * domain.countBucket[iX + iY * domain.pressureBucketNumber[XDIM]] + particle.pressure[iParticle])
				/ (domain.countBucket[iX + iY * domain.pressureBucketNumber[XDIM]] + 1);
			else domain.pressureBucket[iX + iY * domain.pressureBucketNumber[XDIM] + iZ * domain.pressureBucketNumber[XDIM] * domain.pressureBucketNumber[YDIM]]
				= (domain.pressureBucket[iX + iY * domain.pressureBucketNumber[XDIM] + iZ * domain.pressureBucketNumber[XDIM] * domain.pressureBucketNumber[YDIM]] 
					* domain.countBucket[iX + iY * domain.pressureBucketNumber[XDIM] + iZ * domain.pressureBucketNumber[XDIM] * domain.pressureBucketNumber[YDIM]] + particle.pressure[iParticle])
				/ (domain.countBucket[iX + iY * domain.pressureBucketNumber[XDIM] + iZ * domain.pressureBucketNumber[XDIM] * domain.pressureBucketNumber[YDIM]] + 1);
			if (NumberOfDimensions)domain.countBucket[iX + iY * domain.pressureBucketNumber[XDIM]]++;
			else domain.countBucket[iX + iY * domain.pressureBucketNumber[XDIM] + iZ * domain.pressureBucketNumber[XDIM] * domain.pressureBucketNumber[YDIM]];

		}
		

	}

	for (i = 0; i < parameter.OutflowBucketLength; i++) {
		parameter.OutflowAveragePressure += domain.pressureBucket[parameter.OutflowBucketIndex[i]];
	}
	parameter.OutflowAveragePressure /= parameter.OutflowBucketLength;

	for (iParticle = 0; iParticle < particle.totalNumber; iParticle++) {
		if (particle.flagOfBoundaryCondition[iParticle] == GHOST_OR_DUMMY)continue;
		if (particle.type[iParticle] == GHOST)continue;

		particle.bucketPressure[iParticle] =
			(NumberOfDimensions == 2) ? domain.pressureBucket[iX + iY * domain.pressureBucketNumber[XDIM]] - parameter.OutflowAveragePressure
			: domain.pressureBucket[iX + iY * domain.pressureBucketNumber[XDIM] + iZ * domain.pressureBucketNumber[XDIM] * domain.pressureBucketNumber[YDIM]] - parameter.OutflowAveragePressure;

	}

}