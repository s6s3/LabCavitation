#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "define.h"
#include "extern.h"
#include "struct.h"
#include "file.h"
#include "copy.h"
#include "bucket.h"
#include "domain.h"
#include "neigh.h"
#include "density.h"
#include "timer.h"
#include "other.h"
#include "gravity.h"
#include "convection.h"
#include "collision.h"
#include "pressure.h"
#include "gradient.h"
#include "viscosity.h"
#include "sort.h"
#include "object.h"
#include "forcedMotion.h"
#include "memory.h"
#include "nzero.h"
#include "distance.h"
#include "maxmin.h"
#include "init.h"
#include "outflow.h"
#include "inflow.h"
#include "cavitation.h"

//declaration extern.h and struct.h
int     NumberOfDimensions;
int     FlagOfHighSpeedCalculation;
FILE   *FpForLog;


/*--- for cosinTransform.c ---*/
double AverageOfSamplingData;
double FinalTime;
int    FlagOfWindowFunction;
double Phase_current;
double Amplitude_current;
double WaveLength_current;
double AverageWaterLevel_current;

double SolutionOfLeastSquareMethod[20];
double MinimumErrorOfApproximatedEquation;
int    FlagOfLinearApproximation;

structParticle          particle;
structDomain            domain;
structParameter         parameter;
structTimer             timer;
structPhysicalProperty  physicalProperty;
structObject            objectOfRigidBody;
structObject            objectOfWall;
structForcedMotion      forcedMotionOfRigidBody;
structForcedMotion      forcedMotionOfWall;
stack                   ghostStack;


void
INIT_initializeParameters( int argumentCount,char **argumentVector ){

	timer.clockAtStartingTime = clock();

	timer.iTimeStep = 0;
	timer.iTimeStep_copy = -1;

	ghostStack = { 0 };
	STACK_init(&ghostStack);

	FILE_readInputFiles( argumentCount, argumentVector );

	parameter.flagOfBiCG = OFF;

	INFLOW_setInflowLevel();

	INFLOW_setParticleType();

	physicalProperty.saturatedVaporPressure = CAVITATION_getSaturatedVaporPressure(physicalProperty.temperature - 273);

	/*--- Memory of structParticle is allocated in the above FILE_readInputFiles(). ---*/

	COPY_storeInitialParticleProperty();

	DISTANCE_transformUnitOfDistanceFromRatioToMeter();

	MAXMIN_setMaxMinOfRadii();

	TIMER_initializeTimeFunctions();

	NZERO_calculateNZeroAndLambda( particle.averageDistance );
 

	if (parameter.flagOfOutflowBoundaryCondition == ON) {
		if (NumberOfDimensions == 2)parameter.depthOfOutflow = 1.0;
	}

	//INFLOW_setInflowLevel();

	NEIGH_initializeNeighborTable();

	MEMORY_allocateMemoryForCoefficientMatrixOfPressure();

	if(parameter.flagOfForcedMotionOfRigidBody == ON ){
		FORCEDMOTION_initializeForcedMotion(  &forcedMotionOfRigidBody
							,parameter.nameOfSamplingDataFileForForcedMotionOfRigidBody, parameter.typeNumberOfRigidParticle_forForcedMotion );

		/*
		FORCEDMOTION_initializeForcedMotion(  &forcedMotionOfRigidBody
							,parameter.nameOfSamplingDataFileForForcedMotionOfRigidBody );
		*/
		FORCEDMOTION_mainFunctionOfForcedMotionOfRigidBody( &forcedMotionOfRigidBody );
	}

	COPY_setInitialCoodinatesAndVelocities();

	parameter.initTotalOutflow = OUTFLOW_getCurrentNumberOfFluidParticles();
    

	DOMAIN_initializeDomain();

	BUCKET_initializeBucket();

	NEIGH_setNeighborTable( particle.position );//position = 

	DENSITY_calculateParticleNumberDensity( particle.position );

	FILE_countNumberOfParticlesEachType();

	FILE_displayNumberOfParticles();

	fflush(FpForLog);

	FILE_writeProfFile();

	NEIGH_setNeighborTable( particle.position );//position = 

	timer.clockAtBeginningOfMainLoop = clock();
 
	timer.iTimeStep_copy++;

}

