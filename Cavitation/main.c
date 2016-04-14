/*****************************************************************************

 MPS-SW-MAIN
 
 Version:          2.0
 
 Completion date:  February 14, 2008
 
 Copyright:        Kazuya SHIBATA and Seiichi KOSHIZUKA
 hold the copyright of this code.
 (2013-2015 Added Bubble Model for Stirred Tank Analysis by Shogo KAITO)
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include <float.h>

#include "define.h"
#include "extern.h"
#include "struct.h"
#include "file.h"
#include "memory.h"
#include "copy.h"
#include "bucket.h"
#include "domain.h"
#include "distance.h"
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
#include "init.h"
#include "finalize.h"


#include "bubble.h"
#include "buoyancy.h"
#include "outflow.h"
#include "inflow.h"
#include "cavitation.h"

//github test 4/14

int
main( int argumentCount, char **argumentVector ){
	_controlfp(0, _EM_ZERODIVIDE | _EM_UNDERFLOW | _EM_OVERFLOW | _EM_INVALID);//0èúéZåüèo

	INIT_initializeParameters( argumentCount, argumentVector );
    
    for( timer.iTimeStep= 0; timer.iTimeStep < timer.finishTimeStep;  timer.iTimeStep++ ){
		TIMER_setDtAutomatically();
        
        TIMER_putTimeForwardByDt();
        
        TIMER_displayStateOfTimeStep_atAppropriateTime();
		
		INFLOW_setVelocity();

        GRAVITY_calculateGravity();
        
        VISCOSITY_calculateViscosity();
        
		//OUTFLOW_correctOutflowVelocity();

		//OUTFLOW_correctYamakawaOutflowVelocity2();
		OUTFLOW_correctShibataOutflowVelocity();

        CONVECTION_moveParticles( particle.position, particle.velocity );//err?
        
		OUTFLOW_deleteOutflowParticles();
        
        if(parameter.flagOfForcedMotionOfRigidBody == ON ){
            
            FORCEDMOTION_mainFunctionOfForcedMotionOfRigidBody( &forcedMotionOfRigidBody );
            
        }
        
		INFLOW_changeInflowParticles();

        COLLISION_calculateCollisionBetweenParticles();
        
        OTHER_checkThatParticlesAreNotAllGhost();
        
        PRESSURE_calculatePressure();//ERROR!
        
        GRADIENT_correctParticleVelocityAndPositionUsingPressureGradient();
        
        if(parameter.flagOfForcedMotionOfRigidBody == ON ){
            
            FORCEDMOTION_mainFunctionOfForcedMotionOfRigidBody( &forcedMotionOfRigidBody );
            
        }
        
        COPY_updateParticleProperty();
        
        if(parameter.flagOfBubbleCalculation == ON){
            
            //Bubble_calculateBubble();
			CAVITATION_calculateBubble();
        }
        
        FILE_writeCalculationResultInFile();
        
        COLLISION_calculateCollisionBetweenParticles();
        
        NEIGH_setNeighborTable(particle.position);
        
        
        if( YES == TIMER_checkWhetherItIsTimeToFinishProgram() ) break;
        
    }
    
    FINALIZE_finalizeProgram();
    
    return(EXIT_SUCCESS);
    
}



