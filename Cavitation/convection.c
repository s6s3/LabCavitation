#include <stdio.h>
#include <math.h>

#include "define.h"
#include "struct.h"
#include "extern.h"
#include "convection.h"


void
CONVECTION_moveParticles( double **position, double **velocity){
	int iParticle;
	int iDim;
  
	for( iParticle=0; iParticle < particle.totalNumber; iParticle++){
		if( particle.type[iParticle] == GHOST ) continue;

		for( iDim=0; iDim < NumberOfDimensions; iDim++){
			position[iDim][iParticle] +=  (velocity[iDim][iParticle] * timer.dt);
		}

	}

}


