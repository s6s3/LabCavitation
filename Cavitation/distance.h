void
DISTANCE_transformUnitOfDistanceFromRatioToMeter( void );


double
DISTANCE_calculateDistanceBetweenParticles(
										   double        **position
										   ,int            iParticle
										   ,int            jParticle
										   );


double
DISTANCE_calculateSquaredDistance(
								  double       **position
								  ,int            iParticle
								  ,int            jParticle
								  );

double
DISTANCE_calculateDistanceBetweenTwoPositions( double *position1, double *position2 );
