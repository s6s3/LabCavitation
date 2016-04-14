#pragma once
void
CAVITATION_calculateBubble(void);
double
CAVITATION_getBodyPressure();

double
CAVITATION_getSaturatedVaporPressure(double T);

double
CAVITATION_calculateCosine(double distanceIJ, int iParticle, int jParticle);

double
CAVITATION_setInfluenceRadius(void);

void
CAVITATION_setBetaZero(void);

double
CAVITATION_calculateDifferenceRadius(int iParticle, double averagePressure);

double
CAVITATION_calculateVolume(double diameter);

double
CAVITATION_calculateRisingVelocity(double diameter);

void
CAVITATION_calculateBubbleRising(void);

void
CAVITATION_surfaceJudgement(void);

void
CAVITATION_calculateBuoyancy(void);