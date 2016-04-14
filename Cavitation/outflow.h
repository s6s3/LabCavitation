#pragma once


double
OUTFLOW_judgeBounds(double x, double y, double z);

int
OUTFLOW_allowedJudgingFreeSurface(int iParticle);

int
OUTFLOW_getCurrentNumberOfFluidParticles();

double
OUTFLOW_getL(int iParticle, int jParticle);

double
OUTFLOW_getVelocityOfUpstream(int iParticle);

int
OUTFLOW_deleteOutflowParticles();

void
OUTFLOW_correctOutflowVelocity();

void OUTFLOW_correctYamakawaOutflowVelocity();

void OUTFLOW_correctYamakawaOutflowVelocity2();

void OUTFLOW_correctShibataOutflowVelocity();