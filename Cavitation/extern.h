#ifndef INCLUDE_454682
#define INCLUDE_454682

extern int     NumberOfDimensions;
extern int     FlagOfHighSpeedCalculation;
extern FILE   *FpForLog;


/*--- for cosinTransform.c ---*/
extern double AverageOfSamplingData;
extern double FinalTime;
extern int    FlagOfWindowFunction;
extern double Phase_current;
extern double Amplitude_current;
extern double WaveLength_current;
extern double AverageWaterLevel_current;

extern double SolutionOfLeastSquareMethod[20];
extern double MinimumErrorOfApproximatedEquation;
extern int    FlagOfLinearApproximation;



#endif