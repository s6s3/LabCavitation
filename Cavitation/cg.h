#pragma once


int CG_solver(int totalNumber, double** matrix, double* solution, double* sourceTerm, int maxIteration, int minIteration, double error);