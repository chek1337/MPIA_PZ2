#pragma once
#include <chrono>
#define num(i,j) (j+i*n) //Конвертирование двумерной адрессации в одинарную 
const int n = 8;

using namespace std;
using namespace std::chrono;

void VectorGenerator(double* outArray);

void VectorGenerator2(double* outArray);

void VectorOutput(double* x);

double ScalarProduct(double* a, double* b);


void MatrixGenerator(double* outArray);

void MatrixGenerator2(double* outArray);

void MatrixUGenerator(double* outMatrix);

void MatrixMatrixMultiplication(double* matrixA, double* matrixB, double* matrixC);

void MatrixMatrixMultiplicationByLines(double* matrixA, double* matrixB, double* matrixC);

void MatrixVectorMultiplication(double* matrixA, double* vectorX, double* vectorB);

void MatrixOutDense(double* matrixOut);

double MatrixNorm(double* matrixIn);

void MatrixTransposition(double* matrix);

void CalculateY(double* matrixU, double* vectorX, double* vectorB);



void Task1();

void Task2();

void Task3();

void Task7();
