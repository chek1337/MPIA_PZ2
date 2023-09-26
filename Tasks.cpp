#include "stdafx.h"
#include "headers.h"

void Task1()
{
	double* a, * b;
	a = new double[n];
	b = new double[n];
	VectorGenerator(a);
	VectorGenerator(b);

	auto startTime = steady_clock::now();
	double SP = ScalarProduct(a, b);
	auto endTime = steady_clock::now();
	printf("%lf\n", SP);
	printf_s("a*b time %d", duration_cast<milliseconds>(endTime - startTime).count()); // перепроверить спецификатор %d
}

void Task2()
{
	double* matrixA, * matrixB, * matrixC;
	matrixA = new double[n * n];
	matrixB = new double[n * n];
	matrixC = new double[n * n];
	MatrixGenerator(matrixA);
	MatrixGenerator(matrixB);
	MatrixGenerator(matrixC);
	//MatrixOutDense(matrixA);
	//MatrixOutDense(matrixB);
	auto startTime = steady_clock::now();
	MatrixMatrixMultiplication(matrixA, matrixB, matrixC);
	auto endTime = steady_clock::now();
	//MatrixOutDense(matrixC);
	printf_s("AB=C time %d\n", duration_cast<milliseconds>(endTime - startTime).count());
	//printf("Matrix Norm = %lf\n", MatrixNorm(matrixC));

}

void Task3()
{
	double* x, * b, * y;
	b = new double[n];
	x = new double[n];
	y = new double[n];
	double* matrixU;
	matrixU = new double[n * n];
	MatrixUGenerator(matrixU);
	//MatrixOutDense(matrixU);
	VectorGenerator2(x);
	//VectorOutput(x);
	MatrixVectorMultiplication(matrixU, x, b);
	//VectorOutput(b);
	auto startTime = steady_clock::now();
	CalculateY(matrixU, y, b);
	auto endTime = steady_clock::now();
	//VectorOutput(y);
	printf_s("Uy=b time %d", duration_cast<milliseconds>(endTime - startTime).count());
}

void Task7()
{
	double* matrixA, * matrixB, * matrixC;
	matrixA = new double[n * n];
	matrixB = new double[n * n];
	matrixC = new double[n * n];
	MatrixGenerator(matrixA);
	MatrixGenerator(matrixB);
	//MatrixOutDense(matrixA);
	//MatrixOutDense(matrixB);
	MatrixTransposition(matrixB);
	//MatrixOutDense(matrixB);
	auto startTime = steady_clock::now();
	MatrixMatrixMultiplicationByLines(matrixA, matrixB, matrixC);
	auto endTime = steady_clock::now();
	//MatrixOutDense(matrixC);
	printf_s("AB^T=C time %d", duration_cast<milliseconds>(endTime - startTime).count());
}



void VectorGenerator(double* outArray) {
	for (int i = 0; i < n; i++)
	{
		outArray[i] = (i % 10);
	}
}

void VectorGenerator2(double* outArray) {
	for (int i = 0; i < n; i++)
	{
		outArray[i] = (i % 10) + 1;
	}
}


void MatrixGenerator(double* outArray) {
	for (int i = 0; i < n*n; i++)
	{
		outArray[i] = (i % 10);
	}
}

void MatrixGenerator2(double* outArray) {
	for (int i = 0; i < n * n; i++)
	{
		outArray[i] = (i % 10)+1;
	}
}

double ScalarProduct(double* a, double* b) {
	double s = 0;
#pragma omp parallel for reduction(+ : s) num_threads(4)
	for (int i = 0; i < n; i++)
	{
		s += a[i] * b[i];
	}
	return s;
}

void MatrixMatrixMultiplication(double* matrixA, double* matrixB, double* matrixC)
{
	double sum = 0;
#pragma omp parallel for reduction(+ : sum) num_threads(4)
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			sum = 0;
			for (int k = 0; k < n; k++)
			{
				sum += matrixA[num(i, k)] * matrixB[num(k, j)];
			} 
			matrixC[num(i, j)] = sum;
		}	
	}
}

void MatrixVectorMultiplication(double* matrixA, double* vectorX, double* vectorB)
{
	for (int i = 0; i < n; i++)
	{
		double sum = 0;
		for (int j = 0; j < n; j++)
		{
			sum += matrixA[num(i, j)] * vectorX[j];
		}
		vectorB[i] = sum;
	}
}


void MatrixOutDense(double* matrixOut)
{
	for (int i = 0; i < n*n; i++)
	{
		if ((i % n ) == 0)
			cout << endl;
		printf_s("%.1lf\t", matrixOut[i]);
	}
	cout << endl << endl;
}

void VectorOutput(double* x)
{
	for (int i = 0; i < n; i++)
	{
		printf_s("%.1lf\n", x[i]);
	}
	cout << endl;
}

double MatrixNorm(double *matrixIn)
{
	double norm=0;
	for (int i = 0; i < n*n; i++)
	{
		norm += matrixIn[i] * matrixIn[i];
	}
	return sqrt(norm);
}

void MatrixUGenerator(double* outMatrix) {   //Формула для вычисления a(ij) = (i+j) % 10 + 1
	for (int i = 0; i < n; i++)						//Диагональный элемент i'ой строки равен сумме всех внедиагональных + 1
	{
		double sum = 0;
		for (int j = 0; j < i; j++)
		{
			outMatrix[num(i, j)] = (double) 0.0;
		}

		for (int j = i + 1; j < n; j++)
		{
			outMatrix[num(i,j)] = ((i+j) % 10) + 1;
			sum += ((i + j) % 10) + 1;
		}	
		outMatrix[num(i, i)] = sum + 1; 
	}
}

void CalculateY(double *matrixU, double *vectorX, double* vectorB){

	double sum=0, tmp=0;
#pragma omp parallel for reduction(+ : sum, tmp) num_threads(1)
	for (int i = n-1; i >= 0; i--)
	{
		double sum = 0;
		for (int j = n-1; j > i; j--)
		{
			sum += matrixU[num(i, j)]* vectorX[j];
		}
		tmp = (vectorB[i] - sum) / matrixU[num(i, i)];
		vectorX[i] = tmp;
	}
}

void MatrixTransposition(double* matrix)
{
	double tmp = 0;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < i; j++)
		{
			tmp = matrix[num(i, j)];
			matrix[num(i, j)] = matrix[num(j, i)];
			matrix[num(j, i)] = tmp;
		}
	}
}


void MatrixMatrixMultiplicationByLines(double* matrixA, double* matrixB, double* matrixC)
{
	double sum = 0;
#pragma omp parallel for reduction(+ : sum) num_threads(4)
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			sum = 0;
			for (int k = 0; k < n; k++)
			{
				sum += matrixA[num(i, k)] * matrixB[num(j,k)];
			}
			matrixC[num(i, j)] = sum;
		}
	}
}

