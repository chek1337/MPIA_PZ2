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
	MatrixGenerator2(matrixB);
	//MatrixOutDense(matrixA);
	//MatrixOutDense(matrixB);
	auto startTime = steady_clock::now();
	MatrixMatrixMultiplication(matrixA, matrixB, matrixC);
	auto endTime = steady_clock::now();
	//MatrixOutDense(matrixC);
	printf_s("AB=C time %d\n", duration_cast<milliseconds>(endTime - startTime).count());
	//Ne sovsem ponyal nado li vreamya norm'i zaemryat'
	printf("Matrix Norm = %lf\n", MatrixNorm(matrixC)); 

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

	printf_s("\n%f", VectorNorm(x));
}

void Task7()
{
	double* matrixA, * matrixB, * matrixC;
	matrixA = new double[n * n];
	matrixB = new double[n * n];
	matrixC = new double[n * n];
	MatrixGenerator(matrixA);
	MatrixGenerator2(matrixB);
	//MatrixOutDense(matrixA);
	//MatrixOutDense(matrixB);
	MatrixTransposition(matrixB);
	//MatrixOutDense(matrixB);
	auto startTime = steady_clock::now();
	MatrixMatrixMultiplicationByLines(matrixA, matrixB, matrixC);
	auto endTime = steady_clock::now();
	///MatrixOutDense(matrixC);
	printf_s("AB^T=C time %d", duration_cast<milliseconds>(endTime - startTime).count());

	
}


double ScalarProduct(double* a, double* b) {
	double sum = 0;
#pragma omp parallel for reduction(+ : sum) num_threads(4)
	for (int i = 0; i < n; i++)
	{
		sum += a[i] * b[i];
	}
	return sum;
}


void MatrixMatrixMultiplication(double* matrixA, double* matrixB, double* matrixC)
{
	double sum = 0;
#pragma omp parallel for reduction(+ : sum) num_threads(2)
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			sum = 0;
			for (int k = 0; k < n; k++)
			{
				sum += matrixA[i*n + k] * matrixB[k*n+j];
			}
			matrixC[i * n + j] = sum;
		}
	}
}

double MatrixNorm(double* matrixIn)
{
	double norm = 0;
#pragma omp parallel for reduction(+ : norm) num_threads(1)
	for (int i = 0; i < n * n; i++)
	{
		norm += matrixIn[i] * matrixIn[i];
	}
	return sqrt(norm);
}

double VectorNorm(double* a) {
	double norm = 0;
#pragma omp parallel for reduction(+ : norm) num_threads(4)
	for (int i = 0; i < n; i++)
	{
		norm += a[i] * a[i];
	}
	return sqrt(norm);
}


void CalculateY(double* matrixU, double* vectorX, double* vectorB) {

	double sum = 0;
#pragma omp parallel for reduction(+ : sum) num_threads(1)
	for (int i = n - 1; i >= 0; i--)
	{
		double sum = 0;
		for (int j = n - 1; j > i; j--)
		{
			sum += matrixU[i*n+j] * vectorX[j];
		}
		vectorX[i] = (vectorB[i] - sum) / matrixU[i*n + i];
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
				sum += matrixA[i*n + k] * matrixB[j*n + k];
			}
			matrixC[i*n+j] = sum;
		}
	}
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


void MatrixVectorMultiplication(double* matrixA, double* vectorX, double* vectorB)
{
	for (int i = 0; i < n; i++)
	{
		double sum = 0;
		for (int j = 0; j < n; j++)
		{
			sum += matrixA[i*n+j] * vectorX[j];
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


void MatrixUGenerator(double* outMatrix) {   //Формула для вычисления a(ij) = (i+j) % 10 + 1
	for (int i = 0; i < n; i++)						//Диагональный элемент i'ой строки равен сумме всех внедиагональных + 1
	{
		double sum = 0;
		for (int j = 0; j < i; j++)
		{
			outMatrix[i*n+j] = (double) 0.0;
		}

		for (int j = i + 1; j < n; j++)
		{
			outMatrix[i*n+j] = ((i+j) % 10) + 1;
			sum += ((i + j) % 10) + 1;
		}	
		outMatrix[i*n + i] = sum + 1; 
	}
}


void MatrixTransposition(double* matrix)
{
	double tmp = 0;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < i; j++)
		{
			tmp = matrix[i*n+j];
			matrix[i*n+j] = matrix[j*n + i];
			matrix[j*n + i] = tmp;
		}
	}
}


