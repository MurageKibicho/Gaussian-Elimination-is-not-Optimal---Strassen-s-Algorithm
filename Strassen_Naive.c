#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
#include <assert.h>
#define DataType float
//clear && gcc Strassen_Naive.c -lm -o m.o && ./m.o
DataType **SquareMatrix2D_Allocate(int n)
{
	DataType **matrix = malloc(n * sizeof(DataType*));
	for(int i = 0; i < n; i++)
	{
		matrix[i] = calloc(n, sizeof(DataType));
	}
	return matrix;	
}

void SquareMatrix2D_Free(int n, DataType **matrix)
{
	for(int i = 0; i < n; i++)
	{
		free(matrix[i]);
	}
	free(matrix);
}

void SquareMatrix2D_Random(int n, DataType **matrix, DataType bound)
{
	//Generates a random matrix between -bound and bound
	for(int i = 0; i < n; i++)
	{
		for(int j = 0; j < n; j++)
		{
			matrix[i][j] = ((DataType)rand() / RAND_MAX) * 2.0f * bound - bound;
		}
	}
}

void SquareMatrix2D_NaiveMatmul(int n, DataType **C, DataType **A, DataType **B) 
{
	for(int i = 0; i < n; i++)
	{
		for(int j = 0; j < n; j++)
		{
			C[i][j] = 0.0f;
			for(int k = 0; k < n; k++)
			{
				C[i][j] += A[i][k] * B[k][j];
			}
		}
	}
}

void SquareMatrix2D_Add(int n, DataType **C, DataType **A, DataType **B)
{
	for(int i = 0; i < n; i++)
	{
		for(int j = 0; j < n; j++)
		{
			C[i][j] = A[i][j] + B[i][j];
		}
	}
}

void SquareMatrix2D_Sub(int n, DataType **C, DataType **A, DataType **B)
{
	for(int i = 0; i < n; i++)
	{
		for(int j = 0; j < n; j++)
		{
			C[i][j] = A[i][j] - B[i][j];
		}
	}
}

void SquareMatrix2D_Strassen(int n, DataType **C, DataType **A, DataType **B)
{
	if(n <= 64)//New base case
	{
		SquareMatrix2D_NaiveMatmul(n, C, A, B);
		return;
	}
	if(n == 1)//naive base case
	{
		C[0][0] = A[0][0] * B[0][0];
		return;
	}
	
	int k = n / 2;
	//Matrix partitions
	DataType **A11 = SquareMatrix2D_Allocate(k);
	DataType **A12 = SquareMatrix2D_Allocate(k);
	DataType **A21 = SquareMatrix2D_Allocate(k);
	DataType **A22 = SquareMatrix2D_Allocate(k);
	DataType **B11 = SquareMatrix2D_Allocate(k);
	DataType **B12 = SquareMatrix2D_Allocate(k);
	DataType **B21 = SquareMatrix2D_Allocate(k);
	DataType **B22 = SquareMatrix2D_Allocate(k);
	//Intermediate results ie 7 mults
	DataType **M1 = SquareMatrix2D_Allocate(k);
	DataType **M2 = SquareMatrix2D_Allocate(k);
	DataType **M3 = SquareMatrix2D_Allocate(k);
	DataType **M4 = SquareMatrix2D_Allocate(k);
	DataType **M5 = SquareMatrix2D_Allocate(k);
	DataType **M6 = SquareMatrix2D_Allocate(k);
	DataType **M7 = SquareMatrix2D_Allocate(k);
	//Temporary result holders
	DataType **C1 = SquareMatrix2D_Allocate(k);
	DataType **C2 = SquareMatrix2D_Allocate(k);
	
	//Compute constants
	for(int i = 0; i < k; i++)
	{
		for(int j = 0; j < k; j++)
		{
			A11[i][j] = A[i][j];
			A12[i][j] = A[i][j+k];
			A21[i][j] = A[i+k][j];
			A22[i][j] = A[i+k][j+k];

			B11[i][j] = B[i][j];
			B12[i][j] = B[i][j+k];
			B21[i][j] = B[i+k][j];
			B22[i][j] = B[i+k][j+k];
		}
	}
	SquareMatrix2D_Add(k, C1, A11, A22);
	SquareMatrix2D_Add(k, C2, B11, B22);
	SquareMatrix2D_Strassen(k, M1, C1,C2);
	
	SquareMatrix2D_Add(k, C1, A21, A22);
	SquareMatrix2D_Strassen(k, M2, C1, B11);
	
	SquareMatrix2D_Sub(k, C1, B12, B22);
	SquareMatrix2D_Strassen(k, M3, A11, C1);
	
	SquareMatrix2D_Sub(k, C2, B21, B11);
	SquareMatrix2D_Strassen(k, M4, A22, C2);
	
	SquareMatrix2D_Add(k, C1, A11, A12);
	SquareMatrix2D_Strassen(k, M5, C1, B22);
	
	SquareMatrix2D_Sub(k, C1, A21, A11);
	SquareMatrix2D_Add(k, C2, B11, B12);
	SquareMatrix2D_Strassen(k, M6, C1,C2);
	
	SquareMatrix2D_Sub(k, C1, A12, A22);
	SquareMatrix2D_Add(k, C2, B21, B22);
	SquareMatrix2D_Strassen(k, M7, C1,C2);
	
	for(int i = 0; i < k; i++)
	{
		for(int j = 0; j < k; j++)
		{
			C[i][j] = M1[i][j] + M4[i][j] - M5[i][j] + M7[i][j];
			C[i][j+k] =M3[i][j] + M5[i][j];
			C[i+k][j] =M2[i][j] + M4[i][j];
			C[i+k][j+k] =M1[i][j] - M2[i][j] + M3[i][j] + M6[i][j];
		}
	}
	SquareMatrix2D_Free(k,A11);SquareMatrix2D_Free(k,A12);SquareMatrix2D_Free(k,A21);SquareMatrix2D_Free(k,A22);
	SquareMatrix2D_Free(k,B11);SquareMatrix2D_Free(k,B12);SquareMatrix2D_Free(k,B21);SquareMatrix2D_Free(k,B22);
	SquareMatrix2D_Free(k,M1);SquareMatrix2D_Free(k,M2);SquareMatrix2D_Free(k,M3);SquareMatrix2D_Free(k,M4);
	SquareMatrix2D_Free(k,M5);SquareMatrix2D_Free(k,M6);SquareMatrix2D_Free(k,M7);SquareMatrix2D_Free(k,C1);SquareMatrix2D_Free(k,C2);
}

bool SquareMatrix2D_Equal(int n, DataType **A, DataType **B)
{
	DataType epsilon = 1e-1;
	for(int i = 0; i < n; i++)
	{
		for(int j = 0; j < n; j++)
		{
			//printf("%.16f %.16f\n",A[i][j], B[i][j]);
			if(fabsf(A[i][j] - B[i][j]) > epsilon)
			return false;
		}
	}
	return true;
}
bool CheckPowerOfTwo(int n)
{
	return n > 0 && (n & (n - 1)) == 0;
}
void CompareStrassenNaive(int squareMatrixDimension)
{
	assert(CheckPowerOfTwo(squareMatrixDimension) == true);
	DataType **A = SquareMatrix2D_Allocate(squareMatrixDimension);
	DataType **B = SquareMatrix2D_Allocate(squareMatrixDimension);
	DataType **C_naive = SquareMatrix2D_Allocate(squareMatrixDimension);	
	DataType **C_strassen = SquareMatrix2D_Allocate(squareMatrixDimension);
	SquareMatrix2D_Random(squareMatrixDimension, A, 10);
	SquareMatrix2D_Random(squareMatrixDimension, B, 10);
	clock_t start, end;
	double t_naive, t_strassen;

	start = clock();
	SquareMatrix2D_NaiveMatmul(squareMatrixDimension, C_naive, A, B);
	end = clock();
	t_naive = (double)(end - start) / CLOCKS_PER_SEC;

	start = clock();
	SquareMatrix2D_Strassen(squareMatrixDimension, C_strassen, A, B);
	end = clock();
	t_strassen = (double)(end - start) / CLOCKS_PER_SEC;

	bool equalCheck = SquareMatrix2D_Equal(squareMatrixDimension, C_naive, C_strassen);
	printf("N = %d\n", squareMatrixDimension);
	printf("Naive:    %.6f sec\n", t_naive);
	printf("Strassen: %.6f sec\n", t_strassen);
	printf("Equal: %s\n\n", equalCheck ? "YES" : "NO");
	SquareMatrix2D_Free(squareMatrixDimension, A);
	SquareMatrix2D_Free(squareMatrixDimension, B);
	SquareMatrix2D_Free(squareMatrixDimension, C_naive);
	SquareMatrix2D_Free(squareMatrixDimension, C_strassen);
}

int main()
{
	int n = 512;
	CompareStrassenNaive(n);
	return 0;
}
