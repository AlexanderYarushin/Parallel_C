#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "omp.h"

const unsigned N = 3000;

double *createArray(_Bool c){
    double *array = (double*) malloc(sizeof(double) * N);
    for(unsigned i = 0; i < N; ++i){
        if(c) array[i] = rand() / (RAND_MAX * 1.0);

        if(!c) array[i] = 0;
    }
    return array;
}

double **createMatrix(_Bool c){
    double **matrix = (double**) malloc(sizeof(double*) * N);
    for(unsigned i = 0; i < N; ++i){
        matrix[i] = (double*) malloc(sizeof(double) * N);
        for(unsigned j = 0; j < N; ++j){
            if(c) matrix[i][j] = rand() / (RAND_MAX * 1.0);
            if(!c) matrix[i][j] = 0;
        }
    }
    return matrix;
}

double *minus(double* vector1, double* vector2){
    double *result = createArray(0);

    for(unsigned i = 0; i < N; ++i){
        result[i] = vector1[i] - vector2[i];
    }

    return result;
}

double *mul(double** matrix, double* vector){
    double *result = createArray(0);
    double sum = 0;

    for(unsigned i = 0; i < N; ++i){
        sum = 0;
        for(unsigned j = 0; j < N; ++j){
            sum += matrix[i][j] * vector[j];

        }
        result[i] = sum;
    }

    return result;
}

double *gauss(double **a, double *y, int n)
{
    double *x, max;
    int k, index;
    const double eps = 0.00000001;  // точность
    x = malloc(n * sizeof(double));
    k = 0;
    while (k < n)
    {
        // Поиск строки с максимальным a[i][k]
        max = abs(a[k][k]);
        index = k;
        for (int i = k + 1; i < n; i++)
        {
            if (abs(a[i][k]) > max)
            {
                max = abs(a[i][k]);
                index = i;
            }
        }
        // Перестановка строк
        for (int j = 0; j < n; j++)
        {
            double temp = a[k][j];
            a[k][j] = a[index][j];
            a[index][j] = temp;
        }
        double temp = y[k];
        y[k] = y[index];
        y[index] = temp;
        // Нормализация уравнений
        for (int i = k; i < n; i++)
        {
            double temp = a[i][k];
            if (abs(temp) < eps) continue; // для нулевого коэффициента пропустить
            for (int j = 0; j < n; j++)
                a[i][j] = a[i][j] / temp;
            y[i] = y[i] / temp;
            if (i == k)  continue; // уравнение не вычитать само из себя
            for (int j = 0; j < n; j++)
                a[i][j] = a[i][j] - a[k][j];
            y[i] = y[i] - y[k];
        }
        k++;
    }
    // обратная подстановка
    for (k = n - 1; k >= 0; k--)
    {
        x[k] = y[k];
        for (int i = 0; i < k; i++)
            y[i] = y[i] - a[i][k] * x[k];
    }
    return x;
}

double **mulMatrix(double **m1, double **m2){
    double **newMatrix = createMatrix(0);

    for(unsigned i = 0; i < N; ++i){
        for(unsigned j = 0; j < N; ++j){
            for(unsigned k = 0; k < N; ++k){
                newMatrix[i][j] += m1[i][k] * m2[k][j];
                printf("%i\n",k);
            }
        }
    }


    return newMatrix;
}


void logMatrix(double **matrix){
    for(unsigned i = 0; i < N; ++i){
        for(unsigned j = 0; j < N; ++j){
            printf("%0.1lf ", matrix[i][j]);
        }
        if(i != N-1)printf("\n");
    }
}

void logVector(double *vector){
    for(unsigned i = 0; i < N; ++i){
        printf("%0.1lf ", vector[i]);
    }
}


int main()
{
    srand(time(NULL));

    double t1, t2;


    double **a = createMatrix(1);
    double *y = createArray(1);
    double *x;

    double *ax;
    double *axb;

    t1 = omp_get_wtime();
    x = gauss(a, y, N);
    t2 = omp_get_wtime();
    for (int i = 0; i < N; ++i)
       printf("x[%i] %f \n",i,x[i]);

    ax = mul(a,x);
    axb = minus(ax,y);

    double s = 0;

    for(int i = 0 ; i < N; ++i){
        s += pow(axb[i],2);
    }


    printf("%0.7lf \n",sqrt(s));




    return 0;
}
