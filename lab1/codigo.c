#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
int main() {
    unsigned int N;
    int i, *x, *y;
    printf("Introduce el tamano del vector: ");
    scanf(" %u", &N);
    double timeIni, timeFin;
    x = (int*) malloc(N * sizeof(int));
    y = (int*) malloc(N * sizeof(int));
    timeIni = omp_get_wtime();
    // Secuencial
    for (i = 1; i < (int)N; i++) {
        x[i] = y[i-1] * 2;
        y[i] = y[i] + i;
    }
    timeFin = omp_get_wtime();
    printf("Tiempo secuencial = %f milisegundos\n", (timeFin - timeIni) * 1000);
    timeIni = omp_get_wtime();
    // Paralelo

    // Escribe aqui el mismo algoritmo paralelizado
    #pragma omp parallel if (N>20000)
    {
        #pragma omp for
        for (i = 1; i < (int)N; i++)
        {
            y[i] = y[i] + i;
        }
        #pragma omp for
        for (i=1; i < (int)N; i++)
        {
            x[i] = y[i-1] * 2;
        }
    }
    timeFin = omp_get_wtime();
    printf("Tiempo paralelo = %f milisegundos\n", (timeFin - timeIni) * 1000);
    free(x);
    free(y);
    return 0;
}
