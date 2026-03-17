#include <omp.h>
#include <stdio.h>
#include <unistd.h>

void tarea_uno() { sleep(2);// retardo de 2 segundos }
void tarea_dos() { sleep(4);// retardo de 4 segundos }

int main() {
    int tid;
    double timeIni, timeFin;


    #pragma omp parallel private(tid)
    {
        tid = omp_get_thread_num();


        double localStart = omp_get_wtime();

        #pragma omp sections
        {
            #pragma omp section
            {
                printf("(hilo: %d) Ejecutando tarea 1\n", tid);
                tarea_uno();
            }

            #pragma omp section
            {
                printf("(hilo: %d) Ejecutando tarea 2\n", tid);
                tarea_dos();
            }
        }
        // Al salir del bloque sections hay una barrera implícita (sincronización)

        double localEnd = omp_get_wtime();
        printf("(hilo: %d) Tiempo requerido: %f segundos\n", tid, localEnd - localStart);
    }

    return 0;
}
