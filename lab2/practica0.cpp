#include <iostream>  
#include <cstdlib>   
#include <cmath>     
#include <ctime>     // para medir los tiempos
#include <iomanip>   // para la salida
#include <omp.h>



using namespace std;

double get_time() 
{
    struct timespec t;
    clock_gettime(CLOCK_REALTIME, &t);
    return (double)t.tv_sec + ((double)t.tv_nsec) / 1.0e9;
}

int main (int argc, char *argv[])
{

    int D, N;
    bool m_debug=false;
    if (argc < 2) 
    {
        cerr << "Error: Faltan argumentos. 0 para depuracion N D para experimental" << endl;
        return 1;
    }
    if(atoi(argv[1])==0)
    {
        m_debug=true;
        N=6;
        D=4;
        cout<<"MODO DEPURACION"<<endl;
    }
    else
    {
        if (argc < 3) 
        {
                cerr<< "Error: En modo experimental necesitas N y D." << endl;
                return 1;
        }
        cout << "MODO EXPERIMENTAL"<<endl;
        N = atoi(argv[1]); 
        D = atoi(argv[2]);
        m_debug=false;
    }

    double *X   = new double[N * D]; // Matriz Entrada
    double *C   = new double[N * D]; // Matriz Salida


    double *W_K = new double[D * D];
    double *W_Q = new double[D * D];
    double *W_V = new double[D * D];

    double *b_K = new double[D];
    double *b_Q = new double[D];
    double *b_V = new double[D];

    double *K = new double[N * D];
    double *Q = new double[N * D];
    double *V = new double[N * D];

    double *A = new double[N * N]; //matriz atencion



    if(m_debug)//modo depuracion
    {
        //inicializar X
        double X_datos[] = 
            {
                0,  6, 12, 18,
                1,  7, 13, 19,
                2,  8, 14, 20,
                3,  9, 15, 21,
                4, 10, 16, 22,
                5, 11, 17, 23
                
            };
        for(int i=0; i<N*D; i++) X[i] = X_datos[i]; //para rellenarlo

        //inicializar W_K

        double WK_fila[] = {-0.2, -0.1, 0.0, 0.1};

        for(int i=0; i<D; i++) 
            {     
                for(int j=0; j<D; j++) {  
                    W_K[i*D + j] = WK_fila[j];
                }
            }

        //inicializar a W_Q

        double WQ_datos[] = 
            {
                -0.2, -0.2, -0.2, -0.2,
                -0.1, -0.1, -0.1, -0.1,
                0.0,  0.0,  0.0,  0.0,
                0.1,  0.1,  0.1,  0.1
            };
        for(int i=0; i<D*D; i++) 
            W_Q[i] = WQ_datos[i];

        //inicializar W_V

        for(int i=0; i<D*D; i++) 
            W_V[i] = 0.0; 
        for(int i=0; i<D; i++) 
            W_V[i*D + i] = 1.0;

        //inicializar las b 
        for(int i=0; i<D; i++) 
            {
                b_K[i] = -1.0;
                b_Q[i] =  0.1;
                b_V[i] =  0.0;
            }
    }
    else //modo experimental
    {
        srand(time(NULL)); //random
        //inicializacion X
        for(int i=0; i<N*D; i++) 
                X[i] = 10.0 * (double)rand() / RAND_MAX;
        //inicializacion de las W
        for(int i=0; i<D*D; i++) 
            {
                W_K[i] = 0.001 * (((double)rand() / RAND_MAX) - 0.5);
                W_Q[i] = 0.001 * (((double)rand() / RAND_MAX) - 0.5);
                W_V[i] = 0.001 * (((double)rand() / RAND_MAX) - 0.5);
            }
        //inicializacion de las B
        for(int i=0; i<D; i++) 
            {
                b_K[i] = 0.001 * (((double)rand() / RAND_MAX) - 0.5);
                b_Q[i] = 0.001 * (((double)rand() / RAND_MAX) - 0.5);
                b_V[i] = 0.001 * (((double)rand() / RAND_MAX) - 0.5);
            }
    }



    double t_inicio = get_time(); //inicio cronometro
    double inv_sqrt_D = 1.0 / sqrt((double)D);// precalculamos la raiz de D fuera del parallel para que no se ejecute tantas veces como hilos
#pragma omp parallel
{
    //calculo de q,k,v
    #pragma omp for schedule(static)
    for (int n = 0; n < N; n++) 
    {
        for (int d = 0; d < D; d++) {
            
            
            double val_q = b_Q[d];
            double val_k = b_K[d];
            double val_v = b_V[d];


            for (int l = 0; l < D; l++) //multiplicacion de matrices
            {
                double x_val = X[n * D + l];
                val_q += x_val * W_Q[l * D + d];
                val_k += x_val * W_K[l * D + d];
                val_v += x_val * W_V[l * D + d];
            }

            Q[n * D + d] = val_q;
            K[n * D + d] = val_k;
            V[n * D + d] = val_v;
        }
    }


    //calculo de A
    #pragma omp for schedule(static)
    for(int n=0; n<N; n++){
        for(int i=0; i<N; i++){
            double dot_prod = 0.0;
            for (int d = 0; d < D; d++) //mmultiplicacion de matrices
            {
                dot_prod += Q[n * D + d] * K[i * D + d];
            }
            A[n * N + i] = dot_prod * inv_sqrt_D;
        }
    }


    //calculo de alpha
    #pragma omp for schedule(static)
    for (int n = 0; n < N; n++) {

            double suma_fila = 0.0;
            
            for (int i = 0; i < N; i++) {
                double valor_original = A[n * N + i];
                double ex = exp(valor_original);
                
                A[n * N + i] = ex; 
                suma_fila += ex; 
            }


            double inv_suma = 1.0 / suma_fila;
            
            for (int i = 0; i < N; i++) {
                A[n * N + i] *= inv_suma; 
            }
        }


    //calculo de C
    #pragma omp for schedule(static)
    for(int n=0;n < N;n++){
        for(int d=0;d < D;d++)
        {
            double val_salida = 0.0;
            for (int i = 0; i < N; i++) //multiplicacion de matrices
            {
                val_salida += A[n * N + i] * V[i * D + d]; 
            }
            C[n * D + d] = val_salida;
        }
    }
}
    double t_fin = get_time(); // fin de cronometro
    double tiempo_total = t_fin - t_inicio;


    //resultados
    cout << "Tiempo de ejecucion: " << tiempo_total << " s" << endl;
    if(m_debug)
    {
        cout << "Matriz C (Resultado):" << endl;

        cout << fixed << setprecision(1); //para que los numeros salgan formateados
        for (int n = 0; n < N; n++) {
                    for (int d = 0; d < D; d++) {
                        cout << C[n * D + d] << " ";
                    }
                    cout << endl;//al terminar cada fila un salto de linea
                }
    }

    /************************************liberacion de memoria*************************************************/
    delete[] X; delete[] C;
    delete[] W_K; delete[] W_Q; delete[] W_V;
    delete[] b_K; delete[] b_Q; delete[] b_V;
    delete[] K; delete[] Q; delete[] V; delete[] A;
    /*********************************************************************************************************/

    return 0;   
}




