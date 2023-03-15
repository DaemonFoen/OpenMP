#include <iostream>
#include <omp.h>
#include <cmath>

static int N = 16000;


double *makeMatrixA() {
    double *matrix = new double[N * N];

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (i == j) {
                matrix[i * N + j] = 2.0;
            } else {
                matrix[i * N + j] = 1.0;
            }

        }
    }
    return matrix;
}

double *makeVectorB() {
    double *vector = new double[N];
    for (int i = 0; i < N; ++i) {
        vector[i] = (N + 1.0);
    }
    return vector;
}

void printVector(double *vector) {
    std::cout << "[ ";

    for (int i = 0; i < N; ++i) {
        std::cout << vector[i] << " ";
    }
    std::cout << "]";
}

double *makeVectorX() {
    double *vector = new double[N];
    for (int i = 0; i < N; ++i) {
        vector[i] = 0;
    }

    return vector;
}


int main() {
    double *matrixA = makeMatrixA();
    double *vectorB = makeVectorB();
    double *vectorX = makeVectorX();
    double *vectorXBuff = new double[N]();

    double tau = 0.0001;
    double eps = 1e-6;
    double normB = 0;
    double normV = 0;                       // normV = || Ax^n -b ||
    double newNormRatio = 0;

    bool flag = true;


    double startTime = omp_get_wtime();

    //count normB
#pragma omp parallel for schedule(static) reduction(+:normB)
    for (int i = 0; i < N; ++i) {
        normB += vectorB[i] * vectorB[i];
    }

    normB = sqrt(normB);


    while (flag) {
        //count normV
#pragma omp parallel for schedule(static) reduction(+:normV)
        for (int i = 0; i < N; ++i) {
            double valueX = 0;                          //valueX = Ax-b

            const double *m = matrixA + i * N;

            for (int j = 0; j < N; ++j) {
                valueX += m[j] * vectorX[j];            //Ax по i-ому
            }

            valueX -= vectorB[i];                       //Ax-b по i-ому

            vectorXBuff[i] = vectorX[i] - valueX * tau; //x^{n+1}[i]
            normV += valueX * valueX;
        }

        std::swap(vectorX, vectorXBuff);

        normV = sqrt(normV);
        newNormRatio = normV / normB;
        flag = newNormRatio > eps;

        normV = 0;
    }


    double end_time = omp_get_wtime();


    std::cout << "Time passed: " << end_time - startTime << " seconds. Threads used: " << omp_get_max_threads()
              << std::endl;

    delete[] matrixA;
    delete[] vectorX;
    delete[] vectorB;
    delete[] vectorXBuff;

return 0;
}