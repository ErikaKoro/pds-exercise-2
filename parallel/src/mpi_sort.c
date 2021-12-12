#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "quick_select.h"

#define MASTER 0

void findDistance(double *dist, double **points, int dimension, const double *pivot, int pointsPerProc){

    for (int i = 0; i < pointsPerProc; ++i) {
        for (int j = 0; j < dimension; ++j) {
            dist[i] += (points[i][j] - pivot[j]) * (points[i][j] - pivot[j]);
        }
        printf("The distance is %.10f\n", dist[i]);
    }
}

void distributeByMedian(int master, int rank, int dimension, double **holdPoints, int pointsPerProc, int worldSize) {
    double *pivot = (double *) calloc(dimension, sizeof(double));
    double *distance = (double *) calloc(pointsPerProc,sizeof(double)); // for each point hold the squares of the subtractions

    if(rank == master) {
        //printf("The pivot is: ");
        for (int j = 0; j < dimension; j++) {
            pivot[j] = holdPoints[0][j];
            //printf("%.10f ", pivot[j]);
        }
        //putchar('\n');
    }
    //MPI_Barrier(MPI_COMM_WORLD);

    MPI_Bcast(pivot, (int)dimension, MPI_DOUBLE, master, MPI_COMM_WORLD);
    printf(" we exited with flag %d\n", MPI_SUCCESS);

    findDistance(distance, holdPoints, dimension, pivot, pointsPerProc);

    // allocate an array for all the distances coming from the processes
    double *receiver = NULL;

    if (rank == master){
        receiver = (double *) malloc(sizeof(double) * worldSize * pointsPerProc);
    }


    MPI_Gather(distance, (int)pointsPerProc, MPI_DOUBLE, receiver, (int)pointsPerProc, MPI_DOUBLE, MASTER,
               MPI_COMM_WORLD);

    double median;
    if (rank == master){
        median = findMedian(receiver, worldSize * pointsPerProc);


        for(int j = 0; j < worldSize * pointsPerProc; j++) {
            printf("%.10f, ", receiver[j]);
        }
        printf("\n");

        printf("Median: %f", median);
    }

    MPI_Bcast(&median, 1, MPI_DOUBLE, master, MPI_COMM_WORLD);


}


int main(int argc, char **argv) {
    int size;
    int rank;
    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* Read the file back in */
    int64_t numberOfPoints;
    int64_t dimension;
    long int pointsPerProc;
    double **holdThePoints;

    FILE *fh = fopen(argv[1], "rb");
    if (fh != NULL) {
        fread(&dimension, sizeof(int64_t), 1, fh);
        fread(&numberOfPoints, sizeof(int64_t), 1, fh);

        pointsPerProc = numberOfPoints / size;
//        printf("Rank: %d, Points per proc: %ld\n", rank, pointsPerProc);
        holdThePoints = (double **) malloc(pointsPerProc * sizeof(double *));
        for (int i = 0; i < pointsPerProc; i++) {
            holdThePoints[i] = (double *) malloc(dimension * sizeof(double));
        }

        fseek(fh, rank * dimension * pointsPerProc * sizeof(double), SEEK_CUR);

        for (int i = 0; i < pointsPerProc; i++) {
            double *x2 = (double *) malloc(sizeof(double) * dimension); // represents only one point

            fread(x2, sizeof(double), dimension, fh);
            holdThePoints[i] = x2;
        }
        fclose(fh);

        printf("\n");
        printf("Rank: %d\n", rank);
        for (int i = 0; i < pointsPerProc; ++i) {
            for (int j = 0; j < dimension; ++j) {
                printf("%.10f ", holdThePoints[i][j]);
            }
            printf("\n");
        }

        distributeByMedian(0, rank, (int)dimension, holdThePoints, (int)pointsPerProc, size);
    }



    MPI_Finalize();

    return 0;
}


