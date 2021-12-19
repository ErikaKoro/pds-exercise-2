#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "quick_select.h"

/**
 * Finds the exchanges between the processes according to the following algorithm
 *
 *
 * @param counterReceiver The array with the number of points, each process wants to give, that has bees sent to the master
 * @param worldSize
 * @return the array with the infos
 */
int **findExchanges(int *counterReceiver, int worldSize){

    // Allocate memory for the array that will contain with which process should the current process communicate and how many points will they exchange
    int **info = (int **)malloc(worldSize * sizeof(int *));
    for(int i = 0; i < worldSize; i++){
        info[i] = (int *)malloc(2 * worldSize * sizeof (int));
    }

    // This array holds the index for its line in the info array in which we should write the next element
    int *helperIndex = (int *) calloc(worldSize, sizeof (int));

    int nextIndex = worldSize / 2;  // Indicates the index of the second half of the array
    for(int i = 0; i < worldSize/2; i++){
        while(counterReceiver[i] > 0){
            if(counterReceiver[i] == counterReceiver[nextIndex]){
                info[i][helperIndex[i]] = nextIndex;
                helperIndex[i]++;
                info[i][helperIndex[i]] = counterReceiver[i];
                helperIndex[i]++;
                info[nextIndex][helperIndex[nextIndex]] = i;
                helperIndex[nextIndex]++;
                info[nextIndex][helperIndex[nextIndex]] = counterReceiver[i];
                helperIndex[nextIndex]++;
                counterReceiver[i] = 0;
                counterReceiver[nextIndex] = 0;
                nextIndex++;
            }
            else if(counterReceiver[i] < counterReceiver[nextIndex]){
                info[i][helperIndex[i]] = nextIndex;
                helperIndex[i]++;
                info[i][helperIndex[i]] = counterReceiver[i];
                helperIndex[i]++;
                info[nextIndex][helperIndex[nextIndex]] = i;
                helperIndex[nextIndex]++;
                info[nextIndex][helperIndex[nextIndex]] = counterReceiver[i];
                helperIndex[nextIndex]++;
                counterReceiver[nextIndex] -= counterReceiver[i];
                counterReceiver[i] = 0;
            }
            else if(counterReceiver[i] > counterReceiver[nextIndex] && counterReceiver[nextIndex] != 0){
                info[i][helperIndex[i]] = nextIndex;
                helperIndex[i]++;
                info[i][helperIndex[i]] = counterReceiver[nextIndex];
                helperIndex[i]++;
                info[nextIndex][helperIndex[nextIndex]] = i;
                helperIndex[nextIndex]++;
                info[nextIndex][helperIndex[nextIndex]] = counterReceiver[nextIndex];
                helperIndex[nextIndex]++;
                counterReceiver[i] -= counterReceiver[nextIndex];
                counterReceiver[nextIndex] = 0;
                nextIndex++;
            }
            else{
                nextIndex++;
            }
        }
    }
    printf("The array with exchanges is: \n");
    for(int i = 0; i < worldSize; i++){
        printf("The i is %d ", i);
        for (int j = 0; j < helperIndex[i]; ++j) {
            printf("%d ", info[i][j]);
        }
        printf("\n");
    }

    return info;
}

/**
 * Finds the distance of each element, per process, from the pivot
 *
 * @param dist array to hold the results
 * @param points points per process with their coordinates
 * @param dimension the dimension of our space (e.g 3D)
 * @param pivot array with the coordinates of the pivot
 * @param pointsPerProc
 */
void findDistance(int rank, double *dist, double **points, int dimension, const double *pivot, int pointsPerProc){

    for (int i = 0; i < pointsPerProc; ++i) {
        for (int j = 0; j < dimension; ++j) {
            dist[i] += (points[i][j] - pivot[j]) * (points[i][j] - pivot[j]);
        }
        printf("Rank: %d    The distance is %.10f\n", rank, dist[i]);
    }
}

/**
 * Sorts the array with the points per process, holding the ones that each process wants to keep in the first indexes of the array and the other part of the array has the points that it wants to give
 *
 * @param worldSize
 * @param rank
 * @param holdThePoints
 * @param pointsPerProc
 * @param dimension
 * @param distance
 * @param median
 * @return the counter of the points that will remain in the process
 */
int partitionByMedian(int worldSize, int rank, double **holdThePoints, int pointsPerProc, int dimension, double *distance, double median){
    int j = 0;
    int i = 0;
    while(j < pointsPerProc){
        if(distance[j] < median && rank < worldSize/2){
            double temp;
            double *tmp;
            temp = distance[j];
            distance[j] = distance[i];
            distance[i] = temp;

            // In the first indexes of the array we want to keep the points that will remain in the process and in the last the elements that the process wants to give
            tmp = holdThePoints[j];  // so swap the elements
            holdThePoints[j] = holdThePoints[i];
            holdThePoints[i] = tmp;
            i++;
            j++;

        } else if(distance[j] > median && rank >= worldSize/2){
            double temp;
            double *tmp;
            temp = distance[j];
            distance[j] = distance[i];
            distance[i] = temp;

            // In the last indexes of the array we want to keep the points that will be exchanged with other points from other processes
            tmp = holdThePoints[j];
            holdThePoints[j] = holdThePoints[i];
            holdThePoints[i] = tmp;
            j++;
            i++;
        } else {
            j++;
        }
    }

    return i;  // At the end of this function the i index points to the first element of the second part of the array with the points of the process that will be exchanged.
}

/**
 * Chooses the pivot, broadcasts it to the processes, calculates the distances, gathers them to the master, which finds their median
 *
 * @param master master rank
 * @param rank process id
 * @param dimension
 * @param holdPoints the array with the points and their coordinates
 * @param pointsPerProc
 * @param worldSize the number of processes
 * @param communicator "holds" the processes with which we work each time we call recursively the function
 */
void distributeByMedian(double *pivot,int master, int rank, int dimension, double **holdPoints, int pointsPerProc, int worldSize, MPI_Comm communicator) {
    double *distance = (double *) calloc(pointsPerProc,sizeof(double)); // for each point hold the squares of the subtractions

    // Find the distance of each point of the process from the pivot point chosen by the master
    findDistance(rank, distance, holdPoints, dimension, pivot, pointsPerProc);


    //Allocate the master's buffer to hold the distances
    double *receiver = NULL;
    if (rank == master){
        receiver = (double *) malloc(sizeof(double) * worldSize * pointsPerProc);
    }


    //The master gathers from all the processes their points' distances from the pivot
    MPI_Gather(distance, (int)pointsPerProc, MPI_DOUBLE, receiver, (int)pointsPerProc, MPI_DOUBLE, master,
               communicator);

    double median;
    if (rank == master){
        // The master, having the array with all the distances finds their median using the "quickselect" algorithm
        median = findMedian(receiver, worldSize * pointsPerProc);
        printf("\n\n\nDistances master\n");
        for(int j = 0; j < worldSize * pointsPerProc; j++) {
            printf("%.10f, ", receiver[j]);
        }
        printf("\n");

        printf("Median: %f\n", median);
    }


    //Broadcast the median to the processes
    MPI_Bcast(&median, 1, MPI_DOUBLE, master, communicator);


    // Call the function which calculates how many points in each process will remain
    // Sort the array holdPoints by holding in the first indexes the points that will remain in the process and in the last indexes the ones that will be exchanged
    int counter = partitionByMedian(worldSize, rank, holdPoints, pointsPerProc, dimension, distance, median);

    printf("\n\nNow the rank %d has holdThePoints is: ", rank);
    for(int i = 0; i < pointsPerProc; i++){
        for(int j = 0; j < dimension; j++){
            printf("%.10f ", holdPoints[i][j]);
        }
        printf("\n\n");
    }
    int pointsToGive = pointsPerProc - counter;  // The number of points each process wants to give
    printf("Rank: %d    Counter: %d", rank, counter);


    printf("\n\nThe process with rank %d wants to give %d points\n", rank, pointsToGive);

    int *counterReceiver = (int *)malloc(worldSize * sizeof (int));  // Allocate a buffer in which the master will store the pointsToGive from each process

    // Give the number of points each process wants to exchange to the master
    MPI_Gather(&pointsToGive, 1, MPI_INT, counterReceiver, 1, MPI_INT, master, communicator);

    if(rank == master) {
        printf("\n\nThe counter receiver is: ");
        for(int i = 0; i < worldSize; i++) {
            printf(" %d ", counterReceiver[i]);
        }
        printf("\n\n\n");
    }

    if(rank == master) {
        int **exchanges = findExchanges(counterReceiver, worldSize);
//        printf("The array with exchanges is: ");
//        for(int i = 0; i < worldSize; i++){
//            for (int j = 0; j < worldSize; ++j) {
//                printf("%d ", exchanges[i][j]);
//            }
//            printf("\n");
//        }
    }
    else{
        printf("leave me alone\n\n\n");
    }

}


int main(int argc, char **argv) {
    int size;
    int rank;
    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* Read the binary file with our data back in */
    // The first two 64-bit numbers from the file give us the number of our points and the space's dimension
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

        // hold our points with their coordinates read from the file
        holdThePoints = (double **) malloc(pointsPerProc * sizeof(double *));
        for (int i = 0; i < pointsPerProc; i++) {
            holdThePoints[i] = (double *) malloc(dimension * sizeof(double));
        }

        // Seek where is the pointer in our file
        fseek(fh, rank * dimension * pointsPerProc * sizeof(double), SEEK_CUR);

        for (int i = 0; i < pointsPerProc; i++) {
            double *x2 = (double *) malloc(sizeof(double) * dimension); // represents only one point

            fread(x2, sizeof(double), dimension, fh);
            holdThePoints[i] = x2;
        }
        fclose(fh);

        printf("\n\n\n");
        printf("Rank: %d\n", rank);
        for (int i = 0; i < pointsPerProc; ++i) {
            for (int j = 0; j < dimension; ++j) {
                printf("%.10f ", holdThePoints[i][j]);
            }
            printf("\n\n\n");
        }

        double *pivot = (double *) calloc(dimension, sizeof(double));

        if(rank == 0) {
            //printf("The pivot is: ");
            for (int j = 0; j < dimension; j++) {
                pivot[j] = holdThePoints[0][j];
                //printf("%.10f ", pivot[j]);
            }
            //putchar('\n');
        }
        //Broadcast the pivot to the processes
        MPI_Bcast(pivot, (int)dimension, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        distributeByMedian(pivot, 0, rank, (int)dimension, holdThePoints, (int)pointsPerProc, size, MPI_COMM_WORLD);


    }



    MPI_Finalize();

    return 0;
}


