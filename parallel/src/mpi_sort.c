#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>

#define MASTER 0


/**
 * swaps elements
 * @param a element to be swapped
 * @param b element to be swapped
 */
void swap(int* a, int* b){
    int temp = *a;
    *a = *b;
    *b = temp;
}

/**
 * finds the right position of the pivot
 * @param arr the unsorted array
 * @param left index of first element
 * @param right index of last element
 * @return its right position
 */
int partition(int arr[], int left, int right){
    int pivot = arr[right];
    int i = left;

    for(int j = left; j < right; j++){
        // if current element is smaller than the pivot
        if(arr[j] < pivot){
            i++;
            swap(&arr[i], &arr[j]);
        }
    }
    swap(&arr[i], &arr[right]);

    return i;
}


/**
 * picks a random pivot element between left and right and partitions arr[l...r] around the randomly picked element using partition
 * @param arr
 * @param left
 * @param right
 * @return the right posiiton of the element
 */
int randomPartition(int arr[], int left, int right){
    srand(time(NULL));
    int length = right - left + 1;
    int pivot = rand() % length;

    swap(&arr[left + pivot], &arr[right]);

    return partition(arr, left, right);
}






















void distributeByMedian(int rank, long int dimension, double **holdPoints, long int pointsPerProc, int worldSize) {
    double *pivot = (double *) calloc(dimension, sizeof(double));
    double *distance = (double *) calloc(pointsPerProc,sizeof(double)); // for each point hold the squares of the subtractions

    if(rank == 0) {
        //printf("The pivot is: ");
        for (int j = 0; j < dimension; j++) {
            pivot[j] = holdPoints[0][j];
            //printf("%.10f ", pivot[j]);
        }
        putchar('\n');
    }
    //MPI_Barrier(MPI_COMM_WORLD);

    MPI_Bcast(pivot, (int)dimension, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    printf(" we exited with flag %d\n", MPI_SUCCESS);

    //allocate the array of the subtractions
    double **subtraction = (double **) calloc(pointsPerProc,sizeof(double *));
    for (int i = 0; i < pointsPerProc; i++) {
        subtraction[i] = (double *) calloc(dimension, sizeof(double));
    }

    for (int i = 0; i < pointsPerProc; i++) {
        for (int j = 0; j < dimension; j++) {
            subtraction[i][j] = (holdPoints[i][j] - pivot[j]) * (holdPoints[i][j] - pivot[j]);
            distance[i] += subtraction[i][j];
        }
        printf("The distance is %.10f\n", distance[i]);
    }

    // allocate an array for all the coming from the processes
    double *receiver = NULL;

    if (rank == 0){
        receiver = (double *) malloc(sizeof(double) * worldSize * pointsPerProc);
    }


    MPI_Gather(distance, (int)pointsPerProc, MPI_DOUBLE, receiver, (int)pointsPerProc, MPI_DOUBLE, MASTER,
               MPI_COMM_WORLD);

    if (rank == 0){
        for(int j = 0; j < worldSize * pointsPerProc; j++) {
            printf("The receiver is %.10f  ", receiver[j]);
        }
        printf("\n");
    }


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
    }

    distributeByMedian(rank, dimension, holdThePoints, pointsPerProc, size);

    MPI_Finalize();

    return 0;
}


