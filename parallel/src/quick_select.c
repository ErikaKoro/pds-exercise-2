#include <stdio.h>
#include <time.h>
#include <stdlib.h>

// Utility function to swapping of element
void swap(double* a, double* b){
    double temp = *a;
    *a = *b;
    *b = temp;
}

// find the right position of pivot
//retutn its index
int partition(double arr[], int left, int right){
    double pivot = arr[right];

    int i = left;
    int j = left;

    while (j < right){
        if(arr[j] < pivot){
            swap(&arr[i], &arr[j]);
            i++;
        }
        j++;
    }

    swap(&arr[i], &arr[right]);

    return i;
}



//Picks a random pivot element between
//l and r and partitions arr[l..r]
//around the randomly picked element
//using partition()
int randomPartition(double arr[], int left, int right){
    srand(time(NULL));
    int length = right - left + 1;
    int pivot = rand() % length;

    swap(&arr[left + pivot], &arr[right]);

    return partition(arr, left, right);
}

//utility function to find median
void MedianUtil(double arr[], int left, int right, int middle, double *a, double *b){

    if (left <= right) {

        // Find the partition index
        int partitionIndex = randomPartition(arr, left, right);


        if (partitionIndex == middle) {
            *b = arr[partitionIndex];
            if (*a != -1)
                return;
        }


        if (partitionIndex == middle - 1) {
            *a = arr[partitionIndex];
            if (*b != -1)
                return;
        }

        // If partitionIndex >= k then
        // find the index in first half
        // of the arr[]
        if (partitionIndex >= middle){
            return MedianUtil(arr, left, partitionIndex - 1, middle, a, b);
        }
            // If partitionIndex <= k then
            // find the index in second half
            // of the arr[]
        else{
            return MedianUtil(arr, partitionIndex + 1, right, middle, a, b);
        }
    }
    return;
}

// Function to find Median
double findMedian(double *arr, int length){
    double result;

    double a = -1, b = -1;

    // If n is odd
    if (length % 2 == 1) {
        MedianUtil(arr, 0, length - 1, length / 2, &a, &b);
        result = b;
    }
    // If n is even
    else {
        MedianUtil(arr, 0, length - 1, length / 2, &a, &b);
        result = (a + b) / 2;
    }

    // // our length is even(from our data)
    // MedianUtil(arr, 0, length - 1, length/2, &a, &b);
    // result = (a + b) / 2;

    // Print the Median of arr[]
    printf("The median is: %f\n", result);

    return result;

}

