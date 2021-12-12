#include <stdio.h>
#include <time.h>
#include <stdlib.h>

// Utility function to swapping of element
void swap(int* a, int* b){
    int temp = *a;
    *a = *b;
    *b = temp;
}

// find the right position of pivot
//retutn its index
int partition(int arr[], int left, int right){
    int pivot = arr[right];

    int i = left;
    int j = left;

    while (j < right){
        if(arr[j] < pivot){
            swap(&arr[i], &arr[j]);
            i++;
        }
        j++;
    }



    // for(int j = left; j < right; j++){
    //     // if current element is smaller than the pivot
    //     if(arr[j] < pivot){
    //         swap(&arr[i], &arr[j]);
    //         i++;
    //     }
    // }

    swap(&arr[i], &arr[right]);

    return i;
}



//Picks a random pivot element between
//l and r and partitions arr[l..r]
//around the randomly picked element
//using partition()
int randomPartition(int arr[], int left, int right){
    // srand(time(NULL));
    int length = right - left + 1;
    int pivot = rand() % length;

    swap(&arr[left + pivot], &arr[right]);

    return partition(arr, left, right);
}

//utility function to find median
void MedianUtil(int arr[], int left, int right, int middle, int *a, int *b){

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
void findMedian(int arr[], int length){
    double result;

    int a = -1, b = -1;

    // If n is odd
    if (length % 2 == 1) {
        MedianUtil(arr, 0, length - 1, length / 2, &a, &b);
        result = b;
    }

        // If n is even
    else {
        MedianUtil(arr, 0, length - 1, length / 2, &a, &b);
        result = (double)(a + b) / 2;
    }

    // // our length is even(from our data)
    // MedianUtil(arr, 0, length - 1, length/2, &a, &b);
    // result = (a + b) / 2;

    // Print the Median of arr[]
    printf("The median is: %f\n", result);

}

// Driver program to test above methods
int main() {

    int arr[] = {3, 7, 9, 12, 19, 4, 5, 6};


    int n = 8;
    //int par = partition(arr, 0, 6);
    //int ran = randomPartition(arr, 0, 6);
    // printf("The index of partition is: %d\n", ran);
    findMedian(arr, n);
    return 0;
}

