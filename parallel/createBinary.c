#include <stdint.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main(int argc, char **argv){
    if(argc != 4){
        printf("gimme args you stupid\n");
        return -1;
    }
    srand(time(NULL));
    /* Create the file */
    int64_t dimension= atoi(argv[1]);
    int64_t numberOfPoints = atoi(argv[2]);

    FILE *fh = fopen (argv[3], "wb");
    if (fh != NULL) {
        fwrite (&dimension, sizeof (int64_t), 1, fh);
        fwrite (&numberOfPoints, sizeof (int64_t), 1, fh);
        double range = 200.00;
        double div = RAND_MAX / range;
        for (int j = 0; j < numberOfPoints; j++){
            for(int i = 0; i < dimension; i++){
                
                double coordinates = -100 + (rand() / div);
                //printf("%.10f ", coordinates);               
                fwrite(&coordinates, sizeof(double), 1, fh);
            }
            printf("\n");
        }
        fclose (fh);
    }
     printf("\n");
     printf("\n");
     printf("\n");

    /* Read the file back in */
    int64_t numberOfPoints2;
    int64_t dimension2;

    fh = fopen (argv[3], "rb");
    if (fh != NULL) {
         fread(&dimension2, sizeof(int64_t), 1, fh);
         fread(&numberOfPoints2, sizeof(int64_t), 1, fh);

         double *x2 = (double *)malloc(sizeof(double) * dimension2);

         for(int i = 0; i < numberOfPoints; i++){
             fread (x2, sizeof (double), dimension2, fh);

             //for(int j = 0; j < dimension2; j++){
               //  printf("%.10f ", x2[j]);
             //}
             //printf("\n");
        }


        fclose (fh);
    }

    /* Check that it worked */
   // printf ("Value is: %ld\n", x2);
    //printf ("Value is: %ld\n", y2);

    return 0;
}
