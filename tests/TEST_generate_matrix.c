#include <stdio.h> 
#include "LinearAlgebraBasics.h"

int main() {

    int rows = 3;
    int columns = 4;
    
    double *matrix = generate_matrix_double(rows, columns);

    for (int i = 0; i < rows; i++) {
	for (int j = 0; j < columns; j++) {
	    printf("%lf\t", matrix[i * columns + j]);
	}
	printf("\n");
    }

    free(matrix);
    
    return 0;

}
