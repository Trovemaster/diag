# diag
This is a simple wrapper (.f90) for a Lapack diagonilzer DSYEV. It requires an input (ASCII) file containing real matrix elements of a symemtric NxN matrix. Here is an exmaple:

     56
      1       1      2904.93663932
      1       2        13.91569619
      1       3        -3.94951978
      1       4        -0.20100770
      1       5         0.00782046
            
where the first line is the dimension N of the matrix A, while other lines list the matrix indeces and elements. Only one half of the symmetric matrix is required. If only the first line (the dimension) is given, the program will use random numbers to generate a real, symmetric, diagonal dominated matrix. The program will then print the eigenvalues E, and optioaanlly the eigenvectors V and the product V^T M V = D, which should be a diagonal matrix. The diagonalization time is printed out as well. 

