# Matrix_Library

This library simplifies working with matrices and vectors in linear algebra. All the basic functions are implemented here, such as:

* Addition, subtraction, multiplication of matrices and vectors using basic sign operations +, -, *
* Calculate the Hadamard product using the method Adamar(const Matrix& a)
* Trace() 
* Matrix transposition - Transpose()
* Scalar product of vectors - Scalar(const Matrix& a)
* Calculate max element - MaximumRate()
* Find the norm of the vector - EuclideanNorm()
* Calculate the Frobenius norm of the matrix (considered as the root of the sum of squares of all matrix elements) - FrobeniusNorm()
* Calculate the angle between vectors (implemented on the cosine theorem) - Angle(Matrix& a)
* The method Tril(int n, int &change) can bring the matrix to a triangular form
* Calculate the determinant of the matrix - Det()
* Calculate the rank of the matrix - Rank()
* Find the inverse matrix - Inverse()

The method of the main components is fully implemented separately
