#include <vector>

using Matrix = std::vector<std::vector<double>>;

Matrix multiply(Matrix& a, Matrix& b);
Matrix transpose(Matrix& a);
double det(Matrix& a);
Matrix inv(Matrix& a);
void show(Matrix& a);  //  Shows the matrix in a beautiful way

Matrix getData(char* filename);
Matrix getX(Matrix&);
Matrix gety(Matrix& data);

Matrix solve(char* filename);