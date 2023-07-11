#include "hw1.h"

#include <vector>
#include <iostream>
#include <random>
#include <iomanip>

using Matrix = std::vector<std::vector<double>>;

namespace algebra
{
Matrix zeros(size_t n, size_t m)
{
    if(n < 0 || m < 0)
    {
        throw std::logic_error("n and m cannot be negative");
    }
    Matrix result(n, std::vector<double>(m, 0));
    return result;
}

Matrix ones(size_t n, size_t m)
{
    if(n < 0 || m < 0)
    {
        throw std::logic_error("n and m cannot be negatve");
    }
    Matrix result(n, std::vector<double>(m, 1));
    return result;
}

Matrix random(size_t n, size_t m, double min, double max)
{
    if(min > max)
    {
        throw std::logic_error("min cannot be larger than max");
    }
    if(n < 0 || m < 0)
    {
        throw std::logic_error("n and m can not be negative");
    }
    Matrix result(n, std::vector<double>(m, 0));
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(min, max); 
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < m; j++)
        {
            result[i][j] = distribution(generator);
        }
    }
    return result;
}

void show(const Matrix& matrix)
{
    for(int i = 0; i < matrix.size(); i++)
    {
        for(int j = 0; j < matrix[0].size(); j++)
        {
            std::cout << std::left << std::setw(5) << std::setprecision(3) << matrix[i][j];
        }
        std::cout << std::endl;
    }

}

Matrix multiply(const Matrix& matrix, double c)
{
    Matrix result = ones(matrix.size(), matrix[0].size());
    for(int i = 0; i < matrix.size(); i++)
    {
        for(int j = 0; j < matrix[0].size(); j++)
        {
           result[i][j] = c * matrix[i][j];
        }
    }
    return result;
}

Matrix multiply(const Matrix& matrix1, const Matrix& matrix2)
{
    if(matrix1.empty() && matrix2.empty())
    {
        return Matrix();
    }
    if(matrix1.empty() || matrix2.empty())
    {
        throw std::logic_error("one matrix is empty while another is not");
    }
    int m = matrix1.size();
    int n1 = matrix1[0].size();
    int n2 = matrix2.size();
    int p = matrix2[0].size();
    if(n1 != n2)
    {
        throw std::logic_error("n1 must be equal to n2");
    }
    Matrix result = zeros(m, p);
    for(int i = 0; i < m; i++)
    {
        for(int j = 0; j < p; j++)
        {
            for(int k = 0; k < n1; k++)
            {
                result[i][j] += matrix1[i][k] * matrix2[k][j];
            }
        }
    }
    return result;
}

Matrix sum(const Matrix& matrix, double c)
{
    if(matrix.empty())
    {
        return Matrix();
    }
    Matrix result = ones(matrix.size(), matrix[0].size());
    for(int i = 0; i < matrix.size(); i++)
    {
        for(int j = 0; j < matrix[0].size(); j++)
        {
            result[i][j] = matrix[i][j] + c;
        }
    }
    return result;
}

Matrix sum(const Matrix& matrix1, const Matrix& matrix2)
{
    if(matrix1.empty() && matrix2.empty())
    {
        return Matrix();
    }
    if(matrix1.empty() || matrix2.empty())
    {
        throw std::logic_error("error!");
    }
    int m1 = matrix1.size();
    int n1 = matrix1[0].size();
    int m2 = matrix2.size();
    int n2 = matrix2[0].size();
    if(n1 != n2 || m1 != m2)
    {
        throw std::logic_error("the size of two matrices don't match");
    }
    Matrix result = ones(m1, n1);
    for(int i = 0; i < m1; i++)
    {
        for(int j = 0; j < n1; j++)
        {
            result[i][j] = matrix1[i][j] + matrix2[i][j];
        }
    }
    return result;
}

Matrix transpose(const Matrix& matrix)
{
    if(matrix.empty())
    {
        return Matrix();
    }
    int m = matrix.size();
    int n = matrix[0].size();
    Matrix result = zeros(n, m);
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < m; j++)
        {
            result[i][j] = matrix[j][i];
        }
    }
    return result;
}

void modify_swap(Matrix& matrix, size_t r1, size_t r2)
{
    matrix[r1].swap(matrix[r2]);
}

void modify_sum(Matrix& matrix, size_t r1, double c, size_t r2)
{
    int m = matrix.size();
    int n = matrix[0].size();
    for(int j = 0; j < n; j++)
    {
        matrix[r2][j] += c * matrix[r1][j];
    }
    // show(matrix);
    // std::cout << "inside" << std::endl;
}

Matrix minor(const Matrix& matrix, size_t n, size_t m)
{
    if(n < 0 || m < 0)
    {
        throw std::logic_error("error!");
    }
    if(matrix.empty() && n == 0 && m == 0)
    {
        return Matrix();
    }
    else if(matrix.empty())
    {
        throw std::logic_error("error!");
    }
    int row = matrix.size();
    int col = matrix[0].size();
    if(n >= row || m >= col)
    {
        throw std::logic_error("error!");
    }
    Matrix result;
    for(int i = 0; i < row; i++)
    {
        if(i == n)
        {
            continue;
        }
        else
        {
            std::vector<double> new_element;
            for(int j = 0; j < col; j++)
            {
                if(j == m)
                {
                    continue;
                }
                else
                {
                    new_element.push_back(matrix[i][j]);
                }
            }
            result.push_back(new_element);
        }
    }
    return result;
}

double determinant(const Matrix& matrix)
{
    if(matrix.empty())
    {
        return 1;
    }
    int m = matrix.size();
    int n = matrix[0].size();
    if(m != n)
    {
        throw std::logic_error("error!");
    }
    double answer = 1.0;
    Matrix result = matrix;
    show(result);
    std::cout << "1111" << std::endl;
    for(int j = 0; j < n - 1; j++)
    {
        int index = j;
        while(result[index][j] == 0)
        {
            index++;
            if(index == m)
            {
                break;
            }
        }
        if(index >= m)
        {
            continue;
        }
        if(index != j)
        {
            std::cout << index << "  " << j <<std::endl;
            modify_swap(result, index, j);
            std::cout << "OK" << std::endl;
            answer *= -1;
        }
        // show(result);
        // std::cout << "222" << std::endl;
        for(int index2 = j + 1; index2 < m; index2++)
        {
            if(result[index2][j] == 0)
            {
                continue;
            }
            double mul = result[index2][j] / result[j][j];
            //std::cout << "num is " << mul << std::endl;
            modify_sum(result, j, -mul, index2);
        }
        show(result);
        std::cout << "one inning" << std::endl;
    }
    show(result);
    std::cout << "222" << std::endl;
    for(int i = 0; i < m; i++)
    {
        answer *= result[i][i];
    }
    return answer;
}   

Matrix inverse(const Matrix& matrix)
{
    if(matrix.empty())
    {
        return Matrix();
    }
    int m = matrix.size();
    int n = matrix[0].size();
    if(m != n)
    {
        throw std::logic_error("inversion does not exist!");
    }
    if(determinant(matrix) == 0)
    {
        throw std::logic_error("error!");
    }
    Matrix I = zeros(n, n);
    for(int i = 0; i < m; i++)
    {
        for(int j = 0; j < n; j++)
        {
            I[i][j] = determinant(minor(matrix, i, j)) * pow(-1, i + j);
        }
    }
    I = transpose(I);
    I = multiply(I, 1.0 / determinant(matrix));
    return I;
}

Matrix concatenate(const Matrix& matrix1, const Matrix& matrix2, int axis = 0)
{
    if(matrix1.empty() && matrix2.empty())
    {
        return Matrix();
    }
    if(matrix1.empty() || matrix2.empty())
    {
        throw std::logic_error("error!");
    }
    int m1 = matrix1.size();
    int n1 = matrix1[0].size();
    int m2 = matrix2.size();
    int n2 = matrix2[0].size();
    if((axis == 0 && n1 != n2) || (axis == 1 && m1 != m2))
    {
        throw std::logic_error("invalid operation!");
    }
    if(axis != 0 && axis != 1)
    {
        throw std::logic_error("axis is not allowed!");
    }
    Matrix result;
    if(axis == 0)
    {
        for(int i = 0; i < m1; i++)
        {
            result.push_back(matrix1[i]);
        }
        for(int i = 0; i < m2; i++)
        {
            result.push_back(matrix2[i]);
        }
        return result;
    }
    for(int i = 0; i < m1; i++)
    {
        std::vector<double> now;
        for(int j = 0; j < n1; j++)
        {
            now.push_back(matrix1[i][j]);
        }
        for(int j = 0; j < n2; j++)
        {
            now.push_back(matrix2[i][j]);
        }
        result.push_back(now);
    }
    return result;

}

Matrix ero_swap(const Matrix& matrix, size_t r1, size_t r2)
{
    if(matrix.empty())
    {
        throw std::logic_error("error!");
    }
    int m = matrix.size();
    int n = matrix[0].size();
    if(r1 >= m || r2 >= m)
    {
        throw std::logic_error("error!");
    }
    if(r1 == r2)
    {
        return matrix;
    }
    Matrix result;
    for(int i = 0; i < m; i++)
    {
        if(i == r1)
        {
            std::vector<double> now = matrix[r2];
            result.push_back(now);
        }
        else if(i == r2)
        {
            std::vector<double> now = matrix[r1];
            result.push_back(now);
        }
        else
        {
            std::vector<double> now = matrix[i];
            result.push_back(now);
        }
    }
    return result;
}

Matrix ero_multiply(const Matrix& matrix, size_t r, double c)
{
    if(matrix.empty())
    {
        throw std::logic_error("error!");
    }
    int m = matrix.size();
    int n = matrix[0].size();
    if(r >= m)
    {
        throw std::logic_error("error!");
    }
    Matrix result = matrix;
    for(int j = 0; j < n; j++)
    {
        result[r][j] *= c;
    }
    return result;
}

Matrix ero_sum(const Matrix& matrix, size_t r1, double c, size_t r2)
{
    if(matrix.empty())
    {
        throw std::logic_error("error!");
    }
    int m = matrix.size();
    int n = matrix[0].size();
    Matrix result = matrix;
    for(int j = 0; j < n; j++)
    {
        result[r2][j] += c * matrix[r1][j];
    }
    return result;
}

Matrix upper_triangular(const Matrix& matrix)
{
    // show(matrix);
    // std::cout << "111" << std::endl;
    if(matrix.empty())
    {
        return Matrix();
    }
    Matrix result = matrix;
    int m = result.size();
    int n = result[0].size();
    if(m != n)
    {
        throw std::logic_error("error!");
    }
    for(int j = 0; j < n; j++)
    {
        int index = j;
        while(result[index][j] == 0)
        {
            index++;
        }
        if(index == m)
        {
            continue;
        }
        if(index != j)
        {
            modify_swap(result, index, j);
        }
        // show(result);
        // std::cout << "222" << std::endl;
        for(int index2 = j + 1; index2 < m; index2++)
        {
            if(result[index2][j] == 0)
            {
                continue;
            }
            double mul = result[index2][j] / result[j][j];
            std::cout << "num is " << mul << std::endl;
            modify_sum(result, j, -mul, index2);
        }
        // show(result);
        // std::cout << "333" << std::endl;
    }
    return result;
}

}//namespace algebra