#include <vector>
#include <assert.h>>
#include <iostream>
#include <cmath>

using Matrix = std::vector<std::vector<double>>;

Matrix mutiply(Matrix& a, Matrix& b)
{
    int m = a.size();
    int n1 = a[0].size();
    int n2 = b.size();
    int p = b[0].size();
    assert(n1 == n2);
    Matrix result(m, std::vector<double>(p, 0));
    for(int i = 0; i < m; i++)
    {
        for(int j = 0; j < p; j++)
        {
            for(int k = 0; k < n1; k++)
            {
                result[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    return result;
}

Matrix transpose(Matrix& a)
{
    int m = a.size();
    int n = a[0].size();
    Matrix result(m, std::vector<double>(n, 0));
    for(int i = 0; i < m; i++)
    {
        for(int j = 0; j < n; j++)
        {
            result[i][j] = a[j][i];
        }
    }
    return result;
}

double det(Matrix& a)
{
    Matrix b = a;
    int m = a.size();
    int n = a[0].size();
    if(m != n)
    {
        return 0;
    }
    int temp[m];
    int index;
    int total = 1;
    int det = 1;
    for(int i = 0; i < n; i++)
    {
        index = i;
        while(index < m)
        {
            if(a[index][i] != 0)
            {
                break;
            }
            index++;
        }
        if(index == n)
        {
            continue;
        }
        if(index != i)
        {
            det *= std::pow(-1, index - i);
            for(int s = 0; s < n; s++)
            {
                std::swap(a[i][s], a[index][s]);
            }
        }
        for(int r = i + 1; r < m; r++)
        {
            if(a[r][i] != 0)
            {
                int mul = a[r][i] / a[i][i];
                for(int k = 0; k < n; k++)
                {
                    a[r][k] = a[i][i] * a[i][k] - a[i][r] * a[i][k];
                }
            }
            total *= a[i][i];
        }
    }
    for(int i = 0; i < m; i++)
    {
        det *= a[i][i];
    }
    return det / total;
}