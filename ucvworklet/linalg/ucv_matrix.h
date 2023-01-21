#ifndef UCV_MATRIX_H
#define UCV_MATRIX_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <stdbool.h>
#include <cassert>

#include <vtkm/Types.h>

namespace UCVMATH{

typedef struct
{
    int m, n; // m is row, n is column
    double **v;
} mat_t, *mat;

VTKM_EXEC inline mat matrix_new(int m, int n)
{
    mat x = (mat)malloc(sizeof(mat_t));
    x->v = (double **)malloc(sizeof(double *) * m);
    // the init memory are all zero
    x->v[0] = (double *)calloc(sizeof(double), m * n);
    for (int i = 0; i < m; i++)
        x->v[i] = x->v[0] + n * i;
    x->m = m;
    x->n = n;
    return x;
}

VTKM_EXEC inline mat matrix_new_eye(int m, int n)
{
    mat x = (mat)malloc(sizeof(mat_t));
    x->v = (double **)malloc(sizeof(double *) * m);
    x->v[0] = (double *)calloc(sizeof(double), m * n);
    for (int i = 0; i < m; i++)
    {
        x->v[i] = x->v[0] + n * i;
        x->v[i][i] = 1;
    }
    x->m = m;
    x->n = n;
    return x;
}

VTKM_EXEC inline void matrix_delete(mat m)
{
    free(m->v[0]);
    free(m->v);
    free(m);
}

VTKM_EXEC inline void matrix_transpose(mat m)
{
    for (int i = 0; i < m->m; i++)
    {
        for (int j = 0; j < i; j++)
        {
            double t = m->v[i][j];
            m->v[i][j] = m->v[j][i];
            m->v[j][i] = t;
        }
    }
}

// C++ does not support variable-length arrays
/*
inline mat matrix_copy(int n,  int m, double a[m][n])
{
    mat x = matrix_new(m, n);
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            x->v[i][j] = a[i][j];
    return x;
}
*/

VTKM_EXEC inline mat matrix_copy_mat(mat x)
{
    mat y = matrix_new(x->m, x->n);
    for (int i = 0; i < x->m; i++)
        for (int j = 0; j < x->n; j++)
            y->v[i][j] = x->v[i][j];
    return y;
}

VTKM_EXEC inline void matrix_mul_toz(mat x, mat y, mat *z)
{
    if (x->n != y->m)
    {
        printf("error for matrix_mul_toz for dims");
        return;
    }

    for (int i = 0; i < x->m; i++)
        for (int j = 0; j < y->n; j++)
            for (int k = 0; k < x->n; k++)
                (*z)->v[i][j] += x->v[i][k] * y->v[k][j];
}

VTKM_EXEC inline mat matrix_mul(mat x, mat y)
{
    if (x->n != y->m)
        return 0;
    mat r = matrix_new(x->m, y->n);
    for (int i = 0; i < x->m; i++)
        for (int j = 0; j < y->n; j++)
            for (int k = 0; k < x->n; k++)
                r->v[i][j] += x->v[i][k] * y->v[k][j];
    return r;
}

VTKM_EXEC inline mat matrix_minor(mat x, int d)
{
    mat m = matrix_new(x->m, x->n);
    for (int i = 0; i < d; i++)
        m->v[i][i] = 1;
    for (int i = d; i < x->m; i++)
        for (int j = d; j < x->n; j++)
            m->v[i][j] = x->v[i][j];
    return m;
}

/* c = a + b * s */
VTKM_EXEC inline double *vmadd(double a[], double b[], double s, double c[], int n)
{
    for (int i = 0; i < n; i++)
        c[i] = a[i] + s * b[i];
    return c;
}

/* m = I - v v^T */
VTKM_EXEC inline mat vmul(double v[], int n)
{
    mat x = matrix_new(n, n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            x->v[i][j] = -2 * v[i] * v[j];
    for (int i = 0; i < n; i++)
        x->v[i][i] += 1;

    return x;
}

/* ||x|| */
VTKM_EXEC inline double vnorm(double x[], int n)
{
    double sum = 0;
    for (int i = 0; i < n; i++)
        sum += x[i] * x[i];
    return sqrt(sum);
}

/* y = x / d */
VTKM_EXEC inline double *vdiv(double x[], double d, double y[], int n)
{
    for (int i = 0; i < n; i++)
        y[i] = x[i] / d;
    return y;
}

/* take c-th column of m, put in v */
VTKM_EXEC inline double *mcol(mat m, double *v, int c)
{
    for (int i = 0; i < m->m; i++)
        v[i] = m->v[i][c];
    return v;
}

VTKM_EXEC inline void matrix_show(mat m)
{
    for (int i = 0; i < m->m; i++)
    {
        for (int j = 0; j < m->n; j++)
        {
            printf(" %8.3f", m->v[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

VTKM_EXEC inline void householder(mat m, mat *R, mat *Q)
{
    mat q[m->m];
    mat z = m, z1;
    for (int k = 0; k < m->n && k < m->m - 1; k++)
    {
        double e[m->m], x[m->m], a;
        z1 = matrix_minor(z, k);
        if (z != m)
            matrix_delete(z);
        z = z1;

        mcol(z, x, k);
        a = vnorm(x, m->m);
        if (m->v[k][k] > 0)
            a = -a;

        for (int i = 0; i < m->m; i++)
            e[i] = (i == k) ? 1 : 0;

        vmadd(x, e, a, e, m->m);
        vdiv(e, vnorm(e, m->m), e, m->m);
        q[k] = vmul(e, m->m);
        z1 = matrix_mul(q[k], z);
        if (z != m)
            matrix_delete(z);
        z = z1;
    }
    matrix_delete(z);
    *Q = q[0];
    *R = matrix_mul(q[0], m);
    for (int i = 1; i < m->n && i < m->m - 1; i++)
    {
        z1 = matrix_mul(q[i], *Q);
        if (i > 1)
            matrix_delete(*Q);
        *Q = z1;
        matrix_delete(q[i]);
    }
    matrix_delete(q[0]);
    z = matrix_mul(*Q, m);
    matrix_delete(*R);
    *R = z;
    matrix_transpose(*Q);
}

VTKM_EXEC inline bool matrix_is_upper_triangular(mat m, double tol)
{
    // For now, only treat square matricies. */
    assert(m->m == m->n);
    for (int i = 0; i < m->m; i++)
    {
        for (int j = 0; j < i; j++)
        {
            if (fabs(m->v[i][j]) > tol)
            {
                return false;
            }
        }
    }
    return true;
}

// only tested for n*n matrix and it is symetric
VTKM_EXEC inline void eigen_solve_eigenvalues(mat x, double tol, int max_iter, double *eigen_array)
{
    // refer to https://www.andreinc.net/2021/01/25/computing-eigenvalues-and-eigenvectors-using-qr-decomposition
    mat ak = matrix_copy_mat(x);
    mat qq = matrix_new_eye(x->m, x->m);
    int i = 0;
    while (true)
    {
        mat R, Q;
        // the hoiseholder will assign new r and q each time
        householder(ak, &R, &Q);
        mat m = matrix_mul(R, Q);
        // matrix_mul_toz(R, Q, &ak);
        ak = matrix_copy_mat(m);
        // puts("Q");
        // matrix_show(Q);
        // puts("R");
        // matrix_show(R);
        // puts("m");
        // matrix_show(m);

        mat newq = matrix_mul(qq, Q);
        qq = matrix_copy_mat(newq);

        matrix_delete(R);
        matrix_delete(Q);
        matrix_delete(newq);
        i++;
        if (matrix_is_upper_triangular(ak, tol) || i > max_iter)
        {
            // matrix_show(m);
            // printf("iter %d\n", i);
            break;
        }
        // delete m when it is not qualified
        matrix_delete(m);
    }

    // TODO, only consider n*n matrix now
    for (int i = 0; i < ak->m; i++)
    {
        eigen_array[i] = ak->v[i][i];
    }
}

}

#endif