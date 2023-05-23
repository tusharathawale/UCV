
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <assert.h>
#include "./cstm_gsl.hpp"

// for symetric matrix

double in4_0[4][4] = {
    {1, -1, 4, 1},
    {1, 4, -2, 1},
    {1, 4, 2, 1},
    {1, -1, 0, 1}};

double in4_1[4][4] = {
    {0.20, 0.60, 0.40, 0.80},
    {0.60, 1.80, 1.20, 2.40},
    {0.40, 1.20, 0.80, 1.60},
    {0.80, 2.40, 1.60, 3.20}};

double in4_2[4][4] = {
    {1.0, 3.0, 7.0, 8.0},
    {3.0, 2.0, 6.0, 7.0},
    {7.0, 6.0, 5.0, 6.0},
    {8.0, 7.0, 6.0, 5.0}};

double in4_3[4][4] = {
    {0.0, 0.0, 0.0, 0.0},
    {0.0, 20.0, 60.0, 40.0},
    {0.0, 60.0, 180.0, 120.0},
    {0.0, 40.0, 120.0, 80.0}};

int equal_double(double a, double b)
{
    if (fabs(a - b) < 0.0001)
    {
        return 1;
    }
    return 0;
}

int equal_matrix(gsl_matrix *a, gsl_matrix *b)
{
    if (a->size1 != b->size1 || a->size2 != b->size2)
    {
        return 0;
    }

    for (int i = 0; i < a->size1; i++)
    {
        for (int j = 0; j < a->size2; j++)
        {
            double v1 = gsl_matrix_get(a, i, j);
            double v2 = gsl_matrix_get(b, i, j);
            if (equal_double(v1, v2) == 0)
            {
                return 0;
            }
        }
    }
    return 1;
}

void cstm_gsl_matrix_show(gsl_matrix *x)
{
    for (int i = 0; i < x->size1; i++)
    {
        for (int j = 0; j < x->size2; j++)
        {
            if (j != 0)
            {
                printf(",");
            }
            printf("%f", gsl_matrix_get(x, i, j));
        }
        printf("\n");
    }
}

void cstm_gsl_vector_show(gsl_vector *x)
{
    for (int i = 0; i < x->size; i++)
    {
        if (i != 0)
        {
            printf(",");
        }
        printf("%f", gsl_vector_get(x, i));
    }
    printf("\n");
}

void basic_qr_test()
{
    printf("---basic_qr_test\n");

    int msize = 4;
    // mat_t x;
    gsl_matrix *x = gsl_matrix_alloc(msize, msize);

    for (int i = 0; i < msize; i++)
    {
        for (int j = 0; j < msize; j++)
        {
            // x.v[i][j] = in4_0[i][j];
            gsl_matrix_set(x, i, j, in4_0[i][j]);
        }
    }

    puts("original matrix");
    // matrix_show(&x);
    cstm_gsl_matrix_show(x);

    gsl_vector *tau = gsl_vector_alloc(4);
    gsl_linalg_QR_decomp(x, tau);

    puts("tau results");
    cstm_gsl_vector_show(tau);
    puts("x after qr call");
    cstm_gsl_matrix_show(x);
    // it seems that it is unnecessary to compute matirx Q
    // in QR decomposition
    // refer to https://na-inet.jp/na/gslsample/linear_system_qr.html

    /*
    householder(&x, &R, &Q);

    puts("Q");
    matrix_show(&Q);
    puts("R");
    matrix_show(&R);

    //check orthoganal
    bool ifupper = matrix_is_upper_triangular(&R,0.00001);
    assert(ifupper==true);

    // to show their product is the input matrix
    mat_t m = matrix_mul(&Q, &R);
    puts("Q * R");
    matrix_show(&m);

    assert(equal_matrix(&m,&x)==1);
    */

    gsl_matrix_free(x);
    gsl_vector_free(tau);
}

gsl_vector *gsl_matrix_mul_vec_add_vec(gsl_matrix *A, gsl_vector *U, gsl_vector *M)
{
    assert(A->size1 == U->size);
    assert(A->size2 == M->size);

    gsl_vector *AUM = gsl_vector_alloc(4);

    for (int i = 0; i < A->size1; i++)
    {
        // adding vector into the matrix
        for (int j = 0; j < A->size2; j++)
        {
            // AUM.v[i] += (A->v[i][j] * U->v[j]);
            double v1 = gsl_vector_get(AUM, i);
            double v2 = gsl_matrix_get(A, i, j) * gsl_vector_get(U, j);
            gsl_vector_set(AUM, i, v1 + v2);
        }
        // AUM.v[i] += M->v[i];
        double v1 = gsl_vector_get(AUM, i);
        gsl_vector_set(AUM, i, v1 + gsl_vector_get(M, i));
    }

    return AUM;
}

void test_basic_operations()
{

    printf("---test_basic_operations");

    gsl_matrix *a = gsl_matrix_alloc(4, 4);

    gsl_vector *u = gsl_vector_alloc(4);
    gsl_vector *m = gsl_vector_alloc(4);

    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            gsl_matrix_set(a, i, j, i + 1);
        }
        gsl_vector_set(u, i, 1);
        gsl_vector_set(m, i, i + 1);
    }

    gsl_vector *auv = gsl_matrix_mul_vec_add_vec(a, u, m);

    printf("\n");
    gsl_vector_fprintf(stdout, auv, "%g");

    for (int i = 0; i < 4; i++)
    {
        // assert(auv.v[i] == (i + 1) * 5);
        assert(gsl_vector_get(auv, i) == (i + 1) * 5);
    }

    gsl_matrix_free(a);
    gsl_vector_free(u);
    gsl_vector_free(m);
    gsl_vector_free(auv);
}

// refer to http://gnu.ist.utl.pt/software/gsl/manual/html_node/Eigenvalue-and-Eigenvector-Examples.html
void eigen_values_4by4()
{
    printf("--eigen_values_4by4\n");

    gsl_matrix *m = gsl_matrix_alloc(4, 4);

    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            gsl_matrix_set(m, i, j, in4_1[i][j]);
        }
    }

    gsl_vector *eval = gsl_vector_alloc(4);
    gsl_matrix *evec = gsl_matrix_alloc(4, 4);
    gsl_eigen_symmv_workspace *w =
        gsl_eigen_symmv_alloc(4);

    // double result[4]={0};
    // eigen_solve_eigenvalues(&x, 0.0001, 20, result);

    gsl_eigen_symmv(m, eval, evec, w);

    gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_VAL_DESC);

    cstm_gsl_vector_show(eval);

    // the assert only works for debug case
    assert(equal_double(gsl_vector_get(eval, 0), 6.0) == 1);
    assert(equal_double(gsl_vector_get(eval, 1), 0.0) == 1);
    assert(equal_double(gsl_vector_get(eval, 2), 0.0) == 1);
    assert(equal_double(gsl_vector_get(eval, 3), 0.0) == 1);

    gsl_eigen_symmv_free(w);
    gsl_matrix_free(m);
    gsl_vector_free(eval);
    gsl_matrix_free(evec);
}

void eigen_vectors_4by4()
{
    printf("---test eigen_vectors_4by4\n");
    gsl_matrix *m = gsl_matrix_alloc(4, 4);

    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            gsl_matrix_set(m, i, j, in4_3[i][j]);
        }
    }
    puts("original matrix");
    cstm_gsl_matrix_show(m);

    gsl_vector *eval = gsl_vector_alloc(4);
    gsl_matrix *evec = gsl_matrix_alloc(4, 4);
    gsl_eigen_symmv_workspace *w =
        gsl_eigen_symmv_alloc(4);

    gsl_eigen_symmv(m, eval, evec, w);

    gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_VAL_DESC);

    puts("eigen vectors");
    cstm_gsl_vector_show(eval);
    puts("eigen values");
    cstm_gsl_matrix_show(evec);

    assert(equal_double(gsl_vector_get(eval, 0), 280.0) == 1);
    assert(equal_double(gsl_vector_get(eval, 1), 0.0) == 1);
    assert(equal_double(gsl_vector_get(eval, 2), 0.0) == 1);
    assert(equal_double(gsl_vector_get(eval, 3), 0.0) == 1);

    assert(equal_double(gsl_matrix_get(evec, 0, 0), 0.0) == 1);
    assert(equal_double(gsl_matrix_get(evec, 1, 0), -0.2672) == 1);
    assert(equal_double(gsl_matrix_get(evec, 2, 0), -0.8017) == 1);
    assert(equal_double(gsl_matrix_get(evec, 3, 0), -0.5345) == 1);

    gsl_eigen_symmv_free(w);
    gsl_matrix_free(evec);
    gsl_vector_free(eval);
}

// refer to
// https://gist.github.com/bjd2385/7f4685e703f7437e513608f41c65bbd7
// https://www.ibm.com/docs/en/essl/6.3?topic=mos-sgemm-dgemm-cgemm-zgemm-combined-matrix-multiplication-addition-general-matrices-their-transposes-conjugate-transposes
// https://lists.gnu.org/archive/html/help-gsl/2008-11/msg00001.html
void test_invert_4by4matrix()
{
    int dim = 3;

    double a_data[] = {1.0, 0.6, 0.0,
                       0.0, 1.5, 1.0,
                       0.0, 1.0, 1.0};
    // will the matrix view free the memory space automatically?
    // refer to
    // http://gnu.ist.utl.pt/software/gsl/manual/html_node/Matrix-views.html
    // view is just the pointer to original continuous array
    // it does not own its own memory space
    // the actual memory can be either on stack or heap
    gsl_matrix_view x = gsl_matrix_view_array(a_data, dim, dim);

    gsl_matrix *x_copy = gsl_matrix_alloc(dim, dim);

    // deep copy from the src matrix into the dest matrix
    gsl_matrix_memcpy(x_copy, &x.matrix);

    // save original x
    // the value of x is updated after the LU decomposition

    puts("x_copy");

    // matrix_show(&x);
    cstm_gsl_matrix_show(x_copy);

    // mat_t x;
    // gsl_matrix *x = gsl_matrix_alloc(dim, dim);
    // mat_t x_inv;

    // for (int i = 0; i < dim; i++)
    //{
    //     for (int j = 0; j < dim; j++)
    //     {
    // x.v[i][j] = in4_2[i][j];
    //        gsl_matrix_set(&x.matrix, i, j, in4_2[i][j]);
    //    }
    //}

    gsl_permutation *p = gsl_permutation_alloc(dim);
    int s;

    printf("sign %d\n", s);

    // Compute the  inverse of the LU decomposition
    gsl_matrix *x_inv = gsl_matrix_alloc(dim, dim);

    // Compute the LU decomposition of this matrix
    gsl_linalg_LU_decomp(&x.matrix, p, &s);

    gsl_linalg_LU_invert(&x.matrix, p, x_inv);

    gsl_permutation_free(p);

    puts("x_inv");
    cstm_gsl_matrix_show(x_inv);

    gsl_matrix *result = gsl_matrix_alloc(dim, dim);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
                   1.0, x_copy, x_inv,
                   0.0, result);

    // mat_t c1 = matrix_mul(&x, &x_inv);

    puts("x_copy*x_inv");
    cstm_gsl_matrix_show(result);

    for (int i = 0; i < dim; i++)
    {
        for (int j = 0; j < dim; j++)
        {
            if (i == j)
            {
                assert(equal_double(gsl_matrix_get(result, i, j), 1.0) == 1);
            }
            else
            {
                assert(equal_double(gsl_matrix_get(result, i, j), 0.0) == 1);
            }
        }
    }
    gsl_matrix_free(result);
    gsl_matrix_free(x_copy);
}

void test_eigen_vectors_decomposition1()
{
    printf("---test_eigen_vectors_decomposition1\n");

    int dim = 4;
    // mat_t x;
    gsl_matrix *x = gsl_matrix_alloc(dim, dim);
    for (int i = 0; i < dim; i++)
    {
        for (int j = 0; j < dim; j++)
        {
            // x.v[i][j] = in4_3[i][j];
            gsl_matrix_set(x, i, j, in4_3[i][j]);
        }
    }

    //printf("input x\n");
    //cstm_gsl_matrix_show(x);

    // the value of x will be changed during the computation
    // copy another one for reulsts checking
    gsl_matrix *x_copy = gsl_matrix_alloc(4, 4);
    gsl_matrix_memcpy(x_copy, x);

    // the process of eigen vector decomposition
    // 1 eigen value and vector
    gsl_vector *eval = gsl_vector_alloc(4);
    gsl_matrix *evec = gsl_matrix_alloc(4, 4);
    gsl_eigen_symmv_workspace *w =
        gsl_eigen_symmv_alloc(4);

    // double result[4]={0};
    // eigen_solve_eigenvalues(&x, 0.0001, 20, result);
    gsl_eigen_symmv(x, eval, evec, w);

    // 2 create diagonal matrix based on eigen value
    gsl_matrix *dia_matrix = gsl_matrix_alloc(4, 4);
    for (int i = 0; i < 4; i++)
    {
        double sqrtv = sqrt(gsl_vector_get(eval, i));
        gsl_matrix_set(dia_matrix, i, i, sqrtv);
    }

    // 3 matrix multiplication
    //  mat_t A = matrix_mul(&eigen_vectors, &diag);
    gsl_matrix *mul_result = gsl_matrix_alloc(4, 4);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
                   1.0, evec, dia_matrix,
                   0.0, mul_result);

    // cstm_gsl_matrix_show(mul_result);
    gsl_matrix *mul_result_transpose = gsl_matrix_alloc(dim, dim);
    gsl_matrix_memcpy(mul_result_transpose, mul_result);
    gsl_matrix_transpose(mul_result_transpose);

    // evaluate results

    gsl_matrix *mul_result_2 = gsl_matrix_alloc(4, 4);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
                   1.0, mul_result, mul_result_transpose,
                   0.0, mul_result_2);

    printf("original results\n");
    cstm_gsl_matrix_show(x_copy);

    printf("reconstructed results\n");
    cstm_gsl_matrix_show(mul_result_2);

    assert(equal_matrix(mul_result_2, x_copy) == 1);

    gsl_matrix_free(x);
    gsl_matrix_free(dia_matrix);
    gsl_matrix_free(mul_result);
}


void test_eigen_vectors_decomposition2(){
 printf("---test_eigen_vectors_decomposition2\n");

    int dim = 4;
    // mat_t x;
    gsl_matrix *x = gsl_matrix_alloc(dim, dim);
    for (int i = 0; i < dim; i++)
    {
        for (int j = 0; j < dim; j++)
        {
            // x.v[i][j] = in4_3[i][j];
            gsl_matrix_set(x, i, j, in4_3[i][j]);
        }
    }

    gsl_matrix *x_copy = gsl_matrix_alloc(4, 4);
    gsl_matrix_memcpy(x_copy, x);
    
    //allocate new matrix in the function
    gsl_matrix*decomp_result = UCVMATH_CSTM_GSL::gsl_eigen_vector_decomposition(x);

    gsl_matrix *decomp_result_transpose = gsl_matrix_alloc(dim, dim);
    gsl_matrix_memcpy(decomp_result_transpose, decomp_result);
    gsl_matrix_transpose(decomp_result_transpose);

    // evaluate results
    gsl_matrix *decomp_result_check = gsl_matrix_alloc(4, 4);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
                   1.0, decomp_result, decomp_result_transpose,
                   0.0, decomp_result_check);

    printf("original results\n");
    cstm_gsl_matrix_show(x_copy);

    printf("reconstructed results\n");
    cstm_gsl_matrix_show(decomp_result_check);

    assert(equal_matrix(decomp_result_check, x_copy) == 1);

    gsl_matrix_free(x);
    gsl_matrix_free(decomp_result_check);
}

int main()
{
    test_basic_operations();
    basic_qr_test();
    eigen_values_4by4();
    eigen_vectors_4by4();

    test_invert_4by4matrix();
    test_eigen_vectors_decomposition1();
    test_eigen_vectors_decomposition2();

    return 0;
}