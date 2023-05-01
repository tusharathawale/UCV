
#include <gsl/gsl_linalg.h>
#include <assert.h>

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

int equal_matrix(gsl_matrix_view a, gsl_matrix_view b)
{
    return gsl_matrix_equal(&a.matrix, &b.matrix);
}

/*
void eigen_values_4by4()
{
    printf("--eigen_values_4by4\n");

    mat_t x ;
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            x.v[i][j] = in4_1[i][j];
        }
    }

    double result[4]={0};
    eigen_solve_eigenvalues(&x, 0.0001, 20, result);

    for (int i = 0; i < 4; i++)
    {
        printf("%f ", result[i]);
    }

    printf("\n");
    // the assert only works for debug case
    assert(equal_double(result[0], 6.0) == 1);
    assert(equal_double(result[1], 0.0) == 1);
    assert(equal_double(result[2], 0.0) == 1);
    assert(equal_double(result[3], 0.0) == 1);

}

void basic_qr_test()
{
    printf("---basic_qr_test\n");
    mat_t R, Q;
    int msize = 4;
    mat_t x;

    for (int i = 0; i < msize; i++)
    {
        for (int j = 0; j < msize; j++)
        {
            x.v[i][j] = in4_0[i][j];
        }
    }

    puts("original matrix");
    matrix_show(&x);

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
}

void test_invert_4by4matrix()
{
    int dim = 4;
    mat_t x;
    mat_t x_inv;
    for (int i = 0; i < dim; i++)
    {
        for (int j = 0; j < dim; j++)
        {
            x.v[i][j] = in4_2[i][j];
        }
    }

    bool inv_ok = invert4by4matrix(&x, &x_inv);

    puts("x");
    matrix_show(&x);
    puts("x_inv");
    matrix_show(&x_inv);

    mat_t c1 = matrix_mul(&x, &x_inv);

    puts("x*x_inv");
    matrix_show(&c1);

    mat_t c2 = matrix_mul(&x_inv, &x);

    puts("x_inv*x");
    matrix_show(&c2);

    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            if (i == j)
            {
                assert(equal_double(c1.v[i][j], 1.0) == 1);
                assert(equal_double(c1.v[i][j], 1.0) == 1);
            }
            else
            {
                assert(equal_double(c2.v[i][j], 0.0) == 1);
                assert(equal_double(c2.v[i][j], 0.0) == 1);
            }
        }
    }



    assert(inv_ok == true);
}

void eigen_vectors_4by4()
{

    printf("---test eigen_vectors_4by4\n");
    int dim = 4;
    mat_t x;
    for (int i = 0; i < dim; i++)
    {
        for (int j = 0; j < dim; j++)
        {
            x.v[i][j] = in4_3[i][j];
        }
    }

    matrix_show(&x);

    mat_t R, Q;
    householder(&x, &R, &Q);

    puts("Q");
    matrix_show(&Q);
    puts("R");
    matrix_show(&R);

    mat_t m = matrix_mul(&Q, &R);
    puts("Q * R");
    matrix_show(&m);

    // the assert only works for debug case

    double result[4]={0};
    eigen_solve_eigenvalues(&x, 0.0001, 20, result);
    for (int i = 0; i < 4; i++)
    {
        printf(" %f", result[i]);
    }
    printf("\n");
    assert(equal_double(result[0], 0) == 1);
    assert(equal_double(result[1], 280.0) == 1);
    assert(equal_double(result[2], 0.0) == 1);
    assert(equal_double(result[3], 0.0) == 1);

    mat_t eigen_vectors = eigen_solve_eigen_vectors(&x, result, 4, 4, 20);
    matrix_show(&eigen_vectors);
    assert(equal_double(eigen_vectors.v[0][1], 0.0) == 1);
    assert(equal_double(eigen_vectors.v[1][1], 0.2672) == 1);
    assert(equal_double(eigen_vectors.v[2][1], 0.8017) == 1);
    assert(equal_double(eigen_vectors.v[3][1], 0.5345) == 1);

}

void test_eigen_vectors_decomposition()
{
    printf("---test_eigen_vectors_decomposition\n");

    int dim = 4;
    mat_t x;
    for (int i = 0; i < dim; i++)
    {
        for (int j = 0; j < dim; j++)
        {
            x.v[i][j] = in4_3[i][j];
        }
    }

    mat_t A = eigen_vector_decomposition(&x);

    mat_t A_trans =A;
    matrix_transpose(&A_trans);

    mat_t rst = matrix_mul(&A, &A_trans);

    puts("A");
    matrix_show(&A);
    puts("rst");
    matrix_show(&rst);

    assert(equal_matrix(&rst, &x) == 1);

}
*/

gsl_vector* gsl_matrix_mul_vec_add_vec(gsl_matrix *A, gsl_vector *U, gsl_vector *M){
        assert(A->size1 == U->size);
        assert(A->size2 == M->size);

        gsl_vector *AUM = gsl_vector_alloc(4);

        for (int i = 0; i < A->size1; i++)
        {
            // adding vector into the matrix
            for (int j = 0; j < A->size2; j++)
            {
                //AUM.v[i] += (A->v[i][j] * U->v[j]);
                double v1 = gsl_vector_get(AUM,i);
                double v2 = gsl_matrix_get(A,i,j) * gsl_vector_get(U,j);
                gsl_vector_set(AUM,i,v1+v2);

            }
            //AUM.v[i] += M->v[i];
            double v1 = gsl_vector_get(AUM,i);
            gsl_vector_set(AUM,i,v1+gsl_vector_get(M,i));   
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
            gsl_matrix_set(a,i,j,i+1);
        }
        gsl_vector_set(u,i,1);
        gsl_vector_set(m,i,i+1);
    }

    gsl_vector* auv = gsl_matrix_mul_vec_add_vec(a, u, m);
    
    printf("\n");
    gsl_vector_fprintf(stdout, auv, "%g");

    for (int i = 0; i < 4; i++)
    {
        //assert(auv.v[i] == (i + 1) * 5);
        assert(gsl_vector_get(auv,i)==(i + 1) * 5);
    }

    gsl_matrix_free(a);
    gsl_vector_free(u);
    gsl_vector_free(m);
    gsl_vector_free(auv);
}

int main()
{
    test_basic_operations();
    // basic_qr_test();
    // eigen_values_4by4();
    // test_invert_4by4matrix();
    // eigen_vectors_4by4();
    // test_eigen_vectors_decomposition();
    return 0;
}