#ifndef CSTM_GSL
#define CSTM_GSL

#include <vtkm/Types.h>

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>

namespace UCVMATH_CSTM_GSL
{

    VTKM_EXEC inline void cstm_gsl_matrix_show(gsl_matrix *x)
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

    VTKM_EXEC inline void cstm_gsl_vector_show(gsl_vector *x)
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

    VTKM_EXEC inline gsl_matrix* gsl_eigen_vector_decomposition(gsl_matrix *x)
    {
        assert(x->size1 == x->size2);
        int dim = x->size1;

        // the process of eigen vector decomposition
        // 1 eigen value and vector
        gsl_vector *eval = gsl_vector_alloc(dim);
        gsl_matrix *evec = gsl_matrix_alloc(dim, dim);
        gsl_eigen_symmv_workspace *w =
            gsl_eigen_symmv_alloc(4);

        // double result[4]={0};
        // eigen_solve_eigenvalues(&x, 0.0001, 20, result);
        gsl_eigen_symmv(x, eval, evec, w);

        // 2 create diagonal matrix based on eigen value
        gsl_matrix *dia_matrix = gsl_matrix_alloc(dim, dim);
        for (int i = 0; i < 4; i++)
        {
            double sqrtv = sqrt(gsl_vector_get(eval, i));
            gsl_matrix_set(dia_matrix, i, i, sqrtv);
        }

        // 3 matrix multiplication
        //  mat_t A = matrix_mul(&eigen_vectors, &diag);
        gsl_matrix *mul_result = gsl_matrix_alloc(dim, dim);
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
                       1.0, evec, dia_matrix,
                       0.0, mul_result);

        gsl_matrix_free(dia_matrix);

        return mul_result;
    }

    //other wrapper function for the gsl library
    VTKM_EXEC inline gsl_vector* cstm_gsl_vector_alloc(const size_t n){
        return gsl_vector_alloc(n);
    }

    VTKM_EXEC inline void cstm_gsl_vector_free(gsl_vector * v){
        return gsl_vector_free(v);
        
    }

}

#endif
