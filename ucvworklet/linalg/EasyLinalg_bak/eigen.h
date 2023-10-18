// eigen value and eigen vectors

#ifndef EASY_LIALG_EIGEN
#define EASY_LIALG_EIGEN

#include "./basic.h"
#include "./blas.h"

namespace EASYLINALG
{

    template <typename T,
              uint NumRow,
              uint NumCol>
    LIAG_FUNC_MACRO void Householder(
        const Matrix<T, NumRow, NumCol> &A,
        Matrix<T, NumRow, NumCol> &Q,
        Matrix<T, NumRow, NumCol> &R)
    {
        // current code only works for the symetry matrix
        assert(NumRow == NumCol);
        const int msize = NumRow;
        // Matrix<T, NumRow, NumCol> MList[msize];
        Matrix<T, NumRow, NumCol> MTemp;

        // operator both for const and non-const
        // MList[0] = A;
        MTemp = A;

        R = A;
        double ANorm;

        Vec<T, msize> Ve;
        Ve.InitZero();
        Vec<T, msize> Va;
        Va.InitZero();

        for (uint k = 1; k < msize; k++)
        {
            auto HInit = R.GetMinor(k - 1);
            // printf("hhd iter%d\n", k);
            // HInit.Show();

            Va = HInit.GetColum(k - 1);

            // puts("a colum");
            // Va.Show();

            ANorm = Va.GetNorm();

            // printf("a norm %f\n", ANorm);

            if (A[k - 1][k - 1] > 0)
                ANorm = -ANorm;

            for (uint i = 0; i < msize; i++)
                Ve[i] = (i == (k - 1)) ? 1 : 0;

            // puts("ve is");
            // Ve.Show();
            // puts("va is");
            // Va.Show();

            // Ve = ANorm*Ve+Va
            Ve = AXPY(ANorm, Ve, Va);
            // puts("after vmadd v vec");
            // Ve.Show();

            // TODO, only do this when norm is not zero
            if (Ve.GetNorm() > 0.00001)
            {
                Ve.Div(Ve.GetNorm());
            }

            // puts("after vdiv v vec");
            // Ve.Show();

            // I-v*v
            // MList[k] = IMINUSVVT(Ve);
            MTemp = IMINUSVVT(Ve);
            // puts("mat_list[k] H is");
            // MTemp.Show();

            // update Q
            if (k == 1)
            {
                // Q = MList[k];
                Q = MTemp;
            }
            else
            {
                // Q = MatrixMultiply(Q, MList[k]);
                Q = MMMultiply(Q, MTemp);
            }

            // update R
            // puts("before updating R");
            // puts("MList[k]");
            // MTemp.Show();

            // puts("R");
            // R.Show();

            // R = MatrixMultiply(MList[k], R);
            R = MMMultiply(MTemp, R);

            // puts("update R is");
            // R.Show();
        }
    }

    template <typename T>
    LIAG_FUNC_MACRO bool SymmInvertMatrixInner(const Matrix<T, 3, 3> &m, Matrix<T, 3, 3> &inv_m)
    {
        inv_m[0][0] = m[1][1] * m[2][2] - m[1][2] * m[2][1];

        inv_m[1][0] = m[1][2] * m[2][0] - m[1][0] * m[2][2];

        inv_m[2][0] = m[1][0] * m[2][1] - m[1][1] * m[2][0];

        inv_m[0][1] = m[0][2] * m[2][1] - m[0][1] * m[2][2];

        inv_m[1][1] = m[0][0] * m[2][2] - m[0][2] * m[2][0];

        inv_m[2][1] = m[0][1] * m[2][0] - m[0][0] * m[2][1];

        inv_m[0][2] = m[0][1] * m[1][2] - m[0][2] * m[1][1];

        inv_m[1][2] = m[0][2] * m[1][0] - m[0][0] * m[1][2];

        inv_m[2][2] = m[0][0] * m[1][1] - m[0][1] * m[1][0];

        double det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);

        if (det == 0)
        {
            printf("singular matrix\n");
            // matrix_show(m);
            return false;
        }

        det = 1.0 / det;

        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                inv_m[i][j] = inv_m[i][j] * det;
            }
        }

        return true;
    }

    template <typename T>
    LIAG_FUNC_MACRO bool SymmInvertMatrixInner(const Matrix<T, 4, 4> &m, Matrix<T, 4, 4> &inv_m)
    {
        inv_m[0][0] = m[1][1] * m[2][2] * m[3][3] -
                      m[1][1] * m[2][3] * m[3][2] -
                      m[2][1] * m[1][2] * m[3][3] +
                      m[2][1] * m[1][3] * m[3][2] +
                      m[3][1] * m[1][2] * m[2][3] -
                      m[3][1] * m[1][3] * m[2][2];

        inv_m[1][0] = -m[1][0] * m[2][2] * m[3][3] +
                      m[1][0] * m[2][3] * m[3][2] +
                      m[2][0] * m[1][2] * m[3][3] -
                      m[2][0] * m[1][3] * m[3][2] -
                      m[3][0] * m[1][2] * m[2][3] +
                      m[3][0] * m[1][3] * m[2][2];

        inv_m[2][0] = m[1][0] * m[2][1] * m[3][3] -
                      m[1][0] * m[2][3] * m[3][1] -
                      m[2][0] * m[1][1] * m[3][3] +
                      m[2][0] * m[1][3] * m[3][1] +
                      m[3][0] * m[1][1] * m[2][3] -
                      m[3][0] * m[1][3] * m[2][1];

        inv_m[3][0] = -m[1][0] * m[2][1] * m[3][2] +
                      m[1][0] * m[2][2] * m[3][1] +
                      m[2][0] * m[1][1] * m[3][2] -
                      m[2][0] * m[1][2] * m[3][1] -
                      m[3][0] * m[1][1] * m[2][2] +
                      m[3][0] * m[1][2] * m[2][1];

        inv_m[0][1] = -m[0][1] * m[2][2] * m[3][3] +
                      m[0][1] * m[2][3] * m[3][2] +
                      m[2][1] * m[0][2] * m[3][3] -
                      m[2][1] * m[0][3] * m[3][2] -
                      m[3][1] * m[0][2] * m[2][3] +
                      m[3][1] * m[0][3] * m[2][2];

        inv_m[1][1] = m[0][0] * m[2][2] * m[3][3] -
                      m[0][0] * m[2][3] * m[3][2] -
                      m[2][0] * m[0][2] * m[3][3] +
                      m[2][0] * m[0][3] * m[3][2] +
                      m[3][0] * m[0][2] * m[2][3] -
                      m[3][0] * m[0][3] * m[2][2];

        inv_m[2][1] = -m[0][0] * m[2][1] * m[3][3] +
                      m[0][0] * m[2][3] * m[3][1] +
                      m[2][0] * m[0][1] * m[3][3] -
                      m[2][0] * m[0][3] * m[3][1] -
                      m[3][0] * m[0][1] * m[2][3] +
                      m[3][0] * m[0][3] * m[2][1];

        inv_m[3][1] = m[0][0] * m[2][1] * m[3][2] -
                      m[0][0] * m[2][2] * m[3][1] -
                      m[2][0] * m[0][1] * m[3][2] +
                      m[2][0] * m[0][2] * m[3][1] +
                      m[3][0] * m[0][1] * m[2][2] -
                      m[3][0] * m[0][2] * m[2][1];

        inv_m[0][2] = m[0][1] * m[1][2] * m[3][3] -
                      m[0][1] * m[1][3] * m[3][2] -
                      m[1][1] * m[0][2] * m[3][3] +
                      m[1][1] * m[0][3] * m[3][2] +
                      m[3][1] * m[0][2] * m[1][3] -
                      m[3][1] * m[0][3] * m[1][2];

        inv_m[1][2] = -m[0][0] * m[1][2] * m[3][3] +
                      m[0][0] * m[1][3] * m[3][2] +
                      m[1][0] * m[0][2] * m[3][3] -
                      m[1][0] * m[0][3] * m[3][2] -
                      m[3][0] * m[0][2] * m[1][3] +
                      m[3][0] * m[0][3] * m[1][2];

        inv_m[2][2] = m[0][0] * m[1][1] * m[3][3] -
                      m[0][0] * m[1][3] * m[3][1] -
                      m[1][0] * m[0][1] * m[3][3] +
                      m[1][0] * m[0][3] * m[3][1] +
                      m[3][0] * m[0][1] * m[1][3] -
                      m[3][0] * m[0][3] * m[1][1];

        inv_m[3][2] = -m[0][0] * m[1][1] * m[3][2] +
                      m[0][0] * m[1][2] * m[3][1] +
                      m[1][0] * m[0][1] * m[3][2] -
                      m[1][0] * m[0][2] * m[3][1] -
                      m[3][0] * m[0][1] * m[1][2] +
                      m[3][0] * m[0][2] * m[1][1];

        inv_m[0][3] = -m[0][1] * m[1][2] * m[2][3] +
                      m[0][1] * m[1][3] * m[2][2] +
                      m[1][1] * m[0][2] * m[2][3] -
                      m[1][1] * m[0][3] * m[2][2] -
                      m[2][1] * m[0][2] * m[1][3] +
                      m[2][1] * m[0][3] * m[1][2];

        inv_m[1][3] = m[0][0] * m[1][2] * m[2][3] -
                      m[0][0] * m[1][3] * m[2][2] -
                      m[1][0] * m[0][2] * m[2][3] +
                      m[1][0] * m[0][3] * m[2][2] +
                      m[2][0] * m[0][2] * m[1][3] -
                      m[2][0] * m[0][3] * m[1][2];

        inv_m[2][3] = -m[0][0] * m[1][1] * m[2][3] +
                      m[0][0] * m[1][3] * m[2][1] +
                      m[1][0] * m[0][1] * m[2][3] -
                      m[1][0] * m[0][3] * m[2][1] -
                      m[2][0] * m[0][1] * m[1][3] +
                      m[2][0] * m[0][3] * m[1][1];

        inv_m[3][3] = m[0][0] * m[1][1] * m[2][2] -
                      m[0][0] * m[1][2] * m[2][1] -
                      m[1][0] * m[0][1] * m[2][2] +
                      m[1][0] * m[0][2] * m[2][1] +
                      m[2][0] * m[0][1] * m[1][2] -
                      m[2][0] * m[0][2] * m[1][1];

        double det = m[0][0] * inv_m[0][0] +
                     m[0][1] * inv_m[1][0] +
                     m[0][2] * inv_m[2][0] +
                     m[0][3] * inv_m[3][0];

        if (det == 0)
        {
            printf("singular matrix for 4by4\n");
            // matrix_show(m);
            return false;
        }

        det = 1.0 / det;

        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                inv_m[i][j] = inv_m[i][j] * det;
            }
        }

        return true;
    }

    // refer to
    // https://www.cs.upc.edu/~jordicf/Teaching/programming/pdf4/MATH03_Gaussian-4slides.pdf
    template <typename T, uint Size>
    LIAG_FUNC_MACRO Vec<T, Size> SymmBackSubstitution(const Matrix<T, Size, Size> &A, const Vec<T, Size> &b)
    {
        Vec<double, Size> x(0.0);
        bool ifUpper = A.IsUpperTriangular(0.00001);
        // TODO, print sth if the A is not upper triangular
        assert(ifUpper == true);

        for (int i = Size - 1; i >= 0; --i)
        {
            // be careful about the last iteration here
            // it might hit the case where i =-1
            // in that case, there is error
            if (i < 0)
            {
                break;
            }
            double s = 0;
            // start with one element after i
            // end with the last element

            for (uint j = (i + 1); j < Size; ++j)
            {
                s = s + A[i][j] * x[j];
            }

            // at the last iteration, i is goes to -1
            x[i] = (b[i] - s) / A[i][i];
        }
        return x;
    }

    // for small matrix, we use the direact inverse to solve the linear equation
    // for large matrix, we need to use other method such as qr things to solve

    // TODO, checking singularity for inverting the matrix
    // refer to https://inst.eecs.berkeley.edu/~ee127/sp21/livebook/l_lineqs_solving.html
    // computing the A_inv by solving A*A_inv=I computing each colum of A_inv each time
    template <typename T, uint Size>
    LIAG_FUNC_MACRO bool SymmInvertMatrixInner(const Matrix<T, Size, Size> &m, Matrix<T, Size, Size> &inv_m)
    {
        // TODO return false for singular matrix

        // QR decomposition of matrix
        Matrix<T, Size, Size> Q;
        Matrix<T, Size, Size> R;
        Householder(m, Q, R);

        // puts("inv q r");
        // Q.Show();
        // R.Show();

        // Q*R*x=b
        // R*x = Q_t*b
        // b is each colum of the I matrix
        for (uint j = 0; j < Size; j++)
        {
            Vec<double, Size> b(0.0);
            b[j] = 1.0;

            // compute the transpose of Q
            Matrix<T, Size, Size> Q_t = Q;
            Q_t.Transpose();

            // puts("Q_t");
            // Q_t.Show();

            // compute Qt*b
            Vec<T, Size> Q_tb = MMVultiply(Q_t, b);
            // puts("Q_tb");
            // Q_t.Show();

            // using back substition to solve R*x=Q_t*b
            Vec<T, Size> x = SymmBackSubstitution(R, Q_tb);
            // printf("colm %d \n for x\n", j);
            // x.Show();

            // putting associating colum into the jth colum of inv_m matrix
            for (uint i = 0; i < Size; i++)
            {
                inv_m[i][j] = x[i];
            }
        }
        // inv_m.Show();
        return true;
    }

    template <typename T,
              uint Size>
    LIAG_FUNC_MACRO bool SymmInvertMatrix(const Matrix<T, Size, Size> &A, Matrix<T, Size, Size> &AInv)
    {

        // it seems there is issue
        // if we implement different function according to the value of Size
        SymmInvertMatrixInner(A, AInv);
        return true;
    }

    template <typename T, uint Size>
    LIAG_FUNC_MACRO void SymmEigenValuesShift(const Matrix<T, Size, Size> &A, double tol, int maxIter, Vec<T, Size> &eigenValues)
    {
        // refer to https://www.andreinc.net/2021/01/25/computing-eigenvalues-and-eigenvectors-using-qr-decomposition
        Matrix<T, Size, Size> Ak = A;

        // mat_t qq = matrix_new_eye(x->m, x->m);
        Matrix<T, Size, Size> QQ;
        QQ.InitEye();

        Matrix<T, Size, Size> Eye;
        Eye.InitEye();

        Matrix<T, Size, Size> R, Q;
        int actualIter=0;
        for (int i = 0; i < maxIter; ++i)
        {
            actualIter=i;
            // s_k is the last item of the first diagonal
            T s = Ak[Size - 1][Size - 1];
            Matrix<T, Size, Size> smulI = MSCALE(s, Eye);
            Householder(Ak - smulI, Q, R);
            // printf("iter i %d \nR\n", i);
            // R.Show();
            // puts("Q");
            // Q.Show();

            // addinng smulI back
            Ak = MMMultiply(R, Q) + smulI;
            QQ = MMMultiply(QQ, Q);
            
            //Ak.Show();
            if (Ak.IsUpperTriangular(tol))
            {
                break;
            }
        }

        std::cout << "actual iter is " << actualIter << std::endl;

        for (int i = 0; i < Size; i++)
        {
            if (fabs(Ak[i][i] - 0) < tol)
            {
                eigenValues[i] = 0.0;
            }
            else
            {
                eigenValues[i] = Ak[i][i];
            }
        }
    }

    template <typename T, uint Size>
    LIAG_FUNC_MACRO void SymmEigenValues(const Matrix<T, Size, Size> &A, double tol, int maxIter, Vec<T, Size> &eigenValues)
    {
        // refer to https://www.andreinc.net/2021/01/25/computing-eigenvalues-and-eigenvectors-using-qr-decomposition
        // mat_t ak = *x;
        Matrix<T, Size, Size> Ak = A;

        // mat_t qq = matrix_new_eye(x->m, x->m);
        Matrix<T, Size, Size> QQ;
        QQ.InitEye();

        Matrix<T, Size, Size> R, Q;
        for (int i = 0; i < maxIter; ++i)
        {

            Householder(Ak, Q, R);
            // printf("iter i %d \nR\n", i);
            // R.Show();
            // puts("Q");
            // Q.Show();

            Ak = MMMultiply(R, Q);
            QQ = MMMultiply(QQ, Q);

            if (Ak.IsUpperTriangular(tol))
            {
                break;
            }
        }

        for (int i = 0; i < Size; i++)
        {
            if (fabs(Ak[i][i] - 0) < tol)
            {
                eigenValues[i] = 0.0;
            }
            else
            {
                eigenValues[i] = Ak[i][i];
            }
        }
    }

    // using inverse itertaion to compute the eigen vectors
    // https://en.wikipedia.org/wiki/Inverse_iteration
    template <typename T, uint Size>
    LIAG_FUNC_MACRO Matrix<T, Size, Size> SymmEigenVectors(const Matrix<T, Size, Size> &A, const Vec<T, Size> &eigenValues, int maxIter)
    {

        Matrix<T, Size, Size> eigenVectors;
        Matrix<T, Size, Size> A_minus_lambda_i;
        Matrix<T, Size, Size> A_minus_lambda_inv;
        // go through each eigen values
        for (uint j = 0; j < Size; j++)
        {
            if (fabs(eigenValues[j] - 0.0) < 0.00001)
            {
                // the associated eigen vector is all zero
                continue;
            }

            double lambda = eigenValues[j] + 0.000001;
            for (uint ii = 0; ii < Size; ii++)
            {
                for (uint jj = 0; jj < Size; jj++)
                {
                    if (ii == jj)
                    {
                        A_minus_lambda_i[ii][jj] = A[ii][jj] - lambda;
                    }
                    else
                    {
                        A_minus_lambda_i[ii][jj] = A[ii][jj];
                    }
                }
            }

            // puts("a-lambda*I");:1
            // A_minus_lambda_i.Show();

            bool invertok = SymmInvertMatrix(A_minus_lambda_i, A_minus_lambda_inv);
            assert(invertok == true);

            // printf("eigen vec %f a-lambda*I inv\n", eigenValues[j]);
            // A_minus_lambda_inv.Show();

            Vec<T, Size> bPrev;
            Vec<T, Size> bCurr;
            // init as 1
            for (int i = 0; i < Size; i++)
            {
                bPrev[i] = 1.0;
            }
            uint iterNum = 0;
            while (iterNum < maxIter)
            {
                bCurr = MMVultiply(A_minus_lambda_inv, bPrev);
                // puts("check mmv");
                // bCurr.Show();

                bCurr.Normalize();
                // puts("check norm");
                // bCurr.Show();

                if (bCurr.Equal(bPrev))
                {
                    // results convert
                    break;
                }
                bPrev = bCurr;
                iterNum++;
            }

            // assign to the matrix
            for (uint i = 0; i < Size; i++)
            {
                eigenVectors[i][j] = bCurr[i];
            }
        }
        return eigenVectors;
    }

    // The input matrix such as the covariance matrix
    // is supposed to be a symmetric positive definite matrix with all
    // eigen values larger than 0
    template <typename T, uint Size>
    LIAG_FUNC_MACRO Matrix<T, Size, Size> SymmEigenDecomposition(const Matrix<T, Size, Size> &A, double tol, int maxIter, bool useShift = false)
    {
        // solve eigen values
        Vec<T, Size> eigenValues;
        eigenValues.InitZero();
        if (useShift)
        {
            // the values shift converge faster than the naive eigen values case
            // but it is error prone for some test cases
            SymmEigenValuesShift(A, tol, maxIter, eigenValues);
        }
        else
        {
            SymmEigenValues(A, tol, maxIter, eigenValues);
        }

        //std::cout << "eigen values" << std::endl;
        //eigenValues.Show();

        // solve eigen vectors
        Matrix<T, Size, Size> eigenVactors;
        eigenVactors = SymmEigenVectors(A, eigenValues, maxIter);

        // create diagonal matrix
        Matrix<T, Size, Size> diag;
        for (uint i = 0; i < Size; i++)
        {

            if (eigenValues[i] < 0)
            {
                if (fabs(eigenValues[i]) < 0.0002)
                {
                    // make sure all value is >0 and we can compute sqrt for it
                    // there are some numerical errors for computing the eigen values
                    eigenValues[i] = -eigenValues[i];
                }
                else
                {
                    // debug use
                    printf("eigen values are\n");
                    for (int i = 0; i < Size; i++)
                    {
                        printf(" %f ", eigenValues[i]);
                    }
                    printf("the eigen value is supposed to be >=0\n");
                    assert(false);
                }
            }
            diag[i][i] = sqrt(eigenValues[i]);
        }
        Matrix<T, Size, Size> decompA;
        decompA = MMMultiply(eigenVactors, diag);
        return decompA;
    }
}
#endif
