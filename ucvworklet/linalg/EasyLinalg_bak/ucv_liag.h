#ifndef UCV_LIAG
#define UCV_LIAG

#include <iostream>
#include <assert.h>
#include <cmath>
#include <type_traits>

namespace UCVLIAG
{
    // this can be customized in order to run on different devices
    #define LIAG_FUNC_MACRO VTKM_EXEC_CONT

    template <typename T, uint Size>
    class Vec
    {
    public:
        static const int NUM_COMPONENTS = Size;
        T Components[NUM_COMPONENTS];

        LIAG_FUNC_MACRO Vec(){};

        // getting element by [] operator
        // return a reference, so the value can be updated further
        LIAG_FUNC_MACRO T &operator[](uint index)
        {
            assert(index >= 0);
            assert(index < this->NUM_COMPONENTS);
            return this->Components[index];
        }

        // if not adding this, there are some error such as
        /*
        /Users/zw1/Documents/cworkspace/src/EasyLinalg/StaticMemTemplate/mathTemplate.h:141:42: error: no viable overloaded operator[] for type 'const Vec<Vec<double, NUM_COLUMNS>, NUM_ROWS>'
                    this->Components[i][j] = t.Components[i][j];
                                             ^~~~~~~~~~~~ ~
        /Users/zw1/Documents/cworkspace/src/EasyLinalg/StaticMemTemplate/mathTemplate.h:343:15: note: in instantiation of member function 'Matrix<double, 3, 3>::operator=' requested here
                Q = MatrixMultiply(Q, MList[k]);
                  ^
        /Users/zw1/Documents/cworkspace/src/EasyLinalg/StaticMemTemplate/test.cpp:163:5: note: in instantiation of function template specialization 'Householder<double, 3U, 3U>' requested here
        Householder(x,Q,R);
        ^
        /Users/zw1/Documents/cworkspace/src/EasyLinalg/StaticMemTemplate/test.cpp:187:5: note: in instantiation of function template specialization 'testQRInner<3U>' requested here
        testQRInner<3>();
        ^
        /Users/zw1/Documents/cworkspace/src/EasyLinalg/StaticMemTemplate/mathTemplate.h:18:8: note: candidate function not viable: 'this' argument has type 'const Vec<Vec<double, NUM_COLUMNS>, NUM_ROWS>', but method is not marked const
        T &operator[](const int index)

        The failure happens because Non-const functions can only be called by non-const objects.
        However, when a function is declared as const, it can be called on any type of object.

        still curious about results here, some cases, the object that call the [] is a const object
        */
        LIAG_FUNC_MACRO const T &operator[](uint index) const
        {
            assert(index >= 0);
            assert(index < this->NUM_COMPONENTS);
            return this->Components[index];
        }

        LIAG_FUNC_MACRO Vec &operator=(const Vec &t)
        {
            assert(this->NUM_COMPONENTS == t.NUM_COMPONENTS);
            for (uint i = 0; i < this->NUM_COMPONENTS; i++)
            {
                this->Components[i] = t.Components[i];
            }
            return *this;
        }

        LIAG_FUNC_MACRO Vec &operator=(Vec &t)
        {
            assert(this->NUM_COMPONENTS == t.NUM_COMPONENTS);
            for (uint i = 0; i < this->NUM_COMPONENTS; i++)
            {
                this->Components[i] = t.Components[i];
            }
            return *this;
        }

        // this is only used for single vector
        // the vector can be in inserted form
        // so we do not set the init value of vector as zero
        LIAG_FUNC_MACRO void InitZero()
        {
            for (int i = 0; i < this->NUM_COMPONENTS; i++)
            {
                this->Components[i] = 0;
            }
        }

        LIAG_FUNC_MACRO void Show()
        {
            for (int i = 0; i < this->NUM_COMPONENTS; i++)
            {
                printf("%lf ", 1.0 * this->Components[i]);
            }
            printf("\n");
        }

        // ||x||
        LIAG_FUNC_MACRO double GetNorm()
        {
            double sum = 0;
            for (uint i = 0; i < this->NUM_COMPONENTS; i++)
                sum = sum + this->Components[i] * this->Components[i];
            return sqrt(sum);
        }

        // y = y/d
        LIAG_FUNC_MACRO void Div(double d)
        {
            // if the d is a value close to zero
            // just set results to 0
            assert((d - 0.0) > 0.00001);
            for (uint i = 0; i < Size; i++)
            {
                this->Components[i] = static_cast<T>(1.0 * this->Components[i] / d);
            }
            return;
        }

        LIAG_FUNC_MACRO bool Equal(const Vec &t)
        {
            if (this->NUM_COMPONENTS != t.NUM_COMPONENTS)
            {
                return false;
            }

            for (uint i = 0; i < this->NUM_COMPONENTS; i++)
            {
                if (fabs(this->Components[i] - t.Components[i]) > 0.0001)
                {
                    return false;
                }
            }
            return true;
        }

        LIAG_FUNC_MACRO void Normalize()
        {
            double sum = 0;
            for (uint i = 0; i < this->NUM_COMPONENTS; i++)
            {
                sum += this->Components[i] * this->Components[i];
            }
            assert(fabs(sum - 0.0) > 0.000001);
            for (uint i = 0; i < this->NUM_COMPONENTS; i++)
            {
                this->Components[i] = this->Components[i] / sqrt(sum);
            }
        }
    };

    template <typename T, uint RowSize, uint ColSize>
    class Matrix
    {
    public:
        static constexpr int NUM_ROWS = RowSize;
        static constexpr int NUM_COLUMNS = ColSize;
        Vec<Vec<T, NUM_COLUMNS>, NUM_ROWS> Components;

        LIAG_FUNC_MACRO Matrix()
        {
            for (uint i = 0; i < this->NUM_ROWS; i++)
            {
                for (uint j = 0; j < this->NUM_COLUMNS; j++)
                {
                    this->Components[i][j] = 0;
                }
            }
        };

        LIAG_FUNC_MACRO bool IsUpperTriangular(double tol)
        {
            // For now, only treat square matricies.
            if (this->NUM_COLUMNS != this->NUM_ROWS)
            {
                return false;
            }

            for (uint i = 0; i < this->NUM_COLUMNS; i++)
            {
                for (uint j = 0; j < i; j++)
                {
                    if (fabs(this->Components[i][j]) > tol)
                    {
                        return false;
                    }
                }
            }
            return true;
        }

        LIAG_FUNC_MACRO bool IsEqual(const Matrix &t)
        {
            if (this->NUM_ROWS != t.NUM_ROWS || this->NUM_COLUMNS != t.NUM_COLUMNS)
            {
                return false;
            }
            for (uint i = 0; i < this->NUM_ROWS; i++)
            {
                for (uint j = 0; j < this->NUM_COLUMNS; j++)
                {
                    if (fabs(this->Components[i][j] - t.Components[i][j]) > 0.0001)
                    {
                        return false;
                    }
                }
            }
            return true;
        }

        LIAG_FUNC_MACRO Vec<T, NUM_COLUMNS> &operator[](uint rowIndex)
        {
            assert(rowIndex >= 0);
            assert(rowIndex < NUM_ROWS);
            return this->Components[rowIndex];
        }

        LIAG_FUNC_MACRO const Vec<T, NUM_COLUMNS> &operator[](uint rowIndex) const
        {
            assert(rowIndex >= 0);
            assert(rowIndex < NUM_ROWS);
            return this->Components[rowIndex];
        }

        LIAG_FUNC_MACRO Matrix &operator=(Matrix &t)
        {
            assert(this->NUM_ROWS == t.NUM_ROWS);
            assert(this->NUM_COLUMNS == t.NUM_COLUMNS);

            for (uint i = 0; i < this->NUM_ROWS; i++)
            {
                for (uint j = 0; j < this->NUM_COLUMNS; j++)
                {
                    this->Components[i][j] = t.Components[i][j];
                }
            }
            return *this;
        }

        LIAG_FUNC_MACRO Matrix &operator=(const Matrix &t)
        {
            assert(this->NUM_ROWS == t.NUM_ROWS);
            assert(this->NUM_COLUMNS == t.NUM_COLUMNS);

            for (uint i = 0; i < this->NUM_ROWS; i++)
            {
                for (uint j = 0; j < this->NUM_COLUMNS; j++)
                {
                    this->Components[i][j] = t.Components[i][j];
                }
            }
            return *this;
        }

        LIAG_FUNC_MACRO void InitEye()
        {
            assert(this->NUM_ROWS == this->NUM_COLUMNS);
            for (int i = 0; i < this->NUM_ROWS; i++)
            {
                for (int j = 0; j < this->NUM_COLUMNS; j++)
                {
                    if (i == j)
                    {
                        this->Components[i][j] = 1;
                    }
                    else
                    {
                        this->Components[i][j] = 0;
                    }
                }
            }
        }

        LIAG_FUNC_MACRO void Show()
        {
            for (int i = 0; i < this->NUM_ROWS; i++)
            {
                for (int j = 0; j < this->NUM_COLUMNS; j++)
                {
                    printf(" %8.7f", 1.0 * this->Components[i][j]);
                }
                printf("\n");
            }
            printf("\n");
        }

        LIAG_FUNC_MACRO void Transpose()
        {
            for (int i = 0; i < this->NUM_ROWS; i++)
            {
                for (int j = 0; j < i; j++)
                {
                    T t = this->Components[i][j];
                    this->Components[i][j] = this->Components[j][i];
                    this->Components[j][i] = t;
                }
            }
        }

        // only works for the symmetric matrix
        LIAG_FUNC_MACRO Matrix GetMinor(uint d)
        {
            assert(this->NUM_ROWS == this->NUM_COLUMNS);
            Matrix minor;
            for (uint i = 0; i < d; i++)
            {
                minor.Components[i][i] = 1;
            }

            // start form the dth position
            // make it same with the original data

            for (uint i = d; i < NUM_ROWS; i++)
            {
                for (uint j = d; j < NUM_COLUMNS; j++)
                {
                    minor.Components[i][j] = this->Components[i][j];
                }
            }

            return minor;
        }

        // take c-th column of matrix, put in vector
        LIAG_FUNC_MACRO Vec<T, NUM_COLUMNS> GetColum(uint c)
        {
            Vec<T, NUM_COLUMNS> vec;
            for (uint i = 0; i < this->NUM_ROWS; i++)
                vec[i] = this->Components[i][c];
            return vec;
        }
    };

    // m = I - 2* v v^T
    template <typename T, uint Size>
    LIAG_FUNC_MACRO Matrix<T, Size, Size> IMINUSVVT(Vec<T, Size> &v)
    {
        // this member become the static const member?
        Matrix<T, Size, Size> Mat;
        for (uint i = 0; i < Size; i++)
            for (uint j = 0; j < Size; j++)
                Mat[i][j] = -2 * v[i] * v[j];
        // diagnal elements add 1
        for (uint i = 0; i < Size; i++)
            Mat[i][i] = Mat[i][i] + 1;

        return Mat;
    }

    template <typename T,
              uint NumRow,
              uint NumCol,
              uint NumInternal>
    LIAG_FUNC_MACRO Matrix<T, NumRow, NumCol> MMMultiply(
        Matrix<T, NumRow, NumInternal> &leftFactor,
        Matrix<T, NumInternal, NumCol> &rightFactor)
    {
        Matrix<T, NumRow, NumCol> result;

        for (uint i = 0; i < NumRow; i++)
        {
            for (uint j = 0; j < NumCol; j++)
            {
                // init specific position of z as zero
                result[i][j] = 0;
                for (uint k = 0; k < NumInternal; k++)
                {
                    result[i][j] = result[i][j] + leftFactor[i][k] * rightFactor[k][j];
                }
            }
        }
        return result;
    }

    template <typename T,
              uint NumRow,
              uint NumCol>
    LIAG_FUNC_MACRO Vec<T, NumRow> MMVultiply(const Matrix<T, NumRow, NumCol> &inputM, const Vec<T, NumCol> inputV)
    {
        Vec<T, NumRow> result;
        result.InitZero();
        for (int i = 0; i < NumRow; i++)
        {
            for (int j = 0; j < NumCol; j++)
            {
                // for each row
                result[i] += inputM[i][j] * inputV[j];
            }
        }
        return result;
    }

    // matrix multiply vector plus vector
    // A*U+M
    template <typename T,
              uint NumRow,
              uint NumCol>
    LIAG_FUNC_MACRO Vec<T, NumRow> MMVPV(
        Matrix<T, NumRow, NumCol> &A,
        Vec<T, NumCol> &U,
        Vec<T, NumRow> &M)
    {
        Vec<T, NumRow> AUM;
        AUM.InitZero();
        // Be careful about it, init the vector to zero
        for (int i = 0; i < NumRow; i++)
        {
            // adding vector into the matrix
            for (int j = 0; j < NumCol; j++)
            {
                AUM[i] = AUM[i] + (A[i][j] * U[j]);
            }
            AUM[i] = AUM[i] + M[i];
        }

        return AUM;
    }

    // results = a*x + b
    template <typename T, uint Size>
    LIAG_FUNC_MACRO Vec<T, Size> AXPY(double a, Vec<T, Size> &x, Vec<T, Size> &b)
    {
        Vec<T, Size> result;
        for (uint i = 0; i < Size; i++)
        {
            result[i] = a * x[i] + b[i];
        }
        return result;
    }

    // A is matrix, x and y are vectors
    //  results =alpha*A*x+beta*y
    template <typename T, uint Row, uint Col>
    LIAG_FUNC_MACRO Vec<T, Row> DGEMV(const T alpha, const Matrix<T, Row, Col> &A, Vec<T, Col> &x, const T beta, const Vec<T, Row> &y)
    {
        Vec<T, Row> result;
        result.InitZero();
        for (uint i = 0; i < Row; i++)
        {
            for (uint j = 0; j < Col; j++)
            {
                result[i] += A[i][j] * x[j];
            }
            result[i] = result[i] * alpha;
            result[i] = result[i] + beta * y[i];
        }
        return result;
    }

    template <typename T,
              uint NumRow,
              uint NumCol>
    LIAG_FUNC_MACRO void Householder(
        Matrix<T, NumRow, NumCol> &A,
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

    template <typename T, uint Size>
    LIAG_FUNC_MACRO bool SymmInvertMatrixInnerGeneral(const Matrix<T, Size, Size> &m, Matrix<T, Size, Size> &inv_m)
    {
        printf("TODO");
        return true;
    }

    template <typename T,
              uint Size>
    LIAG_FUNC_MACRO bool SymmInvertMatrix(const Matrix<T, Size, Size> &A, Matrix<T, Size, Size> &AInv)
    {

        if (Size == 3)
        {
            SymmInvertMatrixInner(A, AInv);
            return true;
        }
        else if (Size == 4)
        {
            SymmInvertMatrixInner(A, AInv);
            return true;
        }
        else
        {
            // TODO adding size = 2
            SymmInvertMatrixInnerGeneral(A, AInv);
        }
        return false;
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

            // puts("a-lambda*I");
            // A_minus_lambda_i.Show();

            bool invertok = SymmInvertMatrix(A_minus_lambda_i, A_minus_lambda_inv);
            assert(invertok == true);

            // puts("a-lambda*I inv");
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

    template <typename T, uint Size>
    LIAG_FUNC_MACRO Matrix<T, Size, Size> SymmEigenDecomposition(const Matrix<T, Size, Size> &A, double tol, int maxIter)
    {

        // solve eigen values
        Vec<T, Size> eigenValues;
        eigenValues.InitZero();
        SymmEigenValues(A, tol, 20, eigenValues);

        // eigenValues.Show();

        // solve eigen vectors
        Matrix<T, Size, Size> eigenVactors;
        eigenVactors = SymmEigenVectors(A, eigenValues, maxIter);

        // create diagonal matrix
        Matrix<T, Size, Size> diag;
        for (uint i = 0; i < Size; i++)
        {
            diag[i][i] = sqrt(eigenValues[i]);
        }

        Matrix<T, Size, Size> decompA;
        decompA = MMMultiply(eigenVactors, diag);
        return decompA;
    }

}

#endif