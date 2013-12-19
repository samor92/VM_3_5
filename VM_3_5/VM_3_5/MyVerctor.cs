using System;
using System.Collections.Generic;

namespace VM_3_5
{
    class MyVerctor
    {
        public MyVerctor()
        {
        }

        // Генерация рандомного вектора
        public float[] GenerateRandomVector(int vectorSize, float limit1, float limit2)
        {
            const int max = 1000;
            float[] result = new float[vectorSize];
            Random rnd = new Random();
            for (int i = 0; i < vectorSize - 1; i++)
                result[i] = limit1 + rnd.Next(0, max)*(limit2 - limit1)/max;
            return result;
        }

        // Скалярное произведение векторов A B размерности vectorSize
        public float ScalarProduct(int vectorSize, float[] vectorA, float[] vectorB)
        {
            float result = 0;
            for (int i = 0; i < vectorSize - 1; i++)
                result = result + vectorA[i] * vectorB[i];
            return result;
        }

        // Нормализация вектора
        public float[] NormalizeVector(int vectorSize, float[] vectorA)
        {
            float[] result = new float[vectorSize];
            float tmp = (float)Math.Sqrt(ScalarProduct(vectorSize, vectorA, vectorA));
            for (int i = 0; i < vectorSize - 1; i++)
                result[i] = vectorA[i] / tmp;
            return result;
        }

        // Произведение двух квадратных матриц порядка n
        public float[,] MatrixMultiplication(int matrixOrder, float[,] matrixA, float[,] matrixB)
        {
            float[,] result = null;
            for (int i = 0; i < matrixOrder - 1; i++)
            {
                for (int j = 0; j < matrixOrder - 1; j++)
                {
                    result[i, j] = 0;
                    for (int k = 0; k < matrixOrder - 1; k++)
                        result[i, j] = result[i, j] + matrixA[i, k]*matrixB[k, j];
                }
            }
            return result;
        }

        public float[,] CreateSimmetricMatrix(int matrixOrder, int diap,  ref float lMin, ref float[] gMin)
        {
            int r1 = -10;
            int r2 = 10;
            float[,] H = new float[matrixOrder, matrixOrder];
            float[,] tmpA = new float[matrixOrder, matrixOrder];
            float[] omega = new float[matrixOrder];
            float[] lambda = new float[matrixOrder];
            
            omega = GenerateRandomVector(matrixOrder, r1, r2);
            omega = NormalizeVector(matrixOrder,omega);
            lambda = GenerateRandomVector(matrixOrder, -diap, diap);
            for (int i = 0; i < matrixOrder - 1; i++)
                lambda[i] = (float)(1 + 0.05*i);

            int item = 0;
            for (int i=0; i < matrixOrder-1; i++)
                if (Math.Abs(lambda[i]) < Math.Abs(lambda[item]))
                    item = i;
            lMin = lambda[item];

            // Строим матрицу Хаусхолдера
            for (int i = 0; i < matrixOrder - 1; i++)
            {
                for (int j = 0; j < matrixOrder - 1; j++)
                {
                    H[i, j] = -2*omega[i]*omega[j];
                    if (i == j)
                        H[i, j] = 1 + H[i, j];
                }
            }

            // Быстрое умножение H на диагональную матрицу лямбда
            for(int j=0; j < matrixOrder; j++)
                for (int i = 0; i < matrixOrder - 1; i++)
                    tmpA[i, j] = H[i, j]*lambda[j];

            return (MatrixMultiplication(matrixOrder, tmpA, H));
        }

        public void LUDecomposition(int matrixOrder, float[,] A, ref float[,] L, ref float[,] U)
        {
            for (int k = 0; k < matrixOrder - 1; k++)
            {
                for (int i = k; i < matrixOrder - 1; i++)
                {
                    float sum = 0;
                    for (int s = 0; s < k - 1; s++)
                        sum = sum + L[i, s]*U[s, k];
                    L[i, k] = A[i, k] - sum;
                }

                for (int i = k + 1; i < matrixOrder - 1; i++)
                {
                    float sum = 0;
                    for (int s = 0; s < k - 1; s++)
                        sum = sum + L[k, s]*U[s, i];
                    U[k, i] = (A[k, i] - sum)/L[k, k];
                }

                U[k, k] = 1;
                //Обнуление
                for (int i = k + 1; i < matrixOrder - 1; i++)
                {
                    L[k, i] = 0;
                    U[i, k] = 0;
                }
            }
        }

        public float[] SolveSlae(int matrixOrder, float[,] L, float[,] U, float[] b)
        {
            float[] result = new float[matrixOrder];
            float[] y = new float[matrixOrder];

            //Прямой ход
            for (int i = 0; i < matrixOrder - 1; i++)
            {
                float sum = 0;
                for (int j = 0; j < matrixOrder - 1; j++)
                    sum = sum + L[i, j]*y[i];
                y[i] = (b[i] - sum)/L[i, i];
            }

            //Обратный ход
            for (int i = matrixOrder - 1; i > 0; i--)
            {
                float sum = 0;
                for (int j = i + 1; i < matrixOrder - 1; j++)
                    sum = sum + U[i, j]*result[j];
                result[i] = y[i] - sum;
            }

            return result;
        }

        //Вычисляет nu(t)*A*nu , t - транспонирование
        public float CalcSigma(int matrixOrder, float[] nu, float[,] A)
        {
            float result = 0;
            for (int j = 0; j < matrixOrder - 1; j++)
            {
                float sum = 0;
                for (int i = 0; i < matrixOrder - 1; i++)
                    sum = sum + nu[i]*A[i, j];
                result = result + sum*nu[j];
            }
            return result;
        }

        // Вычисление углома между векторами
        public float AngleCos(int vectorSize, float[] vectorA, float[] vectorB)
        {
            return ScalarProduct(vectorSize, vectorA, vectorB)/
                   ((float)Math.Sqrt(ScalarProduct(vectorSize, vectorA, vectorA))*
                    (float) Math.Sqrt(ScalarProduct(vectorSize, vectorB, vectorB)));
        }

        public Dictionary<float, float[]> FindMinEigenvalue(int matrixOrder, float[,] A, float eps, int maxIter, ref int i)
        {
            int r1 = -10;
            int r2 = 10;

            float[] x = new float[matrixOrder];
            float[] nu = new float[matrixOrder];
            float[] nu_pred = new float[matrixOrder];

            float[,] L = null;
            float[,] U = null;

            LUDecomposition(matrixOrder, A, ref L, ref U);

            x = GenerateRandomVector(matrixOrder, r1, r2);
            float sigma = float.PositiveInfinity;
            float sigm_pred = 0;
            int k = 0;
            do
            {
                nu_pred = nu;
                nu = x;
                nu = NormalizeVector(matrixOrder, nu);
                x = SolveSlae(matrixOrder, L, U, nu);
                sigm_pred = sigma;
                sigma = CalcSigma(matrixOrder, nu, A);
                k++;

            } while ((Math.Abs(sigma-sigm_pred)<eps && Math.Abs(i-AngleCos(matrixOrder,nu_pred,nu)) < eps) || (k == maxIter ));

            Dictionary<float, float[]> result = new Dictionary<float, float[]>();
            result.Add(sigma,nu);

            return result;
        }

    }
}
