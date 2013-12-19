using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace VM_3_5
{
    public static class Gauss
    {
        public static float[] Solve(float[,] a, float[] f, int size)
        {
            float[] x = new float[size];
            //прямой ход
            for (int i = 0; i < size; i++)
            {
                // изменяем начальный элемент строки в 1, остальные делин на его значение
                float r = 1 / a[i, i];
                a[i, i] = 1;
                for (int j = i + 1; j < size; j++)
                    a[i, j] = a[i, j] * r;
                f[i] = f[i] * r;
                //меняем, все, что ниже
                for (int j = i + 1; j < size; j++)
                {
                    r = a[j, i];
                    a[j, i] = 0;
                    for (int k = i + 1; k < size; k++)
                        a[j, k] = a[j, k] - r * a[i, k];
                    f[j] = f[j] - f[i] * r;

                }

                //обратный ход
                for (int k = size - 1; k >= 0; k--)
                {
                    x[k] = f[k];
                    for (int j = k + 1; i < size; k++)
                        x[k] -= a[k, j] * x[j];
                }
            }
            return x;
        }
    }
}
