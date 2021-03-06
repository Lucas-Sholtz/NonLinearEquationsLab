using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics;
using MathNet.Numerics.Differentiation;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;

namespace NonLinearEquationsLab
{
    static class NLE
    {
        public const int FSIZE = 2;
        public const double eps = 0.01;
        public static Func<Double[], double>[] GetTwoFunctions()
        {
            //Func<Double[], double> f1 = p => p[0] - 0.5 * Math.Sin((p[0] - p[1]) / 2);
            //Func<Double[], double> f2 = p => p[1] - 0.5 * Math.Cos((p[0] + p[1]) / 2);
            Func<Double[], double> f1 = p => Math.Pow(p[0], 2) / Math.Pow(p[1], 2) - Math.Cos(p[1]) - 2;
            Func<Double[], double> f2 = p => Math.Pow(p[0], 2) + Math.Pow(p[1], 2) - 6;
            return new Func<Double[], double>[FSIZE] { f1, f2 };
        }
        public static Func<Double[], double>[] GetNFunctions(int n)
        {
            Func<Double[], double>[] funcs = new Func<Double[], double>[n];
            for (int i = 0; i < funcs.Length; i++)
            {
                Func<Double[], double> f = p =>
                {
                    double value = 0;
                    for (int j = 0; j < p.Length; j++)
                    {
                        if (j != i)
                        {
                            value += Math.Pow(p[j], 2);
                        }
                        else
                        {
                            value += Math.Pow(p[j], 3);
                        }
                    }
                    for (int j = 1; j <= p.Length; j++)
                    {
                        if (j != i + 1)
                        {
                            value -= Math.Pow(j, 2);
                        }
                        else
                        {
                            value -= Math.Pow(j, 3);
                        }
                    }
                    return value;
                };
                funcs[i] = f;
            }
            return funcs;
        }
        public static double SquareCubicFunction(double[] p, int step)
        {
            double value = 0;
            for (int j = 0; j < p.Length; j++)
            {
                if (j != step)
                {
                    value += Math.Pow(p[j], 2);
                }
                else
                {
                    value += Math.Pow(p[j], 3);
                }
            }
            for (int j = 1; j <= p.Length; j++)
            {
                if (j != step + 1)
                {
                    value -= Math.Pow(j, 2);
                }
                else
                {
                    value -= Math.Pow(j, 3);
                }
            }
            return value;
        }
        private static double GetNFunctionValue(int n, double[] p)
        {
            double value = 0;
            for (int i = 0; i < p.Length; i++)
            {
                if (i == n)
                {
                    value += Math.Pow(p[i], 3);
                }
                else
                {
                    value += Math.Pow(p[i], 2);
                }
            }
            for (int i = 0; i < p.Length; i++)
            {
                if (i == n)
                {
                    value -= Math.Pow((i + 1), 3);
                }
                else
                {
                    value -= Math.Pow((i + 1), 2);
                }
            }
            return value;
        }
        private static double GetNDerivativeValue(int n, double[] p, int index)
        {
            // x1^3 x2^2 ... xn^2 = 1 + 4 + ... + n^2
            // x1^2 x2^3 ... xn^2 = 1 + 8 + ... + n^2
            // ......................................
            // x1^2 x2^2 ... xn^3 = 1 + 4 + ... + n^3

            // and their derivatives
            if (n == index)
            {
                return 3 * Math.Pow(p[n], 2);
            }
            else
            {
                return 2 * p[n];
            }
        }
        public static Func<Double[], double>[][] Derivatives(Func<Double[], double>[] f)
        {
            var jacobian = new Func<Double[], double>[f.Length][];

            for (int i = 0; i < f.Length; i++)
            {
                jacobian[i] = new Func<double[], double>[f.Length];
                for (int j = 0; j < f.Length; j++)
                {
                    jacobian[i][j] = Differentiate.FirstPartialDerivativeFunc(f[i], j);
                }
            }

            return jacobian;
        }
        public static void RunModifiedNewtonMethod(Func<Double[], double>[] f)
        {
            var d = Derivatives(f);

            Vector<double> xn = Vector<double>.Build.Dense(f.Length, 1);
            Vector<double> xnp1 = Vector<double>.Build.Dense(f.Length, 1);
            Matrix<double> A = Matrix<double>.Build.Dense(f.Length, f.Length, 0);
            Vector<double> F = Vector<double>.Build.Dense(f.Length, 0);

            do
            {
                xn = xnp1;
                for (int i = 0; i < A.RowCount; i++)
                {
                    for (int j = 0; j < A.ColumnCount; j++)
                    {
                        A[i, j] = d[i][j].Invoke(xn.ToArray());
                    }
                }
                Console.WriteLine("A:");
                Console.WriteLine(A.ToMatrixString());
                var invertA = A.Inverse();
                Console.WriteLine("Inverse A:");
                Console.WriteLine(invertA.ToMatrixString());
                for (int i = 0; i < F.Count; i++)
                {
                    F[i] = f[i].Invoke(xn.ToArray());
                }
                Console.WriteLine("F:");
                Console.WriteLine(F.ToVectorString());
                xnp1 = xn - invertA * F;
                Console.WriteLine("Xn+1:");
                Console.WriteLine(xnp1.ToVectorString());

            } while ((xnp1 - xn).L2Norm() > eps);

            Console.WriteLine($"Answer:\n{xnp1.ToVectorString()}");

            var delta = Vector<double>.Build.Dense(f.Length, 0);
            for (int i = 0; i < delta.Count; i++)
            {
                delta[i] = f[i].Invoke(xnp1.ToArray());
            }

            Console.WriteLine($"Delta vector:\n{xnp1.ToVectorString()}");
            Console.WriteLine($"Delta norm:{delta.L2Norm()}");
        }
        public static void RunModifiedNewtonMethod(int n)
        {
            Vector<double> xn = Vector<double>.Build.Dense(n, 0);
            Vector<double> xnp1 = Vector<double>.Build.Random(n);
            Matrix<double> A = Matrix<double>.Build.Dense(n, n, 0);
            Vector<double> F = Vector<double>.Build.Dense(n, 0);

            do
            {
                xn = xnp1;
                for (int i = 0; i < A.RowCount; i++)
                {
                    for (int j = 0; j < A.ColumnCount; j++)
                    {
                        A[i, j] = GetNDerivativeValue(j, xn.ToArray(), i);
                    }
                }
                Console.WriteLine("A:");
                Console.WriteLine(A.ToMatrixString());
                var invertA = A.Inverse();
                Console.WriteLine("Inverse A:");
                Console.WriteLine(invertA.ToMatrixString());
                for (int i = 0; i < F.Count; i++)
                {
                    F[i] = GetNFunctionValue(i, xn.ToArray());
                }
                Console.WriteLine("F:");
                Console.WriteLine(F.ToVectorString());
                xnp1 = xn - invertA * F;
                Console.WriteLine("Xn+1:");
                Console.WriteLine(xnp1.ToVectorString());

            } while ((xnp1 - xn).L2Norm() > eps);

            Console.WriteLine($"Answer:\n{xnp1.ToVectorString()}");

            var delta = Vector<double>.Build.Dense(n, 0);
            for (int i = 0; i < delta.Count; i++)
            {
                delta[i] = GetNFunctionValue(i, xnp1.ToArray());
            }

            Console.WriteLine($"Delta vector:\n{delta.ToVectorString()}");
            Console.WriteLine($"Delta norm: {delta.L2Norm()}");
        }
        public static void RunNewtonMethod(Func<Double[], double>[] f)
        {
            var d = Derivatives(f);

            Vector<double> x = Vector<double>.Build.Dense(f.Length, 1);
            Matrix<double> A = Matrix<double>.Build.Dense(f.Length, f.Length, 0);
            Vector<double> F = Vector<double>.Build.Dense(f.Length, 0);
            Vector<double> z = Vector<double>.Build.Dense(f.Length, 0);

            do
            {
                for (int i = 0; i < A.RowCount; i++)
                {
                    for (int j = 0; j < A.ColumnCount; j++)
                    {
                        A[i, j] = d[i][j].Invoke(x.ToArray());
                    }
                }
                Console.WriteLine("A:");
                Console.WriteLine(A.ToMatrixString());
                for (int i = 0; i < F.Count; i++)
                {
                    F[i] = f[i].Invoke(x.ToArray());
                }
                Console.WriteLine("F:");
                Console.WriteLine(F.ToVectorString());
                z = A.Solve(F);
                Console.WriteLine("solve A*z=F \nz:");
                Console.WriteLine(z.ToVectorString());
                Console.WriteLine("z norm:");
                Console.WriteLine(z.L2Norm());
                x = x - z;
                Console.WriteLine("new x:");
                Console.WriteLine(x.ToVectorString());

            } while (z.L2Norm() > eps);

            Console.WriteLine($"Answer:\n {x.ToVectorString()}");

            var delta = Vector<double>.Build.Dense(f.Length, 0);
            for (int i = 0; i < F.Count; i++)
            {
                F[i] = f[i].Invoke(z.ToArray());
            }

            Console.WriteLine($"Delta vector:\n {z.ToVectorString()}");
            Console.WriteLine($"Delta norm: {z.L2Norm()}");
        }
        public static void RunNewtonMethod(int n)
        {
            Vector<double> x = Vector<double>.Build.Random(n);
            Matrix<double> A = Matrix<double>.Build.Dense(n, n, 0);
            Vector<double> F = Vector<double>.Build.Dense(n, 0);
            Vector<double> z = Vector<double>.Build.Dense(n, 0);

            do
            {
                for (int i = 0; i < A.RowCount; i++)
                {
                    for (int j = 0; j < A.ColumnCount; j++)
                    {
                        A[i, j] = GetNDerivativeValue(j, x.ToArray(), i);
                    }
                }
                Console.WriteLine("A:");
                Console.WriteLine(A.ToMatrixString());
                for (int i = 0; i < F.Count; i++)
                {
                    F[i] = GetNFunctionValue(i, x.ToArray());
                }
                Console.WriteLine("F:");
                Console.WriteLine(F.ToVectorString());
                z = A.Solve(F);
                Console.WriteLine("solve A*z=F \nz:");
                Console.WriteLine(z.ToVectorString());
                Console.WriteLine("z norm:");
                Console.WriteLine(z.L2Norm());
                x = x - z;
                Console.WriteLine("new x:");
                Console.WriteLine(x.ToVectorString());

            } while (z.L2Norm() > eps);

            Console.WriteLine($"Answer:\n {x.ToVectorString()}");

            var delta = Vector<double>.Build.Dense(n, 0);
            for (int i = 0; i < F.Count; i++)
            {
                F[i] = GetNFunctionValue(i, z.ToArray());
            }

            Console.WriteLine($"Delta vector:\n{z.ToVectorString()}");
            Console.WriteLine($"Delta norm:{z.L2Norm()}");
        }
        public static void RunRelaxationMethod(Func<Double[], double>[] f)
        {
            Vector<double> xk = Vector<double>.Build.Random(f.Length);
            Console.WriteLine("Enter x vector:");
            for (int i = 0; i < xk.Count; i++)
            {
                xk[i] = double.Parse(Console.ReadLine());
            }

            Vector<double> xkp1 = Vector<double>.Build.Random(f.Length);

            Vector<double> F = Vector<double>.Build.Dense(f.Length, 0);
            for (int i = 0; i < F.Count; i++)
            {
                F[i] = f[i].Invoke(xk.ToArray());
            }

            double t = 0.1;

            //первое приближение
            xkp1 = xk - t * F;
            Console.WriteLine("First approximation:");
            Console.WriteLine(xkp1.ToVectorString());
            do
            {
                xk = xkp1;
                for (int i = 0; i < F.Count; i++)
                {
                    F[i] = f[i].Invoke(xk.ToArray());
                }
                Console.WriteLine(xkp1.ToVectorString());
                xkp1 = xk - t * F;
            } while ((xkp1 - xk).L2Norm() > eps);

            Console.WriteLine(xkp1);

            var delta = Vector<double>.Build.Dense(f.Length, 0);
            for (int i = 0; i < F.Count; i++)
            {
                delta[i] = f[i].Invoke(xkp1.ToArray());
            }
            Console.WriteLine($"Delta vector:\n{delta.ToVectorString()}");
            Console.WriteLine($"Delta norm:{delta.L2Norm()}");
        }
        public static void RunRelaxationMethod(int n)
        {
            Vector<double> xk = Vector<double>.Build.Random(n);
            Console.WriteLine("Enter x vector:");
            for (int i = 0; i < xk.Count; i++)
            {
                xk[i] = double.Parse(Console.ReadLine());
            }

            Vector<double> xkp1 = Vector<double>.Build.Random(n);

            Vector<double> F = Vector<double>.Build.Dense(n);
            for (int i = 0; i < F.Count; i++)
            {
                F[i] = GetNFunctionValue(i, xk.ToArray());
            }

            double t = 0.00001;
                
            //первое приближение
            xkp1 = xk - t * F;
            Console.WriteLine("First approximation:");
            Console.WriteLine(xkp1.ToVectorString());
            do
            {
                xk = xkp1;
                for (int i = 0; i < F.Count; i++)
                {
                    F[i] = GetNFunctionValue(i, xk.ToArray());
                }
                Console.WriteLine(xkp1.ToVectorString());
                xkp1 = xk - t * F;
            } while ((xkp1 - xk).L2Norm() > eps);

            Console.WriteLine(xkp1);

            var delta = Vector<double>.Build.Dense(n, 0);
            for (int i = 0; i < F.Count; i++)
            {
                delta[i] = GetNFunctionValue(i, xk.ToArray());
            }
            Console.WriteLine($"Delta vector:\n{delta.ToVectorString()}");
            Console.WriteLine($"Delta norm:{delta.L2Norm()}");
        }
    }
}
