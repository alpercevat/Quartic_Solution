using System;
using System.Collections.Generic;
using System.Numerics;

namespace Sorting
{
    class Program
    {

            public static void ComputeBasicRoot(double D, double E)
            {
                double root = 0;

                root = -E / D;
                Console.WriteLine(" **********Basic Root********** ");
                Console.WriteLine("Root1: " + root);
            }

            public static Complex[] ComputeQuadraticRoot(double C, double D, double E)
            {
                Complex [] quadratic_roots = new Complex [2];

                double root1 = 0;
                double root2 = 0;
                double equation = 0;
                double K;

                equation = D * D - 4 * C * E;

                if (C == 0)
                {
                    //Console.WriteLine("This equation is not a Quadratic Equation");
                }

                else if (equation > 0)
                {
                    //Console.WriteLine("*******Real Roots*******");
                    root1 = (-D + Math.Sqrt(equation)) / (2 * C);
                    root2 = (-D - Math.Sqrt(equation)) / (2 * C);

                    quadratic_roots[0] = root1;
                    quadratic_roots[1] = root2;

                    Console.WriteLine("Root 1: " + root1);
                    Console.WriteLine("Root 2: " + root2);
                }

                else if (equation == 0)
                {
                    //Console.WriteLine("*********Dejenere Real Roots**********");
                    root1 = (-D) / (2 * C);
                    root2 = (-D) / (2 * C);
                    quadratic_roots[0] = root1;
                    quadratic_roots[1] = root2;

                    Console.WriteLine("Root 1: " + root1);
                    Console.WriteLine("Root 2: " + root2);
                }

                else
                {
                    Complex I = new Complex(0, 1);
                    K = Math.Sqrt((-1) * equation);
                    K = (-1) * K;
                    root1 = (-D) / (2 * C);
                    root2 = (K) / (2 * C);

                    Complex root11;
                    Complex root12;

                    root11 = root1 + I * root2;
                    root12 = root1 - I * root2;

                    quadratic_roots[0] = root11;
                    quadratic_roots[1] = root12;
                    //Console.WriteLine("********İmag Roots*************");
                    Console.WriteLine("Root 1: " + root11);
                    Console.WriteLine("Root 1: " + root12);
                }
            //Console.WriteLine("Roots1: " + quadratic_roots[0]);
            //Console.WriteLine("Roots2: " + quadratic_roots[1]);
            return quadratic_roots;


            }

            public static Complex[] ComputeCubicRoot(double B, double C, double D, double E)
            {
                Complex[] roots = new Complex[3];

                double q, r, sum, u0, u, v0, v, x2, x3, y2, y3, k, q1, sigma1, sigma2, sigma3;
                double theta = 0;

                ///" B2.X^3 + B3.X^2 + B4.X + B5 "///   RESOLVENT CUBİC

                // B2 = B    and B3 = C  and B4 = D and B5 = E 
                //B = B/B;   // Oluşturulan " Resolvent Cubic " denkleminin ilk katsayısı 1'e eşit olmalı.  ------->  B2.X^3
                C = C / B;   // " Resolvent Cubic " denkleminin ikinci katsayısı -------> B3.X^2
                D = D / B;   // " Resolvent Cubic " denkleminin üçüncü katsayısı  -------> B4.X
                E = E / B;   // " Resolvent Cubic " denkleminin dördüncü katsayısı  -------> B5.X^0

                ///" B2.X^3 + B3.X^2 + B4.X + B5 "///
                //Console.WriteLine("****Katsayılar****");
                //Console.WriteLine("B: " + B);
                //Console.WriteLine("C: " + C);
                //Console.WriteLine("D: " + D);
                //Console.WriteLine("E: " + E);

                q = (D / 3) - ((C * C) / 9); // SORUNLU YER
                r = ((D * C - 3 * E) / 6) - ((C * C * C) / 27); // SORUNLU YER 

                //Console.WriteLine("q: " + q);
                //Console.WriteLine("r: " + r);

                sum = Math.Pow(r, 2) + Math.Pow(q, 3);
                //Console.WriteLine("sum: " + sum);
                if (sum > 0)
                {
                    u0 = r + Math.Sqrt(sum);
                    //Console.WriteLine("u0: " + u0);
                    if (u0 < 0)
                    {
                        u0 = (-1) * u0;
                        u = Math.Pow(u0, (1.0 / 3.0));
                        u = (-1) * u;
                        //Console.WriteLine("u: " + u);
                    }
                    else
                    {
                        u = Math.Pow(u0, (1.0 / 3.0));
                        //Console.WriteLine("u: " + u);
                    }

                    v0 = r - Math.Sqrt(sum);
                    //Console.WriteLine("v0: " + v0);
                    if (v0 < 0)
                    {
                        v0 = (-1) * v0;
                        v = Math.Pow(v0, (1.0 / 3.0));
                        v = (-1) * v;
                        //Console.WriteLine("v: " + v);

                    }
                    else
                    {
                        v = Math.Pow(v0, (1.0 / 3.0));
                        //Console.WriteLine("v: " + v);
                    }


                    x2 = -((u + v) / 2) - (C / 3);
                    //Console.WriteLine("x2: " + x2);
                    x3 = -((u + v) / 2) - (C / 3);
                    //Console.WriteLine("x3: " + x3);

                    y2 = ((Math.Pow(3, (1.0 / 2.0))) * (u - v)) / 2;
                    //Console.WriteLine("y2: " + y2);
                    y3 = (-1) * y2;

                    double root1_real = u + v - (C / 3);
                    Complex root1 = new Complex(root1_real, 0);
                    Complex root2 = new Complex(x2, +y2);
                    Complex root3 = new Complex(x2, -y2);


                    roots[0] = root1;
                    roots[1] = root2;
                    roots[2] = root3;

                    //Console.WriteLine("*******Cubic Roots**********");
                    Console.WriteLine("Root1: " + root1);
                    Console.WriteLine("Root2: " + root2);
                    Console.WriteLine("Root3: " + root3);
                    return roots;
                }
                else
                {

                    if (q == 0)
                    {
                        theta = 0;
                    }
                    else if (q < 0)
                    {
                        q1 = (-q) * (-q) * (-q);
                        k = Math.Pow(q1, (1.0 / 2.0));
                        theta = Math.Acos(r / k);
                        //Console.WriteLine("theta: " + theta);
                    }
                    sigma1 = theta / 3;
                    sigma2 = sigma1 - ((2 * Math.PI) / 3);
                    sigma3 = sigma1 + ((2 * Math.PI) / 3);

                    //Console.WriteLine("sigma1: " + sigma1);
                    //Console.WriteLine("sigma2: " + sigma2);
                    //Console.WriteLine("sigma3: " + sigma3);

                    double root1_real = 2 * Math.Pow(-q, (1.0 / 2.0)) * Math.Cos(sigma1) - (C / 3);
                    double root2_real = 2 * Math.Pow(-q, (1.0 / 2.0)) * Math.Cos(sigma2) - (C / 3);
                    double root3_real = 2 * Math.Pow(-q, (1.0 / 2.0)) * Math.Cos(sigma3) - (C / 3);


                    Complex root1 = new Complex(root1_real, 0);
                    Complex root2 = new Complex(root2_real, 0);
                    Complex root3 = new Complex(root3_real, 0);

                    roots[0] = root1;
                    roots[1] = root2;
                    roots[2] = root3;


                    //Console.WriteLine("*******Cubic Roots**********");
                    Console.WriteLine("Root1: " + root1);
                    Console.WriteLine("Root2: " + root2);
                    Console.WriteLine("Root3: " + root3);
                    return roots;
                    

            }

            }

            public static Complex[] ComputeQuarticRoot(double A, double B, double C, double D, double E)
            {
                Complex[] quartic_roots = new Complex[4];

                Complex X1;
                Complex X2;
                Complex X3;
                Complex X4;

                double a1 = A / A;
                double b1 = B / A;
                double c1 = C / A;
                double d1 = D / A;
                double e1 = E / A;

                double cc = b1 / 4;
                double f = c1 - ((3 * b1 * b1) / 8);
                double g = d1 + ((b1 * b1 * b1) / 8) - ((b1 * c1) / 2);
                double h = e1 - ((3 * b1 * b1 * b1 * b1) / 256) + ((b1 * b1 * c1) / 16) - ((b1 * d1) / 4);

                double a2 = 1;
                double b2 = f / 2;
                double c2 = ((f * f) - 4 * h) / 16;
                double d2 = (-(g * g)) / 64;

                //Console.WriteLine(" a2: " + a2);
                //Console.WriteLine(" b2: " + b2);
                //Console.WriteLine(" c2: " + c2);
                //Console.WriteLine(" d2: " + d2);


                Complex[] roots = new Complex[2]; //3 elemnanlı bir array yaratıldı
                Complex[] ordered_roots = new Complex[3];

                roots = ComputeCubicRoot(a2, b2, c2, d2);
                //Console.WriteLine(" Roots: " + roots);

                if (roots[1].Imaginary == 0)
                {
                    if (roots[0].Real > roots[1].Real && roots[0].Real > roots[2].Real)
                    {
                        ordered_roots = roots;
                        //Console.WriteLine(" Ordered Roots: " + ordered_roots);
                    }
                    else if (roots[1].Real > roots[0].Real && roots[1].Real > roots[2].Real)
                    {
                        ordered_roots[0] = roots[1];
                        ordered_roots[1] = roots[0];
                        ordered_roots[2] = roots[2];

                        //Console.WriteLine(" Ordered Root1: " + ordered_roots[0]);
                        //Console.WriteLine(" Ordered Root2: " + ordered_roots[1]);
                        //Console.WriteLine(" Ordered Root3: " + ordered_roots[2]);
                    }
                    else
                    {
                        ordered_roots[0] = roots[2];
                        ordered_roots[1] = roots[0];
                        ordered_roots[2] = roots[1];

                        //Console.WriteLine(" Ordered Root1: " + ordered_roots[0]);
                        //Console.WriteLine(" Ordered Root2: " + ordered_roots[1]);
                        //Console.WriteLine(" Ordered Root3: " + ordered_roots[2]);

                    }

                }

                //Complex Y1     = roots[0];
                double Y1_real = roots[0].Real;
                if (Math.Abs(Y1_real) < 1e-15)
                {
                    Y1_real = 0.0;
                }
                //Console.WriteLine(" Y1_Real: " + Y1_real);

                double Y1_img = roots[0].Imaginary;
                if (Math.Abs(Y1_img) < 1e-15)
                {
                    Y1_img = 0.0;
                }

                Complex Y2 = roots[1];
                double Y2_real = roots[1].Real;
                if (Math.Abs(Y2_real) < 1e-15)
                {
                    Y2_real = 0.0;
                }
                //Console.WriteLine(" Y2_Real: " + Y2_real);

                double Y2_img = roots[1].Imaginary;
                if (Math.Abs(Y2_img) < 1e-15)
                {
                    Y2_img = 0.0;
                }
                //Console.WriteLine(" Y2_img: " + Y2_img);

                Complex Y3 = roots[2];
                //Console.WriteLine(" Y3: " + Y3);
                double Y3_real = roots[2].Real;
                if (Math.Abs(Y3_real) < 1e-15)
                {
                    Y3_real = 0.0;
                }
                //Console.WriteLine(" Y3_Real: " + Y3_real);

                double Y3_img = roots[2].Imaginary;
                if (Math.Abs(Y3_img) < 1e-15)
                {
                    Y3_img = 0.0;
                }
                //Console.WriteLine(" Y3_img: " + Y3_img);


                double o;
                double T;
                T = Y2_real * Y3_real + Y2_img * Y2_img;

                if (T < 0)
                {
                    T = (-1) * T;
                    T = Math.Sqrt(T);
                    T = (-1) * T;
                }

                else
                {
                    T = Math.Sqrt(T);
                }

                if (g > 0)
                {
                    o = 2 * T;
                    //Console.WriteLine(" o: " + o);
                }
                else
                {
                    o = -2 * T;
                    //Console.WriteLine(" o: " + o);
                }


                double X1_real;
                double X1_imag;
                double X2_real;
                double X2_imag;
                double X3_real;
                double X3_imag;
                double X4_real;
                double X4_imag;

                X1 = Complex.Sqrt(Y1_real) + Complex.Sqrt(Y2_real + Y3_real - o) - cc;
                X1_real = X1.Real;
                X1_imag = X1.Imaginary;

                X2 = Complex.Sqrt(Y1_real) - Complex.Sqrt(Y2_real + Y3_real - o) - cc;
                X2_real = X2.Real;
                X2_imag = X2.Imaginary;

                X3 = -Complex.Sqrt(Y1_real) + Complex.Sqrt(Y2_real + Y3_real + o) - cc;
                X3_real = X3.Real;
                X3_imag = X3.Imaginary;

                X4 = -Complex.Sqrt(Y1_real) - Complex.Sqrt(Y2_real + Y3_real + o) - cc;
                X4_real = X4.Real;
                X4_imag = X4.Imaginary;


                if (Math.Abs(X1_real) < 1e-15)
                {
                    X1_real = 0.0;
                }
                if (Math.Abs(X1_imag) < 1e-15)
                {
                    X1_imag = 0.0;
                }

                if (Math.Abs(X2_real) < 1e-15)
                {
                    X2_real = 0.0;
                }
                if (Math.Abs(X2_imag) < 1e-15)
                {
                    X2_imag = 0.0;
                }

                if (Math.Abs(X3_real) < 1e-15)
                {
                    X3_real = 0.0;
                }
                if (Math.Abs(X3_imag) < 1e-15)
                {
                    X3_imag = 0.0;
                }

                if (Math.Abs(X4_real) < 1e-15)
                {
                    X4_real = 0.0;
                }
                if (Math.Abs(X4_imag) < 1e-15)
                {
                    X4_imag = 0.0;
                }

                Complex I = new Complex(0, 1);

                X1 = X1_real + I * X1_imag;
                X2 = X2_real + I * X2_imag;
                X3 = X3_real + I * X3_imag;
                X4 = X4_real + I * X4_imag;
                
                //Console.WriteLine(" roots: " + quartic_roots[0]);
                Console.WriteLine(" X1: " + X1);
                Console.WriteLine(" X2: " + X2);
                Console.WriteLine(" X3: " + X3);
                Console.WriteLine(" X4: " + X4);

            return quartic_roots;
                //Console.WriteLine(" X2: " + X2);
                //Console.WriteLine(" X3: " + X3);
                //Console.WriteLine(" X4: " + X4);


            }

            static void Main(string[] args)
            {
                double A;
                double B;
                double C;
                double D;
                double E;


                Console.WriteLine("Polynomial ----->  A.X^4  +  B.X^3  +  C.X^2  +  D.X  +   E = 0");

                Console.Write("A: ");
                A = double.Parse(Console.ReadLine());

                Console.Write("B: ");
                B = double.Parse(Console.ReadLine());

                Console.Write("C: ");
                C = double.Parse(Console.ReadLine());

                Console.Write("D: ");
                D = double.Parse(Console.ReadLine());

                Console.Write("E: ");
                E = double.Parse(Console.ReadLine());

                if (A != 0)
                {
                    ComputeQuarticRoot(A, B, C, D, E);
                }

                else
                {
                    if (B != 0)
                    {
                        ComputeCubicRoot(B, C, D, E);
                    }
                    else
                    {
                        if (C != 0)
                        {
                            ComputeQuadraticRoot(C, D, E);
                        }
                        else
                        {
                            if (D != 0)
                            {
                                ComputeBasicRoot(D, E);
                            }
                            else
                            {
                                Console.WriteLine("No Solution");
                            }
                        }
                    }
                }



            }


        





    }


}

