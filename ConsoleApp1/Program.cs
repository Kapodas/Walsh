using System;
using System.Drawing;
using System.Drawing.Drawing2D;
using System.Drawing.Imaging;
using System.Numerics;

namespace Hadamard_Matrix
{
    class Program
    {
        const double Pixel = 0.004;
        const double Real_Diam_Circle = 50;
        public static int Pix_Diam = Convert.ToInt32(Real_Diam_Circle / Pixel);

        //Генерация матрицы адамара n*n
        public static Bitmap RotateImg(Bitmap bmp, float angle)
        {
            angle = angle % 360;
            if (angle > 180)
                angle -= 360;
            System.Drawing.Imaging.PixelFormat pf = default(System.Drawing.Imaging.PixelFormat);

            float sin = (float)Math.Abs(Math.Sin(angle * Math.PI / 180.0)); // this function takes radians
            float cos = (float)Math.Abs(Math.Cos(angle * Math.PI / 180.0)); // this one too
            float newImgWidth = sin * bmp.Height + cos * bmp.Width;
            float newImgHeight = sin * bmp.Width + cos * bmp.Height;
            float originX = 0f;
            float originY = 0f;
            if (angle > 0)
            {
                if (angle <= 90)
                    originX = sin * bmp.Height;
                else
                {
                    originX = newImgWidth;
                    originY = newImgHeight - sin * bmp.Width;
                }
            }
            else
            {
                if (angle >= -90)
                    originY = sin * bmp.Width;
                else
                {
                    originX = newImgWidth - sin * bmp.Height;
                    originY = newImgHeight;
                }
            }
            Bitmap newImg = new Bitmap((int)newImgWidth, (int)newImgHeight);
            Graphics g = Graphics.FromImage(newImg);
            g.TranslateTransform(originX, originY); // offset the origin to our calculated values
            g.RotateTransform(angle); // set up rotate
            g.InterpolationMode = System.Drawing.Drawing2D.InterpolationMode.HighQualityBilinear;
            g.DrawImageUnscaled(bmp, 0, 0); // draw the image at 0, 0
            g.Dispose();
            return newImg;
        }
        public static void generate(int M)
        {
            int n = (int)Math.Pow(2, M);
            int[,] hadamard = new int[n, n];
            hadamard[0, 0] = 1;
            for (int k = 1; k < n; k += k)
            {
                for (int i = 0; i < k; i++)
                {
                    for (int j = 0; j < k; j++)
                    {
                        hadamard[i + k, j] = hadamard[i, j];
                        hadamard[i, j + k] = hadamard[i, j];
                        hadamard[i + k, j + k] = -hadamard[i, j];
                    }
                }
            }
            int Lam = 1;

            for (int i = 1; i <= n; i++)
            {
                Lam = i * n;
            }

            int[][,] RevHad = new int[Lam][,];
            for (int k = 0; k < Lam; k++)
            {
                RevHad[k] = new int[n, n];
            }
            RevHad[0] = hadamard;
            int L = 0;
            int Kar = 1;

            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    hadamard[i, j] = -hadamard[i, j];

                }

                for (int P = 0; P < n; P++)
                {
                    for (int j = 0; j < n; j++)
                    {
                        RevHad[Kar][P, j] = hadamard[P, j];

                    }
                }
                Kar++;

            }



            for (int j = 0; j < n; j++)
            {
                for (int i = 0; i < n; i++)
                {
                    hadamard[i, j] = -hadamard[i, j];

                }


                for (int P = 0; P < n; P++)
                {
                    for (int i = 0; i < n; i++)
                    {
                        RevHad[Kar][i, P] = hadamard[i, P];

                    }
                }

                Kar++;
            }



            Bitmap bitmap = new Bitmap(n, n);
            for (L = 0; L < Kar; L++)
            {
                for (int i = 0; i < n - 1; i++)
                {
                    for (int j = 0; j < n - 1; j++)
                    {
                        if (RevHad[L][i, j] == 1) { bitmap.SetPixel(i, j, Color.White); }
                        else { bitmap.SetPixel(i, j, Color.Black); }
                    }
                }
                bitmap.Save(@"C:\Users\osepa\source\repos\ConsoleApp1\Matrix\Matrix\" + "Matrix Size "
                            + n.ToString() + "x" + n.ToString() + @"\Matrix" + L.ToString() + ".bmp", ImageFormat.Bmp);
            }
            float ang = 0;

            for (L = 0; L < Kar; L++)
            {
                bitmap = new Bitmap(@"C:\Users\osepa\source\repos\ConsoleApp1\Matrix\Matrix\" + "Matrix Size "
                            + n.ToString() + "x" + n.ToString() + @"\Matrix" + L.ToString() + ".bmp");

                bitmap = RotateImg(bitmap, ang);

                ang = ang + 10;

                bitmap.Save(@"C:\Users\osepa\source\repos\ConsoleApp1\Matrix\Matrix\Matrix" + L.ToString() + ".bmp", ImageFormat.Bmp);
            }
        }

        //Ввод круга где матрицу разместим

        public static void Circle()
        {
            Bitmap Mask;
            Mask = new Bitmap(Pix_Diam + 1, Pix_Diam + 1);
            Pen myPen = new Pen(Color.Black);
            Brush Br = Brushes.Blue;

            for (int i = 0; i < Pix_Diam; i++)
            {
                for (int j = 0; j < Pix_Diam; j++)
                {
                    Mask.SetPixel(i, j, Color.White);
                }
            }

            Graphics MyGrap = Graphics.FromImage(Mask);
            Graphics graphics = Graphics.FromImage(Mask);

            MyGrap.DrawEllipse(myPen, 0, 0, Pix_Diam, Pix_Diam);
            graphics.FillEllipse(Br, 0, 0, Pix_Diam, Pix_Diam);
            myPen = new Pen(Color.Red);
            MyGrap.DrawEllipse(myPen, Pix_Diam / 4, Pix_Diam / 4, Pix_Diam / 2, Pix_Diam / 2);

            Mask.Save(@"C:\Users\1488\source\repos\Hadamard Matrix\Matrix\" + "Mask.bmp", ImageFormat.Bmp);
        }
        //Нанесение матрицы на окружность
        public static void Hadamard_Mask(int M)
        {
            int n = (int)Math.Pow(2, M);
            Bitmap Circle = new Bitmap(@"C:\Users\osepa\source\repos\ConsoleApp1\Matrix\Mask.bmp");
            Bitmap Matrix;
            Bitmap HadMask = Circle;
            Graphics graphics = Graphics.FromImage(HadMask);
            Image Img;
            int L = 0;
            for (int i = 0; i < Pix_Diam; i++)
            {
                for (int j = 0; j < Pix_Diam; j++)
                {
                    if (Circle.GetPixel(i, j).ToArgb() == Color.Red.ToArgb()
                        && Circle.GetPixel(i + n / 2, j).ToArgb() != Color.White.ToArgb()
                        && Circle.GetPixel(i + n / 2, j).ToArgb() != Color.Black.ToArgb()

                        && Circle.GetPixel(i - n / 2, j).ToArgb() != Color.Black.ToArgb()
                        && Circle.GetPixel(i - n / 2, j).ToArgb() != Color.White.ToArgb()

                        && Circle.GetPixel(i, j + n / 2).ToArgb() != Color.Black.ToArgb()
                        && Circle.GetPixel(i, j + n / 2).ToArgb() != Color.White.ToArgb()

                        && Circle.GetPixel(i, j - n / 2).ToArgb() != Color.Black.ToArgb()
                        && Circle.GetPixel(i, j - n / 2).ToArgb() != Color.White.ToArgb()

                        && Circle.GetPixel(i - n / 2, j - n / 2).ToArgb() != Color.Black.ToArgb()
                        && Circle.GetPixel(i - n / 2, j - n / 2).ToArgb() != Color.White.ToArgb()


                        && Circle.GetPixel(i + n / 2, j - n / 2).ToArgb() != Color.Black.ToArgb()
                        && Circle.GetPixel(i + n / 2, j - n / 2).ToArgb() != Color.White.ToArgb()

                        && Circle.GetPixel(i - n / 2, j + n / 2).ToArgb() != Color.Black.ToArgb()
                        && Circle.GetPixel(i - n / 2, j + n / 2).ToArgb() != Color.White.ToArgb()

                        && Circle.GetPixel(i + n / 2, j + n / 2).ToArgb() != Color.Black.ToArgb()
                        && Circle.GetPixel(i + n / 2, j + n / 2).ToArgb() != Color.White.ToArgb()

                        && Circle.GetPixel(i + n, j).ToArgb() != Color.Black.ToArgb()
                        && Circle.GetPixel(i + n, j).ToArgb() != Color.White.ToArgb()

                        && Circle.GetPixel(i - n, j).ToArgb() != Color.Black.ToArgb()
                        && Circle.GetPixel(i - n, j).ToArgb() != Color.White.ToArgb()

                        && Circle.GetPixel(i, j + n).ToArgb() != Color.Black.ToArgb()
                        && Circle.GetPixel(i, j + n).ToArgb() != Color.White.ToArgb()

                        && Circle.GetPixel(i, j - n).ToArgb() != Color.Black.ToArgb()
                        && Circle.GetPixel(i, j - n).ToArgb() != Color.White.ToArgb()

                        && Circle.GetPixel(i - n, j - n).ToArgb() != Color.Black.ToArgb()
                        && Circle.GetPixel(i - n, j - n).ToArgb() != Color.White.ToArgb()

                        && Circle.GetPixel(i + n, j - n).ToArgb() != Color.Black.ToArgb()
                        && Circle.GetPixel(i + n, j - n).ToArgb() != Color.White.ToArgb()

                        && Circle.GetPixel(i - n, j + n).ToArgb() != Color.Black.ToArgb()
                        && Circle.GetPixel(i - n, j + n).ToArgb() != Color.White.ToArgb()

                        && Circle.GetPixel(i + n, j + n).ToArgb() != Color.Black.ToArgb()
                        && Circle.GetPixel(i + n, j + n).ToArgb() != Color.White.ToArgb()


                        && Circle.GetPixel(i + n / 4, j).ToArgb() != Color.White.ToArgb()
                        && Circle.GetPixel(i + n / 4, j).ToArgb() != Color.Black.ToArgb()

                        && Circle.GetPixel(i - n / 4, j).ToArgb() != Color.Black.ToArgb()
                        && Circle.GetPixel(i - n / 4, j).ToArgb() != Color.White.ToArgb()

                        && Circle.GetPixel(i, j + n / 4).ToArgb() != Color.Black.ToArgb()
                        && Circle.GetPixel(i, j + n / 4).ToArgb() != Color.White.ToArgb()

                        && Circle.GetPixel(i, j - n / 4).ToArgb() != Color.Black.ToArgb()
                        && Circle.GetPixel(i, j - n / 4).ToArgb() != Color.White.ToArgb()

                        && Circle.GetPixel(i - n / 4, j - n / 4).ToArgb() != Color.Black.ToArgb()
                        && Circle.GetPixel(i - n / 4, j - n / 4).ToArgb() != Color.White.ToArgb()


                        && Circle.GetPixel(i + n / 4, j - n / 4).ToArgb() != Color.Black.ToArgb()
                        && Circle.GetPixel(i + n / 4, j - n / 4).ToArgb() != Color.White.ToArgb()

                        && Circle.GetPixel(i - n / 4, j + n / 4).ToArgb() != Color.Black.ToArgb()
                        && Circle.GetPixel(i - n / 4, j + n / 4).ToArgb() != Color.White.ToArgb()

                        && Circle.GetPixel(i + n / 4, j + n / 4).ToArgb() != Color.Black.ToArgb()
                        && Circle.GetPixel(i + n / 4, j + n / 4).ToArgb() != Color.White.ToArgb())

                    {
                        Matrix = new Bitmap(@"C:\Users\osepa\source\repos\Hadamard Matrix\Matrix\Matrix\" + "Matrix Size "
                            + n.ToString() + "x" + n.ToString() + @"\Matrix" + L.ToString() + ".bmp");
                        Img = Matrix;

                        graphics.DrawImage(Img, i, j);
                        L++;
                        if (L > 2 * n) { L = 0; }
                    }


                }

            }
            HadMask.Save(@"C:\Users\osepa\source\repos\Hadamard Matrix\Matrix\Mask\Mask " + n.ToString() + " x "
                            + n.ToString() + ".bmp");
        }
        private static double[] linespace(double x1, double x2, int n)
        {
            double step = (x2 - x1) / (n - 1);
            double[] y = new double[n];
            for (int i = 0; i < n; i++)
            {
                y[i] = x1 + step * i;
            }
            return y;
        }
        //Матрицы Уолша
        public static void Walsh(int M)
        {
            int n = (int)Math.Pow(2, M);
            int L;
            Bitmap[] Matrix = new Bitmap[257];
            
            
                Matrix[0] = new Bitmap(@"C:\Users\osepa\source\repos\ConsoleApp1\Matrix\Matrix\" + "Matrix Size "
                                       + n.ToString() + "x" + n.ToString() + @"\Matrix" + 2.ToString() + ".bmp");
            

            string G;
            string GR;
            
            int[] Gray = new int[256];
            for (int i = 0; i < 256; i++)
            {
                int Perehod;
                GR = "";
                G = Convert.ToString(Convert.ToByte(i), 2);
                
                for (int j = G.Length - 1; j >= 0; j--){ GR = GR + (G[j]); }

                Perehod = Convert.ToInt32(GR, 2);

                Gray[i] = Perehod ^ (Perehod >> 1);
                
                
            }
            Bitmap WM = new Bitmap(256, 256);
            
                for (int i = 0; i < 256; i++)
                {
                    for (int j = 0; j < 256; j++)
                    {
                        if (Matrix[0].GetPixel(Gray[i], j).ToArgb() == Color.White.ToArgb()) { WM.SetPixel(i, j, Color.White); }
                        else { WM.SetPixel(i, j, Color.Black); }

                    }
                }
            //WM.SetResolution(1024, 1024);

                WM.Save(@"C:\Users\osepa\source\repos\ConsoleApp1\Matrix\Matrix\WalshMatrix\Matrix" + 1211.ToString() + ".bmp", ImageFormat.Bmp);
            Bitmap mtr = new Bitmap(16, 16);
            Bitmap str = new Bitmap(1024, 1024);
            int p = 0;
            int m = 0;
            L = 0;
                
            while (p < 128)
            {
                while (m < 128)
                {
                    for (int i = p; i < p + 16; i++)
                    {

                        for (int j = m; j < m + 16; j++)
                        {

                            {
                                if (WM.GetPixel(i, j).ToArgb() == Color.White.ToArgb()) { mtr.SetPixel(i - p, j - m, Color.White); }
                                else { mtr.SetPixel(i - p, j - m, Color.Black); }

                            }
                        }
                    }
                    for (int i = 0; i < 16; i++)
                    {
                        for (int j = 0; j < 16; j++)
                        {
                            if (mtr.GetPixel(i, j).ToArgb() == Color.White.ToArgb())
                            {
                                for (int k = 64 * i; k < 64 * (i + 1); k++)
                                {
                                    for (int b = 64 * j; b < 64 * (j + 1); b++)
                                    {
                                        str.SetPixel(k, b, Color.White);
                                    }
                                }
                            }
                            else
                            {
                                for (int k = 64 * i; k < 64 * (i + 1); k++)
                                {
                                    for (int b = 64 * j; b < 64 * (j + 1); b++)
                                    {
                                        str.SetPixel(k, b, Color.Black);
                                    }
                                }
                            }
                        }
                    }
                    str.Save(@"C:\Users\osepa\source\repos\ConsoleApp1\Matrix\Matrix\WalshMatrix\Matrix" + L.ToString() + ".bmp", ImageFormat.Bmp);
                    L++;
                    m = m + 16; 
                }
                
                
                m = 0;
                p = p + 16;
            }
        }
        //Эксперемент
        public static void Exp(int M)
        {
            Bitmap Pict = new Bitmap(1024, 1024);
            Pict = new Bitmap(@"C:\Users\osepa\source\repos\ConsoleApp1\Matrix\Matrix\example2.bmp");
            
            double[,] VecPic = new double[1024, 1024];
            double[,] Walsh = new double[1024, 1024];
            int n = 0;
            double max;
            double min;
            for (int i = 0; i < 1024; i++)
            {
                for (int j = 0; j < 1024; j++)
                {
                    VecPic[i, j] = Convert.ToDouble(Pict.GetPixel(i, j).R) / 255d;
                    n++;

                }
            }
            max = VecPic[0, 0];
            min = VecPic[0, 0];
            for (int i = 0; i < 1024; i++)
            {
                for (int j = 0; j < 1024; j++)
                {
                    if (max < VecPic[i, j])
                    {
                        max = VecPic[i, j];
                    }
                    if (min > VecPic[i, j])
                    {
                        min = VecPic[i, j];
                    }

                }
            }

            for (int i = 0; i < 1024; i++)
            {
                for (int j = 0; j < 1024; j++)
                {

                    if (VecPic[i, j] == max)
                    {
                        VecPic[i, j] = 1;
                    }

                    if (VecPic[i, j] == min)
                    {
                        VecPic[i, j] = -1;
                    }
                    if (VecPic[i, j] < (max + min) / 2)
                    {
                        VecPic[i, j] = -(VecPic[i, j]); //* 2;
                    }
                    if (VecPic[i, j] > (min + max) / 2)
                    {
                        VecPic[i, j] = VecPic[i, j]; // / 2; 
                    }
                }
            }
            
            double[,] PerPic = new double[1024, 1024];
            double[] Sred = new double[256];
            double[,] PerWalsh = new double[1024, 1024];
            double[,] EndPic = new double[1024, 1024];
            for(int i = 0; i < 1024; i++)
            {
                for(int j = 0; j < 1024; j++)
                {
                    EndPic[i, j] = 0;
                }
            }
            Bitmap Matr;
            for (int L = 0; L < 64; L++)
            {
                Sred[L] = 0;
                n = 0;
              //  if (L < 64)
              //  {
                    Matr = new Bitmap(@"C:\Users\osepa\source\repos\ConsoleApp1\Matrix\Matrix\WalshMatrix\Matrix" + L.ToString() + ".bmp");
             //   }
              //  else { Matr = new Bitmap(@"C:\Users\osepa\source\repos\ConsoleApp1\Matrix\Matrix\WalshMatrix\Matrix" + 1212.ToString() + ".bmp"); }
                for (int i = 0; i < 1024; i++)
                {
                    for (int j = 0; j < 1024; j++)
                    {
                        Walsh[i, j] = Convert.ToDouble(Matr.GetPixel(i, j).R)/ 255d; //.ToArgb();
                       if (Walsh[i, j] == 0)
                        {
                            Walsh[i, j] = -1;
                        } 

                    }
                }

                double kr = 0;
                for(int i = 0; i < 1024; i++)
                {
                    for (int j = 0; j < 1024; j++)
                    {
                        
                        PerPic[i, j] = VecPic[i, j] * Walsh[i, j];
                        if (Walsh[i, j] == 1)
                        {
                            Sred[L] = Sred[L] + PerPic[i, j];
                            kr++;
                        }
                    }
                }
                Sred[L] = Sred[L] / (kr);
                for(int i= 0; i < 1024; i++)
                {
                    for(int j = 0; j < 1024; j++)
                    {
                        PerWalsh[i, j] = Sred[L] * Walsh[i, j];
                        EndPic[i, j] = EndPic[i, j] + PerWalsh[i, j];
                    }
                }
                
            }
            max = EndPic[0, 0];
            min = EndPic[0, 0];
            for (int i = 0; i < 1024; i++)
            {
                for (int j = 0; j < 1024; j++)
                {
                    if (max < EndPic[i, j])
                    {
                        max = EndPic[i, j];
                    }
                    if(min > EndPic[i,j])
                    {
                        min = EndPic[i,j];
                    }

                }
            }

            for (int i = 0; i < 1024; i++)
            {
                for (int j = 0; j < 1024; j++)
                {
                    if (EndPic[i, j] < 0 && EndPic[i, j] != min)
                    {
                        EndPic[i, j] = -EndPic[i, j];
                    }
                    if (EndPic[i, j] > 0 && EndPic[i, j] != max)
                    {
                        EndPic[i, j] = EndPic[i, j];
                    }
                }
            }

            max = EndPic[0, 0];
            min = EndPic[0, 0];
            
            for (int i = 0; i < 1024; i++)
            {
                for (int j = 0; j < 1024; j++)
                {
                    if (max < EndPic[i, j])
                    {
                        max = EndPic[i, j];
                    }
                    if (min > EndPic[i, j])
                    {
                        min = EndPic[i, j];
                    }

                }
            }
            for (int i = 0; i < 1024; i++)
            {
                for(int j = 0; j < 1024; j++)
                {
                    
                    if (EndPic[i, j] == max)
                    {
                        EndPic[i, j] = 255;
                    }
                    
                    if (EndPic[i, j] == min)
                    {
                        EndPic[i, j] = 0;
                    }
                    if (EndPic[i, j] < max && EndPic[i, j] > min)
                    {
                        EndPic[i, j] = EndPic[i, j] * 255d / max;
                    }
                }
            }
            
            Bitmap bitmap = new Bitmap(1024, 1024);
            
            for (int i = 0; i < 1024; i++)
            {
                for (int j = 0; j < 1024; j++) 
                {   
                    //if(EndPic[i, j] < 25.5)
                    //{   
                    //    EndPic[i, j] = EndPic[i, j] * 10;
                    //}
                    
                    bitmap.SetPixel(i, j, Color.FromArgb(Convert.ToInt32(EndPic[i, j]), Convert.ToInt32(EndPic[i, j]), Convert.ToInt32(EndPic[i, j])));
                }
            }
            bitmap.Save(@"C:\Users\osepa\source\repos\ConsoleApp1\Matrix\Pic.bmp");
        }
        
        static public void Main()
        {
            
            string n = Console.ReadLine();
            int M = Convert.ToInt32(n);
            //generate(M);
            //Circle();
            //Hadamard_Mask(M);
            //Walsh(M);
            Exp(M);
        }
    }
}