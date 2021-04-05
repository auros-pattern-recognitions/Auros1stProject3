using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;
using static System.Math;
using static System.Console;
using System.Diagnostics;

//
// SiO2, SiN, Si 물성을 활용하여,
// SiN(30nm) 와 SiO2(20nm) 가 반복되어 있는 구조에 대한 alpha, beta 를 계산한다.
//
namespace Auros1stProject3
{
    class Program
    {
        static void Main(string[] args)
        {
            //
            // "SiO2_new2.txt", "SiN_new2.txt", "Si_new2.txt" 파일로부터
            // 파장별 복소굴절률(n, k)을 읽어온다.
            //
            // 2021.03.29 이지원.
            //
            #region SiO2, SiN, Si 물성을 읽어온다.

            // "SiO2_new.txt" 파일 읽기.
            string[] SiO2_new2 = File.ReadAllLines("SiO2_new2.txt");  // SiO2 물성값 저장.(한 줄씩)
            string[] SingleLineData;                                // 한 줄의 스펙트럼 데이터를 임시로 저장할 배열.

            // 데이터의 첫번째 줄은 column 명이다.
            // 이를 제외하고 데이터를 받기 위해 LenData 변수를 선언한다.
            int LenData = SiO2_new2.Length - 1;
            double[] wavelength_SiO2 = new double[LenData];
            double[] n_SiO2 = new double[LenData];
            double[] k_SiO2 = new double[LenData];

            // SiO2_new 에 받은 데이터를 각 컬럼별로 저장한다.
            int StartIndex = 1;
            LenData = SiO2_new2.Length;
            for (int i = StartIndex; i < LenData; i++)
            {
                // tsv 형식의 데이터를 SingleLineData에 저장한다.
                SingleLineData = SiO2_new2[i].Split((char)0x09);  // 0x09 : 수평 탭.

                // 각 컬럼에 해당하는 데이터를 저장한다.
                wavelength_SiO2[i - 1] = Double.Parse(SingleLineData[0]);
                n_SiO2[i - 1] = Double.Parse(SingleLineData[1]);
                k_SiO2[i - 1] = Double.Parse(SingleLineData[2]);
            }


            // "SiN_new.txt" 파일 읽기.
            string[] SiN_new2 = File.ReadAllLines("SiN_new2.txt");  // SiN 물성값 저장.(한 줄씩)

            // 데이터의 첫번째 줄은 column 명이다.
            // 이를 제외하고 데이터를 받기 위해 LenData 변수를 선언한다.
            LenData = SiN_new2.Length - 1;
            double[] wavelength_SiN = new double[LenData];
            double[] n_SiN = new double[LenData];
            double[] k_SiN = new double[LenData];

            // SiN_new 에 받은 데이터를 각 컬럼별로 저장한다.
            LenData = SiN_new2.Length;
            for (int i = StartIndex; i < LenData; i++)
            {
                // tsv 형식의 데이터를 SingleLineData에 저장한다.
                SingleLineData = SiN_new2[i].Split((char)0x09);  // 0x09 : 수평 탭.

                // 각 컬럼에 해당하는 데이터를 저장한다.
                wavelength_SiN[i - 1] = Double.Parse(SingleLineData[0]);
                n_SiN[i - 1] = Double.Parse(SingleLineData[1]);
                k_SiN[i - 1] = Double.Parse(SingleLineData[2]);
            }


            // "Si_new.txt" 파일 읽기.
            string[] Si_new2 = File.ReadAllLines("Si_new2.txt");  // Si 기판 물성값 저장.(한 줄씩)

            // 데이터의 첫번째 줄은 column 명이다.
            // 이를 제외하고 데이터를 받기 위해 LenData 변수를 선언한다.
            LenData = Si_new2.Length - 1;
            double[] wavelength_Si = new double[LenData];
            double[] n_Si = new double[LenData];
            double[] k_Si = new double[LenData];

            // Si_new 에 받은 데이터를 각 컬럼별로 저장한다.
            LenData = Si_new2.Length;
            for (int i = StartIndex; i < LenData; i++)
            {
                // tsv 형식의 데이터를 SingleLineData에 저장한다.
                SingleLineData = Si_new2[i].Split((char)0x09);  // 0x09 : 수평 탭.

                // 각 컬럼에 해당하는 데이터를 저장한다.
                wavelength_Si[i - 1] = Double.Parse(SingleLineData[0]);
                n_Si[i - 1] = Double.Parse(SingleLineData[1]);
                k_Si[i - 1] = Double.Parse(SingleLineData[2]);
            }
            #endregion

            //
            // 매질에 따른 interface matrix, layer matrix 를 구하여
            // scattering matrix 를 통해 alpha, beta 를 계산한다.
            //
            // 2021.03.30 이지원.
            //
            #region scattering matrix 를 통해 alpha, beta 를 계산한다.

            #region 배열, 변수 선언 및 초기화
            // degree, radian 변환 인라인 함수 정의.
            double degree2radian(double angle) => ((angle * (PI)) / 180.0);

            double AOI_air = degree2radian(65.0);   // 입사각. (라디안) 
            Complex N_air = new Complex(1.0, 0);    // 공기의 굴절률.

            // 프레넬 반사계수를 담을 배열.
            Complex rp, rs;
            // 프레넬 투과계수를 담을 배열.
            Complex tp, ts;

            // 박막 두께 초기화.
            const double SiO2_thickness = 20.0,
                         SiN_thickness = 30.0;

            const int LastLayer = 200;     // 박막 층수. (pair = 층수 / 2)
            int LayerNum = LastLayer + 1;   // 박막의 층수. (1 : Si 기판)
            LenData = wavelength_Si.Length;

            Complex[,,] Sp = new Complex[LenData, 2, 2];       // scattering matrixs.
            Complex[,,] Ss = new Complex[LenData, 2, 2];
            Complex[,,] Ip = new Complex[LayerNum + 1, 2, 2];  // interface matrixs.
            Complex[,,] Is = new Complex[LayerNum + 1, 2, 2];
            Complex[,,] L = new Complex[LayerNum, 2, 2];      // layer matrixs.
                                                              // alpha, beta 이론값을 담을 배열 선언.
            double[] alpha_cal = new double[LenData],
                     beta_cal = new double[LenData];

            // Polarizer 오프셋 각.
            double polarizerAngle = degree2radian(45.0);

            Complex[] Rp = new Complex[LenData],
                      Rs = new Complex[LenData];
            #endregion
            ///////////////////////////////////////////
            // 100 번 반복 시 소요 시간 측정.
            /*
            int times = 100;
            long timeSum = 0;
            for (int time = 0; time < times; time++)
            {
                Stopwatch stopwatch = new Stopwatch();
                stopwatch.Start();
                ///////////////////////////////////////////
            */
            // 층 구조 : air -> SiO2 -> SiN -> ... -> SiO2 -> SiN -> Si
            for (int i = 0; i < LenData; i++)
            {
                Complex Sintheta_j, Costheta_j, Sintheta_k, theta_k, Costheta_k;

                Complex PhaseThickness;

                // 파장에 대한 물질의 복소굴절률을 구한다.
                Complex N_SiO2 = new Complex(n_SiO2[i], -k_SiO2[i]);
                Complex N_SiN = new Complex(n_SiN[i], -k_SiN[i]);
                Complex N_Si = new Complex(n_Si[i], -k_Si[i]);

                #region air -> SiO2 경계면에서의 반사계수, 투과계수, 위상두께를 계산한다.

                // air, SiO2 경계면에서의 굴절각을 구한다. (스넬의 법칙)
                Sintheta_j = new Complex(Sin((double)AOI_air), 0);
                Costheta_j = new Complex(Cos((double)AOI_air), 0);
                Sintheta_k = (N_air / N_SiO2) * Sintheta_j;
                // air, SiO2 경계면에서의 굴절각.
                theta_k = Complex.Asin(Sintheta_k);
                Costheta_k = Complex.Cos(theta_k);

                // air, SiO2 경계면에서의 반사계수를 구한다.
                rp = ((N_SiO2 * Costheta_j) - (N_air * Costheta_k)) /
                               ((N_SiO2 * Costheta_j) + (N_air * Costheta_k));

                rs = ((N_air * Costheta_j) - (N_SiO2 * Costheta_k)) /
                               ((N_air * Costheta_j) + (N_SiO2 * Costheta_k));

                // air, SiO2 경계면에서의 투과계수를 구한다.
                tp = (N_air * Costheta_j * 2) /
                               ((N_SiO2 * Costheta_j) + (N_air * Costheta_k));

                ts = (N_air * Costheta_j * 2) /
                               ((N_air * Costheta_j) + (N_SiO2 * Costheta_k));


                // 위상 두께를 구한다.
                PhaseThickness = (SiO2_thickness * Math.PI * 2) * N_SiO2 * Costheta_k /
                                          wavelength_SiO2[i];

                Complex PlusE = Complex.Exp(PhaseThickness * new Complex(0, 1.0));
                Complex MinusE = Complex.Exp(PhaseThickness * new Complex(0, -1.0));

                #endregion
                #region air -> SiO2 에서의 I, L 을 계산한다.
                // p-편광에 대한 interface matrix.
                Ip[0, 0, 0] = 1 / tp;
                Ip[0, 0, 1] = rp / tp;
                Ip[0, 1, 0] = rp / tp;
                Ip[0, 1, 1] = 1 / tp;

                // s-편광에 대한 interface matrix.
                Is[0, 0, 0] = 1 / ts;
                Is[0, 0, 1] = rs / ts;
                Is[0, 1, 0] = rs / ts;
                Is[0, 1, 1] = 1 / ts;

                // Layer matrix.
                L[0, 0, 0] = PlusE;
                L[0, 0, 1] = 0;
                L[0, 1, 0] = 0;
                L[0, 1, 1] = MinusE;


                #endregion

                #region air, Si 를 제외한 각 층마다 I, L 을 계산한다.
                for (int layer = 1; layer < LayerNum; layer++)
                {
                    switch (layer % 2)
                    {
                        // SiN -> SiO2.
                        case 0:
                            {
                                #region SiN -> SiO2 경계면에서의 반사계수, 투과계수, 위상두께를 구한다.

                                // SiN, SiO2 경계면에서의 굴절각을 구한다. (스넬의 법칙)
                                Sintheta_j = Complex.Sin(theta_k);
                                Costheta_j = Complex.Cos(theta_k);
                                Sintheta_k = (N_SiN / N_SiO2) * Sintheta_j;
                                theta_k = Complex.Asin(Sintheta_k);
                                Costheta_k = Complex.Cos(theta_k);

                                // SiN, SiO2 경계면에서의 반사계수를 구한다.
                                rp = ((N_SiO2 * Costheta_j) - (N_SiN * Costheta_k)) /
                                               ((N_SiO2 * Costheta_j) + (N_SiN * Costheta_k));

                                rs = ((N_SiN * Costheta_j) - (N_SiO2 * Costheta_k)) /
                                               ((N_SiN * Costheta_j) + (N_SiO2 * Costheta_k));

                                // SiN, SiO2 경계면에서의 투과계수를 구한다.
                                tp = (N_SiN * Costheta_j * 2) /
                                               ((N_SiO2 * Costheta_j) + (N_SiN * Costheta_k));

                                ts = (N_SiN * Costheta_j * 2) /
                                               ((N_SiN * Costheta_j) + (N_SiO2 * Costheta_k));

                                // SiO2 의 위상 두께를 구한다.
                                PhaseThickness = (SiO2_thickness * Math.PI * 2) * N_SiO2 * Costheta_k /
                                                          wavelength_SiO2[i];

                                PlusE = Complex.Exp(PhaseThickness * new Complex(0, 1.0));
                                MinusE = Complex.Exp(PhaseThickness * new Complex(0, -1.0));

                                #endregion
                                #region SiN -> SiO2 에서의 I, L 을 계산한다.
                                // p-편광에 대한 interface matrix.
                                Ip[layer, 0, 0] = 1 / tp;
                                Ip[layer, 0, 1] = rp / tp;
                                Ip[layer, 1, 0] = rp / tp;
                                Ip[layer, 1, 1] = 1 / tp;

                                // s-편광에 대한 interface matrix.
                                Is[layer, 0, 0] = 1 / ts;
                                Is[layer, 0, 1] = rs / ts;
                                Is[layer, 1, 0] = rs / ts;
                                Is[layer, 1, 1] = 1 / ts;

                                // Layer matrix.
                                L[layer, 0, 0] = PlusE;
                                L[layer, 0, 1] = 0;
                                L[layer, 1, 0] = 0;
                                L[layer, 1, 1] = MinusE;

                                #endregion
                            }
                            break;

                        // SiO2 -> SiN.
                        case 1:
                            {
                                #region SiO2 -> SiN 경계면에서의 반사계수, 투과계수, 위상두께를 구한다.

                                // SiO2, SiN 경계면에서의 굴절각을 구한다. (스넬의 법칙)
                                Sintheta_j = Complex.Sin(theta_k);
                                Costheta_j = Complex.Cos(theta_k);
                                Sintheta_k = (N_SiO2 / N_SiN) * Sintheta_j;
                                theta_k = Complex.Asin(Sintheta_k);
                                Costheta_k = Complex.Cos(theta_k);

                                // SiO2, SiN 경계면에서의 반사계수를 구한다.
                                rp = ((N_SiN * Costheta_j) - (N_SiO2 * Costheta_k)) /
                                               ((N_SiN * Costheta_j) + (N_SiO2 * Costheta_k));

                                rs = ((N_SiO2 * Costheta_j) - (N_SiN * Costheta_k)) /
                                               ((N_SiO2 * Costheta_j) + (N_SiN * Costheta_k));

                                // SiO2, SiN 경계면에서의 투과계수를 구한다.
                                tp = (N_SiO2 * Costheta_j * 2) /
                                               ((N_SiN * Costheta_j) + (N_SiO2 * Costheta_k));

                                ts = (N_SiO2 * Costheta_j * 2) /
                                               ((N_SiO2 * Costheta_j) + (N_SiN * Costheta_k));

                                // SiN 의 위상 두께를 구한다.
                                PhaseThickness = (SiN_thickness * Math.PI * 2) * N_SiN * Costheta_k /
                                                          wavelength_SiO2[i];

                                PlusE = Complex.Exp(PhaseThickness * new Complex(0, 1.0));
                                MinusE = Complex.Exp(PhaseThickness * new Complex(0, -1.0));

                                #endregion
                                #region SiO2 -> SiN 에서의 I, L 을 계산한다.
                                // p-편광에 대한 interface matrix.
                                Ip[layer, 0, 0] = 1 / tp;
                                Ip[layer, 0, 1] = rp / tp;
                                Ip[layer, 1, 0] = rp / tp;
                                Ip[layer, 1, 1] = 1 / tp;

                                // s-편광에 대한 interface matrix.
                                Is[layer, 0, 0] = 1 / ts;
                                Is[layer, 0, 1] = rs / ts;
                                Is[layer, 1, 0] = rs / ts;
                                Is[layer, 1, 1] = 1 / ts;

                                // Layer matrix.
                                L[layer, 0, 0] = PlusE;
                                L[layer, 0, 1] = 0;
                                L[layer, 1, 0] = 0;
                                L[layer, 1, 1] = MinusE;

                                #endregion
                            }
                            break;

                        default:
                            break;
                    }
                }
                #endregion

                #region SiN -> Si 경계면에서의 반사계수, 투과계수를 계산한다.

                // SiN, Si 경계면에서의 굴절각을 구한다. (스넬의 법칙)
                Sintheta_j = Complex.Sin(theta_k);
                Costheta_j = Complex.Cos(theta_k);
                Sintheta_k = (N_SiN / N_Si) * Sintheta_j;
                theta_k = Complex.Asin(Sintheta_k);
                Costheta_k = Complex.Cos(theta_k);

                // SiO2, SiN 경계면에서의 반사계수를 구한다.
                rp = ((N_Si * Costheta_j) - (N_SiN * Costheta_k)) /
                               ((N_Si * Costheta_j) + (N_SiN * Costheta_k));

                rs = ((N_SiN * Costheta_j) - (N_Si * Costheta_k)) /
                               ((N_SiN * Costheta_j) + (N_Si * Costheta_k));

                // SiO2, SiN 경계면에서의 투과계수를 구한다.
                tp = (N_SiN * Costheta_j * 2) /
                               ((N_Si * Costheta_j) + (N_SiN * Costheta_k));

                ts = (N_SiN * Costheta_j * 2) /
                               ((N_SiN * Costheta_j) + (N_Si * Costheta_k));

                #endregion
                #region SiN -> Si 에서의 I 를 계산한다.

                Ip[LayerNum, 0, 0] = 1 / tp;
                Ip[LayerNum, 0, 1] = rp / tp;
                Ip[LayerNum, 1, 0] = rp / tp;
                Ip[LayerNum, 1, 1] = 1 / tp;

                Is[LayerNum, 0, 0] = 1 / ts;
                Is[LayerNum, 0, 1] = rs / ts;
                Is[LayerNum, 1, 0] = rs / ts;
                Is[LayerNum, 1, 1] = 1 / ts;

                #endregion

                #region scattering matrix 를 계산한다.

                #region strassen algorithm

                Sp[i, 0, 0] = new Complex(1, 0);
                Sp[i, 1, 1] = new Complex(1, 0);
                Ss[i, 0, 0] = new Complex(1, 0);
                Ss[i, 1, 1] = new Complex(1, 0);
                Complex[] C = new Complex[10];
                Complex[] P = new Complex[7];
                //Complex[,] tempP = new Complex[2,2];
                //Complex[,] tempS = new Complex[2,2];
                //tempP[0, 0] = new Complex(1, 0);
                //tempP[1, 1] = new Complex(1, 0);
                //tempS[0, 0] = new Complex(1, 0);
                //tempS[1, 1] = new Complex(1, 0);
                for (int layer = 0; layer < LayerNum; layer++)
                {
                    //strassen algorithm
                    //I 곱하기(p편광)
                    C[0] = Ip[layer, 0, 1] - Ip[layer, 1, 1];
                    C[1] = Sp[i, 0, 0] + Sp[i, 0, 1];
                    C[2] = Sp[i, 1, 0] + Sp[i, 1, 1];
                    C[3] = Ip[layer, 1, 0] - Ip[layer, 0, 0];
                    C[4] = Sp[i, 0, 0] + Sp[i, 1, 1];
                    C[5] = Ip[layer, 0, 0] + Ip[layer, 1, 1];
                    C[6] = Sp[i, 0, 1] - Sp[i, 1, 1];
                    C[7] = Ip[layer, 1, 0] + Ip[layer, 1, 1];
                    C[8] = Sp[i, 0, 0] - Sp[i, 1, 0];
                    C[9] = Ip[layer, 0, 0] + Ip[layer, 0, 1];

                    P[0] = Sp[i, 0, 0] * C[0];
                    P[1] = C[1] * Ip[layer, 1, 1];
                    P[2] = C[2] * Ip[layer, 0, 0];
                    P[3] = Sp[i, 1, 1] * C[3];
                    P[4] = C[4] * C[5];
                    P[5] = C[6] * C[7];
                    P[6] = C[8] * C[9];

                    Sp[i, 0, 0] = P[4] + P[3] - P[1] + P[5];
                    Sp[i, 0, 1] = P[0] + P[1];
                    Sp[i, 1, 0] = P[2] + P[3];
                    Sp[i, 1, 1] = P[0] + P[4] - P[2] - P[6];

                    //L곱하기(p편광)
                    C[0] = L[layer, 0, 1] - L[layer, 1, 1];
                    C[1] = Sp[i, 0, 0] + Sp[i, 0, 1];
                    C[2] = Sp[i, 1, 0] + Sp[i, 1, 1];
                    C[3] = L[layer, 1, 0] - L[layer, 0, 0];
                    C[4] = Sp[i, 0, 0] + Sp[i, 1, 1];
                    C[5] = L[layer, 0, 0] + L[layer, 1, 1];
                    C[6] = Sp[i, 0, 1] - Sp[i, 1, 1];
                    C[7] = L[layer, 1, 0] + L[layer, 1, 1];
                    C[8] = Sp[i, 0, 0] - Sp[i, 1, 0];
                    C[9] = L[layer, 0, 0] + L[layer, 0, 1];
                    //
                    P[0] = Sp[i, 0, 0] * C[0];
                    P[1] = C[1] * L[layer, 1, 1];
                    P[2] = C[2] * L[layer, 0, 0];
                    P[3] = Sp[i, 1, 1] * C[3];
                    P[4] = C[4] * C[5];
                    P[5] = C[6] * C[7];
                    P[6] = C[8] * C[9];

                    Sp[i, 0, 0] = P[4] + P[3] - P[1] + P[5];
                    Sp[i, 0, 1] = P[0] + P[1];
                    Sp[i, 1, 0] = P[2] + P[3];
                    Sp[i, 1, 1] = P[0] + P[4] - P[2] - P[6];

                    //I 곱하기(S편광)
                    C[0] = Is[layer, 0, 1] - Is[layer, 1, 1];
                    C[1] = Ss[i, 0, 0] + Ss[i, 0, 1];
                    C[2] = Ss[i, 1, 0] + Ss[i, 1, 1];
                    C[3] = Is[layer, 1, 0] - Is[layer, 0, 0];
                    C[4] = Ss[i, 0, 0] + Ss[i, 1, 1];
                    C[5] = Is[layer, 0, 0] + Is[layer, 1, 1];
                    C[6] = Ss[i, 0, 1] - Ss[i, 1, 1];
                    C[7] = Is[layer, 1, 0] + Is[layer, 1, 1];
                    C[8] = Ss[i, 0, 0] - Ss[i, 1, 0];
                    C[9] = Is[layer, 0, 0] + Is[layer, 0, 1];
                    ////
                    P[0] = Ss[i, 0, 0] * C[0];
                    P[1] = C[1] * Is[layer, 1, 1];
                    P[2] = C[2] * Is[layer, 0, 0];
                    P[3] = Ss[i, 1, 1] * C[3];
                    P[4] = C[4] * C[5];
                    P[5] = C[6] * C[7];
                    P[6] = C[8] * C[9];

                    Ss[i, 0, 0] = P[4] + P[3] - P[1] + P[5];
                    Ss[i, 0, 1] = P[0] + P[1];
                    Ss[i, 1, 0] = P[2] + P[3];
                    Ss[i, 1, 1] = P[0] + P[4] - P[2] - P[6];


                    //L곱하기(P편광, S편광)
                    C[0] = L[layer, 0, 1] - L[layer, 1, 1];
                    C[1] = Ss[i, 0, 0] + Ss[i, 0, 1];
                    C[2] = Ss[i, 1, 0] + Ss[i, 1, 1];
                    C[3] = L[layer, 1, 0] - L[layer, 0, 0];
                    C[4] = Ss[i, 0, 0] + Ss[i, 1, 1];
                    C[5] = L[layer, 0, 0] + L[layer, 1, 1];
                    C[6] = Ss[i, 0, 1] - Ss[i, 1, 1];
                    C[7] = L[layer, 1, 0] + L[layer, 1, 1];
                    C[8] = Ss[i, 0, 0] - Ss[i, 1, 0];
                    C[9] = L[layer, 0, 0] + L[layer, 0, 1];
                    //
                    P[0] = Ss[i, 0, 0] * C[0];
                    P[1] = C[1] * L[layer, 1, 1];
                    P[2] = C[2] * L[layer, 0, 0];
                    P[3] = Ss[i, 1, 1] * C[3];
                    P[4] = C[4] * C[5];
                    P[5] = C[6] * C[7];
                    P[6] = C[8] * C[9];

                    Ss[i, 0, 0] = P[4] + P[3] - P[1] + P[5];
                    Ss[i, 0, 1] = P[0] + P[1];
                    Ss[i, 1, 0] = P[2] + P[3];
                    Ss[i, 1, 1] = P[0] + P[4] - P[2] - P[6];

                }
                //마지막 Interface matrix계산 

                //I 곱하기(P편광)
                C[0] = Ip[LayerNum, 0, 1] - Ip[LayerNum, 1, 1];
                C[1] = Sp[i, 0, 0] + Sp[i, 0, 1];
                C[2] = Sp[i, 1, 0] + Sp[i, 1, 1];
                C[3] = Ip[LayerNum, 1, 0] - Ip[LayerNum, 0, 0];
                C[4] = Sp[i, 0, 0] + Sp[i, 1, 1];
                C[5] = Ip[LayerNum, 0, 0] + Ip[LayerNum, 1, 1];
                C[6] = Sp[i, 0, 1] - Sp[i, 1, 1];
                C[7] = Ip[LayerNum, 1, 0] + Ip[LayerNum, 1, 1];
                C[8] = Sp[i, 0, 0] - Sp[i, 1, 0];
                C[9] = Ip[LayerNum, 0, 0] + Ip[LayerNum, 0, 1];

                P[0] = Sp[i, 0, 0] * C[0];
                P[1] = C[1] * Ip[LayerNum, 1, 1];
                P[2] = C[2] * Ip[LayerNum, 0, 0];
                P[3] = Sp[i, 1, 1] * C[3];
                P[4] = C[4] * C[5];
                P[5] = C[6] * C[7];
                P[6] = C[8] * C[9];

                Sp[i, 0, 0] = P[4] + P[3] - P[1] + P[5];
                Sp[i, 0, 1] = P[0] + P[1];
                Sp[i, 1, 0] = P[2] + P[3];
                Sp[i, 1, 1] = P[0] + P[4] - P[2] - P[6];

                //I 곱하기(S편광)
                C[0] = Is[LayerNum, 0, 1] - Is[LayerNum, 1, 1];
                C[1] = Ss[i, 0, 0] + Ss[i, 0, 1];
                C[2] = Ss[i, 1, 0] + Ss[i, 1, 1];
                C[3] = Is[LayerNum, 1, 0] - Is[LayerNum, 0, 0];
                C[4] = Ss[i, 0, 0] + Ss[i, 1, 1];
                C[5] = Is[LayerNum, 0, 0] + Is[LayerNum, 1, 1];
                C[6] = Ss[i, 0, 1] - Ss[i, 1, 1];
                C[7] = Is[LayerNum, 1, 0] + Is[LayerNum, 1, 1];
                C[8] = Ss[i, 0, 0] - Ss[i, 1, 0];
                C[9] = Is[LayerNum, 0, 0] + Is[LayerNum, 0, 1];

                P[0] = Ss[i, 0, 0] * C[0];
                P[1] = C[1] * Is[LayerNum, 1, 1];
                P[2] = C[2] * Is[LayerNum, 0, 0];
                P[3] = Ss[i, 1, 1] * C[3];
                P[4] = C[4] * C[5];
                P[5] = C[6] * C[7];
                P[6] = C[8] * C[9];

                Ss[i, 0, 0] = P[4] + P[3] - P[1] + P[5];
                Ss[i, 0, 1] = P[0] + P[1];
                Ss[i, 1, 0] = P[2] + P[3];
                Ss[i, 1, 1] = P[0] + P[4] - P[2] - P[6];

                #endregion

                #region 행렬 계산
                /*
                Complex[,] tempP = {
                    {new Complex(1, 0), new Complex(0, 0) },
                    {new Complex(0, 0), new Complex(1, 0) }
                };

                Complex[,] tempS = {
                    {new Complex(1, 0), new Complex(0, 0) },
                    {new Complex(0, 0), new Complex(1, 0) }
                };
                Complex tempP00, tempP01, tempP10, tempP11;
                Complex tempS00, tempS01, tempS10, tempS11;

                // 층의 개수만큼 I, L 을 곱해준다.
                for (int layer = 0; layer < LayerNum; layer++)
                {
                    tempP00 = tempP[0, 0];  // 1
                    tempP01 = tempP[0, 1];  // 0
                    tempP10 = tempP[1, 0];  // 0
                    tempP11 = tempP[1, 1];  // 1

                    // I 곱하기. (p편광)
                    tempP[0, 0] = tempP00 * Ip[layer, 0, 0] + tempP01 * Ip[layer, 1, 0];
                    tempP[0, 1] = tempP00 * Ip[layer, 0, 1] + tempP01 * Ip[layer, 1, 1];
                    tempP[1, 0] = tempP10 * Ip[layer, 0, 0] + tempP11 * Ip[layer, 1, 0];
                    tempP[1, 1] = tempP10 * Ip[layer, 0, 1] + tempP11 * Ip[layer, 1, 1];

                    tempP00 = tempP[0, 0];  // 1
                    tempP01 = tempP[0, 1];  // 0
                    tempP10 = tempP[1, 0];  // 0
                    tempP11 = tempP[1, 1];  // 1

                    // L 곱하기. (p편광)
                    tempP[0, 0] = tempP00 * L[layer, 0, 0] + tempP01 * L[layer, 1, 0];
                    tempP[0, 1] = tempP00 * L[layer, 0, 1] + tempP01 * L[layer, 1, 1];
                    tempP[1, 0] = tempP10 * L[layer, 0, 0] + tempP11 * L[layer, 1, 0];
                    tempP[1, 1] = tempP10 * L[layer, 0, 1] + tempP11 * L[layer, 1, 1];

                    // I 곱하기. (s편광)
                    tempS00 = tempS[0, 0];  // 1
                    tempS01 = tempS[0, 1];  // 0
                    tempS10 = tempS[1, 0];  // 0
                    tempS11 = tempS[1, 1];  // 1

                    tempS[0, 0] = tempS00 * Is[layer, 0, 0] + tempS01 * Is[layer, 1, 0];
                    tempS[0, 1] = tempS00 * Is[layer, 0, 1] + tempS01 * Is[layer, 1, 1];
                    tempS[1, 0] = tempS10 * Is[layer, 0, 0] + tempS11 * Is[layer, 1, 0];
                    tempS[1, 1] = tempS10 * Is[layer, 0, 1] + tempS11 * Is[layer, 1, 1];

                    // L 곱하기. (s편광)
                    tempS00 = tempS[0, 0];  // 1
                    tempS01 = tempS[0, 1];  // 0
                    tempS10 = tempS[1, 0];  // 0
                    tempS11 = tempS[1, 1];  // 1

                    tempS[0, 0] = tempS00 * L[layer, 0, 0] + tempS01 * L[layer, 1, 0];
                    tempS[0, 1] = tempS00 * L[layer, 0, 1] + tempS01 * L[layer, 1, 1];
                    tempS[1, 0] = tempS10 * L[layer, 0, 0] + tempS11 * L[layer, 1, 0];
                    tempS[1, 1] = tempS10 * L[layer, 0, 1] + tempS11 * L[layer, 1, 1];
                }
                // 마지막 interface matrix 를 곱해준다. (Si 기판에 관한 것)
                // I 곱하기. (p편광)
                tempP00 = tempP[0, 0];  // 1
                tempP01 = tempP[0, 1];  // 0
                tempP10 = tempP[1, 0];  // 0
                tempP11 = tempP[1, 1];  // 1

                tempP[0, 0] = tempP00 * Ip[LayerNum, 0, 0] + tempP01 * Ip[LayerNum, 1, 0];
                tempP[0, 1] = tempP00 * Ip[LayerNum, 0, 1] + tempP01 * Ip[LayerNum, 1, 1];
                tempP[1, 0] = tempP10 * Ip[LayerNum, 0, 0] + tempP11 * Ip[LayerNum, 1, 0];
                tempP[1, 1] = tempP10 * Ip[LayerNum, 0, 1] + tempP11 * Ip[LayerNum, 1, 1];

                // I 곱하기. (s편광)
                tempS00 = tempS[0, 0];  // 1
                tempS01 = tempS[0, 1];  // 0
                tempS10 = tempS[1, 0];  // 0
                tempS11 = tempS[1, 1];  // 1

                tempS[0, 0] = tempS00 * Is[LayerNum, 0, 0] + tempS01 * Is[LayerNum, 1, 0];
                tempS[0, 1] = tempS00 * Is[LayerNum, 0, 1] + tempS01 * Is[LayerNum, 1, 1];
                tempS[1, 0] = tempS10 * Is[LayerNum, 0, 0] + tempS11 * Is[LayerNum, 1, 0];
                tempS[1, 1] = tempS10 * Is[LayerNum, 0, 1] + tempS11 * Is[LayerNum, 1, 1];

                // 계산이 완료된 temp 행렬의 각 요소를 대응되는 scattering matrix 의 각 요소에 대입한다.
                Sp[i, 0, 0] = tempP[0, 0];
                Sp[i, 0, 1] = tempP[0, 1];
                Sp[i, 1, 0] = tempP[1, 0];
                Sp[i, 1, 1] = tempP[1, 1];

                Ss[i, 0, 0] = tempS[0, 0];
                Ss[i, 0, 1] = tempS[0, 1];
                Ss[i, 1, 0] = tempS[1, 0];
                Ss[i, 1, 1] = tempS[1, 1];
               */
                #endregion

                #endregion
                #region alpha, beta 를 계산한다.
                Rp[i] = Sp[i, 1, 0] / Sp[i, 0, 0];
                Rs[i] = Ss[i, 1, 0] / Ss[i, 0, 0];

                // 총 반사계수비. (복소반사계수비)
                Complex rho = Rp[i] / Rs[i];

                // Psi, Delta.
                double Psi = Atan(rho.Magnitude);
                double Delta = rho.Phase;


                alpha_cal[i] = (Pow(Tan(Psi), 2.0) - Pow(Tan(polarizerAngle), 2.0)) /
                                       (Pow(Tan(Psi), 2.0) + Pow(Tan(polarizerAngle), 2.0));

                beta_cal[i] = (2.0 * Tan(Psi) * Cos(Delta) * Tan(polarizerAngle)) /
                                       (Pow(Tan(Psi), 2.0) + Pow(Tan(polarizerAngle), 2.0));
                #endregion
            }
            /*
                /////////////////////////////////////////////
                // 소요시간 측정 범위 종료 지점.
                stopwatch.Stop();
                timeSum += stopwatch.ElapsedMilliseconds;
                //WriteLine($"소요 시간: {stopwatch.ElapsedMilliseconds}ms");
            }
            long avgTime = timeSum / times;
            //WriteLine($"avgTime: {avgTime}(ms)");
             /////////////////////////////////////////////
            */

            #region 출력
            for (int i = 0; i < LenData; i++)
            {
                WriteLine(
                    wavelength_Si[i] + "\t"
                    + alpha_cal[i] + "\t"
                    + beta_cal[i]);
            }
            #endregion
            #endregion

            #region 파일 쓰기.
            // 파일 쓰기.
            //using (StreamWriter NewSpectrumOutputFile = new StreamWriter("100pairs.dat"))
            //{
            //    // 컬럼 명 쓰기.
            //    NewSpectrumOutputFile.WriteLine(
            //        "wavelength(nm)"    + "\t"
            //        + "alpha"           + "\t"
            //        + "beta");    // 컬럼명 쓰기.
            //    // WriteLine(Columns);
            //
            //
            //    // 스펙트럼 데이터 쓰기.
            //    for (int i = 0; i < LenData; i++)
            //    {
            //        // tsv 데이터 형식으로 데이터를 쓴다.
            //        NewSpectrumOutputFile.WriteLine(
            //            wavelength_Si[i]    + "\t"
            //            + alpha_cal[i]      + "\t"
            //            + beta_cal[i]);
            //    }
            //}

            #endregion
        }
    }
}
