﻿using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Media;
using static System.Math;


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
            double[] wavelength_SiO2    = new double[LenData];
            double[] n_SiO2             = new double[LenData];
            double[] k_SiO2             = new double[LenData];

            // SiO2_new 에 받은 데이터를 각 컬럼별로 저장한다.
            int StartIndex = 1;
            LenData = SiO2_new2.Length;
            for (int i = StartIndex; i < LenData; i++)
            {
                // tsv 형식의 데이터를 SingleLineData에 저장한다.
                SingleLineData = SiO2_new2[i].Split((char)0x09);  // 0x09 : 수평 탭.

                // 각 컬럼에 해당하는 데이터를 저장한다.
                wavelength_SiO2[i - 1]  = Double.Parse(SingleLineData[0]);
                n_SiO2[i - 1]           = Double.Parse(SingleLineData[1]);
                k_SiO2[i - 1]           = Double.Parse(SingleLineData[2]);
            }


            // "SiN_new.txt" 파일 읽기.
            string[] SiN_new2 = File.ReadAllLines("SiN_new2.txt");  // SiN 물성값 저장.(한 줄씩)

            // 데이터의 첫번째 줄은 column 명이다.
            // 이를 제외하고 데이터를 받기 위해 LenData 변수를 선언한다.
            LenData = SiN_new2.Length - 1;
            double[] wavelength_SiN = new double[LenData];
            double[] n_SiN          = new double[LenData];
            double[] k_SiN          = new double[LenData];

            // SiN_new 에 받은 데이터를 각 컬럼별로 저장한다.
            LenData = SiN_new2.Length;
            for (int i = StartIndex; i < LenData; i++)
            {
                // tsv 형식의 데이터를 SingleLineData에 저장한다.
                SingleLineData = SiN_new2[i].Split((char)0x09);  // 0x09 : 수평 탭.

                // 각 컬럼에 해당하는 데이터를 저장한다.
                wavelength_SiN[i - 1]   = Double.Parse(SingleLineData[0]);
                n_SiN[i - 1]            = Double.Parse(SingleLineData[1]);
                k_SiN[i - 1]            = Double.Parse(SingleLineData[2]);
            }
            

            // "Si_new.txt" 파일 읽기.
            string[] Si_new2 = File.ReadAllLines("Si_new2.txt");  // Si 기판 물성값 저장.(한 줄씩)

            // 데이터의 첫번째 줄은 column 명이다.
            // 이를 제외하고 데이터를 받기 위해 LenData 변수를 선언한다.
            LenData = Si_new2.Length - 1;
            double[] wavelength_Si  = new double[LenData];
            double[] n_Si           = new double[LenData];
            double[] k_Si           = new double[LenData];

            // Si_new 에 받은 데이터를 각 컬럼별로 저장한다.
            LenData = Si_new2.Length;
            for (int i = StartIndex; i < LenData; i++)
            {
                // tsv 형식의 데이터를 SingleLineData에 저장한다.
                SingleLineData = Si_new2[i].Split((char)0x09);  // 0x09 : 수평 탭.

                // 각 컬럼에 해당하는 데이터를 저장한다.
                wavelength_Si[i - 1]    = Double.Parse(SingleLineData[0]);
                n_Si[i - 1]             = Double.Parse(SingleLineData[1]);
                k_Si[i - 1]             = Double.Parse(SingleLineData[2]);
            }
            #endregion

            //
            // 매질에 따른 interface matrix, layer matrix 를 구하여
            // scattering matrix 를 통해 alpha, beta 를 계산한다.
            //
            // 2021.03.30 이지원.
            //
            #region scattering matrix 를 계산한다.
            Complex[,,] S = new Complex[LenData, 2, 2];   
            Complex[,,] I = new Complex[LenData, 2, 2];   // interface matrix.
            Complex[,,] L = new Complex[LenData, 2, 2];   // layer matrix.


            // degree, radian 변환 인라인 함수 정의.
            double degree2radian(double angle) => ((angle * (PI)) / 180.0);

            double AOI_air = degree2radian(65.0);   // 입사각. (라디안) 
            Complex N_air = new Complex(1.0, 0);    // 공기의 굴절률.

            // 프레넬 반사계수를 담을 배열.
            Complex[] rp = new Complex[LenData],
                      rs = new Complex[LenData];
            // 프레넬 투과계수를 담을 배열.
            Complex[] tp = new Complex[LenData],
                      ts = new Complex[LenData];

            // 박막 두께 초기화.
            double SiO2_thickness   = 20.0,
                   SiN_thickness    = 30.0;

            const int LastLayer = 100;      // 박막 층수.
            int LayerNum = LastLayer + 1;   // 박막의 층수. (1 : Si 기판)
            // 층 구조 : air -> SiO2 -> SiN -> ... -> SiO2 -> SiN -> Si
            for (int layer = 0; layer < LayerNum; layer++)
            {
                // layer 에 따른 프레넬 반사, 투과계수를 구한다.
                switch (layer)
                {
                    case 0:
                        {
                            #region air -> SiO2 일 때 프레넬 반사계수, 투과계수를 구한다.
                            // 프레넬 반사계수, 투과계수를 구한다.
                            for (int i = 0; i < LenData; i++)
                            {
                                // 파장에 대한 물질의 복소굴절률을 구한다.
                                Complex N_SiO2 = new Complex(n_SiO2[i], -k_SiO2[i]);

                                // air, SiO2 경계면에서의 굴절각을 구한다. (스넬의 법칙)
                                Complex Sintheta_j = new Complex(Sin(AOI_air), 0);
                                Complex Costheta_j = new Complex(Cos(AOI_air), 0);
                                Complex Sintheta_k = (N_air / N_SiO2) * Sintheta_j;
                                Complex theta_k = Complex.Asin(Sintheta_k);
                                // air, SiO2 경계면에서의 굴절각.
                                Complex Costheta_k = Complex.Cos(theta_k);

                                // air, SiO2 경계면에서의 반사계수를 구한다.
                                rp[i] = ((N_SiO2 * Costheta_j) - (N_air * Costheta_k)) /
                                               ((N_SiO2 * Costheta_j) + (N_air * Costheta_k));

                                rs[i] = ((N_air * Costheta_j) - (N_SiO2 * Costheta_k)) /
                                               ((N_air * Costheta_j) + (N_SiO2 * Costheta_k));

                                // air, SiO2 경계면에서의 투과계수를 구한다.
                                tp[i] = (N_air * Costheta_j * 2.0) /
                                               ((N_SiO2 * Costheta_j) + (N_air * Costheta_k));

                                ts[i] = (N_air * Costheta_j * 2.0) /
                                               ((N_air * Costheta_j) + (N_SiO2 * Costheta_k));
                            }
                            #endregion
                        }
                        break;
                    
                    case LastLayer:
                        {
                            #region SiN -> Si

                            #endregion
                        }
                        break;

                    default:
                        {
                            #region SiO2 -> SiN
                            if (layer % 2 != 0)
                            {

                            }
                            #endregion
                            #region SiN-> SiO2
                            else
                            {

                            }
                            #endregion
                        }
                        break;
                }
                

            }

            #endregion
        }
    }
}
