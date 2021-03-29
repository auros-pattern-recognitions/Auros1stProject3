using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;

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
            // "SiO2_new.txt", "SiN_new.txt", "Si_new.txt" 파일로부터
            // 파장별 복소굴절률(n, k)을 읽어온다.
            //
            // 2021.03.29 이지원.
            //
            #region SiO2, SiN, Si 물성을 읽어온다.

            // "SiO2_new.txt" 파일 읽기.
            string[] SiO2_new = File.ReadAllLines("SiO2_new.txt");  // SiO2 물성값 저장.(한 줄씩)
            string[] SingleLineData;                                // 한 줄의 스펙트럼 데이터를 임시로 저장할 배열.

            // 데이터의 첫번째 줄은 column 명이다.
            // 이를 제외하고 데이터를 받기 위해 LenData 변수를 선언한다.
            int LenData = SiO2_new.Length - 1;
            double[] wavelength_SiO2    = new double[LenData];
            double[] n_SiO2             = new double[LenData];
            double[] k_SiO2             = new double[LenData];

            // SiO2_new 에 받은 데이터를 각 컬럼별로 저장한다.
            int StartIndex = 1;
            LenData = SiO2_new.Length;
            for (int i = StartIndex; i < LenData; i++)
            {
                // tsv 형식의 데이터를 SingleLineData에 저장한다.
                SingleLineData = SiO2_new[i].Split((char)0x09);  // 0x09 : 수평 탭.

                // 각 컬럼에 해당하는 데이터를 저장한다.
                wavelength_SiO2[i - 1]  = Double.Parse(SingleLineData[0]);
                n_SiO2[i - 1]           = Double.Parse(SingleLineData[1]);
                k_SiO2[i - 1]           = Double.Parse(SingleLineData[2]);
            }


            // "SiN_new.txt" 파일 읽기.
            string[] SiN_new = File.ReadAllLines("SiN_new.txt");  // SiN 물성값 저장.(한 줄씩)

            // 데이터의 첫번째 줄은 column 명이다.
            // 이를 제외하고 데이터를 받기 위해 LenData 변수를 선언한다.
            LenData = SiN_new.Length - 1;
            double[] wavelength_SiN = new double[LenData];
            double[] n_SiN          = new double[LenData];
            double[] k_SiN          = new double[LenData];

            // SiN_new 에 받은 데이터를 각 컬럼별로 저장한다.
            LenData = SiN_new.Length;
            for (int i = StartIndex; i < LenData; i++)
            {
                // tsv 형식의 데이터를 SingleLineData에 저장한다.
                SingleLineData = SiN_new[i].Split((char)0x09);  // 0x09 : 수평 탭.

                // 각 컬럼에 해당하는 데이터를 저장한다.
                wavelength_SiN[i - 1]   = Double.Parse(SingleLineData[0]);
                n_SiN[i - 1]            = Double.Parse(SingleLineData[1]);
                k_SiN[i - 1]            = Double.Parse(SingleLineData[2]);
            }
            

            // "Si_new.txt" 파일 읽기.
            string[] Si_new = File.ReadAllLines("Si_new.txt");  // Si 기판 물성값 저장.(한 줄씩)

            // 데이터의 첫번째 줄은 column 명이다.
            // 이를 제외하고 데이터를 받기 위해 LenData 변수를 선언한다.
            LenData = Si_new.Length - 1;
            double[] wavelength_Si  = new double[LenData];
            double[] n_Si           = new double[LenData];
            double[] k_Si           = new double[LenData];

            // Si_new 에 받은 데이터를 각 컬럼별로 저장한다.
            LenData = Si_new.Length;
            for (int i = StartIndex; i < LenData; i++)
            {
                // tsv 형식의 데이터를 SingleLineData에 저장한다.
                SingleLineData = Si_new[i].Split((char)0x09);  // 0x09 : 수평 탭.

                // 각 컬럼에 해당하는 데이터를 저장한다.
                wavelength_Si[i - 1]    = Double.Parse(SingleLineData[0]);
                n_Si[i - 1]             = Double.Parse(SingleLineData[1]);
                k_Si[i - 1]             = Double.Parse(SingleLineData[2]);
            }
            #endregion

            #region scattering matrix 를 계산한다.

            int LayerNum = 100;   // 박막의 층수.

            Complex[,] S = new Complex[2, 2];   // scattering matrix.
            Complex[,] I = new Complex[2, 2];   // interface matrix.
            Complex[,] L = new Complex[2, 2];   // Layer matrix.

            #endregion
        }
    }
}
