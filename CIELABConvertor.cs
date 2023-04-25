using System;
using System.Drawing;

namespace nQuant.Master
{
	internal class CIELABConvertor
	{
		private static readonly double XYZ_WHITE_REFERENCE_X = 95.047;
		private static readonly double XYZ_WHITE_REFERENCE_Y = 100;
		private static readonly double XYZ_WHITE_REFERENCE_Z = 108.883;
		private static readonly double XYZ_EPSILON = 0.008856;
		private static readonly double XYZ_KAPPA = 903.3;

		internal struct Lab {
			internal double alpha, A, B, L;
		}

		private static float pivotXyzComponent(double component)
		{
			return component > XYZ_EPSILON
					? (float) Math.Cbrt(component)
					: (float)((XYZ_KAPPA * component + 16) / 116.0);
		}
		
		private static double gammaToLinear(int channel)
		{
			var c = channel / 255.0;
			return c < 0.04045 ? c / 12.92 : Math.Pow((c + 0.055) / 1.055, 2.4);
		}

		internal static Lab RGB2LAB(Color c1)
		{
			var sr = gammaToLinear(c1.R);
			var sg = gammaToLinear(c1.G);
			var sb = gammaToLinear(c1.B);
			var x = pivotXyzComponent(100 * (sr * 0.4124 + sg * 0.3576 + sb * 0.1805) / XYZ_WHITE_REFERENCE_X);
			var y = pivotXyzComponent(100 * (sr * 0.2126 + sg * 0.7152 + sb * 0.0722) / XYZ_WHITE_REFERENCE_Y);
			var z = pivotXyzComponent(100 * (sr * 0.0193 + sg * 0.1192 + sb * 0.9505) / XYZ_WHITE_REFERENCE_Z);

			Lab lab = new()
			{
				alpha = c1.A,
				L = Math.Max(0, 116 * y - 16),
				A = 500 * (x - y),
				B = 200 * (y - z)
			};
			return lab;
		}

		internal static Color LAB2RGB(Lab lab)
		{
			var fy = (lab.L + 16.0) / 116.0;
			var fx = lab.A / 500 + fy;
			var fz = fy - lab.B / 200.0;
			var tmp = fx * fx * fx;
			var xr = tmp > XYZ_EPSILON ? tmp : (116.0 * fx - 16) / XYZ_KAPPA;
			var yr = lab.L > XYZ_KAPPA * XYZ_EPSILON ? fy * fy * fy : lab.L / XYZ_KAPPA;
			tmp = fz * fz * fz;
			var zr = tmp > XYZ_EPSILON ? tmp : (116.0 * fz - 16) / XYZ_KAPPA;
			var x = xr * XYZ_WHITE_REFERENCE_X;
			var y = yr * XYZ_WHITE_REFERENCE_Y;
			var z = zr * XYZ_WHITE_REFERENCE_Z;

			var alpha = Math.Clamp((int)lab.alpha, Byte.MinValue, Byte.MaxValue);
			double r = (x * 3.2406 + y * -1.5372 + z * -0.4986) / 100.0;
			double g = (x * -0.9689 + y * 1.8758 + z * 0.0415) / 100.0;
			double b = (x * 0.0557 + y * -0.2040 + z * 1.0570) / 100.0;
			r = r > 0.0031308 ? 1.055 * Math.Pow(r, 1 / 2.4) - 0.055 : 12.92 * r;
			g = g > 0.0031308 ? 1.055 * Math.Pow(g, 1 / 2.4) - 0.055 : 12.92 * g;
			b = b > 0.0031308 ? 1.055 * Math.Pow(b, 1 / 2.4) - 0.055 : 12.92 * b;
			
			return Color.FromArgb(alpha, Math.Clamp((int)(r * Byte.MaxValue), Byte.MinValue, Byte.MaxValue), Math.Clamp((int)(g * Byte.MaxValue), Byte.MinValue, Byte.MaxValue), Math.Clamp((int)(b * Byte.MaxValue), Byte.MinValue, Byte.MaxValue));
		}

		/*******************************************************************************
		* Conversions.
		******************************************************************************/

		private static float Deg2Rad(float deg)
		{
			return (float) (deg * (Math.PI / 180.0));
		}

		internal static float L_prime_div_k_L_S_L(Lab lab1, Lab lab2)
		{
			var k_L = 1.0f;
			var deltaLPrime = lab2.L - lab1.L;	
			var barLPrime = (lab1.L + lab2.L) / 2.0;
			var S_L = 1 + ((0.015 * Math.Pow(barLPrime - 50.0, 2.0)) / Math.Sqrt(20 + Math.Pow(barLPrime - 50.0, 2.0)));
			return (float)(deltaLPrime / (k_L * S_L));
		}

		internal static float C_prime_div_k_L_S_L(Lab lab1, Lab lab2, out float a1Prime, out float a2Prime, out float CPrime1, out float CPrime2)
		{
			var k_C = 1.0f;
			var pow25To7 = 6103515625f; /* pow(25, 7) */
			var C1 = Math.Sqrt((lab1.A * lab1.A) + (lab1.B * lab1.B));
			var C2 = Math.Sqrt((lab2.A * lab2.A) + (lab2.B * lab2.B));
			var barC = (C1 + C2) / 2.0;
			var G = 0.5 * (1 - Math.Sqrt(Math.Pow(barC, 7) / (Math.Pow(barC, 7) + pow25To7)));
			a1Prime = (float)((1f + G) * lab1.A);
			a2Prime = (float)((1f + G) * lab2.A);

			CPrime1 = (float)(Math.Sqrt((a1Prime * a1Prime) + (lab1.B * lab1.B)));
			CPrime2 = (float)(Math.Sqrt((a2Prime * a2Prime) + (lab2.B * lab2.B)));
			var deltaCPrime = CPrime2 - CPrime1;
			var barCPrime = (CPrime1 + CPrime2) / 2.0;

			var S_C = 1f + (0.045 * barCPrime);
			return (float)(deltaCPrime / (k_C * S_C));
		}

		internal static float H_prime_div_k_L_S_L(Lab lab1, Lab lab2, float a1Prime, float a2Prime, float CPrime1, float CPrime2, out float barCPrime, out float barhPrime)
		{
			var k_H = 1.0;
			var deg360InRad = Deg2Rad(360f);
			var deg180InRad = Deg2Rad(180f);
			var CPrimeProduct = CPrime1 * CPrime2;
			float hPrime1;
			if (lab1.B == 0.0 && a1Prime == 0.0)
				hPrime1 = 0.0f;
			else {
				hPrime1 = (float)Math.Atan2(lab1.B, a1Prime);
				/*
				* This must be converted to a hue angle in degrees between 0
				* and 360 by addition of 2π to negative hue angles.
				*/
				if (hPrime1 < 0)
					hPrime1 += deg360InRad;
			}
			float hPrime2;
			if (lab2.B == 0.0 && a2Prime == 0.0)
				hPrime2 = 0.0f;
			else {
				hPrime2 = (float)Math.Atan2(lab2.B, a2Prime);
				/*
				* This must be converted to a hue angle in degrees between 0
				* and 360 by addition of 2π to negative hue angles.
				*/
				if (hPrime2 < 0)
					hPrime2 += deg360InRad;
			}
			float deltahPrime;
			if (CPrimeProduct == 0.0)
				deltahPrime = 0;
			else {
				/* Avoid the Math.abs() call */
				deltahPrime = hPrime2 - hPrime1;
				if (deltahPrime < -deg180InRad)
					deltahPrime += deg360InRad;
				else if (deltahPrime > deg180InRad)
					deltahPrime -= deg360InRad;
			}

			var deltaHPrime = 2f * Math.Sqrt(CPrimeProduct) * Math.Sin(deltahPrime / 2f);
			var hPrimeSum = hPrime1 + hPrime2;
			if ((CPrime1 * CPrime2) == 0.0) {
				barhPrime = hPrimeSum;
			}
			else {
				if (Math.Abs(hPrime1 - hPrime2) <= deg180InRad)
					barhPrime = (hPrimeSum / 2f);
				else {
					if (hPrimeSum < deg360InRad)
						barhPrime = (hPrimeSum + deg360InRad) / 2f;
					else
						barhPrime = (hPrimeSum - deg360InRad) / 2f;
				}
			}

			barCPrime = (CPrime1 + CPrime2) / 2f;
			var T = 1f - (0.17 * Math.Cos(barhPrime - Deg2Rad(30f))) +
				(0.24 * Math.Cos(2.0 * barhPrime)) +
				(0.32 * Math.Cos((3.0 * barhPrime) + Deg2Rad(6f))) -
				(0.20 * Math.Cos((4.0 * barhPrime) - Deg2Rad(63f)));
			var S_H = 1f + (0.015 * barCPrime * T);
			return (float)(deltaHPrime / (k_H * S_H));
		}

		internal static float R_T(float barCPrime, float barhPrime, float C_prime_div_k_L_S_L, float H_prime_div_k_L_S_L)
		{
			var pow25To7 = 6103515625f; /* Math.Pow(25, 7) */
			var deltaTheta = Deg2Rad(30f) * Math.Exp(-Math.Pow((barhPrime - Deg2Rad(275f)) / Deg2Rad(25f), 2.0));
			var R_C = 2.0 * Math.Sqrt(Math.Pow(barCPrime, 7.0) / (Math.Pow(barCPrime, 7.0) + pow25To7));
			var R_T = (-Math.Sin(2f * deltaTheta)) * R_C;
			return (float) (R_T * C_prime_div_k_L_S_L * H_prime_div_k_L_S_L);
		}

		/* From the paper "The CIEDE2000 Color-Difference Formula: Implementation Notes, */
		/* Supplementary Test Data, and Mathematical Observations", by */
		/* Gaurav Sharma, Wencheng Wu and Edul N. Dalal, */
		/* Color Res. Appl., vol. 30, no. 1, pp. 21-30, Feb. 2005. */
		/* Return the CIEDE2000 Delta E color difference measure squared, for two Lab values */
		internal static float CIEDE2000(Lab lab1, Lab lab2)
		{
			var deltaL_prime_div_k_L_S_L = L_prime_div_k_L_S_L(lab1, lab2);
			var deltaC_prime_div_k_L_S_L = C_prime_div_k_L_S_L(lab1, lab2, out float a1Prime, out float a2Prime, out float CPrime1, out float CPrime2);
			var deltaH_prime_div_k_L_S_L = H_prime_div_k_L_S_L(lab1, lab2, a1Prime, a2Prime, CPrime1, CPrime2, out float barCPrime, out float barhPrime);
			var deltaR_T = R_T(barCPrime, barhPrime, deltaC_prime_div_k_L_S_L, deltaH_prime_div_k_L_S_L);
			return
				(float) (Math.Pow(deltaL_prime_div_k_L_S_L, 2.0) +
				Math.Pow(deltaC_prime_div_k_L_S_L, 2.0) +
				Math.Pow(deltaH_prime_div_k_L_S_L, 2.0) +
				deltaR_T);
		}
		
		internal static double Y_Diff(Color c1, Color c2)
		{
			Func<Color, double> color2Y = c => {
				var sr = gammaToLinear(c.R);
				var sg = gammaToLinear(c.G);
				var sb = gammaToLinear(c.B);
				return sr * 0.2126 + sg * 0.7152 + sb * 0.0722;
			};
		
			var y = color2Y(c1);
			var y2 = color2Y(c2);
            return Math.Abs(y2 - y) / XYZ_WHITE_REFERENCE_Y;
		}
	}
}
