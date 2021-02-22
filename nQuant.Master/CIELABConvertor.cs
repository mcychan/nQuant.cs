using System;
using System.Collections.Generic;
using System.Drawing;
using System.Text;

namespace PnnQuant
{
    internal class CIELABConvertor
    {
		internal struct Lab {
            internal double alpha, A, B, L;
	    }

        internal static Lab RGB2LAB(Color c1)
	    {
			float r = c1.R / 255.0f, g = c1.G / 255.0f, b = c1.B / 255.0f;
			float x, y, z;

		    r = (r > 0.04045) ? (float) Math.Pow((r + 0.055) / 1.055, 2.4) : r / 12.92f;
		    g = (g > 0.04045) ? (float) Math.Pow((g + 0.055) / 1.055, 2.4) : g / 12.92f;
		    b = (b > 0.04045) ? (float) Math.Pow((b + 0.055) / 1.055, 2.4) : b / 12.92f;

		    x = (float)(r * 0.4124 + g * 0.3576 + b * 0.1805) / 0.95047f;
		    y = (float)(r * 0.2126 + g * 0.7152 + b * 0.0722) / 1.00000f;
		    z = (float)(r * 0.0193 + g * 0.1192 + b * 0.9505) / 1.08883f;

		    x = (x > 0.008856) ? (float) Math.Cbrt(x) : (float)(7.787 * x) + 16.0f / 116.0f;
		    y = (y > 0.008856) ? (float) Math.Cbrt(y) : (float)(7.787 * y) + 16.0f / 116.0f;
		    z = (z > 0.008856) ? (float) Math.Cbrt(z) : (float)(7.787 * z) + 16.0f / 116.0f;

            Lab lab = new Lab
            {
                alpha = c1.A,
                L = (116 * y) - 16,
                A = 500 * (x - y),
                B = 200 * (y - z)
            };
            return lab;
	    }

        internal static Color LAB2RGB(Lab lab)
        {
			var y = (float)(lab.L + 16) / 116;
			var x = (float)lab.A / 500 + y;
			var z = y - (float)lab.B / 200;
			float r, g, b;

		    x = (float)(0.95047 * ((x * x * x > 0.008856) ? x * x * x : (x - 16.0 / 116.0) / 7.787));
		    y = (float)(1.00000 * ((y * y * y > 0.008856) ? y * y * y : (y - 16.0 / 116.0) / 7.787));
		    z = (float)(1.08883 * ((z * z * z > 0.008856) ? z * z * z : (z - 16.0 / 116.0) / 7.787));

		    r = x *  3.2406f + y * -1.5372f + z * -0.4986f;
		    g = x * -0.9689f + y *  1.8758f + z *  0.0415f;
		    b = x *  0.0557f + y * -0.2040f + z *  1.0570f;

		    r = (r > 0.0031308) ? (float)(1.055 * Math.Pow(r, 1.0 / 2.4) - 0.055) : 12.92f * r;
		    g = (g > 0.0031308) ? (float)(1.055 * Math.Pow(g, 1.0 / 2.4) - 0.055) : 12.92f * g;
		    b = (b > 0.0031308) ? (float)(1.055 * Math.Pow(b, 1.0 / 2.4) - 0.055) : 12.92f * b;

            var alpha = Math.Clamp((int) lab.alpha, Byte.MinValue, Byte.MaxValue);
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
    }
}
