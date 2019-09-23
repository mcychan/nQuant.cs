using System;
using System.Collections.Generic;
using System.Drawing;
using System.Text;

namespace PnnQuant
{
    internal class CIELABConvertor
    {
        internal struct Lab {
            internal float alpha, A, B, L;
	    }

        internal static Lab RGB2LAB(Color c1)
	    {
		    double r = c1.R / 255.0, g = c1.G / 255.0, b = c1.B / 255.0;
		    double x, y, z;

		    r = (r > 0.04045) ? Math.Pow((r + 0.055) / 1.055, 2.4) : r / 12.92;
		    g = (g > 0.04045) ? Math.Pow((g + 0.055) / 1.055, 2.4) : g / 12.92;
		    b = (b > 0.04045) ? Math.Pow((b + 0.055) / 1.055, 2.4) : b / 12.92;

		    x = (r * 0.4124 + g * 0.3576 + b * 0.1805) / 0.95047;
		    y = (r * 0.2126 + g * 0.7152 + b * 0.0722) / 1.00000;
		    z = (r * 0.0193 + g * 0.1192 + b * 0.9505) / 1.08883;

		    x = (x > 0.008856) ? Math.Pow(x, (1.0 / 3.0)) : (7.787 * x) + 16.0 / 116.0;
		    y = (y > 0.008856) ? Math.Pow(y, (1.0 / 3.0)) : (7.787 * y) + 16.0 / 116.0;
		    z = (z > 0.008856) ? Math.Pow(z, (1.0 / 3.0)) : (7.787 * z) + 16.0 / 116.0;

            Lab lab = new Lab
            {
                alpha = c1.A,
                L = (float)((116 * y) - 16),
                A = (float)(500 * (x - y)),
                B = (float)(200 * (y - z))
            };
            return lab;
	    }

        internal static Color LAB2RGB(Lab lab)
        {
		    double y = (lab.L + 16) / 116;
		    double x = lab.A / 500 + y;
		    double z = y - lab.B / 200;
		    double r, g, b;

		    x = 0.95047 * ((x * x * x > 0.008856) ? x * x * x : (x - 16.0 / 116.0) / 7.787);
		    y = 1.00000 * ((y * y * y > 0.008856) ? y * y * y : (y - 16.0 / 116.0) / 7.787);
		    z = 1.08883 * ((z * z * z > 0.008856) ? z * z * z : (z - 16.0 / 116.0) / 7.787);

		    r = x *  3.2406 + y * -1.5372 + z * -0.4986;
		    g = x * -0.9689 + y *  1.8758 + z *  0.0415;
		    b = x *  0.0557 + y * -0.2040 + z *  1.0570;

		    r = (r > 0.0031308) ? (1.055 * Math.Pow(r, 1.0 / 2.4) - 0.055) : 12.92 * r;
		    g = (g > 0.0031308) ? (1.055 * Math.Pow(g, 1.0 / 2.4) - 0.055) : 12.92 * g;
		    b = (b > 0.0031308) ? (1.055 * Math.Pow(b, 1.0 / 2.4) - 0.055) : 12.92 * b;

            var alpha = Math.Clamp((int) lab.alpha, Byte.MinValue, Byte.MaxValue);
            return Color.FromArgb(alpha, (int)(Math.Max(0, Math.Min(1, r)) * Byte.MaxValue), (int)(Math.Max(0, Math.Min(1, g)) * Byte.MaxValue), (int)(Math.Max(0, Math.Min(1, b)) * Byte.MaxValue));
	    }

	    /*******************************************************************************
	    * Conversions.
	    ******************************************************************************/

	    private static double deg2Rad(double deg)
	    {
		    return (deg * (Math.PI / 180.0));
	    }

	    internal static double L_prime_div_k_L_S_L(Lab lab1, Lab lab2)
	    {
		    const double k_L = 1.0;
		    double deltaLPrime = lab2.L - lab1.L;	
		    double barLPrime = (lab1.L + lab2.L) / 2.0;
		    double S_L = 1 + ((0.015 * Math.Pow(barLPrime - 50.0, 2.0)) / Math.Sqrt(20 + Math.Pow(barLPrime - 50.0, 2.0)));
		    return deltaLPrime / (k_L * S_L);
	    }

        internal static double C_prime_div_k_L_S_L(Lab lab1, Lab lab2, out double a1Prime, out double a2Prime, out double CPrime1, out double CPrime2)
	    {
		    const double k_C = 1.0;
		    const double pow25To7 = 6103515625.0; /* pow(25, 7) */
		    double C1 = Math.Sqrt((lab1.A * lab1.A) + (lab1.B * lab1.B));
		    double C2 = Math.Sqrt((lab2.A * lab2.A) + (lab2.B * lab2.B));
		    double barC = (C1 + C2) / 2.0;
		    double G = 0.5 * (1 - Math.Sqrt(Math.Pow(barC, 7) / (Math.Pow(barC, 7) + pow25To7)));
		    a1Prime = (1.0 + G) * lab1.A;
		    a2Prime = (1.0 + G) * lab2.A;

		    CPrime1 = Math.Sqrt((a1Prime * a1Prime) + (lab1.B * lab1.B));
		    CPrime2 = Math.Sqrt((a2Prime * a2Prime) + (lab2.B * lab2.B));
		    double deltaCPrime = CPrime2 - CPrime1;
		    double barCPrime = (CPrime1 + CPrime2) / 2.0;
		
		    double S_C = 1 + (0.045 * barCPrime);
		    return deltaCPrime / (k_C * S_C);
	    }

        internal static double H_prime_div_k_L_S_L(Lab lab1, Lab lab2, double a1Prime, double a2Prime, double CPrime1, double CPrime2, out double barCPrime, out double barhPrime)
	    {
		    const double k_H = 1.0;
		    double deg360InRad = deg2Rad(360.0);
		    double deg180InRad = deg2Rad(180.0);
		    double CPrimeProduct = CPrime1 * CPrime2;
		    double hPrime1;
		    if (lab1.B == 0.0 && a1Prime == 0.0)
			    hPrime1 = 0.0;
		    else {
			    hPrime1 = Math.Atan2(lab1.B, a1Prime);
			    /*
			    * This must be converted to a hue angle in degrees between 0
			    * and 360 by addition of 2􏰏 to negative hue angles.
			    */
			    if (hPrime1 < 0)
				    hPrime1 += deg360InRad;
		    }
		    double hPrime2;
		    if (lab2.B == 0.0 && a2Prime == 0.0)
			    hPrime2 = 0.0;
		    else {
			    hPrime2 = Math.Atan2(lab2.B, a2Prime);
			    /*
			    * This must be converted to a hue angle in degrees between 0
			    * and 360 by addition of 2􏰏 to negative hue angles.
			    */
			    if (hPrime2 < 0)
				    hPrime2 += deg360InRad;
		    }
		    double deltahPrime;
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

		    double deltaHPrime = 2.0 * Math.Sqrt(CPrimeProduct) * Math.Sin(deltahPrime / 2.0);
		    double hPrimeSum = hPrime1 + hPrime2;
		    if ((CPrime1 * CPrime2) == 0.0) {
			    barhPrime = hPrimeSum;
		    }
		    else {
			    if (Math.Abs(hPrime1 - hPrime2) <= deg180InRad)
				    barhPrime = (hPrimeSum / 2.0);
			    else {
				    if (hPrimeSum < deg360InRad)
					    barhPrime = (hPrimeSum + deg360InRad) / 2.0;
				    else
					    barhPrime = (hPrimeSum - deg360InRad) / 2.0;
			    }
		    }

		    barCPrime = (CPrime1 + CPrime2) / 2.0;
		    double T = 1.0 - (0.17 * Math.Cos(barhPrime - deg2Rad(30.0))) +
			    (0.24 * Math.Cos(2.0 * barhPrime)) +
			    (0.32 * Math.Cos((3.0 * barhPrime) + deg2Rad(6.0))) -
			    (0.20 * Math.Cos((4.0 * barhPrime) - deg2Rad(63.0)));
		    double S_H = 1 + (0.015 * barCPrime * T);
		    return deltaHPrime / (k_H * S_H);
	    }

        internal static double R_T(double barCPrime, double barhPrime, double C_prime_div_k_L_S_L, double H_prime_div_k_L_S_L)
	    {
		    const double pow25To7 = 6103515625.0; /* Math.Pow(25, 7) */
		    double deltaTheta = deg2Rad(30.0) * Math.Exp(-Math.Pow((barhPrime - deg2Rad(275.0)) / deg2Rad(25.0), 2.0));
		    double R_C = 2.0 * Math.Sqrt(Math.Pow(barCPrime, 7.0) / (Math.Pow(barCPrime, 7.0) + pow25To7));
		    double R_T = (-Math.Sin(2.0 * deltaTheta)) * R_C;
		    return R_T * C_prime_div_k_L_S_L * H_prime_div_k_L_S_L;
	    }

	    /* From the paper "The CIEDE2000 Color-Difference Formula: Implementation Notes, */
	    /* Supplementary Test Data, and Mathematical Observations", by */
	    /* Gaurav Sharma, Wencheng Wu and Edul N. Dalal, */
	    /* Color Res. Appl., vol. 30, no. 1, pp. 21-30, Feb. 2005. */
	    /* Return the CIEDE2000 Delta E color difference measure squared, for two Lab values */
        internal static double CIEDE2000(Lab lab1, Lab lab2)
	    {
		    double deltaL_prime_div_k_L_S_L = L_prime_div_k_L_S_L(lab1, lab2);
		    double a1Prime, a2Prime, CPrime1, CPrime2;
		    double deltaC_prime_div_k_L_S_L = C_prime_div_k_L_S_L(lab1, lab2, out a1Prime, out a2Prime, out CPrime1, out CPrime2);
		    double barCPrime, barhPrime;
            double deltaH_prime_div_k_L_S_L = H_prime_div_k_L_S_L(lab1, lab2, a1Prime, a2Prime, CPrime1, CPrime2, out barCPrime, out barhPrime);
		    double deltaR_T = R_T(barCPrime, barhPrime, deltaC_prime_div_k_L_S_L, deltaH_prime_div_k_L_S_L);
		    return
			    Math.Pow(deltaL_prime_div_k_L_S_L, 2.0) +
			    Math.Pow(deltaC_prime_div_k_L_S_L, 2.0) +
			    Math.Pow(deltaH_prime_div_k_L_S_L, 2.0) +
			    deltaR_T;
	    }
    }
}
