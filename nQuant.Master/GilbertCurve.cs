using System;
using System.Collections.Generic;
using System.Drawing;

/* Generalized Hilbert ("gilbert") space-filling curve for rectangular domains of arbitrary (non-power of two) sizes.
Copyright (c) 2021 - 2026 Miller Cy Chan
* A general rectangle with a known orientation is split into three regions ("up", "right", "down"), for which the function calls itself recursively, until a trivial path can be produced. */

namespace nQuant.Master
{
	class GilbertCurve
	{
		internal sealed class ErrorBox
		{
			internal double yDiff { get; set; }

			private readonly float[] p;
			internal ErrorBox()
			{
				p = new float[4];
			}

			internal ErrorBox(Color c)
			{
				p = new float[] {
					c.R,
					c.G,
					c.B,
					c.A
				};
			}

			internal float this[int i]
			{
				get { return p[i]; }
				set { p[i] = value; }
			}

			internal int Length
			{
				get { return p.Length; }
			}
		}

		private byte ditherMax, DITHER_MAX;
		private float beta;
		private float[] weights;
		private readonly bool dither, m_hasAlpha, sortedByYDiff;
		private readonly int width, height;
		private readonly double weight;
		private readonly int[] pixels;
		private readonly Color[] palette;
		private readonly int[] qPixels;
		private readonly Ditherable ditherable;
		private readonly float[] saliencies;
		private List<ErrorBox> errorq;

		private readonly int margin, thresold;
		private const float BLOCK_SIZE = 343f;

		private GilbertCurve(int width, int height, int[] pixels, Color[] palette, int[] qPixels, Ditherable ditherable, float[] saliencies, double weight, bool dither)
		{
			this.width = width;
			this.height = height;
			this.pixels = pixels;
			this.palette = palette;
			this.qPixels = qPixels;
			this.ditherable = ditherable;
			this.saliencies = saliencies;
			this.dither = dither;
			this.m_hasAlpha = weight < 0;

			this.weight = Math.Abs(weight);
			margin = weight < .0025 ? 12 : weight < .004 ? 8 : 6;
			sortedByYDiff = saliencies != null && palette.Length >= 128 && weight >= .02 && (!m_hasAlpha || weight < .18);
			errorq = new();
			beta = palette.Length > 4 ? (float) (.6f - .00625f * palette.Length) : 1;
			if (palette.Length > 4) {
				var boundary = .005 - .0000625 * palette.Length;
				beta = (float) (weight > boundary ? Math.Max(.25, beta - palette.Length * weight) : Math.Min(1.5, beta + palette.Length * weight));
				if (palette.Length > 16 && palette.Length <= 32 && weight < .003)
					beta += .075f;
				else if (weight < .0015 || (palette.Length > 32 && palette.Length < 256))
					beta += .1f;
				if (palette.Length >= 64 && (weight > .012 && weight < .0125) || (weight > .025 && weight < .03))
					beta += .05f;
				else if (palette.Length > 32 && palette.Length < 64 && weight < .015)
					beta = .55f;
				else if (palette.Length > 16 && palette.Length <= 32 && weight < .005)
					beta += (float)(.05 + weight * palette.Length);
            }
			else
				beta *= .95f;

			if (palette.Length > 64 || (palette.Length > 4 && weight > .02))
				beta *= .4f;
			if (palette.Length > 64 && weight < .02)
				beta = .18f;

			DITHER_MAX = (byte)(weight < .015 ? (weight > .0025) ? 25 : 16 : 9);
			if (weight > .99)
			{
				beta = (float)weight;
				DITHER_MAX = 25;
			}

			var edge = m_hasAlpha ? 1 : Math.Exp(weight) + .25;
			var deviation = !m_hasAlpha && weight > .002 ? .25 : 1;
			ditherMax = (m_hasAlpha || DITHER_MAX > 9) ? (byte) BitmapUtilities.Sqr(Math.Sqrt(DITHER_MAX) + edge * deviation) : (byte)(DITHER_MAX * 1.5);
			int density = palette.Length > 16 ? 3200 : 1500;
			if (palette.Length / weight > 5000 && (weight > .045 || (weight > .01 && palette.Length < 64)))
				ditherMax = (byte) BitmapUtilities.Sqr(5 + edge);
			else if (weight < .03 && palette.Length / weight < density && palette.Length >= 16 && palette.Length < 256)
				ditherMax = (byte) BitmapUtilities.Sqr(5 + edge);
			thresold = DITHER_MAX > 9 ? -112 : -64;
			weights = new float[0];
		}


		private static float NormalDistribution(float x, float peak)
		{
			const float mean = .5f, stdDev = .1f;

			// Calculate the probability density function (PDF)
			double exponent = -Math.Pow(x - mean, 2) / (2 * Math.Pow(stdDev, 2));
			double pdf = (1 / (stdDev * Math.Sqrt(2 * Math.PI))) * Math.Exp(exponent);
			double maxPdf = 1 / (stdDev * Math.Sqrt(2 * Math.PI)); // Peak at x = mean
			double scaledPdf = (pdf / maxPdf) * peak;
			return (float) Math.Max(0.0, Math.Min(peak, scaledPdf));
		}


		private int DitherPixel(int x, int y, Color c2, float beta)
		{
			int bidx = x + y * width;
			Color pixel = Color.FromArgb(pixels[bidx]);
			int r_pix = c2.R;
			int g_pix = c2.G;
			int b_pix = c2.B;
			int a_pix = c2.A;

			var qPixelIndex = qPixels[bidx];
			var strength = 1 / 3f;
			int acceptedDiff = Math.Max(2, palette.Length - margin);
			if (palette.Length <= 4 && saliencies[bidx] > .2f && saliencies[bidx] < .25f)
				c2 = BlueNoise.Diffuse(pixel, palette[qPixelIndex], beta * 2 / saliencies[bidx], strength, x, y);
			else if (palette.Length <= 4 || CIELABConvertor.Y_Diff(pixel, c2) < (2 * acceptedDiff)) {
				if (palette.Length > 64)
				{
					var kappa = saliencies[bidx] < .6f ? beta * .15f / saliencies[bidx] : beta * .4f / saliencies[bidx];
					c2 = BlueNoise.Diffuse(pixel, palette[qPixelIndex], kappa, strength, x, y);
				}
				else if(palette.Length > 16 && weight < .005)
					c2 = BlueNoise.Diffuse(pixel, palette[qPixelIndex], beta * NormalDistribution(saliencies[bidx], .5f) + beta, strength, x, y);
				else
					c2 = BlueNoise.Diffuse(pixel, palette[qPixelIndex], beta * .5f / saliencies[bidx], strength, x, y);
			}

			var gamma = (palette.Length <= 32 && weight < .01 && weight > .007) ? 1 - beta : beta;
			if (palette.Length > 4 && CIELABConvertor.Y_Diff(pixel, c2) > (gamma * acceptedDiff))
			{
				if (margin > 6 || gamma > beta)
				{
					var kappa = saliencies[bidx] < .4f ? beta * .4f * saliencies[bidx] : beta * .4f / saliencies[bidx];
					var c1 = Color.FromArgb(a_pix, r_pix, g_pix, b_pix);
					if (palette.Length > 32)
						kappa = beta * NormalDistribution(saliencies[bidx], 2f);
					else
					{
						if (weight >= .0015 && saliencies[bidx] < .6)
							c1 = pixel;
						if (weight < .005 && saliencies[bidx] < .6)
							kappa = beta * NormalDistribution(saliencies[bidx], weight < .0008 ? 2.5f : 1.75f);
						else if (palette.Length >= 32 || CIELABConvertor.Y_Diff(c1, c2) > (gamma * Math.PI * acceptedDiff))
						{
							var ub = 1 - palette.Length / 320.0;
							if (saliencies[bidx] > .15 && saliencies[bidx] < ub)
								kappa = beta * (!sortedByYDiff && weight < .0025 ? .55f : .5f) / saliencies[bidx];
							else
								kappa = beta * NormalDistribution(saliencies[bidx], weight < .0025 ? 1.82f : 2f);
						}
					}

					c2 = BlueNoise.Diffuse(c1, palette[qPixelIndex], kappa, strength, x, y);
				}
				else if (palette.Length <= 32 && weight >= .004)
					c2 = BlueNoise.Diffuse(c2, palette[qPixelIndex], beta * NormalDistribution(saliencies[bidx], .25f), strength, x, y);
				else
					c2 = Color.FromArgb(a_pix, r_pix, g_pix, b_pix);
			}

			if (DITHER_MAX < 16 && palette.Length > 4 && saliencies[bidx] < .6f && CIELABConvertor.Y_Diff(pixel, c2) > margin - 1)
				c2 = Color.FromArgb(a_pix, r_pix, g_pix, b_pix);
			if (palette.Length > 32 && saliencies[bidx] > .95)
			{
				var kappa = beta * Math.Max(.05f, .75f - palette.Length / 128f) * saliencies[bidx];
				c2 = BlueNoise.Diffuse(pixel, palette[qPixelIndex], kappa, strength, x, y);
			}

			return ditherable.DitherColorIndex(palette, c2.ToArgb(), bidx);
		}

		private void DiffusePixel(int x, int y)
		{
			int bidx = x + y * width;
			Color pixel = Color.FromArgb(pixels[bidx]);
			var error = new ErrorBox(pixel);
			int i = sortedByYDiff ? weights.Length - 1 : 0;
			float maxErr = DITHER_MAX - 1;
			foreach (var eb in errorq)
			{
				if(i < 0 || i >= weights.Length)
					break;

				for (int j = 0; j < eb.Length; ++j)
				{
					error[j] += eb[j] * weights[i];
					if(error[j] > maxErr)
						maxErr = error[j];
				}
				i += sortedByYDiff ? -1 : 1;
			}

			int r_pix = (int)Math.Min(Byte.MaxValue, Math.Max(error[0], 0.0));
			int g_pix = (int)Math.Min(Byte.MaxValue, Math.Max(error[1], 0.0));
			int b_pix = (int)Math.Min(Byte.MaxValue, Math.Max(error[2], 0.0));
			int a_pix = (int)Math.Min(Byte.MaxValue, Math.Max(error[3], 0.0));

			Color c2 = Color.FromArgb(a_pix, r_pix, g_pix, b_pix);
			if (saliencies != null && dither && !sortedByYDiff && (!m_hasAlpha || pixel.A < a_pix))
			{
				if ((palette.Length >= 256 && saliencies[bidx] > .99f) || (m_hasAlpha && (pixel.A - a_pix) < (.5 * margin)))
                    qPixels[bidx] = ditherable.DitherColorIndex(palette, c2.ToArgb(), bidx);
				else
					qPixels[bidx] = DitherPixel(x, y, c2, beta);
			}
			else if (palette.Length <= 32 && a_pix > 0xF0)
			{
				qPixels[bidx] = ditherable.DitherColorIndex(palette, c2.ToArgb(), bidx);

				int acceptedDiff = Math.Max(2, palette.Length - margin);
				if (saliencies != null && (CIELABConvertor.Y_Diff(pixel, c2) > acceptedDiff || CIELABConvertor.U_Diff(pixel, c2) > (2 * acceptedDiff)))
				{
					var strength = 1 / 3f;
					c2 = BlueNoise.Diffuse(pixel, palette[qPixels[bidx]], 1 / saliencies[bidx], strength, x, y);
					qPixels[bidx] = ditherable.DitherColorIndex(palette, c2.ToArgb(), bidx);
				}
			}
			else
				qPixels[bidx] = ditherable.DitherColorIndex(palette, c2.ToArgb(), bidx);

			if(errorq.Count >= DITHER_MAX)
				errorq.RemoveAt(0);
			else if (errorq.Count > 0)
				InitWeights(errorq.Count);

			c2 = palette[qPixels[bidx]];
			if (palette.Length > 256)
				qPixels[bidx] = (short)ditherable.GetColorIndex(c2.ToArgb());

			error[0] = r_pix - c2.R;
			error[1] = g_pix - c2.G;
			error[2] = b_pix - c2.B;
			error[3] = a_pix - c2.A;

			var denoise = palette.Length > 2;
			var diffuse = BlueNoise.TELL_BLUE_NOISE[bidx & 4095] > thresold;
			error.yDiff = sortedByYDiff ? CIELABConvertor.Y_Diff(pixel, c2) : 1;
			var illusion = !diffuse && BlueNoise.TELL_BLUE_NOISE[(int)(error.yDiff * 4096) & 4095] > thresold;

			var unaccepted = false;
			var errLength = denoise ? error.Length - 1 : 0;
			for (int j = 0; j < errLength; ++j)
			{
				if (Math.Abs(error[j]) >= ditherMax)
				{
					if (sortedByYDiff && saliencies != null)
						unaccepted = true;

					if (diffuse)
						error[j] = (float)Math.Tanh(error[j] / maxErr * 20) * (ditherMax - 1);
					else if(illusion)
						error[j] = (float)(error[j] / maxErr * error.yDiff) * (ditherMax - 1);
					else
						error[j] /= (float)(1 + Math.Sqrt(ditherMax));
				}

				if (sortedByYDiff && saliencies == null && Math.Abs(error[j]) >= DITHER_MAX)
					unaccepted = true;
			}

			if (unaccepted) {
				if (saliencies != null)
					qPixels[bidx] = DitherPixel(x, y, c2, beta);
				else if (CIELABConvertor.Y_Diff(pixel, c2) > 3 && CIELABConvertor.U_Diff(pixel, c2) > 3)
					qPixels[bidx] = DitherPixel(x, y, c2, 1.25f);

				if (palette.Length > 256) {
					c2 = palette[qPixels[bidx]];
					qPixels[bidx] = (short)ditherable.GetColorIndex(c2.ToArgb());
				}
			}

			errorq.Add(error);
			if (sortedByYDiff)
				errorq.Sort((o1, o2) => o2.yDiff.CompareTo(o1.yDiff));
		}

		private void Generate2d(int x, int y, int ax, int ay, int bx, int by) {
			int w = Math.Abs(ax + ay);
			int h = Math.Abs(bx + by);
			int dax = Math.Sign(ax);
			int day = Math.Sign(ay);
			int dbx = Math.Sign(bx);
			int dby = Math.Sign(by);

			if (h == 1) {
				for (int i = 0; i < w; ++i){
					DiffusePixel(x, y);
					x += dax;
					y += day;
				}
				return;
			}

			if (w == 1) {
				for (int i = 0; i < h; ++i){
					DiffusePixel(x, y);
					x += dbx;
					y += dby;
				}
				return;
			}

			int ax2 = ax / 2;
			int ay2 = ay / 2;
			int bx2 = bx / 2;
			int by2 = by / 2;

			int w2 = Math.Abs(ax2 + ay2);
			int h2 = Math.Abs(bx2 + by2);

			if (2 * w > 3 * h) {
				if ((w2 % 2) != 0 && w > 2) {
					ax2 += dax;
					ay2 += day;
				}
				Generate2d(x, y, ax2, ay2, bx, by);
				Generate2d(x + ax2, y + ay2, ax - ax2, ay - ay2, bx, by);
				return;
			}

			if ((h2 % 2) != 0 && h > 2) {
				bx2 += dbx;
				by2 += dby;
			}

			Generate2d(x, y, bx2, by2, ax2, ay2);
			Generate2d(x + bx2, y + by2, ax, ay, bx - bx2, by - by2);
			Generate2d(x + (ax - dax) + (bx2 - dbx), y + (ay - day) + (by2 - dby), -bx2, -by2, -(ax - ax2), -(ay - ay2));
		}

		private void InitWeights(int size) {
			/* Dithers all pixels of the image in sequence using
			 * the Gilbert path, and distributes the error in
			 * a sequence of pixels size.
			 */
			float weightRatio = (float) Math.Pow(BLOCK_SIZE + 1f, 1f / (size - 1f));
			float weight = 1f, sumweight = 0f;
			weights = new float[size];
			for (int c = 0; c < size; ++c)
			{
				errorq.Add(new ErrorBox());
				sumweight += (weights[size - c - 1] = 1.0f / weight);
				weight *= weightRatio;
			}
			if (sortedByYDiff)
				errorq.Sort((o1, o2) => o2.yDiff.CompareTo(o1.yDiff));

			weight = 0f; /* Normalize */
			for (int c = 0; c < size; ++c)
				weight += (weights[c] /= sumweight);
			weights[0] += 1f - weight;
		}

		private void Run()
		{
			if(!sortedByYDiff)
				InitWeights(DITHER_MAX);
			
			if (width >= height)
				Generate2d(0, 0, width, 0, 0, height);
			else
				Generate2d(0, 0, 0, height, width, 0);
		}

		public static int[] Dither(int width, int height, int[] pixels, Color[] palette, Ditherable ditherable, float[] saliencies = null, double weight = 1.0, bool dither = true)
		{
			var qPixels = new int[pixels.Length];
			new GilbertCurve(width, height, pixels, palette, qPixels, ditherable, saliencies, weight, dither).Run();
			return qPixels;
		}
	}
}
