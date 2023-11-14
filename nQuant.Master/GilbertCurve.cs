using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;

/* Generalized Hilbert ("gilbert") space-filling curve for rectangular domains of arbitrary (non-power of two) sizes.
Copyright (c) 2021 - 2023 Miller Cy Chan
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
		private float[] weights;
		private readonly bool sortedByYDiff;
		private readonly int width;
		private readonly int height;
		private readonly int[] pixels;
		private readonly Color[] palette;
		private readonly int[] qPixels;
		private readonly Ditherable ditherable;
		private readonly float[] saliencies;
		private List<ErrorBox> errorq;
		private readonly int[] lookup;

		private readonly int margin, thresold;
		private const float BLOCK_SIZE = 343f;

		private GilbertCurve(int width, int height, int[] pixels, Color[] palette, int[] qPixels, Ditherable ditherable, float[] saliencies, double weight)
		{
			this.width = width;
			this.height = height;
			this.pixels = pixels;
			this.palette = palette;
			this.qPixels = qPixels;
			this.ditherable = ditherable;
			this.saliencies = saliencies;
			var hasAlpha = weight < 0;

			errorq = new();
			weight = Math.Abs(weight);
			margin = weight < .003 ? 12 : 6;
			sortedByYDiff = !hasAlpha && palette.Length >= 128 && weight >= .04;
			DITHER_MAX = (byte)(weight < .01 ? (weight > .0025) ? 25 : 16 : 9);
			var edge = hasAlpha ? 1 : Math.Exp(weight) + .25;
			ditherMax = (hasAlpha || DITHER_MAX > 9) ? (byte) BitmapUtilities.Sqr(Math.Sqrt(DITHER_MAX) + edge) : DITHER_MAX;
			if (palette.Length / weight > 5000 && (weight > .045 || (weight > .01 && palette.Length <= 64)))
				ditherMax = (byte) BitmapUtilities.Sqr(5 + edge);
			else if (palette.Length / weight < 3200 && palette.Length >= 16 && palette.Length < 256)
				ditherMax = (byte) BitmapUtilities.Sqr(5 + edge);
			thresold = DITHER_MAX > 9 ? -112 : -64;
			weights = new float[0];
			lookup = new int[65536];
		}

		private void DitherPixel(int x, int y)
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
			if (palette.Length <= 32 && a_pix > 0xF0)
			{
				int offset = ditherable.GetColorIndex(c2.ToArgb());
				if (lookup[offset] == 0)
					lookup[offset] = ditherable.DitherColorIndex(palette, c2.ToArgb(), bidx) + 1;
				qPixels[bidx] = lookup[offset] - 1;
				
				if(saliencies != null && CIELABConvertor.Y_Diff(pixel, c2) > palette.Length - margin) {
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

			var errLength = denoise ? error.Length - 1 : 0;
			for (int j = 0; j < errLength; ++j)
			{
				if (Math.Abs(error[j]) >= ditherMax)
				{
					if (diffuse)
						error[j] = (float)Math.Tanh(error[j] / maxErr * 20) * (ditherMax - 1);
					else if(illusion)
						error[j] = (float)(error[j] / maxErr * error.yDiff) * (ditherMax - 1);
					else
						error[j] /= (float)(1 + Math.Sqrt(ditherMax));
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
					DitherPixel(x, y);
					x += dax;
					y += day;
				}
				return;
			}

			if (w == 1) {
				for (int i = 0; i < h; ++i){
					DitherPixel(x, y);
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

		public static int[] Dither(int width, int height, int[] pixels, Color[] palette, Ditherable ditherable, float[] saliencies = null, double weight = 1.0)
		{
			var qPixels = new int[pixels.Length];
			new GilbertCurve(width, height, pixels, palette, qPixels, ditherable, saliencies, weight).Run();
			return qPixels;
		}
	}
}
