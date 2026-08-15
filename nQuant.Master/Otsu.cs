/* Otsu's Image Segmentation Method
  Copyright (C) 2009 Tolga Birdal
  Copyright (c) 2018-2026 Miller Cy Chan
*/

using nQuant.Master;
using System;
using System.Collections.Generic;
using System.Drawing;
using System.Drawing.Imaging;
using System.Linq;

namespace OtsuThreshold
{
	public class Otsu : Ditherable
	{
		protected byte alphaThreshold = 0xF;
		protected bool hasSemiTransparency = false;

		protected int m_transparentPixelIndex = -1;
		protected Color m_transparentColor = Color.Transparent;

		protected float[] saliencies;
		internal Dictionary<int, CIELABConvertor.Lab> pixelMap = new();
		protected readonly Dictionary<int, ushort> nearestMap = new();

		// function is used to compute the q values in the equation
		private static float Px(int init, int end, int[] hist)
		{
			int sum = 0;
			for (int i = init; i <= end; ++i)
				sum += hist[i];

			return sum;
		}

		// function is used to compute the mean values in the equation (mu)
		private static float Mx(int init, int end, int[] hist)
		{
			int sum = 0;
			for (int i = init; i <= end; ++i)
				sum += i * hist[i];

			return sum;
		}

		// finds the maximum element in a vector
		private static short FindMax(float[] vec, int n)
		{
			float maxVec = 0;
			short idx = 0;

			for (short i = 1; i < n - 1; ++i)
			{
				if (vec[i] > maxVec)
				{
					maxVec = vec[i];
					idx = i;
				}
			}
			return idx;
		}

		// simply computes the image histogram
		private void GetHistogram(int[] pixels, int[] hist)
		{
			foreach(var pixel in pixels)
			{
				var c = Color.FromArgb(pixel);
				if (c.A <= alphaThreshold)
					continue;

				hist[c.R]++;
				hist[c.G]++;
				hist[c.B]++;
			}
		}

		private short GetOtsuThreshold(int[] pixels)
		{
			var vet = new float[256];
			var hist = new int[256];

			GetHistogram(pixels, hist);

			// loop through all possible t values and maximize between class variance
			for (int k = 1; k != Byte.MaxValue; ++k)
			{
				float p1 = Px(0, k, hist);
				float p2 = Px(k + 1, Byte.MaxValue, hist);
				float p12 = p1 * p2;
				if (p12 == 0)
					p12 = 1;
				float diff = (Mx(0, k, hist) * p2) - (Mx(k + 1, Byte.MaxValue, hist) * p1);
				vet[k] = diff * diff / p12;
			}

			return FindMax(vet, 256);
		}

		private void Threshold(int[] pixels, int[] dest, short thresh, float weight = 1f)
		{
			var maxThresh = (byte)thresh;
			if (thresh >= 200)
			{
				weight = .78f;
				maxThresh = (byte)(thresh * weight);
				thresh = 200;
			}

			var minThresh = (byte)(thresh * (m_transparentPixelIndex >= 0 ? .9f : weight));
			var shadow = m_transparentPixelIndex >= 0 ? 3.5 : 3;
			for (int i = 0; i < pixels.Length; ++i)
			{
				var c = Color.FromArgb(pixels[i]);
				if (c.A < alphaThreshold && c.R + c.G + c.B > maxThresh * 3)
					dest[i] = Color.FromArgb(c.A, Byte.MaxValue, Byte.MaxValue, Byte.MaxValue).ToArgb();
				else if (c.R + c.G + c.B < minThresh * shadow)
					dest[i] = Color.FromArgb(c.A, 0, 0, 0).ToArgb();
			}
		}

		private int[] CannyFilter(int width, int[] pixelsGray, double lowerThreshold, double higherThreshold, bool dither) {
			int height = pixelsGray.Length / width;
			int area = width * height;

			var pixelsCanny = Enumerable.Repeat(Color.White.ToArgb(), area).ToArray();

			var gx = new int[3, 3]{{-1, 0, 1}, {-2, 0, 2}, {-1, 0, 1}};
			var gy = new int[3, 3]{{-1, -2, -1}, {0, 0, 0}, {1, 2, 1}};
			var G = new double[area];
			var theta = new int[area];
			var largestG = 0.0;

			// Sobel Gradient Calculation
			for (int i = 1; i < height - 1; ++i) {
				for (int j = 1; j < width - 1; ++j) {
					// find gx and gy for each pixel
					var gxValue = 0.0;
					var gyValue = 0.0;
					for (int x = -1; x <= 1; ++x) {
						for (int y = -1; y <= 1; ++y) {
							var c = Color.FromArgb(pixelsGray[(i + x) * width +  j + y]);
							gxValue += gx[1 - x, 1 - y] * c.G;
							gyValue += gy[1 - x, 1 - y] * c.G;
						}
					}

					int center = i * width + j;
					// calculate G and theta
					G[center] = Math.Sqrt(Math.Pow(gxValue, 2) + Math.Pow(gyValue, 2));
					var atanResult = Math.Atan2(gyValue, gxValue) * 180.0 / Math.PI;
					theta[center] = (int)(180.0 + atanResult);

					if (G[center] > largestG)
						largestG = G[center];

					// to the nearest 45 degrees
					theta[center] = (int) Math.Round(theta[center] / 45.0) * 45;
				}
			}

			// Non-Maximum Suppression (NMS) -> Ensures 1-pixel thin line candidates
			var suppressedG = new double[area];
			for (int i = 1; i < height - 1; ++i) {
				for (int j = 1; j < width - 1; ++j) {
					int center = i * width + j;
					var currentG = G[center];

					if (theta[center] == 0 || theta[center] == 180 || theta[center] == 360) {
						if (currentG >= G[center - 1] && currentG >= G[center + 1]) suppressedG[center] = currentG;
					} else if (theta[center] == 45 || theta[center] == 225) {
						if (currentG >= G[center + width + 1] && currentG >= G[center - width - 1]) suppressedG[center] = currentG;
					} else if (theta[center] == 90 || theta[center] == 270) {
						if (currentG >= G[center + width] && currentG >= G[center - width]) suppressedG[center] = currentG;
					} else { // 135 or 315
						if (currentG >= G[center + width - 1] && currentG >= G[center - width + 1]) suppressedG[center] = currentG;
					}
				}
			}

			// Saliency-Driven Hysteresis Thresholding
			// A temporary tracking grid to mark confirmed clean edges (0 = background, 255 = edge)
			var edges = new int[area];

			for (int i = 1; i < height - 1; ++i) {
				for (int j = 1; j < width - 1; ++j) {
					int center = i * width + j;
					// Scale thresholds locally using the saliency value of the pixel
					// Highly salient areas get responsive sensitivity adjustment
					var saliency = (saliencies != null) ? saliencies[center] : 1.0f;
					var factor = 1.0 - (saliency * 0.6);
					if (factor < 0.15)
						factor = 0.15; // Ensures a 15% threshold floor always remains

					var localMin = lowerThreshold * largestG * factor; 
					var localMax = higherThreshold * largestG * factor;

					if (suppressedG[center] >= localMax) {
						edges[center] = 255; // Strong edge
					} else if (suppressedG[center] >= localMin) {
						// Hysteresis check: link weak edges to strong neighbors
						var connected = false;
						for (int x = -1; x <= 1 && !connected; ++x) {
							for (int y = -1; y <= 1; y++) {
								if (suppressedG[center + x * width + y] >= localMax) {
									connected = true;
									break;
								}
							}
						}
						if (connected) {
							edges[center] = 255;
						}
					}
				}
			}
			
			var dilatedEdges = new int[area];
			Array.Copy(edges, 0, dilatedEdges, 0, area);

			// Render to Out-Pixel Array
			for (int i = 0; i < area; ++i) {
				var c = Color.FromArgb(pixelsGray[i]);
				if (dilatedEdges[i] == 255) {
					// Draw clean black edge lines
					pixelsCanny[i] = Color.FromArgb(c.A, 0, 0, 0).ToArgb();
				} else {
					// Keep background white (or pass through original image depending on needs)
					pixelsCanny[i] = Color.FromArgb(c.A, 255, 255, 255).ToArgb();
				}
			}

			return pixelsCanny;
		}

		public ushort DitherColorIndex(Color[] palette, int pixel, int pos)
		{
			if (nearestMap.TryGetValue(pixel, out var k))
				return k;

			var c = Color.FromArgb(pixel);
			if (c.A <= alphaThreshold)
				return 0;

			double mindist = 1e100;
			for (int i = 0; i < palette.Length; ++i)
			{
				var c2 = palette[i];
				var curdist = BitmapUtilities.Sqr(c2.A - c.A);
				if (curdist > mindist)
					continue;

				curdist += BitmapUtilities.Sqr(c2.R - c.R);
				if (curdist > mindist)
					continue;

				curdist += BitmapUtilities.Sqr(c2.G - c.G);
				if (curdist > mindist)
					continue;

				curdist += BitmapUtilities.Sqr(c2.B - c.B);
				if (curdist > mindist)
					continue;

				mindist = curdist;
				k = (ushort)i;
			}
			nearestMap[pixel] = k;
			return k;
		}

		internal void GetLab(int argb, out CIELABConvertor.Lab lab1)
		{
			if (!pixelMap.TryGetValue(argb, out lab1))
			{
				lab1 = CIELABConvertor.RGB2LAB(Color.FromArgb(argb));
				pixelMap[argb] = lab1;
			}
		}

		public int GetColorIndex(int argb)
		{
			return BitmapUtilities.GetARGBIndex(argb, hasSemiTransparency, m_transparentPixelIndex > -1);
		}

		public Bitmap ConvertToGrayScale(Bitmap srcimg)
		{
			var iWidth = srcimg.Width;
			var iHeight = srcimg.Height;

			var pixelFormat = srcimg.PixelFormat;
			var bitDepth = Image.GetPixelFormatSize(srcimg.PixelFormat);
			if (bitDepth != 32 && bitDepth != 24)
				pixelFormat = PixelFormat.Format32bppArgb;

			var sourceImg = srcimg.Clone(new Rectangle(0, 0, iWidth, iHeight), pixelFormat);
			var data = sourceImg.LockBits(new Rectangle(0, 0, iWidth, iHeight), ImageLockMode.WriteOnly, sourceImg.PixelFormat);
			bitDepth = Image.GetPixelFormatSize(sourceImg.PixelFormat);
			var DJ = (byte)(bitDepth >> 3);

			unsafe
			{
				var ptr = (byte*)data.Scan0;

				float min1 = Byte.MaxValue;
				float max1 = .0f;
				int remain = data.Stride - iWidth * DJ;

				for (int i = 0; i < iHeight; ++i)
				{
					for (int j = 0; j < iWidth; ++j)
					{
						if (DJ > 3 && ptr[3] <= alphaThreshold) {
							ptr += DJ;
							continue;
						}
						
						if (min1 > ptr[1])
							min1 = ptr[1];

						if (max1 < ptr[1])
							max1 = ptr[1];
						ptr += DJ;
					}
					ptr += remain;
				}

				ptr = (byte*)data.Scan0;

				for (int i = 0; i < iHeight; ++i)
				{
					for (int j = 0; j < iWidth; ++j)
					{
						ptr[0] = ptr[1] = ptr[2] = (byte)((ptr[1] - min1) * (Byte.MaxValue / (max1 - min1)));
						ptr += DJ;
					}
					ptr += remain;
				}
			}

			sourceImg.UnlockBits(data);
			return sourceImg;
		}

		private void ConvertToGrayScale(int[] pixels, int[] dest)
		{
			float min1 = Byte.MaxValue;
			float max1 = .0f;

			foreach (var pixel in pixels)
			{
				int alfa = (pixel >> 24) & 0xff;
				if (alfa <= alphaThreshold)
					continue;

				int green = (pixel >> 8) & 0xff;
				if (min1 > green)
					min1 = green;

				if (max1 < green)
					max1 = green;
			}

			saliencies = new float[pixels.Length];
			var saliencyBase = .1f;
			for (int i = 0; i < pixels.Length; ++i)
			{
				var pixel = pixels[i];
				GetLab(pixel, out var lab1);
				saliencies[i] = (float) (saliencyBase + (1 - saliencyBase) * lab1.L / 100f * lab1.alpha / 255f);

				int alfa = (pixel >> 24) & 0xff;
				if (alfa <= alphaThreshold)
					continue;

				int green = (pixel >> 8) & 0xff;
				var grey = (int)((green - min1) * (Byte.MaxValue / (max1 - min1)));
				dest[i] = Color.FromArgb(alfa, grey, grey, grey).ToArgb();
			}
		}
		

		public Bitmap ConvertGrayScaleToBinary(Bitmap srcimg, bool isGrayscale = false, bool dither = true)
		{
			int bitmapWidth = srcimg.Width;
			int bitmapHeight = srcimg.Height;

			var pixels = new int[bitmapWidth * bitmapHeight];
			int semiTransCount = 0;
			if (!BitmapUtilities.GrabPixels(srcimg, pixels, ref semiTransCount, ref m_transparentColor, ref m_transparentPixelIndex, alphaThreshold))
				return srcimg;
			hasSemiTransparency = semiTransCount > 0;

			var pixelsGray = (int[]) pixels.Clone();
			if (!isGrayscale)
				ConvertToGrayScale(pixels, pixelsGray);

			var otsuThreshold = GetOtsuThreshold(pixelsGray);
			double lowerThreshold = .03, higherThreshold = .1;
			if(!dither) {
				lowerThreshold = otsuThreshold / 3.0;
				higherThreshold = otsuThreshold;
			}
			pixels = CannyFilter(bitmapWidth, pixelsGray, lowerThreshold, higherThreshold, dither);
			Threshold(pixelsGray, pixels, otsuThreshold);

			var dest = new Bitmap(bitmapWidth, bitmapHeight, PixelFormat.Format1bppIndexed);
			var palettes = dest.Palette.Entries;
			if (m_transparentPixelIndex >= 0)
			{
				palettes[0] = m_transparentColor;
				palettes[1] = Color.Black;
			}
			else
			{
				palettes[0] = Color.Black;
				palettes[1] = Color.White;
			}

			var qPixels = GilbertCurve.Dither(bitmapWidth, bitmapHeight, pixels, palettes, this);
			if (m_transparentPixelIndex >= 0)
			{
				var k = qPixels[m_transparentPixelIndex];
				if (palettes[k] != m_transparentColor)
					BitmapUtilities.Swap(ref palettes[0], ref palettes[1]);
			}

			nearestMap.Clear();
			return BitmapUtilities.ProcessImagePixels(dest, palettes, qPixels, m_transparentPixelIndex >= 0);
		}

	}
}
