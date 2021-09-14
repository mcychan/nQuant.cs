/* Otsu's Image Segmentation Method
  Copyright (C) 2009 Tolga Birdal
  Copyright (c) 2018 Miller Cy Chan
*/

using nQuant.Master;
using System;
using System.Collections.Generic;
using System.Drawing;
using System.Drawing.Imaging;

namespace OtsuThreshold
{
	public class Otsu
	{
		protected byte alphaThreshold = 0;
		protected bool hasSemiTransparency = false;

		protected int m_transparentPixelIndex = -1;
		protected Color m_transparentColor = Color.Transparent;
		protected readonly Dictionary<int, ushort> nearestMap = new Dictionary<int, ushort>();

		// function is used to compute the q values in the equation
		private static float Px(int init, int end, int[] hist)
		{
			int sum = 0;
			int i;
			for (i = init; i <= end; i++)
				sum += hist[i];

			return (float)sum;
		}

		// function is used to compute the mean values in the equation (mu)
		private static float Mx(int init, int end, int[] hist)
		{
			int sum = 0;
			int i;
			for (i = init; i <= end; i++)
				sum += i * hist[i];

			return (float)sum;
		}

		// finds the maximum element in a vector
		private static short findMax(float[] vec, int n)
		{
			float maxVec = 0;
			short idx = 0;
			short i;

			for (i = 1; i < n - 1; i++)
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
		private static unsafe void getHistogram(byte* p, int w, int h, int ws, byte DJ, int[] hist)
		{
			for (uint i = 0; i < h; i++)
			{
				for (uint j = 0; j < w * DJ; j += DJ)
				{
					var index = i * ws + j;
					hist[p[index]]++;
				}
			}
		}

		private static short getOtsuThreshold(Bitmap bitmap)
		{
			var bitmapWidth = bitmap.Width;
			var bitmapHeight = bitmap.Height;
			var bitDepth = Image.GetPixelFormatSize(bitmap.PixelFormat);
			var DJ = (byte)(bitDepth >> 3);

			var vet = new float[256];
			var hist = new int[256];

			var data = bitmap.LockBits(new Rectangle(0, 0, bitmapWidth, bitmapHeight), ImageLockMode.ReadOnly, bitmap.PixelFormat);

			unsafe
			{
				var pRowSource = (byte*)data.Scan0;
				getHistogram(pRowSource, bitmapWidth, bitmapHeight, data.Stride, DJ, hist);

				// loop through all possible t values and maximize between class variance
				for (int k = 1; k != Byte.MaxValue; k++)
				{
					float p1 = Px(0, k, hist);
					float p2 = Px(k + 1, Byte.MaxValue, hist);
					float p12 = p1 * p2;
					if (p12 == 0)
						p12 = 1;
					float diff = (Mx(0, k, hist) * p2) - (Mx(k + 1, Byte.MaxValue, hist) * p1);
					vet[k] = diff * diff / p12;
				}

				bitmap.UnlockBits(data);
			}
			return findMax(vet, 256);
		}

		private static bool threshold(Bitmap bitmap, short thresh)
		{
			var bitmapWidth = bitmap.Width;
			var bitmapHeight = bitmap.Height;

			var pixelFormat = thresh < 200 ? PixelFormat.Format24bppRgb : PixelFormat.Format32bppArgb;
			var data = bitmap.LockBits(new Rectangle(0, 0, bitmapWidth, bitmapHeight), ImageLockMode.ReadOnly, pixelFormat);

			var bitDepth = Image.GetPixelFormatSize(pixelFormat);
			var DJ = (byte)(bitDepth >> 3);

			unsafe
			{
				var pRowDest = (byte*)data.Scan0;

				for (int i = 0; i < bitmapHeight; ++i)
				{
					var ptr = &pRowDest[i * data.Stride];
					for (int j = 0; j < bitmapWidth * DJ; j += DJ)
					{
						ptr[j] = (byte)((ptr[j] > (byte)thresh) ? Byte.MaxValue : 0);
						ptr[j + 1] = (byte)((ptr[j + 1] > (byte)thresh) ? Byte.MaxValue : 0);
						ptr[j + 2] = (byte)((ptr[j + 2] > (byte)thresh) ? Byte.MaxValue : 0);
					}
				}
			}

			bitmap.UnlockBits(data);
			return true;
		}

		protected ushort NearestColorIndex(Color[] palette, int nMaxColors, int pixel)
		{
			if (nearestMap.TryGetValue(pixel, out var k))
				return k;

			var c = Color.FromArgb(pixel);
			if (c.A <= alphaThreshold)
				return 0;

			double mindist = 1e100;
			for (int i = 0; i < nMaxColors; ++i)
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

		protected bool GrabPixels(Bitmap source, int[] pixels)
		{
			return BitmapUtilities.GrabPixels(source, pixels, ref hasSemiTransparency, ref m_transparentColor, ref m_transparentPixelIndex);
		}

		public static Bitmap ConvertToGrayScale(Bitmap srcimg)
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
						byte grey = Math.Min(ptr[0], ptr[1]);
						grey = Math.Min(ptr[1], ptr[2]);
						if (min1 > grey)
							min1 = grey;

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
						byte grey = Math.Min(ptr[0], ptr[1]);
						grey = Math.Min(ptr[1], ptr[2]);
						ptr[0] = ptr[1] = ptr[2] = (byte)((grey - min1) * (Byte.MaxValue / (max1 - min1)));
						ptr += DJ;
					}
					ptr += remain;
				}
			}

			sourceImg.UnlockBits(data);
			return sourceImg;
		}
		protected int GetColorIndex(int argb)
		{
			return BitmapUtilities.GetARGBIndex(argb, hasSemiTransparency, m_transparentPixelIndex > -1);
		}

		public Bitmap ConvertGrayScaleToBinary(Bitmap srcimg)
		{
			var sourceImg = ConvertToGrayScale(srcimg);
			var otsuThreshold = getOtsuThreshold(sourceImg);

			if (!threshold(sourceImg, otsuThreshold))
				return srcimg;

			int bitmapWidth = sourceImg.Width;
			int bitmapHeight = sourceImg.Height;

			var pixels = new int[bitmapWidth * bitmapHeight];
			if (!GrabPixels(sourceImg, pixels))
				return sourceImg;

			var dest = new Bitmap(bitmapWidth, bitmapHeight, PixelFormat.Format1bppIndexed);
			var palettes = dest.Palette.Entries;
			if (m_transparentPixelIndex >= 0)
			{
				palettes[0] = Color.Transparent;
				palettes[1] = Color.Black;
			}
			else
			{
				palettes[0] = Color.Black;
				palettes[1] = Color.White;
			}

			var qPixels = GilbertCurve.Dither(bitmapWidth, bitmapHeight, pixels, palettes, NearestColorIndex, GetColorIndex);

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
