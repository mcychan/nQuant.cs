using System;
using System.Collections.Generic;
using System.Drawing;
using System.Drawing.Imaging;
using System.Numerics.Tensors;
using System.Threading.Tasks;

/* Graph regularized harmonic mean vector quantization algorithm with CIELAB color space
Copyright (c) 2026 Miller Cy Chan */

namespace nQuant.Master
{
	public class GRHVQuantizer : Ditherable
	{
		private int _k;
		private readonly float _gamma;
		private readonly float _alpha;
		private readonly float _w;

		protected byte alphaThreshold = 0xF;
		protected bool dither = true, hasSemiTransparency = false;
		protected int m_transparentPixelIndex = -1;
		protected Color m_transparentColor = Color.Transparent;
		protected Color[] m_palette;
		protected readonly Random rand = new();
		protected readonly Dictionary<int, ushort[]> closestMap = new();
		protected readonly Dictionary<int, ushort> nearestMap = new();

		protected double PR = 0.299, PG = 0.587, PB = 0.114, PA = .3333;
		protected double mDivn = 1;
		protected static readonly float[,] coeffs = new float[,] {
			{0.299f, 0.587f, 0.114f},
			{-0.14713f, -0.28886f, 0.436f},
			{0.615f, -0.51499f, -0.10001f}
		};
		
		protected float[] saliencies;
		private Dictionary<int, CIELABConvertor.Lab> pixelMap = new();
		private static readonly double TRANS_RATE = 1 - (512 + 101) / 768.0;

		public GRHVQuantizer(float alpha = 0.02f, float w = 0.5f, float gamma = 1.0f)
		{
			_alpha = alpha;
			_w = Math.Min(1.0f, Math.Max(0.0f, w));
			_gamma = gamma;
		}

		internal bool HasAlpha
		{
			get => m_transparentPixelIndex >= 0;
		}

		public int GetColorIndex(int argb)
		{
			return BitmapUtilities.GetARGBIndex(argb, hasSemiTransparency, HasAlpha);
		}

		internal void GetLab(int argb, out CIELABConvertor.Lab lab1)
		{
			if (!pixelMap.TryGetValue(argb, out lab1))
			{
				lab1 = CIELABConvertor.RGB2LAB(Color.FromArgb(argb));
				pixelMap[argb] = lab1;
			}
		}

		private static float DistanceSquared(ReadOnlySpan<float> a, ReadOnlySpan<float> b)
		{
			return (float) Math.Pow(TensorPrimitives.Distance(a, b), 2);
		}
		
		/// <summary>
		/// Builds a 4-neighborhood spatial graph Laplacian L = D - A using SIMD-accelerated TensorPrimitives.
		/// </summary>
		internal float[] BuildGridLaplacian(float[] X, int width, int height, int N)
		{
			// L dimensions: N x N
			var L = new float[N * N]; // Column-major flat array

			double sigma = 0.1;
			double invTwoSigmaSq = -1.0 / (2.0 * sigma * sigma);

			// 4-Neighbor directional offsets
			ReadOnlySpan<int> dx = stackalloc int[4] { 0, 0, -1, 1 };
			ReadOnlySpan<int> dy = stackalloc int[4] { -1, 1, 0, 0 };

			int D = HasAlpha ? 4 : 3;

			for (int y = 0; y < height; ++y)
			{
				for (int x = 0; x < width; ++x)
				{
					int u = y * width + x;

					var uLab = X.AsSpan(u * D, D);
					var degreeSum = 0.0f;

					for (int i = 0; i < 4; ++i)
					{
						int nx = x + dx[i];
						int ny = y + dy[i];

						if (nx >= 0 && nx < width && ny >= 0 && ny < height)
						{
							int v = ny * width + nx;

							var vLab = X.AsSpan(v * D, D);

							var colorDistSq = DistanceSquared(uLab, vLab);
							var w = (float) Math.Exp(colorDistSq * invTwoSigmaSq);

							// Column-major indexing: L[u, v] -> u + v * N
							L[u + v * N] = -w;
							degreeSum += w;
						}
					}

					// Diagonal degree matrix D[u, u]
					L[u + u * N] = degreeSum;
				}
			}

			return L;
		}

		protected float[] Optimize(float[] xData, float[] lData, int n, int d)
		{
			// 1. Matrix multiplications (LX, L2X)
			var lxData = MultiplyMatricesColumnMajor(lData, n, n, xData, n, d);
			var l2xData = MultiplyMatricesColumnMajor(lData, n, n, lxData, n, d);

			// Convert F_cont to ROW-MAJOR layout for optimal SIMD cache performance: [N x D]
			var fDataRowMajor = new float[n * d];
			var alphaSq = _alpha * _alpha;

			for (int i = 0; i < n; ++i)
			{
				for (int col = 0; col < d; ++col)
				{
					// Source is col-major: i + col * n -> Dest row-major: i * d + col
					var xVal = xData[i + col * n];
					var lxVal = lxData[i + col * n];
					var l2xVal = l2xData[i + col * n];
					fDataRowMajor[i * d + col] = xVal - (_alpha * lxVal) + (alphaSq * l2xVal);
				}
			}

			// Row-Major Adjacency buffer A = D - L and static total affinities
			var aDataRowMajor = new double[n * n];
			var totalAffinities = new double[n];

			for (int i = 0; i < n; ++i)
			{
				var sumAffinity = 0.0;
				for (int j = 0; j < n; ++j)
				{
					if (i != j)
					{
						var w = Math.Max(0.0, -lData[i + j * n]);
						aDataRowMajor[i * n + j] = w;
						sumAffinity += w;
					}
				}
				totalAffinities[i] = sumAffinity;
			}

			// 2. Centroid Initialization (Row-Major Codebook [K x D])
			var cData = new float[_k * d];
			var labels = new int[n];
			var rand = new Random(42);

			// Copy first centroid
			int firstIdx = rand.Next(n);
			fDataRowMajor.AsSpan(firstIdx * d, d).CopyTo(cData.AsSpan(0, d));

			var minDistances = new double[n];
			Array.Fill(minDistances, double.MaxValue);

			for (int c = 1; c < _k; ++c)
			{
				var totalDist = 0.0;
				var prevCentroid = cData.AsSpan((c - 1) * d, d);

				for (int i = 0; i < n; ++i)
				{
					var fRow = fDataRowMajor.AsSpan(i * d, d);
					double dist = DistanceSquared(fRow, prevCentroid);
					if (dist < minDistances[i])
						minDistances[i] = dist;
					totalDist += minDistances[i];
				}

				var target = rand.NextDouble() * totalDist;
				var currentSum = 0.0;
				int chosenIdx = n - 1;

				for (int i = 0; i < n; ++i)
				{
					currentSum += minDistances[i];
					if (currentSum >= target) { chosenIdx = i; break; }
				}

				fDataRowMajor.AsSpan(chosenIdx * d, d).CopyTo(cData.AsSpan(c * d, d));
			}

			// Precompute squared norms of F
			var fNorms = new double[n];
			// Initial label assignments
			for (int i = 0; i < n; ++i)
			{
				var fRow = fDataRowMajor.AsSpan(i * d, d);
				int bestCluster = 0;
				var minSqDist = double.MaxValue;

				for (int g = 0; g < _k; ++g)
				{
					var sqDist = DistanceSquared(fRow, cData.AsSpan(g * d, d));
					if (sqDist < minSqDist)
					{
						minSqDist = sqDist;
						bestCluster = g;
					}
				}
				labels[i] = bestCluster;

				fNorms[i] = TensorPrimitives.SumOfSquares(fRow);
			}

			var centroidNorms = new double[_k];
			var harmonicWeights = new double[_k];
			var newLabels = new int[n];
			int maxIterations = Math.Max(10, (int)(8 * Math.Log2(_k)));
			int iter = 0;

			// Convergence threshold: stop early if < 0.1% of pixels change clusters
			int convergenceThreshold = (int)(n * 0.005);

			// 3. High-Speed Optimization Loop
			while (iter < maxIterations)
			{
				iter++;

				// Fast O(N^2) total cluster objective calculation without dense matrix operations
				Parallel.For(0, _k, j =>
				{
					var obj = 0.0;
					for (int i = 0; i < n; ++i)
					{
						if (labels[i] != j)
							continue;
						for (int col = 0; col < n; ++col)
						{
							if (labels[col] == j)
								obj += lData[i + col * n]; // Column-major lookup for L
						}
					}
					var val = Math.Max(obj, 1e-8);
					harmonicWeights[j] = 1.0 / (val * val + _gamma);
				});

				// Cache centroid squared norms using TensorPrimitives
				Parallel.For(0, _k, g =>
				{
					centroidNorms[g] = TensorPrimitives.SumOfSquares(cData.AsSpan(g * d, d));
				});

				int changesCount = 0;

				// Parallel node assignment evaluation
				Parallel.For(0, n, i =>
				{
					int currentCluster = labels[i];
					int bestCluster = currentCluster;
					var minCost = double.MaxValue;
					
					var fRow = fDataRowMajor.AsSpan(i * d, d);
					var aRow = aDataRowMajor.AsSpan(i * n, n);
					var totalAff = totalAffinities[i];

					for (int g = 0; g < _k; ++g)
					{
						var laplacianCost = 0.0;
						if (totalAff > 0)
						{
							var alignedAffinity = 0.0;
							for (int j = 0; j < n; ++j)
							{
								if (labels[j] == g)
									alignedAffinity += aRow[j];
							}
							var misalignedAffinity = totalAff - alignedAffinity;
							laplacianCost = (misalignedAffinity / totalAff) * harmonicWeights[g];
						}

						// Vectorized Dot Product via TensorPrimitives
						var dot = TensorPrimitives.Dot(fRow, cData.AsSpan(g * d, d));
						var distance = Math.Sqrt(Math.Max(0.0, fNorms[i] + centroidNorms[g] - 2.0 * dot));

						var combinedCost = (_w * laplacianCost) + ((1.0 - _w) * distance);

						if (combinedCost < minCost)
						{
							minCost = combinedCost;
							bestCluster = g;
						}
					}

					newLabels[i] = bestCluster;
					if (bestCluster != currentCluster)
						System.Threading.Interlocked.Increment(ref changesCount);
				});

				if (changesCount > 0)
				{
					Array.Copy(newLabels, labels, n);
					UpdateCodebook(fDataRowMajor, labels, cData, n, d, _k);
				}

				if (changesCount <= convergenceThreshold)
				{
					break;
				}
			}

			return cData;
		}

		private static void UpdateCodebook(float[] fDataRow, int[] labels, float[] cData, int n, int d, int k)
		{
			Array.Clear(cData, 0, cData.Length);
			int[] counts = new int[k];

			for (int i = 0; i < n; ++i)
			{
				int label = labels[i];
				counts[label]++;
				
				var centroidRow = cData.AsSpan(label * d, d);
				var fRow = fDataRow.AsSpan(i * d, d);
				
				TensorPrimitives.Add(centroidRow, fRow, centroidRow);
			}

			for (int g = 0; g < k; ++g)
			{
				if (counts[g] > 0)
				{
					var centroidRow = cData.AsSpan(g * d, d);
					TensorPrimitives.Divide(centroidRow, counts[g], centroidRow);
				}
			}
		}

		private static float[] MultiplyMatricesColumnMajor(float[] A, int rowsA, int colsA, float[] B, int rowsB, int colsB)
		{
			var C = new float[rowsA * colsB];
			Span<float> temp = stackalloc float[rowsA];

			for (int j = 0; j < colsB; ++j)
			{
				int cColOffset = j * rowsA;
				int bColOffset = j * rowsB;

				for (int k = 0; k < colsA; ++k)
				{
					var bVal = B[bColOffset + k];
					if (bVal == 0.0)
						continue;

					int aColOffset = k * rowsA;

					var aCol = A.AsSpan(aColOffset, rowsA);
                    var cCol = C.AsSpan(cColOffset, rowsA);

					// temp = aCol * bVal
					TensorPrimitives.Multiply(aCol, bVal, temp);
					
					// cCol = cCol + temp
					TensorPrimitives.Add(cCol, temp, cCol);
				}
			}
			return C;
		}

		internal void Grhvquan(Bitmap source, int[] pixels, ref Color[] palettes, ref int nMaxColors)
		{
			saliencies = new float[pixels.Length];
			var saliencyBase = .1f;

			/* Build histogram */
			for (int i = 0; i < pixels.Length; ++i)
			{
				var pixel = pixels[i];
				var c = Color.FromArgb(pixel);
				if (c.A <= alphaThreshold)
					c = m_transparentColor;

				GetLab(pixel, out var lab1);

				if (saliencies != null)
					saliencies[i] = (float) (saliencyBase + (1 - saliencyBase) * lab1.L / 100f * lab1.alpha / 255f);
			}

			mDivn = Math.Min(0.9, nMaxColors * 1.0 / pixelMap.Count);

			if (pixelMap.Count <= nMaxColors)
			{
				/* Fill palette */
				palettes = new Color[pixelMap.Count];
				int i = 0;
				foreach (var pixel in pixelMap.Keys)
				{
					var c = Color.FromArgb(pixel);
					palettes[i++] = c;

					if (i > 1 && c.A == 0)
						BitmapUtilities.Swap(ref palettes[i - 1], ref palettes[0]);
				}
				nMaxColors = i;
				Console.WriteLine("Maximum number of colors: " + palettes.Length);
				return;
			}

			int sampleSize = 64;

            using Bitmap resizedBitmap = new Bitmap(source, new Size(sampleSize, sampleSize));

			int N = sampleSize * sampleSize;
			int d = HasAlpha ? 4 : 3; // CIELAB alpha channels
			var X = new float[N * d]; // Column-major 1D flat array [N x 3]

			for (int y = 0; y < sampleSize; ++y)
			{
				for (int x = 0; x < sampleSize; ++x)
				{
					int idx = y * sampleSize + x;
					var pixel = resizedBitmap.GetPixel(x, y);
					GetLab(pixel.ToArgb(), out var lab1);
					X[idx + 0 * N] = (float) lab1.L;
					X[idx + 1 * N] = (float) lab1.A;
					X[idx + 2 * N] = (float) lab1.B;
					if (HasAlpha)
						X[idx + 3 * N] = pixel.A * 1.0f;
				}
			}

			// Build 4-neighbor spatial-color Laplacian L = D - A [N x N]
			var L = BuildGridLaplacian(X, sampleSize, sampleSize, N);

			var codebook = Optimize(X, L, N, d);
			/* Fill palette */
			for (int k = 0; k < nMaxColors; ++k)
			{
				int offset = k * d;
				var lab1 = new CIELABConvertor.Lab
				{
					alpha = (hasSemiTransparency || HasAlpha) ? Math.Round(codebook[offset + 3]) : Byte.MaxValue,
					L = codebook[offset], A = codebook[offset + 1], B = codebook[offset + 2]
				};
				palettes[k] = CIELABConvertor.LAB2RGB(lab1);
            }
		}

		internal ushort NearestColorIndex(Color[] palette, int pixel, int pos)
		{
			var nMaxColors = palette.Length;
			int offset = GetColorIndex(pixel);
			if (nearestMap.TryGetValue(offset, out var k))
				return k;

			var c = Color.FromArgb(pixel);
			if (c.A <= alphaThreshold)
				c = m_transparentColor;

			if (palette.Length > 2 && HasAlpha && c.A > alphaThreshold)
				k = 1;

			double mindist = 1e100;
			GetLab(pixel, out var lab1);

			for (int i = k; i < nMaxColors; ++i)
			{
				var c2 = palette[i];

				var curdist = hasSemiTransparency ? BitmapUtilities.Sqr(c2.A - c.A) * TRANS_RATE : 0;
				if (curdist > mindist)
					continue;

				GetLab(c2.ToArgb(), out var lab2);
				if (nMaxColors <= 4)
				{
					curdist = BitmapUtilities.Sqr(c2.R - c.R) + BitmapUtilities.Sqr(c2.G - c.G) + BitmapUtilities.Sqr(c2.B - c.B);
					if(hasSemiTransparency)
						curdist += BitmapUtilities.Sqr(c2.A - c.A);
				}
				else if (hasSemiTransparency || nMaxColors < 16)
				{
					curdist += BitmapUtilities.Sqr(lab2.L - lab1.L);
					if (curdist > mindist)
						continue;
				
					curdist += BitmapUtilities.Sqr(lab2.A - lab1.A);
					if (curdist > mindist)
						continue;
				
					curdist += BitmapUtilities.Sqr(lab2.B - lab1.B);
				}
				else if (nMaxColors > 32)
				{
					curdist += Math.Abs(lab2.L - lab1.L);
					if (curdist > mindist)
						continue;

					curdist += Math.Sqrt(BitmapUtilities.Sqr(lab2.A - lab1.A) + BitmapUtilities.Sqr(lab2.B - lab1.B));
				}
				else
				{
					var deltaL_prime_div_k_L_S_L = CIELABConvertor.L_prime_div_k_L_S_L(lab1, lab2);
					curdist += BitmapUtilities.Sqr(deltaL_prime_div_k_L_S_L);
					if (curdist > mindist)
						continue;

					var deltaC_prime_div_k_L_S_L = CIELABConvertor.C_prime_div_k_L_S_L(lab1, lab2, out var a1Prime, out var a2Prime, out var CPrime1, out var CPrime2);
					curdist += BitmapUtilities.Sqr(deltaC_prime_div_k_L_S_L);
					if (curdist > mindist)
						continue;

					var deltaH_prime_div_k_L_S_L = CIELABConvertor.H_prime_div_k_L_S_L(lab1, lab2, a1Prime, a2Prime, CPrime1, CPrime2, out var barCPrime, out var barhPrime);
					curdist += BitmapUtilities.Sqr(deltaH_prime_div_k_L_S_L);
					if (curdist > mindist)
						continue;

					curdist += CIELABConvertor.R_T(barCPrime, barhPrime, deltaC_prime_div_k_L_S_L, deltaH_prime_div_k_L_S_L);
				}

				if (curdist > mindist)
					continue;
				mindist = curdist;
				k = (ushort)i;
			}
			nearestMap[offset] = k;
			return k;
		}
		
		protected ushort ClosestColorIndex(Color[] palette, int pixel, int pos)
		{
			var c = Color.FromArgb(pixel);
			if (c.A <= alphaThreshold)
				return NearestColorIndex(palette, pixel, pos);

			var nMaxColors = palette.Length;
			int offset = GetColorIndex(pixel);
			if (!closestMap.TryGetValue(offset, out var closest))
			{
				closest = new ushort[4];
				closest[2] = closest[3] = ushort.MaxValue;

				for (ushort k = 0; k < nMaxColors; ++k)
				{
					var c2 = palette[k];
					var err = PR * BitmapUtilities.Sqr(c.R - c2.R);
					if (err >= closest[3])
						continue;

					err += PG * BitmapUtilities.Sqr(c.G - c2.G);
					if (err >= closest[3])
						continue;

					err += PB * BitmapUtilities.Sqr(c.B - c2.B);
					if (err >= closest[3])
						continue;

					if (hasSemiTransparency)
						err += PA * BitmapUtilities.Sqr(c.A - c2.A);

					if (err < closest[2])
					{
						closest[1] = closest[0];
						closest[3] = closest[2];
						closest[0] = k;
						closest[2] = (ushort) err;
					}
					else if (err < closest[3])
					{
						closest[1] = k;
						closest[3] = (ushort) err;
					}
				}

				if (closest[3] == ushort.MaxValue)
					closest[1] = closest[0];

				closestMap[offset] = closest;
			}

			int idx = 1;
			if (closest[2] == 0 || (rand.Next(closest[3] + closest[2])) <= closest[3])
				idx = 0;

			var MAX_ERR = palette.Length;
			if (closest[idx + 2] >= MAX_ERR || closest[idx] == 0 || palette[closest[idx]].A < c.A)
				return NearestColorIndex(palette, pixel, pos);
			return closest[idx];
		}

		public ushort DitherColorIndex(Color[] palette, int pixel, int pos)
		{
			var nMaxColors = palette.Length;
			if (nMaxColors <= 4)
				return NearestColorIndex(palette, pixel, pos);
			return ClosestColorIndex(palette, pixel, pos);
		}

		internal void Clear()
		{
			m_palette = null;
			saliencies = null;
			closestMap.Clear();
			nearestMap.Clear();
		}

		protected int[] Dither(int[] pixels, Color[] palettes, int width, int height, int frameIndex, bool dither)
		{
			this.dither = dither;
			if (hasSemiTransparency)
				mDivn *= -1;

			var qPixels = GilbertCurve.Dither(width, height, pixels, palettes, this, saliencies, mDivn, frameIndex, dither);

			if (!dither && palettes.Length > 32)
			{
				var delta = BitmapUtilities.Sqr(palettes.Length) / pixelMap.Count;
				var mDivn = delta > 0.023 ? 1.0f : (float)(36.921 * delta + 0.906);
				BlueNoise.Dither(width, height, pixels, palettes, this, qPixels, mDivn);
			}

			return qPixels;
		}

		protected bool IsValidFormat(PixelFormat pixelFormat, int nMaxColors)
		{
			if (pixelFormat == PixelFormat.Undefined)
				return false;

			int bitDepth = Image.GetPixelFormatSize(pixelFormat);
			return Math.Pow(2, bitDepth) >= nMaxColors;
		}

		internal int[] GrabPixels(Bitmap source, int nMaxColors, ref bool hasSemiTransparency)
		{
			var bitmapWidth = source.Width;
			var bitmapHeight = source.Height;
			var pixels = new int[bitmapWidth * bitmapHeight];
			int semiTransCount = 0;
			if (!BitmapUtilities.GrabPixels(source, pixels, ref semiTransCount, ref m_transparentColor, ref m_transparentPixelIndex, alphaThreshold, nMaxColors))
				return null;
			this.hasSemiTransparency = hasSemiTransparency = semiTransCount > 0;
			return pixels;
		}

		public Bitmap QuantizeImage(Bitmap source, PixelFormat pixelFormat, int nMaxColors, int frameIndex, bool dither)
		{
			_k = nMaxColors;

            if (nMaxColors <= 32)
				PR = PG = PB = PA = 1;
			else
			{
				PR = coeffs[0, 0]; PG = coeffs[0, 1]; PB = coeffs[0, 2];
			}

			if (!IsValidFormat(pixelFormat, nMaxColors))
			{
				if (nMaxColors > 256)
					pixelFormat = HasAlpha ? PixelFormat.Format16bppArgb1555 : PixelFormat.Format16bppRgb565;
				else
					pixelFormat = (nMaxColors > 16) ? PixelFormat.Format8bppIndexed : (nMaxColors > 2) ? PixelFormat.Format4bppIndexed : PixelFormat.Format1bppIndexed;
			}

			var bitmapWidth = source.Width;
			var bitmapHeight = source.Height;

			var dest = new Bitmap(bitmapWidth, bitmapHeight, pixelFormat);
			var pixels = GrabPixels(source, nMaxColors, ref hasSemiTransparency);
			if (pixels == null)
				return dest;

			if (m_palette == null) {
				var palettes = dest.Palette.Entries;
				if (palettes.Length != nMaxColors)
					palettes = new Color[nMaxColors];

				if (nMaxColors > 2) {
					Grhvquan(source, pixels, ref palettes, ref nMaxColors);
				}
				else
				{
					if (HasAlpha)
					{
						palettes[0] = m_transparentColor;
						palettes[1] = Color.Black;
					}
					else
					{
						palettes[0] = Color.Black;
						palettes[1] = Color.White;
					}
				}
				m_palette = palettes;
			}

			var qPixels = Dither(pixels, m_palette, bitmapWidth, bitmapHeight, frameIndex, dither);

			if (nMaxColors > 256)
				return BitmapUtilities.ProcessImagePixels(dest, qPixels, hasSemiTransparency, m_transparentPixelIndex);

			return BitmapUtilities.ProcessImagePixels(dest, m_palette, qPixels, HasAlpha);
		}

		public Bitmap QuantizeImage(Bitmap source, PixelFormat pixelFormat, int nMaxColors, bool dither)
		{
			return QuantizeImage(source, pixelFormat, nMaxColors, 0, dither);
		}
	}
}