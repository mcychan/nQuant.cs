using nQuant.Master;
using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using System.Numerics.Tensors;
using System.Threading.Tasks;

/* Density-weighted pairwise nearest neighbor K-harmonic means algorithm with CIELAB color space advanced version
Copyright (c) 2018-2026 Miller Cy Chan
* error measure; time used is proportional to number of bins squared - WJ */

namespace PnnQuant
{
	public class DwPnnKhmLABQuantizer : PnnLABQuantizer
	{
		public DwPnnKhmLABQuantizer()
		{
		}

		/// <summary>
		/// Refines the initial PNN palette using Density-Weighted K-harmonic means 
		/// optimized with System.Numerics.Tensors.
		/// </summary>
		internal void RefinePaletteWithDwKhMeans(CIELABConvertor.Lab[] uniqueLabs, float[] weights, ref Color[] palettes, int maxIterations = 10)
		{
			int k = palettes.Length;
			var centroids = new CIELABConvertor.Lab[k];

			// Initialize centroids from the initial palette
			for (int i = 0; i < k; ++i)
			{
				GetLab(palettes[i].ToArgb(), out centroids[i]);
			}

			var newL = new float[k];
			var newA = new float[k];
			var newB = new float[k];
			var newAlpha = new float[k];
			var weightSum = new float[k];

			const float p = 2.0f; // K-harmonic means power parameter

			for (int iter = 0; iter < maxIterations; ++iter)
			{
				Array.Clear(weightSum, 0, weightSum.Length);
				Array.Clear(newL, 0, newL.Length);
				Array.Clear(newA, 0, newA.Length);
				Array.Clear(newB, 0, newB.Length);
				Array.Clear(newAlpha, 0, newAlpha.Length);

				// K-harmonic means update phase with density weighting
				Parallel.For(0, uniqueLabs.Length, i =>
				{
					var lab = uniqueLabs[i];
					Span<float> p1 = stackalloc float[3] { lab.L, lab.A, lab.B };
					Span<float> p2 = stackalloc float[3];

					Span<float> distances = stackalloc float[k];
					float sumInvDistP = 0f;

					for (int j = 0; j < k; j++)
					{
						p2[0] = centroids[j].L;
						p2[1] = centroids[j].A;
						p2[2] = centroids[j].B;

						var dist = TensorPrimitives.Distance(p1, p2);
						if (dist < 1e-5f)
							dist = 1e-5f; // Prevent division by zero

						var invDistP = (float)Math.Pow(dist, -p);
						distances[j] = invDistP;
						sumInvDistP += invDistP;
					}

					// Accumulate weighted contributions across all clusters using KHM soft membership
					var w = weights[i];
					for (int j = 0; j < k; j++)
					{
						// KHM membership weight factor: d(x, c_j)^(-p-2) / (sum d(x, c_l)^(-p))^2
						float dist = TensorPrimitives.Distance(p1, stackalloc float[3] { centroids[j].L, centroids[j].A, centroids[j].B });
						if (dist < 1e-5f)
							dist = 1e-5f;

						var khmWeight = w * (float)(Math.Pow(dist, -p - 2.0) / Math.Pow(sumInvDistP, 2.0));

						lock (newL) // Thread-safe accumulation for parallel loop
						{
							newL[j] += lab.L * khmWeight;
							newA[j] += lab.A * khmWeight;
							newB[j] += lab.B * khmWeight;
							newAlpha[j] += (float)lab.alpha * khmWeight;
							weightSum[j] += khmWeight;
						}
					}
				});

				var changed = false;
				for (int j = 0; j < k; j++)
				{
					if (weightSum[j] > 1e-5f)
					{
						var updatedLab = new CIELABConvertor.Lab
						{
							L = newL[j] / weightSum[j],
							A = newA[j] / weightSum[j],
							B = newB[j] / weightSum[j],
							alpha = newAlpha[j] / weightSum[j]
						};

						centroids[j] = updatedLab;
						palettes[j] = CIELABConvertor.LAB2RGB(updatedLab);
						changed = true;
					}
				}

				if (!changed)
					break;
			}
		}

		internal override void Pnnquan(int[] pixels, ref Color[] palettes, ref int nMaxColors)
		{
			var bins = Getbins(pixels, ref palettes, out var maxbins, ref nMaxColors);

			// Build synchronized uniqueLabs and weights arrays from the non-empty bins
			var uniqueLabsList = new List<CIELABConvertor.Lab>();
			var finalWeightsList = new List<float>();

			for (int i = 0; i < maxbins; ++i)
			{
				var bin = bins[i];
				if (bin == null)
					continue;

				uniqueLabsList.Add(new CIELABConvertor.Lab
				{
					alpha = bin.ac,
					L = bin.Lc,
					A = bin.Ac,
					B = bin.Bc
				});
				finalWeightsList.Add(bin.cnt * (1.0f + (float)Math.Abs(bin.Lc / 100.0f)));
			}

			int maxIterations = nMaxColors <= 16 ? 5 : (nMaxColors <= 64 ? 10 : 18);
			RefinePaletteWithDwKhMeans(uniqueLabsList.ToArray(), finalWeightsList.ToArray(), ref palettes, maxIterations);
		}

	}
}
