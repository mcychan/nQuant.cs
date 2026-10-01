using nQuant.Master;
using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using System.Numerics.Tensors;
using System.Threading.Tasks;

/* Density-weighted pairwise nearest neighbor K-means algorithm with CIELAB color space advanced version
Copyright (c) 2018-2026 Miller Cy Chan
* error measure; time used is proportional to number of bins squared - WJ */

namespace PnnQuant
{
	public class DwPnnKmeansLABQuantizer : PnnLABQuantizer
	{
		public DwPnnKmeansLABQuantizer()
		{
		}

		/// <summary>
		/// Refines the initial PNN palette using Density-Weighted K-Means 
		/// optimized with System.Numerics.Tensors.
		/// </summary>
		internal void RefinePaletteWithDwKMeans(CIELABConvertor.Lab[] uniqueLabs, float[] weights, ref Color[] palettes, int maxIterations = 10)
		{
			int k = palettes.Length;
			var centroids = new CIELABConvertor.Lab[k];

			// Initialize centroids from the initial palette
			for (int i = 0; i < k; ++i)
			{
				GetLab(palettes[i].ToArgb(), out centroids[i]);
			}

			var assignments = new int[uniqueLabs.Length];
			var newL = new float[k];
			var newA = new float[k];
			var newB = new float[k];
			var newAlpha = new float[k];
			var weightSum = new float[k];

			for (int iter = 0; iter < maxIterations; iter++)
			{
				Array.Clear(weightSum, 0, weightSum.Length);
				Array.Clear(newL, 0, newL.Length);
				Array.Clear(newA, 0, newA.Length);
				Array.Clear(newB, 0, newB.Length);
				Array.Clear(newAlpha, 0, newAlpha.Length);

				// Assignment phase with density weighting
				Parallel.For(0, uniqueLabs.Length, i =>
				{
					var lab = uniqueLabs[i];
					var minCost = float.MaxValue;
					int bestCluster = 0;

					// Use stackalloc to avoid allocation overhead and ensure safe span references
					Span<float> p1 = stackalloc float[3];
					p1[0] = lab.L;
					p1[1] = lab.A;
					p1[2] = lab.B;

					Span<float> p2 = stackalloc float[3];

					for (int j = 0; j < k; j++)
					{
						p2[0] = centroids[j].L;
						p2[1] = centroids[j].A;
						p2[2] = centroids[j].B;

						// Use TensorPrimitives for high-speed Euclidean distance calculation in LAB space
						float dist = TensorPrimitives.Distance(p1, p2);

						// Apply density weight w_i to penalize/reward transitions
						float cost = dist * weights[i];

						if (cost < minCost)
						{
							minCost = cost;
							bestCluster = j;
						}
					}
					assignments[i] = bestCluster;
				});

				// Update phase (Weighted Centroids)
				for (int i = 0; i < uniqueLabs.Length; ++i)
				{
					int cIdx = assignments[i];
					float w = weights[i];

					newL[cIdx] += uniqueLabs[i].L * w;
					newA[cIdx] += uniqueLabs[i].A * w;
					newB[cIdx] += uniqueLabs[i].B * w;
					newAlpha[cIdx] += (float)uniqueLabs[i].alpha * w;
					weightSum[cIdx] += w;
				}

				bool changed = false;
				for (int j = 0; j < k; j++)
				{
					if (weightSum[j] > 0f)
					{
						var updatedLab = new CIELABConvertor.Lab
						{
							L = newL[j] / weightSum[j],
							A = newA[j] / weightSum[j],
							B = newB[j] / weightSum[j],
							alpha = newAlpha[j] / weightSum[j]
						};

						// Check convergence movement threshold if needed
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
			RefinePaletteWithDwKMeans(uniqueLabsList.ToArray(), finalWeightsList.ToArray(), ref palettes, maxIterations);
		}

	}
}
