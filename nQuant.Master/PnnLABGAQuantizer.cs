using nQuant.Master;
using nQuant.Master.Ga;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Drawing;
using System.IO;
using System.Linq;
using System.Text;

/* Fast pairwise nearest neighbor based genetic algorithm with CIELAB color space
Copyright (c) 2023 Miller Cy Chan
* error measure; time used is proportional to number of bins squared - WJ */

namespace PnnQuant
{
	public class PnnLABGAQuantizer : Chromosome<PnnLABGAQuantizer>, IDisposable
	{
		// boolean variable to ensure dispose
		// method executes only once
		private bool disposedValue;

		private double _fitness = double.NegativeInfinity;
		public double ratioX = 0, ratioY = 0;
		private double[] _objectives;
		private PnnLABQuantizer m_pq;
		public double[] ConvertedObjectives { get; private set; }
		public double[] Objectives { get => _objectives; }

		private readonly Random _random;

		private static int _dp = 1, _nMaxColors = 256;
		private static double minRatio = 0, maxRatio = 1.0;
		private static List<int[]> m_pixelsList;
		private static List<int> _bitmapWidths;

		private static readonly ConcurrentDictionary<string, double[]> _fitnessMap = new();

		public PnnLABGAQuantizer(PnnLABQuantizer pq, List<Bitmap> sources, int nMaxColors)
		{
			// increment value when criteria violation occurs
			_objectives = new double[4];
			_random = new Random(sources[0].Width * sources[0].Height);
			m_pq = new PnnLABQuantizer(pq);
			if (pq.IsGA)
				return;

			_nMaxColors = nMaxColors;
			_bitmapWidths = new List<int>();
			m_pixelsList = new List<int[]>();
			
			var hasSemiTransparency = false;
			foreach(var source in sources) {
				_bitmapWidths.Add(source.Width);
				m_pixelsList.Add(m_pq.GrabPixels(source, _nMaxColors, ref hasSemiTransparency));
			}
			minRatio = (hasSemiTransparency || nMaxColors < 64) ? .0111 : .85;
			maxRatio = Math.Min(1.0, nMaxColors / ((nMaxColors < 64) ? 400.0 : 50.0));
			_dp = maxRatio < .1 ? 10000 : 100;
		}

		private PnnLABGAQuantizer(PnnLABQuantizer pq, List<int[]> pixelsList, List<int> bitmapWidths, int nMaxColors)
		{
			m_pq = new PnnLABQuantizer(pq);
			m_pixelsList = pixelsList;
			_bitmapWidths = bitmapWidths;
			_random = new Random(m_pixelsList[0].Length);
			_nMaxColors = nMaxColors;
		}

		private string RatioKey
		{
			get
			{
				var sb = new StringBuilder();
				sb.Append((int) (ratioX * _dp));
				var difference = Math.Abs(ratioX - ratioY);
				if (difference <= 0.0000001)
					return sb.ToString();

				sb.Append(';').Append((int) (ratioY * _dp * 100));
				return sb.ToString();
			}
		}

		private void CalculateError(double[] errors)
		{
			var maxError = maxRatio < .1 ? .5 : .0625;
			if(m_pq.HasAlpha)
				maxError = 1;

			double fitness = 0;
			int length = m_pixelsList.Select(pixels => pixels.Length).Sum();
			for (int i = 0; i < errors.Length; ++i)
				errors[i] /= maxError * length;

			for (int i = 0; i < errors.Length; ++i)
			{
				if (i > 0)
					errors[i] /= 2.55;
				fitness -= errors[i];
			}

			_objectives = errors;
			_fitness = fitness;
		}


		private void CalculateFitness()
		{
			var ratioKey = RatioKey;
			if (_fitnessMap.TryGetValue(ratioKey, out _objectives))
			{
				_fitness = -1f * _objectives.Sum();
				return;
			}

			_objectives = new double[4];
			m_pq.SetRatio(ratioX, ratioY);
			var palette = new Color[_nMaxColors];
			m_pq.Pnnquan(m_pixelsList[0], ref palette, ref _nMaxColors);

			int threshold = maxRatio < .1 ? -64 : -112;
			var errors = new double[_objectives.Length];
			m_pixelsList.ForEach(pixels => {
				for (int i = 0; i < pixels.Length; ++i)
				{
					if (BlueNoise.TELL_BLUE_NOISE[i & 4095] > threshold)
						continue;

					m_pq.GetLab(pixels[i], out CIELABConvertor.Lab lab1);
					var qPixelIndex = m_pq.NearestColorIndex(palette, pixels[i], i);
					m_pq.GetLab(palette[qPixelIndex].ToArgb(), out CIELABConvertor.Lab lab2);

					if (m_pq.HasAlpha)
					{
						errors[0] += BitmapUtilities.Sqr(lab2.L - lab1.L);
						errors[1] += BitmapUtilities.Sqr(lab2.A - lab1.A);
						errors[2] += BitmapUtilities.Sqr(lab2.B - lab1.B);
						errors[3] += BitmapUtilities.Sqr(lab2.alpha - lab1.alpha) / Math.Exp(1.5);
					}
					else
					{
						errors[0] += Math.Abs(lab2.L - lab1.L);
						errors[1] += Math.Sqrt(BitmapUtilities.Sqr(lab2.A - lab1.A) + BitmapUtilities.Sqr(lab2.B - lab1.B));
					}
				}
			});
			CalculateError(errors);
			_fitnessMap[ratioKey] = _objectives;
		}

		public List<Bitmap> QuantizeImage(bool dither)
		{
			m_pq.SetRatio(ratioX, ratioY);
			var palette = new Color[_nMaxColors];
			m_pq.Pnnquan(m_pixelsList[0], ref palette, ref _nMaxColors);

			var bitmaps = m_pixelsList.Select((pixels, i) => m_pq.QuantizeImage(pixels, _bitmapWidths[i], _nMaxColors, dither));
			m_pq.Clear();
			return bitmaps.ToList();
		}

		protected virtual void Dispose(bool disposing)
		{
			// check if already disposed
			if (!disposedValue)
			{
				if (disposing)
				{
					// free managed objects here
					m_pixelsList = null;
					_fitnessMap.Clear();
				}

				// free unmanaged objects here

				// set the bool value to true
				disposedValue = true;
			}
		}

		// The consumer object can call
		// the below dispose method
		public void Dispose()
		{
			// Invoke the above virtual
			// dispose(bool disposing) method
			Dispose(disposing: true);

			// Notify the garbage collector
			// about the cleaning event
			GC.SuppressFinalize(this);
		}

		private double Randrange(double min, double max)
		{
			return min + _random.NextDouble() * (max - min);
		}

		public double[] Ratios
		{
			get => new double[] { ratioX, ratioY };
		}

		public void SetRatio(double ratioX, double ratioY)
		{
			var difference = Math.Abs(ratioX - ratioY);
			if (difference <= minRatio)
				ratioY = ratioX;
			this.ratioX = Math.Min(Math.Max(ratioX, minRatio), maxRatio);
			this.ratioY = Math.Min(Math.Max(ratioY, minRatio), maxRatio);
		}

		public float Fitness
		{
			get => (float) _fitness;
		}
		
		private double RotateLeft(double u, double v, double delta = 0.0) {
			var theta = Math.PI * Randrange(minRatio, maxRatio) / Math.Exp(delta);
			var result = u * Math.Sin(theta) + v * Math.Cos(theta);
			if(result <= minRatio || result >= maxRatio)
				result = RotateLeft(u, v, delta + .5);
			return result;
		}
		
		private double RotateRight(double u, double v, double delta = 0.0) {
			var theta = Math.PI * Randrange(minRatio, maxRatio) / Math.Exp(delta);
			var result = u * Math.Cos(theta) - v * Math.Sin(theta);
			if(result <= minRatio || result >= maxRatio)
				result = RotateRight(u, v, delta + .5);
			return result;
		}

		public PnnLABGAQuantizer Crossover(PnnLABGAQuantizer mother, int numberOfCrossoverPoints, float crossoverProbability)
		{
			var child = MakeNewFromPrototype();
			if (_random.Next(100) <= crossoverProbability)
				return child;

			var ratioX = RotateRight(this.ratioX, mother.Ratios[1]);
			var ratioY = RotateLeft(this.ratioY, mother.Ratios[0]);
			child.SetRatio(ratioX, ratioY);
			child.CalculateFitness();
			return child;
		}

		private double BoxMuller(double value) {
			var r1 = Randrange(minRatio, maxRatio);
			return Math.Sqrt(-2 * Math.Log(value)) * Math.Cos(2 * Math.PI * r1);
		}

		public void Mutation(int mutationSize, float mutationProbability)
		{
			// check probability of mutation operation
			if (_random.Next(100) > mutationProbability)
				return;

			var ratioX = this.ratioX;
			var ratioY = this.ratioY;
			if (_random.NextDouble() > .5)
				ratioX = BoxMuller(ratioX);
			else
				ratioY = BoxMuller(ratioY);

			SetRatio(ratioX, ratioY);
			CalculateFitness();
		}

		public void ResizeConvertedObjectives(int numObj)
		{
			ConvertedObjectives = new double[numObj];
		}

		public PnnLABGAQuantizer MakeNewFromPrototype()
		{
			var child = new PnnLABGAQuantizer(m_pq, m_pixelsList, _bitmapWidths, _nMaxColors);
			var minRatio2 = 2 * minRatio;
			if (minRatio2 > 1)
				minRatio2 = 0;
			var ratioX = Randrange(minRatio, maxRatio);
			var ratioY = ratioX < minRatio2 ? Randrange(minRatio, maxRatio) : ratioX;
			child.SetRatio(ratioX, ratioY);
			child.CalculateFitness();
			return child;
		}

		public Random Random
		{
			get => _random;
		}

		public String Result
		{
			get
			{
				var difference = Math.Abs(ratioX - ratioY);
				if (difference <= 0.0000001)
					return ratioX.ToString("0.######");
				return ratioX.ToString("0.######") + ", " + ratioY.ToString("0.######");
			}
		}

		public bool HasAlpha
		{
            get => m_pq.HasAlpha;
		}

		public static int MaxColors
		{
			get => _nMaxColors;
		}

		public bool Dominates(PnnLABGAQuantizer other)
		{
			var better = false;
			for (int f = 0; f < Objectives.Length; ++f)
			{
				if (Objectives[f] > other.Objectives[f])
					return false;

				if (Objectives[f] < other.Objectives[f])
					better = true;
			}
			return better;
		}

	}
}
