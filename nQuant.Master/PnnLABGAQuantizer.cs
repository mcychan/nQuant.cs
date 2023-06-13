using nQuant.Master;
using nQuant.Master.Ga;
using System;
using System.Collections.Concurrent;
using System.Drawing;
using System.Linq;

/* Fast pairwise nearest neighbor based genetic algorithm with CIELAB color space advanced version
Copyright (c) 2023 Miller Cy Chan
* error measure; time used is proportional to number of bins squared - WJ */

namespace PnnQuant
{
        public class PnnLABGAQuantizer : Chromosome<PnnLABGAQuantizer>, IDisposable
        {
		// boolean variable to ensure dispose
		// method executes only once
		private bool disposedValue;

		private float _fitness = float.NegativeInfinity;
		public double ratio = 0;
		private double[] _objectives;
		private PnnLABQuantizer m_pq;
		public double[] ConvertedObjectives { get; private set; }
		public double[] Objectives { get => _objectives; }

		private readonly Random _random;
		
		private static int _dp = 1, _nMaxColors = 256;
		private static double minRatio = 0, maxRatio = 1.0;
		private static int[] m_pixels;
		private static int _bitmapWidth;

		private static readonly ConcurrentDictionary<short, double[]> fitnessMap = new();
		private static readonly ConcurrentDictionary<short, Color[]> paletteMap = new();

		public PnnLABGAQuantizer(PnnLABQuantizer pq, Bitmap source, int nMaxColors) {
			// increment value when criteria violation occurs
			_objectives = new double[4];
			_bitmapWidth = source.Width;
			_random = new Random(_bitmapWidth * source.Height);
            
			m_pq = new PnnLABQuantizer(pq);
			if(pq.IsGA)
				return;

			_nMaxColors = nMaxColors;

			var hasSemiTransparency = false;
			m_pixels = m_pq.GrabPixels(source, _nMaxColors, ref hasSemiTransparency);
			minRatio = (hasSemiTransparency || nMaxColors < 64) ? 0 : .9;
			maxRatio = Math.Min(1.0, nMaxColors / ((nMaxColors < 64) ? 500.0 : 50.0));
			_dp = maxRatio < .1 ? 1000 : 10;
			if(hasSemiTransparency)
				maxRatio = .1;
		}

		private PnnLABGAQuantizer(PnnLABQuantizer pq, int[] pixels, int bitmapWidth, int nMaxColors)
		{
			m_pq = pq;
			m_pixels = pixels;
			_bitmapWidth = bitmapWidth;
			_random = new Random(pixels.Length);
			_nMaxColors = nMaxColors;
		}

		private short RatioKey
		{
			get
			{
				return (short)(ratio * _dp);
			}
		}

		private void CalculateFitness() {
			var ratioKey = RatioKey;
			if (fitnessMap.TryGetValue(ratioKey, out _objectives)) {
				_fitness = -1f * (float) _objectives.Sum();
				m_pq.Palette = paletteMap[ratioKey];
				return;
			}

			_objectives = new double[4];
			m_pq.Ratio = Ratio;
			var palette = new Color[_nMaxColors];
			m_pq.Pnnquan(m_pixels, ref palette, ref _nMaxColors);
			m_pq.Palette = palette;

			var errors = new double[_objectives.Length];
			for (int i = 0; i < m_pixels.Length; ++i) {
				if(BlueNoise.RAW_BLUE_NOISE[i & 4095] > -112)
					continue;

				m_pq.GetLab(m_pixels[i], out CIELABConvertor.Lab lab1);
				var qPixelIndex = m_pq.NearestColorIndex(palette, m_pixels[i], i);
				m_pq.GetLab(palette[qPixelIndex].ToArgb(), out CIELABConvertor.Lab lab2);

				if (m_pq.HasAlpha) {
					errors[0] += BitmapUtilities.Sqr(lab2.L - lab1.L);
					errors[1] += BitmapUtilities.Sqr(lab2.A - lab1.A);
					errors[2] += BitmapUtilities.Sqr(lab2.B - lab1.B);
					errors[3] += BitmapUtilities.Sqr(lab2.alpha - lab1.alpha) / Math.Exp(1.5);
				}
				else {
					errors[0] += Math.Abs(lab2.L - lab1.L);
					errors[1] += Math.Sqrt(BitmapUtilities.Sqr(lab2.A - lab1.A) + BitmapUtilities.Sqr(lab2.B - lab1.B));
				}
			}
			_objectives = errors;
			_fitness = -1f * (float) _objectives.Sum();
			fitnessMap[ratioKey] = _objectives;
			paletteMap[ratioKey] = palette;
		}
		
		public Bitmap QuantizeImage(bool dither) {
			var ratioKey = RatioKey;
			m_pq.Ratio = ratio;
			if (!paletteMap.TryGetValue(ratioKey, out Color[] palette))
				m_pq.Pnnquan(m_pixels, ref palette, ref _nMaxColors);
			m_pq.Palette = palette;
			return m_pq.QuantizeImage(m_pixels, _bitmapWidth, _nMaxColors, dither);
		}

		protected virtual void Dispose(bool disposing) {
			// check if already disposed
			if (!disposedValue) {
				if (disposing) {
					// free managed objects here
					m_pixels = null;
					fitnessMap.Clear();
					paletteMap.Clear();
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

		public double Ratio
		{
			get => ratio;
			set => ratio = Math.Min(Math.Max(value, minRatio), maxRatio);
		}

		public float Fitness {
			get => _fitness;
		}

		public PnnLABGAQuantizer Crossover(PnnLABGAQuantizer mother, int numberOfCrossoverPoints, float crossoverProbability) {
			var child = MakeNewFromPrototype();
			if (_random.Next(100) <= crossoverProbability)
				return child;
			
			double ratio = (Ratio + mother.Ratio) * .5;
			child.Ratio = ratio;
			child.CalculateFitness();
			return child;
		}

		public void Mutation(int mutationSize, float mutationProbability) {
			// check probability of mutation operation
			if (_random.Next(100) > mutationProbability)
				return;
			
			ratio = Randrange(minRatio, maxRatio);
			CalculateFitness();
		}

		public void ResizeConvertedObjectives(int numObj) {
			ConvertedObjectives = new double[numObj];
		}

		public PnnLABGAQuantizer MakeNewFromPrototype() {
			var child = new PnnLABGAQuantizer(m_pq, m_pixels, _bitmapWidth, _nMaxColors);
			child.Ratio = Randrange(minRatio, maxRatio);
			child.CalculateFitness();
			return child;
		}

		public Random Random {
			get => _random;
		}

		public static int MaxColors {
			get => _nMaxColors;
		}	

	}
}
