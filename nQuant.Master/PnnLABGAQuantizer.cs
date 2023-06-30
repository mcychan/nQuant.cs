using nQuant.Master;
using nQuant.Master.Ga;
using System;
using System.Collections.Concurrent;
using System.Drawing;
using System.Linq;
using System.Text;

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
        public double ratioX = 0, ratioY = 0;
        private double[] _objectives;
        private PnnLABQuantizer m_pq;
        public double[] ConvertedObjectives { get; private set; }
        public double[] Objectives { get => _objectives; }

        private readonly Random _random;

        private static int _dp = 1, _nMaxColors = 256;
        private static double minRatio = 0, maxRatio = 1.0;
        private static int[] m_pixels;
        private static int _bitmapWidth;

        private static readonly ConcurrentDictionary<string, double[]> _fitnessMap = new();

        public PnnLABGAQuantizer(PnnLABQuantizer pq, Bitmap source, int nMaxColors)
        {
            // increment value when criteria violation occurs
            _objectives = new double[4];
            _bitmapWidth = source.Width;
            _random = new Random(_bitmapWidth * source.Height);

            m_pq = new PnnLABQuantizer(pq);
            if (pq.IsGA)
                return;

            _nMaxColors = nMaxColors;

            var hasSemiTransparency = false;
            m_pixels = m_pq.GrabPixels(source, _nMaxColors, ref hasSemiTransparency);
            minRatio = (hasSemiTransparency || nMaxColors < 64) ? .01 : .85;
            maxRatio = Math.Min(1.0, nMaxColors / ((nMaxColors < 64) ? 500.0 : 50.0));
            _dp = maxRatio < .1 ? 10000 : 100;
        }

        private PnnLABGAQuantizer(PnnLABQuantizer pq, int[] pixels, int bitmapWidth, int nMaxColors)
        {
            m_pq = new PnnLABQuantizer(pq);
            m_pixels = pixels;
            _bitmapWidth = bitmapWidth;
            _random = new Random(pixels.Length);
            _nMaxColors = nMaxColors;
        }

        private string RatioKey
        {
            get
            {
                var sb = new StringBuilder();
                sb.Append((int) (ratioX * _dp)).Append(";");
                sb.Append((int) (ratioY * _dp));
                return sb.ToString();
            }
        }

        private void CalculateFitness()
        {
            var ratioKey = RatioKey;
            if (_fitnessMap.TryGetValue(ratioKey, out _objectives))
            {
                _fitness = -1f * (float)_objectives.Sum();
                return;
            }

            _objectives = new double[4];
            m_pq.SetRatio(ratioX, ratioY);
            var palette = new Color[_nMaxColors];
            m_pq.Pnnquan(m_pixels, ref palette, ref _nMaxColors);
            m_pq.Palette = palette;

            int threshold = maxRatio < .1 ? -64 : -112;
            var errors = new double[_objectives.Length];
            for (int i = 0; i < m_pixels.Length; ++i)
            {
                if (BlueNoise.RAW_BLUE_NOISE[i & 4095] > threshold)
                    continue;

                m_pq.GetLab(m_pixels[i], out CIELABConvertor.Lab lab1);
                var qPixelIndex = m_pq.NearestColorIndex(palette, m_pixels[i], i);
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
            _objectives = errors;
            _fitness = -1f * (float)_objectives.Sum();
            _fitnessMap[ratioKey] = _objectives;
        }

        public Bitmap QuantizeImage(bool dither)
        {
            m_pq.SetRatio(ratioX, ratioY);
            var palette = new Color[_nMaxColors];
            m_pq.Pnnquan(m_pixels, ref palette, ref _nMaxColors);
            m_pq.Palette = palette;
            return m_pq.QuantizeImage(m_pixels, _bitmapWidth, _nMaxColors, dither);
        }

        protected virtual void Dispose(bool disposing)
        {
            // check if already disposed
            if (!disposedValue)
            {
                if (disposing)
                {
                    // free managed objects here
                    m_pixels = null;
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
            this.ratioX = Math.Min(Math.Max(ratioX, minRatio), maxRatio);
            this.ratioY = Math.Min(Math.Max(ratioY, minRatio), maxRatio);
        }

        public float Fitness
        {
            get => _fitness;
        }

        public PnnLABGAQuantizer Crossover(PnnLABGAQuantizer mother, int numberOfCrossoverPoints, float crossoverProbability)
        {
            var child = MakeNewFromPrototype();
            if (_random.Next(100) <= crossoverProbability)
                return child;

            var ratioX = Math.Sqrt(this.ratioX * mother.Ratios[1]);
            var ratioY = Math.Sqrt(this.ratioY * mother.Ratios[0]);
            child.SetRatio(ratioX, ratioY);
            child.CalculateFitness();
            return child;
        }

        public void Mutation(int mutationSize, float mutationProbability)
        {
            // check probability of mutation operation
            if (_random.Next(100) > mutationProbability)
                return;

            var ratioX = this.ratioX;
            var ratioY = this.ratioY;
            if (_random.NextDouble() > .5)
                ratioX = .5 * (ratioX + Randrange(minRatio, maxRatio));
            else
                ratioY = .5 * (ratioY + Randrange(minRatio, maxRatio));

            SetRatio(ratioX, ratioY);
            CalculateFitness();
        }

        public void ResizeConvertedObjectives(int numObj)
        {
            ConvertedObjectives = new double[numObj];
        }

        public PnnLABGAQuantizer MakeNewFromPrototype()
        {
            var child = new PnnLABGAQuantizer(m_pq, m_pixels, _bitmapWidth, _nMaxColors);
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

        public static int MaxColors
        {
            get => _nMaxColors;
        }

    }
}
