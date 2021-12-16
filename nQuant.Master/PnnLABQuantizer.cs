using nQuant.Master;
using System;
using System.Collections.Generic;
using System.Drawing;

namespace PnnQuant
{
    public class PnnLABQuantizer : PnnQuantizer
    {
        private double ratio = 1.0;
        private readonly Dictionary<int, CIELABConvertor.Lab> pixelMap = new();
        private sealed class Pnnbin
        {
            internal float ac, Lc, Ac, Bc;
            internal float cnt;
            internal int nn, fw, bk, tm, mtm;
            internal float err;
        }
        private void GetLab(int argb, out CIELABConvertor.Lab lab1)
        {
            if (!pixelMap.TryGetValue(argb, out lab1))
            {
                lab1 = CIELABConvertor.RGB2LAB(Color.FromArgb(argb));
                pixelMap[argb] = lab1;
            }
        }
        private void Find_nn(Pnnbin[] bins, int idx, bool texicab)
        {
            int nn = 0;
            var err = 1e100;

            var bin1 = bins[idx];
            var n1 = bin1.cnt;
            var lab1 = new CIELABConvertor.Lab
            {
                alpha = bin1.ac, L = bin1.Lc, A = bin1.Ac, B = bin1.Bc
            };
            for (int i = bin1.fw; i != 0; i = bins[i].fw)
            {
                float n2 = bins[i].cnt;
                float nerr2 = (n1 * n2) / (n1 + n2);
                if (nerr2 >= err)
                    continue;

                var lab2 = new CIELABConvertor.Lab
                {
                    alpha = bins[i].ac, L = bins[i].Lc, A = bins[i].Ac, B = bins[i].Bc
                };
                double alphaDiff = hasSemiTransparency ? (lab2.alpha - lab1.alpha) / Math.Exp(1.5) : 0;
                double nerr = nerr2 * BitmapUtilities.Sqr(alphaDiff);
                if (nerr >= err)
                    continue;

                if (hasSemiTransparency || !texicab)
                {
                    nerr += (1 - ratio) * nerr2 * BitmapUtilities.Sqr(lab2.L - lab1.L);
                    if (nerr >= err)
                        continue;

                    nerr += (1 - ratio) * nerr2 * BitmapUtilities.Sqr(lab2.A - lab1.A);
                    if (nerr >= err)
                        continue;

                    nerr += (1 - ratio) * nerr2 * BitmapUtilities.Sqr(lab2.B - lab1.B);
                }
                else
                {
                    nerr += (1 - ratio) * nerr2 * Math.Abs(lab2.L - lab1.L);
                    if (nerr >= err)
                        continue;

                    nerr += (1 - ratio) * nerr2 * Math.Sqrt(BitmapUtilities.Sqr(lab2.A - lab1.A) + BitmapUtilities.Sqr(lab2.B - lab1.B));
                }

                if (nerr >= err)
                    continue;

                var deltaL_prime_div_k_L_S_L = CIELABConvertor.L_prime_div_k_L_S_L(lab1, lab2);
                nerr += ratio * nerr2 * BitmapUtilities.Sqr(deltaL_prime_div_k_L_S_L);
                if (nerr >= err)
                    continue;

                var deltaC_prime_div_k_L_S_L = CIELABConvertor.C_prime_div_k_L_S_L(lab1, lab2, out var a1Prime, out var a2Prime, out var CPrime1, out var CPrime2);
                nerr += ratio * nerr2 * BitmapUtilities.Sqr(deltaC_prime_div_k_L_S_L);
                if (nerr >= err)
                    continue;

                var deltaH_prime_div_k_L_S_L = CIELABConvertor.H_prime_div_k_L_S_L(lab1, lab2, a1Prime, a2Prime, CPrime1, CPrime2, out var barCPrime, out var barhPrime);
                nerr += ratio * nerr2 * BitmapUtilities.Sqr(deltaH_prime_div_k_L_S_L);
                if (nerr >= err)
                    continue;

                nerr += ratio * nerr2 * CIELABConvertor.R_T(barCPrime, barhPrime, deltaC_prime_div_k_L_S_L, deltaH_prime_div_k_L_S_L);
                if (nerr >= err)
                    continue;

                err = (float)nerr;
                nn = i;
            }
            bin1.err = (float) err;
            bin1.nn = nn;
        }
        
        protected override QuanFn GetQuanFn(int nMaxColors, short quan_rt) {
            if (quan_rt > 0) {
                if (quan_rt > 1)
                    return cnt => (int) Math.Pow(cnt, 0.75);
                if (nMaxColors< 64)
                    return cnt => (int) Math.Sqrt(cnt);
                return cnt => (float) Math.Sqrt(cnt);
            }
            return cnt => cnt;
        }
        protected override void Pnnquan(int[] pixels, ref Color[] palettes, ref int nMaxColors, short quan_rt)
        {
            bool noBias = m_transparentPixelIndex >= 0 || hasSemiTransparency || nMaxColors < 64;
            if (noBias)
                PR = PG = PB = 1.0;
            else if (pixels.Length < BitmapUtilities.Sqr(512))
            {
                PR = 0.299; PG = 0.587; PB = 0.114;
            }

            var bins = new Pnnbin[ushort.MaxValue + 1];

            /* Build histogram */
            foreach (var pixel in pixels)
            {
                // !!! Can throw gamma correction in here, but what to do about perceptual
                // !!! nonuniformity then?
                int index = BitmapUtilities.GetARGBIndex(pixel, hasSemiTransparency, m_transparentPixelIndex > -1);
                GetLab(pixel, out var lab1);
                if (bins[index] == null)
                    bins[index] = new Pnnbin();
                bins[index].ac += (float)lab1.alpha;
                bins[index].Lc += (float)lab1.L;
                bins[index].Ac += (float)lab1.A;
                bins[index].Bc += (float)lab1.B;
                bins[index].cnt += 1.0f;
            }

            /* Cluster nonempty bins at one end of array */
            int maxbins = 0;
            for (int i = 0; i < bins.Length; ++i)
            {
                if (bins[i] == null)
                    continue;

                var d = 1.0f / bins[i].cnt;
                bins[i].ac *= d;
                bins[i].Lc *= d;
                bins[i].Ac *= d;
                bins[i].Bc *= d;

                bins[maxbins++] = bins[i];
            }

            var proportional = BitmapUtilities.Sqr(nMaxColors) / maxbins;
            if ((m_transparentPixelIndex > -1 || hasSemiTransparency) && nMaxColors < 32)
                quan_rt = -1;
            
            var weight = nMaxColors * 1.0 / maxbins;
            if (weight > .0015 && weight < .002)
                quan_rt = 2;

            var quanFn = GetQuanFn(nMaxColors, quan_rt);

            int j = 0;
            for (; j < maxbins - 1; ++j)
            {
                bins[j].fw = j + 1;
                bins[j + 1].bk = j;

                bins[j].cnt = quanFn(bins[j].cnt);                
            }
            bins[j].cnt = quanFn(bins[j].cnt);

            var texicab = proportional > .025;
            int h, l, l2;
            if (quan_rt != 0 && nMaxColors < 64)
            {
                if (proportional > .018 && proportional < .022)
                    ratio = Math.Min(1.0, proportional + weight * Math.Exp(3.872));
                else if (proportional > .1)
                    ratio = Math.Min(1.0, 1.0 - weight);
                else if (proportional > .03)
                    ratio = Math.Min(1.0, weight * Math.Exp(2.28));
                else
                    ratio = Math.Min(1.0, proportional - weight * Math.Exp(1.997));
            }
            else if (nMaxColors > 256)
                ratio = Math.Min(m_transparentPixelIndex >= 0 ? 0.0 : 1.0, 1 - 1.0 / proportional);
            else
                ratio = Math.Min(m_transparentPixelIndex >= 0 ? 0.0 : 1.0, 0.14 * Math.Exp(4.681 * proportional));

            if (quan_rt < 0)
                ratio = Math.Min(m_transparentPixelIndex >= 0 ? 0.0 : 1.0, weight * Math.Exp(1.997));

            /* Initialize nearest neighbors and build heap of them */
            var heap = new int[bins.Length + 1];
            for (int i = 0; i < maxbins; ++i)
            {
                Find_nn(bins, i, texicab);
                /* Push slot on heap */
                var err = bins[i].err;

                for (l = ++heap[0]; l > 1; l = l2)
                {
                    l2 = l >> 1;
                    if (bins[h = heap[l2]].err <= err)
                        break;
                    heap[l] = h;
                }
                heap[l] = i;
            }

            if (quan_rt > 0 && nMaxColors < 64 && (proportional < .023 || proportional > .05) && proportional < .1)
                ratio = Math.Min(1.0, proportional - weight * Math.Exp(2.347));

            /* Merge bins which increase error the least */
            int extbins = maxbins - nMaxColors;
            for (int i = 0; i < extbins;)
            {
                Pnnbin tb;
                /* Use heap to find which bins to merge */
                for (; ; )
                {
                    int b1 = heap[1];
                    tb = bins[b1]; /* One with least error */
                    /* Is stored error up to date? */
                    if ((tb.tm >= tb.mtm) && (bins[tb.nn].mtm <= tb.tm))
                        break;
                    if (tb.mtm == ushort.MaxValue) /* Deleted node */
                        b1 = heap[1] = heap[heap[0]--];
                    else /* Too old error value */
                    {
                        Find_nn(bins, b1, texicab && proportional < 1);
                        tb.tm = i;
                    }
                    /* Push slot down */
                    var err = bins[b1].err;
                    for (l = 1; (l2 = l + l) <= heap[0]; l = l2)
                    {
                        if ((l2 < heap[0]) && (bins[heap[l2]].err > bins[heap[l2 + 1]].err))
                            ++l2;
                        if (err <= bins[h = heap[l2]].err)
                            break;
                        heap[l] = h;
                    }
                    heap[l] = b1;
                }

                /* Do a merge */
                var nb = bins[tb.nn];
                var n1 = tb.cnt;
                var n2 = nb.cnt;
                var d = 1.0f / (float)(n1 + n2);
                tb.ac = d * (n1 * tb.ac + n2 * nb.ac);
                tb.Lc = d * (n1 * tb.Lc + n2 * nb.Lc);
                tb.Ac = d * (n1 * tb.Ac + n2 * nb.Ac);
                tb.Bc = d * (n1 * tb.Bc + n2 * nb.Bc);
                tb.cnt += nb.cnt;
                tb.mtm = ++i;

                /* Unchain deleted bin */
                bins[nb.bk].fw = nb.fw;
                bins[nb.fw].bk = nb.bk;
                nb.mtm = ushort.MaxValue;
            }

            /* Fill palette */
            if (extbins < 0)
                palettes = new Color[maxbins];

            int k = 0;
            for (int i = 0; ; ++k)
            {
                var lab1 = new CIELABConvertor.Lab
                {
                    alpha = Math.Round(bins[i].ac), L = bins[i].Lc, A = bins[i].Ac, B = bins[i].Bc
                };
                palettes[k] = CIELABConvertor.LAB2RGB(lab1);
                if (m_transparentPixelIndex >= 0 && lab1.alpha == 0)
                {
                    BitmapUtilities.Swap(ref palettes[0], ref palettes[k]);
                    palettes[0] = m_transparentColor;
                }

                if ((i = bins[i].fw) == 0)
                    break;
            }
	    
	        if (k < nMaxColors - 1)
            {
                nMaxColors = k + 1;
                Console.WriteLine("Maximum number of colors: " + palettes.Length);
            }
        }

        protected override ushort NearestColorIndex(Color[] palette, int pixel)
        {
            if (nearestMap.TryGetValue(pixel, out var k))
                return k;

            var c = Color.FromArgb(pixel);
            if (c.A <= alphaThreshold)
                return 0;

            double mindist = int.MaxValue;
            var nMaxColors = palette.Length;
            GetLab(pixel, out var lab1);

            for (int i = 0; i < nMaxColors; ++i)
            {
                var c2 = palette[i];

                var curdist = hasSemiTransparency ? Math.Abs(c2.A - c.A) / Math.Exp(0.75) : 0;
                if (curdist > mindist)
                    continue;

                GetLab(c2.ToArgb(), out var lab2);
                if (nMaxColors <= 4)
                {
                    curdist = BitmapUtilities.Sqr(c2.R - c.R) + BitmapUtilities.Sqr(c2.G - c.G) + BitmapUtilities.Sqr(c2.B - c.B);
                    if(hasSemiTransparency)
                        curdist += BitmapUtilities.Sqr(c2.A - c.A);
                }
                else if (nMaxColors > 32 || hasSemiTransparency)
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
            nearestMap[pixel] = k;
            return k;
        }
        protected override ushort ClosestColorIndex(Color[] palette, int pixel)
        {
            ushort k = 0;
	        var c = Color.FromArgb(pixel);
	        if (c.A <= alphaThreshold)
                return 0;
		
            if (!closestMap.TryGetValue(pixel, out var closest))
            {
                closest = new ushort[4];
                closest[2] = closest[3] = ushort.MaxValue;

                var nMaxColors = palette.Length;
                for (; k < nMaxColors; ++k)
                {
                    var c2 = palette[k];
                    var err = PR * BitmapUtilities.Sqr(c2.R - c.R) + PG * BitmapUtilities.Sqr(c2.G - c.G) + PB * BitmapUtilities.Sqr(c2.B - c.B);  
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
		    
                closestMap[pixel] = closest;
            }

            var MAX_ERR = palette.Length;
            if (closest[2] == 0 || (rand.Next(short.MaxValue) % (closest[3] + closest[2])) <= closest[3]) {
                if (closest[2] > MAX_ERR)
                    return NearestColorIndex(palette, pixel);
                return closest[0];
            }
            if (closest[3] > MAX_ERR)
                return NearestColorIndex(palette, pixel);
            return closest[1];
        }

        public override ushort DitherColorIndex(Color[] palette, int nMaxColors, int pixel)
        {
            if(hasSemiTransparency)
                return NearestColorIndex(palette, pixel);
            return ClosestColorIndex(palette, pixel);
        }

        protected override int[] Dither(int[] pixels, Color[] palettes, int nMaxColors, int width, int height, bool dither)
        {
            int[] qPixels;
            if (hasSemiTransparency)
                qPixels = GilbertCurve.Dither(width, height, pixels, palettes, this, 1.75f);
            else if (nMaxColors < 64 && nMaxColors > 32)
                qPixels = BitmapUtilities.Quantize_image(width, height, pixels, palettes, this, hasSemiTransparency, dither);
            else if (nMaxColors <= 32)
                qPixels = GilbertCurve.Dither(width, height, pixels, palettes, this, 1.5f);
            else
                qPixels = GilbertCurve.Dither(width, height, pixels, palettes, this);

            if (!dither)
            {
                var delta = BitmapUtilities.Sqr(nMaxColors) / pixelMap.Count;
                var weight = delta > 0.023 ? 1.0f : (float)(36.921 * delta + 0.906);
                return BlueNoise.Dither(width, height, pixels, palettes, this, qPixels, weight);
            }

            pixelMap.Clear();
            return qPixels;
        }

    }
}
