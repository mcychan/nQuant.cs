using System;
using System.Collections.Generic;
using System.Drawing;

namespace PnnQuant
{
    public class PnnLABQuantizer : PnnQuantizer
    {
        private double PR = .299, PG = .587, PB = .114;
        private double ratio = 1.0;
        private readonly Dictionary<int, CIELABConvertor.Lab> pixelMap = new Dictionary<int, CIELABConvertor.Lab>();
        private sealed class Pnnbin
        {
            internal float ac, Lc, Ac, Bc;
            internal int cnt;
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
        private void Find_nn(Pnnbin[] bins, int idx)
        {
            int nn = 0;
            var err = 1e100;

            var bin1 = bins[idx];
            int n1 = bin1.cnt;
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
                double alphaDiff = hasSemiTransparency ? Math.Abs(lab2.alpha - lab1.alpha) : 0;
                double nerr = nerr2 * Sqr(alphaDiff) / Math.Exp(1.5);
                if (nerr >= err)
                    continue;

                nerr += (1 - ratio) * nerr2 * Sqr(lab2.L - lab1.L);
                if (nerr >= err)
                    continue;

                nerr += (1 - ratio) * nerr2 * Sqr(lab2.A - lab1.A);
                if (nerr >= err)
                    continue;

                nerr += (1 - ratio) * nerr2 * Sqr(lab2.B - lab1.B);

                if (nerr >= err)
                    continue;

                var deltaL_prime_div_k_L_S_L = CIELABConvertor.L_prime_div_k_L_S_L(lab1, lab2);
                nerr += ratio * nerr2 * Sqr(deltaL_prime_div_k_L_S_L);
                if (nerr >= err)
                    continue;

                var deltaC_prime_div_k_L_S_L = CIELABConvertor.C_prime_div_k_L_S_L(lab1, lab2, out var a1Prime, out var a2Prime, out var CPrime1, out var CPrime2);
                nerr += ratio * nerr2 * Sqr(deltaC_prime_div_k_L_S_L);
                if (nerr >= err)
                    continue;

                var deltaH_prime_div_k_L_S_L = CIELABConvertor.H_prime_div_k_L_S_L(lab1, lab2, a1Prime, a2Prime, CPrime1, CPrime2, out var barCPrime, out var barhPrime);
                nerr += ratio * nerr2 * Sqr(deltaH_prime_div_k_L_S_L);
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
        protected override void Pnnquan(int[] pixels, Color[] palettes, int nMaxColors, short quan_sqrt)
        {
            if (hasSemiTransparency)
                PR = PG = PB = 1.0;

            var bins = new Pnnbin[65536];

            /* Build histogram */
            foreach (var pixel in pixels)
            {
                // !!! Can throw gamma correction in here, but what to do about perceptual
                // !!! nonuniformity then?
                var c = Color.FromArgb(pixel);

                int index = GetARGBIndex(pixel, hasSemiTransparency);
                GetLab(pixel, out var lab1);
                if (bins[index] == null)
                    bins[index] = new Pnnbin();
                bins[index].ac += (float)lab1.alpha;
                bins[index].Lc += (float)lab1.L;
                bins[index].Ac += (float)lab1.A;
                bins[index].Bc += (float)lab1.B;
                bins[index].cnt++;
            }

            /* Cluster nonempty bins at one end of array */
            int maxbins = 0;            
            
            for (int i = 0; i < bins.Length; ++i)
            {
                if (bins[i] == null)
                    continue;

                var d = 1.0f / (float)bins[i].cnt;
                bins[i].ac *= d;
                bins[i].Lc *= d;
                bins[i].Ac *= d;
                bins[i].Bc *= d;

                bins[maxbins++] = bins[i];
            }

            var proportional = Sqr(nMaxColors) / maxbins;
            if (nMaxColors < 16 || (hasSemiTransparency && nMaxColors < 32))
                quan_sqrt = -1;
            else if ((proportional < .022 || proportional > .5) && nMaxColors < 64)
                quan_sqrt = 0;

            if (quan_sqrt > 0)
                bins[0].cnt = (int)Math.Sqrt(bins[0].cnt);
            for (int i = 0; i < maxbins - 1; ++i)
            {
                bins[i].fw = i + 1;
                bins[i + 1].bk = i;

                if (quan_sqrt > 0)
                    bins[i + 1].cnt = (int) Math.Sqrt(bins[i + 1].cnt);
            }            

            int h, l, l2;
            if (quan_sqrt != 0 && nMaxColors < 64)
                ratio = Math.Min(1.0, proportional - nMaxColors * Math.Exp(4.172) / pixelMap.Count);
            else if (quan_sqrt > 0)
                ratio = Math.Min(1.0, Math.Pow(nMaxColors, 1.05) / pixelMap.Count);
            else
                ratio = Math.Min(1.0, Math.Pow(nMaxColors, 2.07) / maxbins);

            if (quan_sqrt < 0)
            {
                ratio += 0.45;
                ratio = Math.Min(1.0, ratio);
            }

            /* Initialize nearest neighbors and build heap of them */
            var heap = new int[bins.Length + 1];
            for (int i = 0; i < maxbins; ++i)
            {
                Find_nn(bins, i);
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

            if (quan_sqrt > 0 && nMaxColors < 64)
                ratio = Math.Min(1.0, proportional - nMaxColors * Math.Exp(4.12) / pixelMap.Count);

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
                    if (tb.mtm == 0xFFFF) /* Deleted node */
                        b1 = heap[1] = heap[heap[0]--];
                    else /* Too old error value */
                    {
                        Find_nn(bins, b1);
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
                nb.mtm = 0xFFFF;
            }

            /* Fill palette */
            int k = 0;
            for (int i = 0; ; ++k)
            {
                var lab1 = new CIELABConvertor.Lab
                {
                    alpha = Math.Round(bins[i].ac), L = bins[i].Lc, A = bins[i].Ac, B = bins[i].Bc
                };
                palettes[k] = CIELABConvertor.LAB2RGB(lab1);
                if (m_transparentPixelIndex >= 0 && palettes[k] == m_transparentColor)
                    Swap(ref palettes[0], ref palettes[k]);

                if ((i = bins[i].fw) == 0)
                    break;
            }
        }

        protected override ushort NearestColorIndex(Color[] palette, int nMaxColors, int argb)
        {
            if (nearestMap.TryGetValue(argb, out var k))
                return k;

            var c = Color.FromArgb(argb);
            if (c.A <= alphaThreshold)
                return 0;

            double mindist = ushort.MaxValue;
            GetLab(argb, out var lab1);

            for (int i = 0; i < nMaxColors; ++i)
            {
                var c2 = palette[i];

                var curdist = Sqr(c2.A - c.A);
                if (curdist > mindist)
                    continue;

                GetLab(c2.ToArgb(), out var lab2);
                if (nMaxColors > 32 || hasSemiTransparency)
                {
                    curdist += PR * Sqr(c2.R - c.R);
                    if (curdist > mindist)
                        continue;

                    curdist += PG * Sqr(c2.G - c.G);
                    if (curdist > mindist)
                        continue;

                    curdist += PB * Sqr(c2.B - c.B);
                    if (PB < 1)
                    {
                        if (curdist > mindist)
                            continue;

                        curdist += Sqr(lab2.B - lab1.B) / 2.0;
                    }
                }
                else
                {  
                    var deltaL_prime_div_k_L_S_L = CIELABConvertor.L_prime_div_k_L_S_L(lab1, lab2);
                    curdist += Sqr(deltaL_prime_div_k_L_S_L);
                    if (curdist > mindist)
                        continue;

                    var deltaC_prime_div_k_L_S_L = CIELABConvertor.C_prime_div_k_L_S_L(lab1, lab2, out var a1Prime, out var a2Prime, out var CPrime1, out var CPrime2);
                    curdist += Sqr(deltaC_prime_div_k_L_S_L);
                    if (curdist > mindist)
                        continue;

                    var deltaH_prime_div_k_L_S_L = CIELABConvertor.H_prime_div_k_L_S_L(lab1, lab2, a1Prime, a2Prime, CPrime1, CPrime2, out var barCPrime, out var barhPrime);
                    curdist += Sqr(deltaH_prime_div_k_L_S_L);
                    if (curdist > mindist)
                        continue;

                    curdist += CIELABConvertor.R_T(barCPrime, barhPrime, deltaC_prime_div_k_L_S_L, deltaH_prime_div_k_L_S_L);
                }

                if (curdist > mindist)
                    continue;
                mindist = curdist;
                k = (ushort)i;
            }
            nearestMap[argb] = k;
            return k;
        }
        protected override ushort ClosestColorIndex(Color[] palette, int nMaxColors, int pixel)
        {
            ushort k = 0;
            if (!closestMap.TryGetValue(pixel, out var closest))
            {
                closest = new ushort[5];
                closest[2] = closest[3] = ushort.MaxValue;
                GetLab(pixel, out var lab1);

                for (; k < nMaxColors; ++k)
                {
                    var c2 = palette[k];
                    GetLab(c2.ToArgb(), out var lab2);

                    closest[4] = (ushort) (Sqr(lab2.L - lab1.L) + Sqr(lab2.A - lab1.A) + Sqr(lab2.B - lab1.B));
                    if (closest[4] < closest[2])
                    {
                        closest[1] = closest[0];
                        closest[3] = closest[2];
                        closest[0] = (ushort)k;
                        closest[2] = closest[4];
                    }
                    else if (closest[4] < closest[3])
                    {
                        closest[1] = (ushort)k;
                        closest[3] = closest[4];
                    }
                }

                if (closest[3] == ushort.MaxValue)
                    closest[2] = 0;
            }

            if (closest[2] == 0 || (rand.Next(short.MaxValue) % (closest[3] + closest[2])) <= closest[3])
                k = closest[0];
            else
                k = closest[1];

            closestMap[pixel] = closest;
            return k;
        }

        protected override int[] Quantize_image(int[] pixels, Color[] palette, int nMaxColors, int width, int height, bool dither)
        {
            var qPixels = new int[width * height];
            int pixelIndex = 0;
            if (dither)
            {
                const short DJ = 4;
                const short BLOCK_SIZE = 256;
                const short DITHER_MAX = 16;
                int err_len = (width + 2) * DJ;
                var clamp = new short[DJ * BLOCK_SIZE];
                var limtb = new short[2 * BLOCK_SIZE];

                for (short i = 0; i < BLOCK_SIZE; ++i)
                {
                    clamp[i] = 0;
                    clamp[i + BLOCK_SIZE] = i;
                    clamp[i + BLOCK_SIZE * 2] = Byte.MaxValue;
                    clamp[i + BLOCK_SIZE * 3] = Byte.MaxValue;

                    limtb[i] = -DITHER_MAX;
                    limtb[i + BLOCK_SIZE] = DITHER_MAX;
                }
                for (short i = -DITHER_MAX; i <= DITHER_MAX; ++i)
                    limtb[i + BLOCK_SIZE] = i;

                bool noBias = hasSemiTransparency || nMaxColors < 64;
                int dir = 1;
                var row0 = new int[err_len];
                var row1 = new int[err_len];
                for (int i = 0; i < height; ++i)
                {
                    if (dir < 0)
                        pixelIndex += width - 1;

                    int cursor0 = DJ, cursor1 = width * DJ;
                    row1[cursor1] = row1[cursor1 + 1] = row1[cursor1 + 2] = row1[cursor1 + 3] = 0;
                    for (int j = 0; j < width; ++j)
                    {
                        var c = Color.FromArgb(pixels[pixelIndex]);
                        var ditherPixel = CalcDitherPixel(c, clamp, row0, cursor0, noBias);
                        int r_pix = ditherPixel[0];
                        int g_pix = ditherPixel[1];
                        int b_pix = ditherPixel[2];
                        int a_pix = ditherPixel[3];

                        var c1 = Color.FromArgb(a_pix, r_pix, g_pix, b_pix);
                        if (c.A == 0 && a_pix > 0)
                            qPixels[pixelIndex] = 0;
                        else
                            qPixels[pixelIndex] = noBias ? NearestColorIndex(palette, nMaxColors, c1.ToArgb()) : ClosestColorIndex(palette, nMaxColors, c1.ToArgb());

                        var c2 = palette[qPixels[pixelIndex]];
                        if (nMaxColors > 256)
                            qPixels[pixelIndex] = hasSemiTransparency ? c2.ToArgb() : GetARGBIndex(c2.ToArgb(), false);

                        r_pix = limtb[r_pix - c2.R + BLOCK_SIZE];
                        g_pix = limtb[g_pix - c2.G + BLOCK_SIZE];
                        b_pix = limtb[b_pix - c2.B + BLOCK_SIZE];
                        a_pix = limtb[a_pix - c2.A + BLOCK_SIZE];

                        int k = r_pix * 2;
                        row1[cursor1 - DJ] = r_pix;
                        row1[cursor1 + DJ] += (r_pix += k);
                        row1[cursor1] += (r_pix += k);
                        row0[cursor0 + DJ] += (r_pix += k);

                        k = g_pix * 2;
                        row1[cursor1 + 1 - DJ] = g_pix;
                        row1[cursor1 + 1 + DJ] += (g_pix += k);
                        row1[cursor1 + 1] += (g_pix += k);
                        row0[cursor0 + 1 + DJ] += (g_pix += k);

                        k = b_pix * 2;
                        row1[cursor1 + 2 - DJ] = b_pix;
                        row1[cursor1 + 2 + DJ] += (b_pix += k);
                        row1[cursor1 + 2] += (b_pix += k);
                        row0[cursor0 + 2 + DJ] += (b_pix += k);

                        k = a_pix * 2;
                        row1[cursor1 + 3 - DJ] = a_pix;
                        row1[cursor1 + 3 + DJ] += (a_pix += k);
                        row1[cursor1 + 3] += (a_pix += k);
                        row0[cursor0 + 3 + DJ] += (a_pix += k);

                        cursor0 += DJ;
                        cursor1 -= DJ;
                        pixelIndex += dir;
                    }
                    if ((i % 2) == 1)
                        pixelIndex += width + 1;

                    dir *= -1;
                    Swap(ref row0, ref row1);
                }
                return qPixels;
            }

            if (m_transparentPixelIndex >= 0 || nMaxColors < 64)
            {
                for (int i = 0; i < qPixels.Length; ++i)
                    qPixels[i] = NearestColorIndex(palette, nMaxColors, pixels[i]);
            }
            else
            {
                for (int i = 0; i < qPixels.Length; ++i)
                    qPixels[i] = ClosestColorIndex(palette, nMaxColors, pixels[i]);
            }

            return qPixels;
        }

    }
}