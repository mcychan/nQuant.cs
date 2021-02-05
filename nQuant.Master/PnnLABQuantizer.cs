using System;
using System.Collections.Generic;
using System.Drawing;
using System.Drawing.Imaging;
using System.Threading.Tasks;

namespace PnnQuant
{
    public class PnnLABQuantizer : PnnQuantizer
    {
        private double PR = .2126, PG = .7152, PB = .0722;
        private double ratio = 1.0;
        private Dictionary<int, CIELABConvertor.Lab> pixelMap = new Dictionary<int, CIELABConvertor.Lab>();
        private sealed class Pnnbin
        {
            internal double ac, Lc, Ac, Bc;
            internal int cnt;
            internal int nn, fw, bk, tm, mtm;
            internal double err;
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
            double err = 1e100;

            var bin1 = bins[idx];
            int n1 = bin1.cnt;
            CIELABConvertor.Lab lab1;
            lab1.alpha = bin1.ac; lab1.L = bin1.Lc; lab1.A = bin1.Ac; lab1.B = bin1.Bc;
            for (int i = bin1.fw; i != 0; i = bins[i].fw)
            {
                double n2 = bins[i].cnt;
                double nerr2 = (n1 * n2) / (n1 + n2);
                if (nerr2 >= err)
                    continue;

                CIELABConvertor.Lab lab2;
                lab2.alpha = bins[i].ac; lab2.L = bins[i].Lc; lab2.A = bins[i].Ac; lab2.B = bins[i].Bc;
                double alphaDiff = lab2.alpha - lab1.alpha;
                double nerr = nerr2 * Sqr(alphaDiff) * alphaDiff / 3.0;
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

                double deltaL_prime_div_k_L_S_L = CIELABConvertor.L_prime_div_k_L_S_L(lab1, lab2);
                nerr += ratio * nerr2 * Sqr(deltaL_prime_div_k_L_S_L);
                if (nerr >= err)
                    continue;

                double deltaC_prime_div_k_L_S_L = CIELABConvertor.C_prime_div_k_L_S_L(lab1, lab2, out double a1Prime, out double a2Prime, out double CPrime1, out double CPrime2);
                nerr += ratio * nerr2 * Sqr(deltaC_prime_div_k_L_S_L);
                if (nerr >= err)
                    continue;

                double deltaH_prime_div_k_L_S_L = CIELABConvertor.H_prime_div_k_L_S_L(lab1, lab2, a1Prime, a2Prime, CPrime1, CPrime2, out double barCPrime, out double barhPrime);
                nerr += ratio * nerr2 * Sqr(deltaH_prime_div_k_L_S_L);
                if (nerr >= err)
                    continue;

                nerr += ratio * nerr2 * CIELABConvertor.R_T(barCPrime, barhPrime, deltaC_prime_div_k_L_S_L, deltaH_prime_div_k_L_S_L);
                if (nerr >= err)
                    continue;

                err = nerr;
                nn = i;
            }
            bin1.err = err;
            bin1.nn = nn;
        }
        private void Pnnquan(int[] pixels, Color[] palettes, int nMaxColors, bool quan_sqrt)
        {
            var bins = new Pnnbin[65536];

            /* Build histogram */
            for (int i = 0; i < pixels.Length; ++i)
            {
                // !!! Can throw gamma correction in here, but what to do about perceptual
                // !!! nonuniformity then?
                int index = GetARGBIndex(pixels[i], hasSemiTransparency);
                GetLab(pixels[i], out CIELABConvertor.Lab lab1);
                if (bins[index] == null)
                    bins[index] = new Pnnbin();
                bins[index].ac += lab1.alpha;
                bins[index].Lc += lab1.L;
                bins[index].Ac += lab1.A;
                bins[index].Bc += lab1.B;
                bins[index].cnt++;
            }

            /* Cluster nonempty bins at one end of array */
            int maxbins = 0;
            var heap = new int[65537];
            for (int i = 0; i < bins.Length; ++i)
            {
                if (bins[i] == null)
                    continue;

                double d = 1.0 / (double) bins[i].cnt;
                bins[i].ac *= d;
                bins[i].Lc *= d;
                bins[i].Ac *= d;
                bins[i].Bc *= d;

                if (quan_sqrt)
                    bins[i].cnt = (int)Math.Sqrt(bins[i].cnt);
                bins[maxbins++] = bins[i];
            }

            for (int i = 0; i < maxbins - 1; i++)
            {
                bins[i].fw = (i + 1);
                bins[i + 1].bk = i;
            }
            // !!! Already zeroed out by calloc()
            //	bins[0].bk = bins[i].fw = 0;

            int h, l, l2;
            ratio = 0.0;
            /* Initialize nearest neighbors and build heap of them */
            for (int i = 0; i < maxbins; i++)
            {
                Find_nn(bins, i);
                /* Push slot on heap */
                double err = bins[i].err;

                for (l = ++heap[0]; l > 1; l = l2)
                {
                    l2 = l >> 1;
                    if (bins[h = heap[l2]].err <= err)
                        break;
                    heap[l] = h;
                }
                heap[l] = i;
            }

			ratio = Math.Min(1.0, Math.Pow(nMaxColors, 2.25) / pixelMap.Count);
            /* Merge bins which increase error the least */
            int extbins = maxbins - nMaxColors;
            for (int i = 0; i < extbins;)
            {
                Pnnbin tb = null;
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
                    double err = bins[b1].err;
                    for (l = 1; (l2 = l + l) <= heap[0]; l = l2)
                    {
                        if ((l2 < heap[0]) && (bins[heap[l2]].err > bins[heap[l2 + 1]].err))
                            l2++;
                        if (err <= bins[h = heap[l2]].err)
                            break;
                        heap[l] = h;
                    }
                    heap[l] = b1;
                }

                /* Do a merge */
                var nb = bins[tb.nn];
                double n1 = tb.cnt;
                double n2 = nb.cnt;
                double d = 1.0 / (n1 + n2);
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
                CIELABConvertor.Lab lab1;
                lab1.alpha = (int)Math.Round(bins[i].ac);
                lab1.L = bins[i].Lc; lab1.A = bins[i].Ac; lab1.B = bins[i].Bc;
                palettes[k] = CIELABConvertor.LAB2RGB(lab1);
                if (m_transparentPixelIndex >= 0 && palettes[k] == m_transparentColor)
                    Swap(ref palettes[0], ref palettes[k]);

                if ((i = bins[i].fw) == 0)
                    break;
            }
        }

        protected override ushort NearestColorIndex(Color[] palette, int nMaxColors, int argb)
        {
            ushort k = 0;
            Color c = Color.FromArgb(argb);

            double mindist = 1e100;
            GetLab(argb, out CIELABConvertor.Lab lab1);

            for (int i = 0; i < nMaxColors; ++i)
            {
                Color c2 = palette[i];

                double curdist = Sqr(c2.A - c.A);
                if (curdist > mindist)
                    continue;

                if (nMaxColors > 32)
                {
                    curdist += PR * Sqr(c2.R - c.R);
                    if (curdist > mindist)
                        continue;

                    curdist += PG * Sqr(c2.G - c.G);
                    if (curdist > mindist)
                        continue;

                    curdist += PB * Sqr(c2.B - c.B);
                }
                else
                {
                    GetLab(c2.ToArgb(), out CIELABConvertor.Lab lab2);

                    double deltaL_prime_div_k_L_S_L = CIELABConvertor.L_prime_div_k_L_S_L(lab1, lab2);
                    curdist += Sqr(deltaL_prime_div_k_L_S_L);
                    if (curdist > mindist)
                        continue;

                    double deltaC_prime_div_k_L_S_L = CIELABConvertor.C_prime_div_k_L_S_L(lab1, lab2, out double a1Prime, out double a2Prime, out double CPrime1, out double CPrime2);
                    curdist += Sqr(deltaC_prime_div_k_L_S_L);
                    if (curdist > mindist)
                        continue;

                    double deltaH_prime_div_k_L_S_L = CIELABConvertor.H_prime_div_k_L_S_L(lab1, lab2, a1Prime, a2Prime, CPrime1, CPrime2, out double barCPrime, out double barhPrime);
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
            return k;
        }
        protected override ushort ClosestColorIndex(Color[] palette, int nMaxColors, int pixel)
        {
            ushort k = 0;
            if (!closestMap.TryGetValue(pixel, out ushort[] closest))
            {
                closest = new ushort[5];
                closest[2] = closest[3] = ushort.MaxValue;
                GetLab(pixel, out CIELABConvertor.Lab lab1);

                for (; k < nMaxColors; k++)
                {
                    Color c2 = palette[k];
                    GetLab(c2.ToArgb(), out CIELABConvertor.Lab lab2);

                    closest[4] = (ushort)(Sqr(lab2.alpha - lab1.alpha) + CIELABConvertor.CIEDE2000(lab2, lab1));
                    //closest[4] = (short) (Math.abs(lab2.alpha - lab1.alpha) + Math.abs(lab2.L - lab1.L) + Math.abs(lab2.A - lab1.A) + Math.abs(lab2.B - lab1.B));
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

        public override Bitmap QuantizeImage(Bitmap source, PixelFormat pixelFormat, int nMaxColors, bool dither)
        {
            int bitmapWidth = source.Width;
            int bitmapHeight = source.Height;

            Bitmap dest = new Bitmap(bitmapWidth, bitmapHeight, pixelFormat);
            if (!IsValidFormat(pixelFormat, nMaxColors))
                return dest;

            var pixels = new int[bitmapWidth * bitmapHeight];
            if (!GrabPixels(source, pixels))
                return dest;
            
            var palette = dest.Palette;
            var palettes = palette.Entries;
            if (palettes.Length != nMaxColors)
                palettes = new Color[nMaxColors];
            if (nMaxColors > 256)
                dither = true;

            if (hasSemiTransparency)
                PR = PG = PB = 1;

            bool quan_sqrt = true;
            if (nMaxColors > 2)
                Pnnquan(pixels, palettes, nMaxColors, quan_sqrt);
            else
            {
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
            }

            var qPixels = Quantize_image(pixels, palettes, nMaxColors, bitmapWidth, bitmapHeight, dither);
            if (m_transparentPixelIndex >= 0)
            {
                var k = qPixels[m_transparentPixelIndex];
                if (nMaxColors > 2)
                    palettes[k] = m_transparentColor;
                else if (palettes[k] != m_transparentColor)
                    Swap(ref palettes[0], ref palettes[1]);
            }
            pixelMap.Clear();
            closestMap.Clear();

            if (nMaxColors > 256)
                return ProcessImagePixels(dest, qPixels, hasSemiTransparency, m_transparentPixelIndex);

            for (int i = 0; i < palettes.Length; ++i)
                palette.Entries[i] = palettes[i];
            return ProcessImagePixels(dest, palette, qPixels);
        }

    }
}