﻿using System;
using System.Collections.Generic;
using System.Drawing;
using System.Drawing.Imaging;

namespace PnnQuant
{
    public class PnnLABQuantizer : PnnQuantizer
    {
        private Dictionary<Color, CIELABConvertor.Lab> pixelMap = new Dictionary<Color, CIELABConvertor.Lab>();	

	    private class Pnnbin {
		    public double ac = 0, Lc = 0, Ac = 0, Bc = 0, err = 0;
		    public int cnt = 0;
		    public int nn, fw, bk, tm, mtm;
	    }
	
	    private CIELABConvertor.Lab getLab(Color c)
	    {
		    CIELABConvertor.Lab lab1;
		    if (!pixelMap.TryGetValue(c, out lab1)) {
			    lab1 = CIELABConvertor.RGB2LAB(c);
			    pixelMap[c] = lab1;
		    }
		    return lab1;
	    }

	    private void find_nn(Pnnbin[] bins, int idx)
	    {
		    int nn = 0;
		    double err = 1e100;

		    Pnnbin bin1 = bins[idx];
		    int n1 = bin1.cnt;
		    double wa = bin1.ac;
		    double wL = bin1.Lc;
		    double wA = bin1.Ac;
		    double wB = bin1.Bc;
		    for (int i = bin1.fw; i != 0; i = bins[i].fw) {
			    double nerr = Math.Pow((bins[i].ac - wa), 2) + Math.Pow((bins[i].Lc - wL), 2) + Math.Pow((bins[i].Ac - wA), 2) + Math.Pow((bins[i].Bc - wB), 2);
			    double n2 = bins[i].cnt;
			    nerr *= (n1 * n2) / (n1 + n2);
			    if (nerr >= err)
				    continue;
			    err = nerr;
			    nn = i;
		    }
		    bin1.err = err;
		    bin1.nn = nn;
	    }

	    private int pnnquan(Color[] pixels, Pnnbin[] bins, ColorPalette palette, int nMaxColors)
	    {
		    int[] heap = new int[65537];
		    double err, n1, n2;
		    int l, l2, h, b1, maxbins, extbins;

		    /* Build histogram */
		    foreach (var c in pixels) {
			    // !!! Can throw gamma correction in here, but what to do about perceptual
			    // !!! nonuniformity then?
			    int index = getARGBIndex(c);
			    CIELABConvertor.Lab lab1 = getLab(c);
			    if(bins[index] == null)
				    bins[index] = new Pnnbin();
			    Pnnbin tb = bins[index];
			    tb.ac += c.A;
			    tb.Lc += lab1.L;
			    tb.Ac += lab1.A;
			    tb.Bc += lab1.B;
			    tb.cnt++;
		    }

		    /* Cluster nonempty bins at one end of array */
		    maxbins = 0;

		    for (int i = 0; i < 65536; ++i) {
			    if (bins[i] == null)
				    continue;

			    double d = 1.0 / (double)bins[i].cnt;
			    bins[i].ac *= d;
			    bins[i].Lc *= d;
			    bins[i].Ac *= d;
			    bins[i].Bc *= d;
			    bins[maxbins++] = bins[i];
		    }

		    for (int i = 0; i < maxbins - 1; i++) {
			    bins[i].fw = (i + 1);
			    bins[i + 1].bk = i;
		    }
		    // !!! Already zeroed out by calloc()
		    //	bins[0].bk = bins[i].fw = 0;

		    /* Initialize nearest neighbors and build heap of them */
		    for (int i = 0; i < maxbins; i++) {
			    find_nn(bins, i);
			    /* Push slot on heap */
			    err = bins[i].err;
			    for (l = ++heap[0]; l > 1; l = l2) {
				    l2 = l >> 1;
				    if (bins[h = heap[l2]].err <= err)
					    break;
				    heap[l] = h;
			    }
			    heap[l] = i;
		    }

		    /* Merge bins which increase error the least */
		    extbins = maxbins - nMaxColors;
		    for (int i = 0; i < extbins; ) {
			    /* Use heap to find which bins to merge */
			    for (;;) {
				    Pnnbin tb = bins[b1 = heap[1]]; /* One with least error */
											   /* Is stored error up to date? */
				    if ((tb.tm >= tb.mtm) && (bins[tb.nn].mtm <= tb.tm))
					    break;
				    if (tb.mtm == 0xFFFF) /* Deleted node */
					    b1 = heap[1] = heap[heap[0]--];
				    else /* Too old error value */
				    {
					    find_nn(bins, b1);
					    tb.tm = i;
				    }
				    /* Push slot down */
				    err = bins[b1].err;
				    for (l = 1; (l2 = l + l) <= heap[0]; l = l2) {
					    if ((l2 < heap[0]) && (bins[heap[l2]].err > bins[heap[l2 + 1]].err))
						    l2++;
					    if (err <= bins[h = heap[l2]].err)
						    break;
					    heap[l] = h;
				    }
				    heap[l] = b1;
			    }

			    /* Do a merge */
			    Pnnbin tb1 = bins[b1];
			    Pnnbin nb = bins[tb1.nn];
			    n1 = tb1.cnt;
			    n2 = nb.cnt;
			    double d = 1.0 / (n1 + n2);
			    tb1.ac = d * (n1 * tb1.ac + n2 * nb.ac);
			    tb1.Lc = d * (n1 * tb1.Lc + n2 * nb.Lc);
			    tb1.Ac = d * (n1 * tb1.Ac + n2 * nb.Ac);
			    tb1.Bc = d * (n1 * tb1.Bc + n2 * nb.Bc);
			    tb1.cnt += nb.cnt;
			    tb1.mtm = ++i;

			    /* Unchain deleted bin */
			    bins[nb.bk].fw = nb.fw;
			    bins[nb.fw].bk = nb.bk;
			    nb.mtm = 0xFFFF;
		    }

		    /* Fill palette */
		    short k = 0;
		    for (int i = 0;; ++k) {
			    CIELABConvertor.Lab lab1 = new CIELABConvertor.Lab();
			    lab1.alpha = (int) Math.Round(bins[i].ac);
			    lab1.L = bins[i].Lc; lab1.A = bins[i].Ac; lab1.B = bins[i].Bc;
			    palette.Entries[k] = CIELABConvertor.LAB2RGB(lab1);
			    if (hasTransparency && palette.Entries[k] == m_transparentColor) {
                    Color temp = palette.Entries[0];
				    palette.Entries[0] = palette.Entries[k];
                    palette.Entries[k] = temp;
                }

			    if ((i = bins[i].fw) == 0)
				    break;
		    }

		    return 0;
	    }

	    private int nearestColorIndex(ColorPalette palette, Color c)
	    {
		    int k = 0;
		    double mindist = int.MaxValue;
		    CIELABConvertor.Lab lab1 = getLab(c);

            var nMaxColors = palette.Entries.Length;
            for (int i = 0; i < nMaxColors; ++i) {
                Color c2 = palette.Entries[i];
			    CIELABConvertor.Lab lab2 = getLab(c2);
			
			    double curdist = Math.Pow(c2.A - c.A, 2.0);
			    if (curdist > mindist)
				    continue;

                if (nMaxColors < 256) {
				    double deltaL_prime_div_k_L_S_L = CIELABConvertor.L_prime_div_k_L_S_L(lab1, lab2);
				    curdist += Math.Pow(deltaL_prime_div_k_L_S_L, 2.0);
				    if (curdist > mindist)
					    continue;

				    double a1Prime, a2Prime, CPrime1, CPrime2;
				    double deltaC_prime_div_k_L_S_L = CIELABConvertor.C_prime_div_k_L_S_L(lab1, lab2, out a1Prime, out a2Prime, out CPrime1, out CPrime2);
				    curdist += Math.Pow(deltaC_prime_div_k_L_S_L, 2.0);
				    if (curdist > mindist)
					    continue;

				    double barCPrime, barhPrime;
				    double deltaH_prime_div_k_L_S_L = CIELABConvertor.H_prime_div_k_L_S_L(lab1, lab2, a1Prime, a2Prime, CPrime1, CPrime2, out barCPrime, out barhPrime);
				    curdist += Math.Pow(deltaH_prime_div_k_L_S_L, 2.0);
				    if (curdist > mindist)
					    continue;

				    curdist += CIELABConvertor.R_T(barCPrime, barhPrime, deltaC_prime_div_k_L_S_L, deltaH_prime_div_k_L_S_L);
				    if (curdist > mindist)
					    continue;
			    }
			    else {
				    curdist += Math.Pow(lab2.L - lab1.L, 2.0);
				    if (curdist > mindist)
					    continue;

				    curdist += Math.Pow(lab2.A - lab1.A, 2.0);
				    if (curdist > mindist)
					    continue;

				    curdist += Math.Pow(lab2.B - lab1.B, 2.0);
				    if (curdist > mindist)
					    continue;
			    }

			    mindist = curdist;
			    k = i;
		    }
		    return k;
	    }

        private int closestColorIndex(ColorPalette palette, Color c)
	    {
            int k = 0;
            ushort[] closest;
            if (!closestMap.TryGetValue(c, out closest))
            {
                closest = new ushort[5];
                closest[2] = closest[3] = ushort.MaxValue;
                CIELABConvertor.Lab lab1 = getLab(c);

                var nMaxColors = palette.Entries.Length;
                for (; k < nMaxColors; k++)
                {
                    Color c2 = palette.Entries[k];
				    CIELABConvertor.Lab lab2 = getLab(c2);

                    closest[4] = (ushort)(Math.Pow(lab2.alpha - lab1.alpha, 2) + CIELABConvertor.CIEDE2000(lab2, lab1));
                    //closest[4] = (ushort) (Math.abs(lab2.alpha - lab1.alpha) + Math.abs(lab2.L - lab1.L) + Math.abs(lab2.A - lab1.A) + Math.abs(lab2.B - lab1.B));
				if (closest[4] < closest[2]) {
					closest[1] = closest[0];
					closest[3] = closest[2];
                    closest[0] = (ushort)k;
					closest[2] = closest[4];
				}
				else if (closest[4] < closest[3]) {
                    closest[1] = (ushort) k;
					closest[3] = closest[4];
				}
			}

                if (closest[3] == ushort.MaxValue)
                    closest[2] = 0;
            }

            Random rand = new Random();
            if (closest[2] == 0 || (rand.Next(short.MaxValue) % (closest[3] + closest[2])) <= closest[3])
                k = closest[0];
            else
                k = closest[1];

            closestMap[c] = closest;
            return k;
	    }

        private bool quantize_image(Color[] pixels, ColorPalette palette, int[] qPixels, int width, int height, bool dither)
        {
            var nMaxColors = palette.Entries.Length;

            int pixelIndex = 0;
            if (dither)
            {
                bool odd_scanline = false;
                int[] row0, row1;
                int a_pix, r_pix, g_pix, b_pix, dir, k;
                const int DJ = 4;
                const int DITHER_MAX = 20;
                int err_len = (width + 2) * DJ;
                int[] clamp = new int[DJ * 256];
                int[] limtb = new int[512];
                int[] erowerr = new int[err_len];
                int[] orowerr = new int[err_len];
                int[] lookup = new int[65536];

                for (int i = 0; i < 256; i++)
                {
                    clamp[i] = 0;
                    clamp[i + 256] = (short)i;
                    clamp[i + 512] = Byte.MaxValue;
                    clamp[i + 768] = Byte.MaxValue;

                    limtb[i] = -DITHER_MAX;
                    limtb[i + 256] = DITHER_MAX;
                }
                for (int i = -DITHER_MAX; i <= DITHER_MAX; i++)
                    limtb[i + 256] = i;

                for (int i = 0; i < height; i++)
                {
                    if (odd_scanline)
                    {
                        dir = -1;
                        pixelIndex += (int)(width - 1);
                        row0 = orowerr;
                        row1 = erowerr;
                    }
                    else
                    {
                        dir = 1;
                        row0 = erowerr;
                        row1 = orowerr;
                    }

                    int cursor0 = DJ, cursor1 = width * DJ;
                    row1[cursor1] = row1[cursor1 + 1] = row1[cursor1 + 2] = row1[cursor1 + 3] = 0;
                    for (int j = 0; j < width; j++)
                    {
                        Color c = pixels[pixelIndex];
                        r_pix = clamp[((row0[cursor0] + 0x1008) >> 4) + c.R];
                        g_pix = clamp[((row0[cursor0 + 1] + 0x1008) >> 4) + c.G];
                        b_pix = clamp[((row0[cursor0 + 2] + 0x1008) >> 4) + c.B];
                        a_pix = clamp[((row0[cursor0 + 3] + 0x1008) >> 4) + c.A];

                        Color c1 = Color.FromArgb(a_pix, r_pix, g_pix, b_pix);
                        int offset = getARGBIndex(c1);
                        if (lookup[offset] == 0)
                            lookup[offset] = nearestColorIndex(palette, c1) + 1;
                        qPixels[pixelIndex] = lookup[offset] - 1;

                        Color c2 = palette.Entries[qPixels[pixelIndex]];
                        if (c2.A < Byte.MaxValue && c.A == Byte.MaxValue)
                        {
                            lookup[offset] = nearestColorIndex(palette, pixels[pixelIndex]) + 1;
                            qPixels[pixelIndex] = lookup[offset] - 1;
                        }

                        r_pix = limtb[r_pix - c2.R + 256];
                        g_pix = limtb[g_pix - c2.G + 256];
                        b_pix = limtb[b_pix - c2.B + 256];
                        a_pix = limtb[a_pix - c2.A + 256];

                        k = r_pix * 2;
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

                    odd_scanline = !odd_scanline;
                }
                return true;
            }

            if (hasTransparency || nMaxColors < 256)
            {
                for (int i = 0; i < qPixels.Length; i++)
                    qPixels[i] = nearestColorIndex(palette, pixels[i]);
            }
            else
            {
                for (int i = 0; i < qPixels.Length; i++)
                    qPixels[i] = closestColorIndex(palette, pixels[i]);
            }

            return true;
        }

        public override bool QuantizeImage(Bitmap source, Bitmap dest, int nMaxColors, bool dither)
        {
            int bitDepth = Image.GetPixelFormatSize(source.PixelFormat);
            int bitmapWidth = source.Width;
            int bitmapHeight = source.Height;

            hasTransparency = hasSemiTransparency = false;
            int pixelIndex = 0;
            var pixels = new Color[bitmapWidth * bitmapHeight];
            if (bitDepth <= 16)
            {
                for (int y = 0; y < bitmapHeight; y++)
                {
                    for (int x = 0; x < bitmapWidth; x++)
                    {
                        Color color = source.GetPixel(x, y);
                        if (color.A < Byte.MaxValue)
                        {
                            hasSemiTransparency = true;
                            if (color.A == 0)
                            {
                                hasTransparency = true;
                                m_transparentColor = color;
                            }
                        }
                        pixels[pixelIndex++] = color;
                    }
                }
            }

            // Lock bits on 3x8 source bitmap
            else
            {
                BitmapData data = source.LockBits(new Rectangle(0, 0, bitmapWidth, bitmapHeight), ImageLockMode.ReadOnly, source.PixelFormat);
                int strideSource;

                unsafe
                {
                    var pRowSource = (byte*)data.Scan0;

                    // Compensate for possible negative stride
                    if (data.Stride > 0)
                        strideSource = data.Stride;
                    else
                    {
                        pRowSource += bitmapHeight * data.Stride;
                        strideSource = -data.Stride;
                    }

                    // First loop: gather color information
                    for (uint y = 0; y < bitmapHeight; y++)
                    {
                        var pPixelSource = pRowSource;
                        // For each row...
                        for (uint x = 0; x < bitmapWidth; x++)
                        {	// ...for each pixel...
                            byte pixelBlue = *pPixelSource++;
                            byte pixelGreen = *pPixelSource++;
                            byte pixelRed = *pPixelSource++;
                            byte pixelAlpha = bitDepth < 32 ? Byte.MaxValue : *pPixelSource++;

                            var argb = Color.FromArgb(pixelAlpha, pixelRed, pixelGreen, pixelBlue);
                            if (pixelAlpha < Byte.MaxValue)
                            {
                                hasSemiTransparency = true;
                                if (pixelAlpha == 0)
                                {
                                    hasTransparency = true;
                                    m_transparentColor = argb;
                                }
                            }
                            pixels[pixelIndex++] = argb;
                        }
                        pRowSource += strideSource;
                    }
                }

                source.UnlockBits(data);
            }             

            var qPixels = new int[bitmapWidth * bitmapHeight];
            if (nMaxColors > 256)
            {
                hasSemiTransparency = false;
                quantize_image(pixels, qPixels, bitmapWidth, bitmapHeight);
                return ProcessImagePixels(dest, qPixels);
            }

            var bins = new Pnnbin[65536];
            var palette = dest.Palette;
            if (nMaxColors > 2)
                pnnquan(pixels, bins, palette, nMaxColors);
            else
            {
                if (hasTransparency)
                {
                    palette.Entries[0] = Color.Transparent;
                    palette.Entries[1] = Color.Black;
                }
                else
                {
                    palette.Entries[0] = Color.Black;
                    palette.Entries[1] = Color.White;
                }
            }

            quantize_image(pixels, palette, qPixels, bitmapWidth, bitmapHeight, dither);
            pixelMap.Clear();
            closestMap.Clear();

            return ProcessImagePixels(dest, palette, qPixels);
        }

    }
}