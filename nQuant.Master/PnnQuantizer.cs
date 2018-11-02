using System;
using System.Collections.Generic;
using System.Drawing;
using System.Drawing.Imaging;

/* Fast pairwise nearest neighbor based algorithm for multilevel thresholding
Copyright (C) 2004-2016 Mark Tyler and Dmitry Groshev
Copyright (c) 2018 Miller Cy Chan
* error measure; time used is proportional to number of bins squared - WJ */

namespace PnnQuant
{
    public class PnnQuantizer
    {
	    protected bool hasTransparency = false, hasSemiTransparency = false;
        protected Color m_transparentColor = Color.Transparent;
        protected Random rand = new Random();
        protected Dictionary<int, ushort[]> closestMap = new Dictionary<int, ushort[]>();

        private sealed class Pnnbin
        {
		    public float ac, rc, gc, bc;            
		    public int cnt;
            public int nn, fw, bk, tm, mtm;
            public double err;
	    }

        protected int getARGBIndex(int argb)
	    {
            Color c = Color.FromArgb(argb);
		    if(hasSemiTransparency)
			    return (c.A & 0xF0) << 8 | (c.R & 0xF0) << 4 | (c.G & 0xF0) | (c.B >> 4);
		    if (hasTransparency)
			    return (c.A & 0x80) << 8 | (c.R & 0xF8) << 7 | (c.G & 0xF8) << 2 | (c.B >> 3);
		    return (c.R & 0xF8) << 8 | (c.G & 0xFC) << 3 | (c.B >> 3);
	    }

	    private void find_nn(Pnnbin[] bins, int idx)
	    {
		    int nn = 0;
		    double err = 1e100;

		    var bin1 = bins[idx];
		    var n1 = bin1.cnt;
		    var wa = bin1.ac;
		    var wr = bin1.rc;
		    var wg = bin1.gc;
		    var wb = bin1.bc;
		    for (int i = bin1.fw; i != 0; i = bins[i].fw) {
			    double nerr, n2;

			    nerr = Math.Pow((bins[i].ac - wa), 2) + Math.Pow((bins[i].rc - wr), 2) + Math.Pow((bins[i].gc - wg), 2) + Math.Pow((bins[i].bc - wb), 2);
			    n2 = bins[i].cnt;
			    nerr *= (n1 * n2) / (n1 + n2);
			    if (nerr >= err)
				    continue;
			    err = nerr;
			    nn = i;
		    }
		    bin1.err = err;
		    bin1.nn = nn;
	    }

        private int pnnquan(int[] pixels, Pnnbin[] bins, ColorPalette palette, int nMaxColors, bool quan_sqrt)
	    {
		    var heap = new int[65537];
		    double err;
		    int l, l2, h, b1, maxbins, extbins;

		    /* Build histogram */
            for(int i=0; i<pixels.Length; ++i)
            {
			    // !!! Can throw gamma correction in here, but what to do about perceptual
			    // !!! nonuniformity then?
                Color c = Color.FromArgb(pixels[i]);
                int index = getARGBIndex(pixels[i]);
                if (bins[index] == null)
                    bins[index] = new Pnnbin();
			    var tb = bins[index];
			    tb.ac += c.A;
			    tb.rc += c.R;
			    tb.gc += c.G;
			    tb.bc += c.B;
			    tb.cnt++;
		    }

            /* Cluster nonempty bins at one end of array */
            maxbins = 0;

            for (int i = 0; i < heap.Length - 1; ++i)
            {
                if (bins[i] == null)
                    continue;

                float d = 1.0f / (float) bins[i].cnt;
                bins[i].ac *= d;
                bins[i].rc *= d;
                bins[i].gc *= d;
                bins[i].bc *= d;
                if (quan_sqrt)
                    bins[i].cnt = (int) Math.Sqrt(bins[i].cnt);
                bins[maxbins++] = bins[i];
            }

            for (int i = 0; i < maxbins - 1; i++)
            {
                bins[i].fw = (i + 1);
                bins[i + 1].bk = i;
            }
            // !!! Already zeroed out by calloc()
            //	bins[0].bk = bins[i].fw = 0;

            /* Initialize nearest neighbors and build heap of them */
            for (int i = 0; i < maxbins; i++)
            {
                find_nn(bins, i);
                /* Push slot on heap */
                err = bins[i].err;
                for (l = ++heap[0]; l > 1; l = l2)
                {
                    l2 = l >> 1;
                    if (bins[h = heap[l2]].err <= err)
                        break;
                    heap[l] = h;
                }
                heap[l] = i;
            }

            /* Merge bins which increase error the least */
            extbins = maxbins - nMaxColors;
            for (int i = 0; i < extbins; )
            {
                /* Use heap to find which bins to merge */
                for (; ; )
                {
                    var tb = bins[b1 = heap[1]]; /* One with least error */
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
                var tb1 = bins[b1];
                var nb = bins[tb1.nn];
                float n1 = tb1.cnt;
                float n2 = nb.cnt;
                float d = 1.0f / (n1 + n2);
                tb1.ac = d * (n1 * tb1.ac + n2 * nb.ac);
                tb1.rc = d * (n1 * tb1.rc + n2 * nb.rc);
                tb1.gc = d * (n1 * tb1.gc + n2 * nb.gc);
                tb1.bc = d * (n1 * tb1.bc + n2 * nb.bc);
                tb1.cnt += nb.cnt;
                tb1.mtm = ++i;

                /* Unchain deleted bin */
                bins[nb.bk].fw = nb.fw;
                bins[nb.fw].bk = nb.bk;
                nb.mtm = 0xFFFF;
            }

		    /* Fill palette */
		    int k = 0;
		    for (int i = 0;; ++k) {
			    var alpha = (int) Math.Round(bins[i].ac);
			    palette.Entries[k] = Color.FromArgb(alpha, (int) Math.Round(bins[i].rc), (int) Math.Round(bins[i].gc), (int) Math.Round(bins[i].bc));
			    if (hasTransparency && palette.Entries[k] == m_transparentColor) {
                    Color temp = palette.Entries[0];
				    palette.Entries[0] = palette.Entries[k];
                    palette.Entries[k] = temp;
                }

			    if ((i = bins[i].fw) == 0)
				    break;
		    }

		    return k;
	    }

        private short nearestColorIndex(ColorPalette palette, int nMaxColors, int[] squares3, int pixel)
	    {
            short k = 0;
            Color c = Color.FromArgb(pixel);

		    int curdist, mindist = int.MaxValue;
            for (short i = 0; i < nMaxColors; i++)
            {
			    Color c2 = palette.Entries[i];
			    int adist = Math.Abs(c2.A - c.A);
			    curdist = squares3[adist];
			    if (curdist > mindist)
				    continue;

			    int rdist = Math.Abs(c2.R - c.R);
			    curdist += squares3[rdist];
			    if (curdist > mindist)
				    continue;

			    int gdist = Math.Abs(c2.G - c.G);
			    curdist += squares3[gdist];
			    if (curdist > mindist)
				    continue;

			    int bdist = Math.Abs(c2.B - c.B);
			    curdist += squares3[bdist];
			    if (curdist > mindist)
				    continue;

			    mindist = curdist;
			    k = i;
		    }
		    return k;
	    }

        private short closestColorIndex(ColorPalette palette, int nMaxColors, int[] squares3, int pixel)
	    {
            short k = 0;
            Color c = Color.FromArgb(pixel);

            ushort[] closest;
            if (!closestMap.TryGetValue(pixel, out closest))
            {
                closest = new ushort[5];
                closest[2] = closest[3] = ushort.MaxValue;

			    for (; k < nMaxColors; k++) {
				    Color c2 = palette.Entries[k];
                    closest[4] = (ushort)(Math.Abs(c.A - c2.A) + Math.Abs(c.R - c2.R) + Math.Abs(c.G - c2.G) + Math.Abs(c.B - c2.B));
				    if (closest[4] < closest[2]) {
					    closest[1] = closest[0];
					    closest[3] = closest[2];
                        closest[0] = (ushort)k;
					    closest[2] = closest[4];
				    }
				    else if (closest[4] < closest[3]) {
                        closest[1] = (ushort)k;
					    closest[3] = closest[4];
				    }
			    }

                if (closest[3] == ushort.MaxValue)
				    closest[2] = 0;
		    }

		    if (closest[2] == 0 || (rand.Next(short.MaxValue) % (closest[3] + closest[2])) <= closest[3])
			    k = (short) closest[0];
		    else
			    k = (short) closest[1];

            closestMap[pixel] = closest;
		    return k;
	    }

        private bool quantize_image(int[] pixels, ColorPalette palette, int nMaxColors, short[] qPixels, int width, int height, bool dither)
	    {
		    var sqr_tbl = new int[Byte.MaxValue + Byte.MaxValue + 1];

		    for (int i = (-Byte.MaxValue); i <= Byte.MaxValue; i++)
			    sqr_tbl[i + Byte.MaxValue] = i * i;

		    var squares3 = new int[sqr_tbl.Length - Byte.MaxValue];
		    for (int i = 0; i < squares3.Length; i++)
			    squares3[i] = sqr_tbl[i + Byte.MaxValue];

		    int pixelIndex = 0;
		    if (dither) {
			    bool odd_scanline = false;
                short[] row0, row1;
			    int a_pix, r_pix, g_pix, b_pix, dir, k;
			    const int DJ = 4;
			    const int DITHER_MAX = 20;
			    int err_len = (width + 2) * DJ;
			    int[] clamp = new int[DJ * 256];
			    int[] limtb = new int[512];
                short[] erowerr = new short[err_len];
                short[] orowerr = new short[err_len];
                short[] lookup = new short[65536];

			    for (int i = 0; i < 256; i++) {
				    clamp[i] = 0;
				    clamp[i + 256] = (short) i;
				    clamp[i + 512] = Byte.MaxValue;
				    clamp[i + 768] = Byte.MaxValue;

				    limtb[i] = -DITHER_MAX;
				    limtb[i + 256] = DITHER_MAX;
			    }
			    for (int i = -DITHER_MAX; i <= DITHER_MAX; i++)
				    limtb[i + 256] = i;

                for (int i = 0; i < height; i++)
                {
				    if (odd_scanline) {
					    dir = -1;
					    pixelIndex += (int) (width - 1);
					    row0 = orowerr;
					    row1 = erowerr;
				    }
				    else {
					    dir = 1;
					    row0 = erowerr;
					    row1 = orowerr;
				    }
				
				    int cursor0 = DJ, cursor1 = width * DJ;
				    row1[cursor1] = row1[cursor1 + 1] = row1[cursor1 + 2] = row1[cursor1 + 3] = 0;
				    for (int j = 0; j < width; j++) {
                        Color c = Color.FromArgb(pixels[pixelIndex]);
					    r_pix = clamp[((row0[cursor0] + 0x1008) >> 4) + c.R];
					    g_pix = clamp[((row0[cursor0 + 1] + 0x1008) >> 4) + c.G];
					    b_pix = clamp[((row0[cursor0 + 2] + 0x1008) >> 4) + c.B];
					    a_pix = clamp[((row0[cursor0 + 3] + 0x1008) >> 4) + c.A];

					    Color c1 = Color.FromArgb(a_pix, r_pix, g_pix, b_pix);
					    int offset = getARGBIndex(c1.ToArgb());
					    if (lookup[offset] == 0)
                            lookup[offset] = (short) (nearestColorIndex(palette, nMaxColors, squares3, c1.ToArgb()) + 1);
					    qPixels[pixelIndex] = (short) (lookup[offset] - 1);

					    Color c2 = palette.Entries[qPixels[pixelIndex]];

					    r_pix = limtb[r_pix - c2.R + 256];
					    g_pix = limtb[g_pix - c2.G + 256];
					    b_pix = limtb[b_pix - c2.B + 256];
					    a_pix = limtb[a_pix - c2.A + 256];

					    k = r_pix * 2;
					    row1[cursor1 - DJ] = (short) r_pix;
                        row1[cursor1 + DJ] += (short)(r_pix += k);
                        row1[cursor1] += (short) (r_pix += k);
                        row0[cursor0 + DJ] += (short) (r_pix += k);

					    k = g_pix * 2;
                        row1[cursor1 + 1 - DJ] = (short) g_pix;
                        row1[cursor1 + 1 + DJ] += (short) (g_pix += k);
                        row1[cursor1 + 1] += (short) (g_pix += k);
                        row0[cursor0 + 1 + DJ] += (short) (g_pix += k);

					    k = b_pix * 2;
                        row1[cursor1 + 2 - DJ] = (short) b_pix;
                        row1[cursor1 + 2 + DJ] += (short) (b_pix += k);
                        row1[cursor1 + 2] += (short) (b_pix += k);
                        row0[cursor0 + 2 + DJ] += (short) (b_pix += k);

					    k = a_pix * 2;
                        row1[cursor1 + 3 - DJ] = (short) a_pix;
                        row1[cursor1 + 3 + DJ] += (short) (a_pix += k);
                        row1[cursor1 + 3] += (short) (a_pix += k);
                        row0[cursor0 + 3 + DJ] += (short) (a_pix += k);

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
                    qPixels[i] = nearestColorIndex(palette, nMaxColors, squares3, pixels[i]);
		    }
		    else {
			    for (int i = 0; i < qPixels.Length; i++)
                    qPixels[i] = closestColorIndex(palette, nMaxColors, squares3, pixels[i]);
		    }

		    return true;
	    }

        protected bool quantize_image(int[] pixels, short[] qPixels, int width, int height)
	    {
		    int pixelIndex = 0;
		    bool odd_scanline = false;
		    short[] row0, row1;
		    int a_pix, r_pix, g_pix, b_pix, dir, k;
		    const int DJ = 4;
		    const int DITHER_MAX = 20;
		    int err_len = (width + 2) * DJ;
		    int[] clamp = new int[DJ * 256];
		    int[] limtb = new int[512];
            short[] erowerr = new short[err_len];
            short[] orowerr = new short[err_len];
		    int[] lookup = new int[65536];

		    for (int i = 0; i < 256; i++) {
			    clamp[i] = 0;
			    clamp[i + 256] = (short) i;
			    clamp[i + 512] = Byte.MaxValue;
			    clamp[i + 768] = Byte.MaxValue;

			    limtb[i] = -DITHER_MAX;
			    limtb[i + 256] = DITHER_MAX;
		    }
		    for (int i = -DITHER_MAX; i <= DITHER_MAX; i++)
			    limtb[i + 256] = i;

		    for (int i = 0; i < height; i++) {
			    if (odd_scanline) {
				    dir = -1;
				    pixelIndex += (int)(width - 1);
				    row0 = orowerr;
				    row1 = erowerr;
			    }
			    else {
				    dir = 1;
				    row0 = erowerr;
				    row1 = orowerr;
			    }
			
			    int cursor0 = DJ, cursor1 = (width * DJ);
			    row1[cursor1] = row1[cursor1 + 1] = row1[cursor1 + 2] = row1[cursor1 + 3] = 0;
			    for (int j = 0; j < width; j++) {
				    Color c = Color.FromArgb(pixels[pixelIndex]);

				    r_pix = clamp[((row0[cursor0] + 0x1008) >> 4) + c.R];
				    g_pix = clamp[((row0[cursor0 + 1] + 0x1008) >> 4) + c.G];
				    b_pix = clamp[((row0[cursor0 + 2] + 0x1008) >> 4) + c.B];
				    a_pix = clamp[((row0[cursor0 + 3] + 0x1008) >> 4) + c.A];

				    Color c1 = Color.FromArgb(a_pix, r_pix, g_pix, b_pix);
				    int offset = getARGBIndex(c1.ToArgb());
				    if (lookup[offset] == 0) {
                        Color rgba1 = Color.FromArgb(c1.A, (c1.R & 0xF8), (c1.G & 0xFC), (c1.B & 0xF8));
					    if (hasSemiTransparency)
						    rgba1 = Color.FromArgb((c1.A & 0xF0), (c1.R & 0xF0), (c1.G & 0xF0), (c1.B & 0xF0));
					    else if (hasTransparency)
                            rgba1 = Color.FromArgb((c1.A & 0x80), (c1.R & 0xF8), (c1.G & 0xF8), (c1.B & 0xF8));
					    lookup[offset] = rgba1.ToArgb();
				    }

                    qPixels[pixelIndex] = (short) offset;
				    Color c2 = Color.FromArgb(lookup[offset]);				

				    r_pix = limtb[r_pix - c2.R + 256];
				    g_pix = limtb[g_pix - c2.G + 256];
				    b_pix = limtb[b_pix - c2.B + 256];
				    a_pix = limtb[a_pix - c2.A + 256];

				    k = r_pix * 2;
				    row1[cursor1 - DJ] = (short) r_pix;
                    row1[cursor1 + DJ] += (short) (r_pix += k);
                    row1[cursor1] += (short) (r_pix += k);
                    row0[cursor0 + DJ] += (short) (r_pix += k);

				    k = g_pix * 2;
                    row1[cursor1 + 1 - DJ] = (short) g_pix;
                    row1[cursor1 + 1 + DJ] += (short) (g_pix += k);
                    row1[cursor1 + 1] += (short) (g_pix += k);
                    row0[cursor0 + 1 + DJ] += (short) (g_pix += k);

				    k = b_pix * 2;
                    row1[cursor1 + 2 - DJ] = (short) b_pix;
                    row1[cursor1 + 2 + DJ] += (short) (b_pix += k);
                    row1[cursor1 + 2] += (short) (b_pix += k);
                    row0[cursor0 + 2 + DJ] += (short) (b_pix += k);

				    k = a_pix * 2;
                    row1[cursor1 + 3 - DJ] = (short) a_pix;
                    row1[cursor1 + 3 + DJ] += (short) (a_pix += k);
                    row1[cursor1 + 3] += (short) (a_pix += k);
                    row0[cursor0 + 3 + DJ] += (short) (a_pix += k);

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

        protected bool ProcessImagePixels(Bitmap dest, ColorPalette palette, short[] qPixels)
        {
            dest.Palette = palette;

            int bpp = Image.GetPixelFormatSize(dest.PixelFormat);
            int w = dest.Width;
            int h = dest.Height;

            BitmapData targetData = dest.LockBits(new Rectangle(0, 0, w, h), ImageLockMode.WriteOnly, dest.PixelFormat);

            int pixelIndex = 0;
            int strideDest;

            unsafe
            {
                var pRowDest = (byte*)targetData.Scan0;

                // Compensate for possible negative stride
                if (targetData.Stride > 0)
                    strideDest = targetData.Stride;
                else
                {
                    pRowDest += h * targetData.Stride;
                    strideDest = -targetData.Stride;
                }

                // Second loop: fill indexed bitmap
                for (int y = 0; y < h; y++)
                {	// For each row...
                    for (int x = 0; x < w; x++)
                    {	// ...for each pixel...
                        byte nibbles = 0;
                        byte index = (byte)qPixels[pixelIndex++];

                        switch (bpp)
                        {
                            case 8:
                                pRowDest[x] = index;
                                break;
                            case 4:
                                // First pixel is the high nibble. From and To indices are 0..16
                                nibbles = pRowDest[x / 2];
                                if ((x & 1) == 0)
                                {
                                    nibbles &= 0x0F;
                                    nibbles |= (byte)(index << 4);
                                }
                                else
                                {
                                    nibbles &= 0xF0;
                                    nibbles |= index;
                                }

                                pRowDest[x / 2] = nibbles;
                                break;
                            case 1:
                                // First pixel is MSB. From and To are 0 or 1.
                                int pos = x / 8;
                                byte mask = (byte)(128 >> (x & 7));
                                if (index == 0)
                                    pRowDest[pos] &= (byte)~mask;
                                else
                                    pRowDest[pos] |= mask;
                                break;
                        }
                    }

                    pRowDest += strideDest;
                }
            }

            dest.UnlockBits(targetData);
            return true;
        }

        protected bool ProcessImagePixels(Bitmap dest, short[] qPixels)
        {
            int bpp = Image.GetPixelFormatSize(dest.PixelFormat);
            if (bpp < 16)
                return false;

            int w = dest.Width;
            int h = dest.Height;

            BitmapData targetData = dest.LockBits(new Rectangle(0, 0, w, h), ImageLockMode.WriteOnly, dest.PixelFormat);

            int pixelIndex = 0;
            int strideDest;

            unsafe
            {
                var pRowDest = (byte*)targetData.Scan0;

                // Compensate for possible negative stride
                if (targetData.Stride > 0)
                    strideDest = targetData.Stride;
                else
                {
                    pRowDest += h * targetData.Stride;
                    strideDest = -targetData.Stride;
                }

                // Second loop: fill indexed bitmap
                for (int y = 0; y < h; y++)
                {
                    // For each row...
                    for (int x = 0; x < w * 2; )
                    {
                        ushort argb = (ushort) qPixels[pixelIndex++];
                        pRowDest[x++] = (byte)(argb & 0xFF);
                        pRowDest[x++] = (byte)(argb >> 8);
                    }
                    pRowDest += strideDest;
                }
            }

            dest.UnlockBits(targetData);
            return true;
        }

        protected bool IsValidFormat(PixelFormat pixelFormat, int nMaxColors)
        {
            int bitDepth = Image.GetPixelFormatSize(pixelFormat);
            return Math.Pow(2, bitDepth) >= nMaxColors;
        }

	    public virtual bool QuantizeImage(Bitmap source, Bitmap dest, int nMaxColors, bool dither)
	    {
		    int bitDepth = Image.GetPixelFormatSize(source.PixelFormat);
            if (!IsValidFormat(dest.PixelFormat, nMaxColors))
                return false;

		    int bitmapWidth = source.Width;
		    int bitmapHeight = source.Height;

		    hasTransparency = hasSemiTransparency = false;
		    int pixelIndex = 0;
		    var pixels = new int[bitmapWidth * bitmapHeight];
		    if (bitDepth <= 16) {
			    for (int y = 0; y < bitmapHeight; y++) {
				    for (int x = 0; x < bitmapWidth; x++) {
					    Color color = source.GetPixel(x, y);
					    if (color.A < Byte.MaxValue) {
						    hasSemiTransparency = true;
						    if (color.A == 0) {
							    hasTransparency = true;
							    m_transparentColor = color;
						    }
					    }
					    pixels[pixelIndex++] = color.ToArgb();
				    }
			    }                
		    }

		    // Lock bits on 3x8 source bitmap
		    else {
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
                            pixels[pixelIndex++] = argb.ToArgb();
                        }
                        pRowSource += strideSource;
                    }
                }			    

			    source.UnlockBits(data);
		    }

            var qPixels = new short[bitmapWidth * bitmapHeight];
		    if (nMaxColors > 256) {
                hasSemiTransparency = false;
			    quantize_image(pixels, qPixels, bitmapWidth, bitmapHeight);
			    return ProcessImagePixels(dest, qPixels);
		    }

		    var bins = new Pnnbin[65536];
		    var palette = dest.Palette;
		    bool quan_sqrt = nMaxColors > Byte.MaxValue;
		    if (nMaxColors > 2)
                nMaxColors = pnnquan(pixels, bins, palette, nMaxColors, quan_sqrt);
		    else {
                if (hasTransparency)
                {
				    palette.Entries[0] = Color.Transparent;
				    palette.Entries[1] = Color.Black;
			    }
			    else {
				    palette.Entries[0] = Color.Black;
				    palette.Entries[1] = Color.White;
			    }
		    }

            quantize_image(pixels, palette, nMaxColors, qPixels, bitmapWidth, bitmapHeight, dither);
		    closestMap.Clear();

		    return ProcessImagePixels(dest, palette, qPixels);
	    }
    }

}
