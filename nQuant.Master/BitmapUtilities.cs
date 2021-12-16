using System;
using System.Drawing;
using System.Drawing.Imaging;
using System.Runtime.Serialization;

namespace nQuant.Master
{
    class BitmapUtilities
    {
        protected const int PropertyTagTypeByte = 1;
        public const int PropertyTagIndexTransparent = 0x5104;        

        public static void Swap<T>(ref T x, ref T y)
        {
            T t = y;
            y = x;
            x = t;
        }
        public static double Sqr(double value)
        {
            return value * value;
        }
        public static int GetARGB1555(int argb)
        {
            var c = Color.FromArgb(argb);
            return (c.A & 0x80) << 8 | (c.R & 0xF8) << 7 | (c.G & 0xF8) << 2 | (c.B >> 3);
        }

        public static int GetARGBIndex(int argb, bool hasSemiTransparency, bool hasTransparency)
        {
            var c = Color.FromArgb(argb);
            if (hasSemiTransparency)
                return (c.A & 0xF0) << 8 | (c.R & 0xF0) << 4 | (c.G & 0xF0) | (c.B >> 4);
            if (hasTransparency)
                return GetARGB1555(argb);
            return (c.R & 0xF8) << 8 | (c.G & 0xFC) << 3 | (c.B >> 3);
        }

        protected static int[] CalcDitherPixel(Color c, short[] clamp, int[] rowerr, int cursor, bool noBias)
        {
            var ditherPixel = new int[4];
            if (noBias)
            {
                ditherPixel[0] = clamp[((rowerr[cursor] + 0x1008) >> 4) + c.R];
                ditherPixel[1] = clamp[((rowerr[cursor + 1] + 0x1008) >> 4) + c.G];
                ditherPixel[2] = clamp[((rowerr[cursor + 2] + 0x1008) >> 4) + c.B];
                ditherPixel[3] = clamp[((rowerr[cursor + 3] + 0x1008) >> 4) + c.A];
            }
            else
            {
                ditherPixel[0] = clamp[((rowerr[cursor] + 0x2010) >> 5) + c.R];
                ditherPixel[1] = clamp[((rowerr[cursor + 1] + 0x1008) >> 4) + c.G];
                ditherPixel[2] = clamp[((rowerr[cursor + 2] + 0x2010) >> 5) + c.B];
                ditherPixel[3] = c.A;
            }
            return ditherPixel;
        }
        public static int[] Quantize_image(int width, int height, int[] pixels, Color[] palette, Ditherable ditherable, bool hasSemiTransparency, bool dither)
        {
            var qPixels = new int[width * height];
            int nMaxColors = palette.Length;
            int pixelIndex = 0;
            if (dither)
            {
                const short DJ = 4;
                const short BLOCK_SIZE = 256;
                const short DITHER_MAX = 20;
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
                {
                    limtb[i + BLOCK_SIZE] = i;
                    if (nMaxColors > 16 && i % 4 == 3)
                        limtb[i + BLOCK_SIZE] = 0;
                }

                bool noBias = hasSemiTransparency || nMaxColors < 64;
                int dir = 1;
                var row0 = new int[err_len];
                var row1 = new int[err_len];
                var lookup = new int[65536];
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
                        if (noBias && a_pix > 0xF0)
                        {
                            int offset = ditherable.GetColorIndex(c1.ToArgb());
                            if (lookup[offset] == 0)
                                lookup[offset] = (c.A == 0) ? 1 : ditherable.DitherColorIndex(palette, nMaxColors, c1.ToArgb()) + 1;
                            qPixels[pixelIndex] = lookup[offset] - 1;
                        }
                        else
                            qPixels[pixelIndex] = (c.A == 0) ? 0 : ditherable.DitherColorIndex(palette, nMaxColors, c1.ToArgb());

                        var c2 = palette[qPixels[pixelIndex]];
                        if (nMaxColors > 256)
                            qPixels[pixelIndex] = hasSemiTransparency ? c2.ToArgb() : ditherable.GetColorIndex(c2.ToArgb());

                        r_pix = limtb[r_pix - c2.R + BLOCK_SIZE];
                        g_pix = limtb[g_pix - c2.G + BLOCK_SIZE];
                        b_pix = limtb[b_pix - c2.B + BLOCK_SIZE];
                        a_pix = limtb[a_pix - c2.A + BLOCK_SIZE];

                        int k = r_pix * 2;
                        row1[cursor1 - DJ] = r_pix;
                        row1[cursor1 + DJ] += (r_pix += k);
                        row1[cursor1] += (r_pix += k);
                        row0[cursor0 + DJ] += (r_pix + k);

                        k = g_pix * 2;
                        row1[cursor1 + 1 - DJ] = g_pix;
                        row1[cursor1 + 1 + DJ] += (g_pix += k);
                        row1[cursor1 + 1] += (g_pix += k);
                        row0[cursor0 + 1 + DJ] += (g_pix + k);

                        k = b_pix * 2;
                        row1[cursor1 + 2 - DJ] = b_pix;
                        row1[cursor1 + 2 + DJ] += (b_pix += k);
                        row1[cursor1 + 2] += (b_pix += k);
                        row0[cursor0 + 2 + DJ] += (b_pix + k);

                        k = a_pix * 2;
                        row1[cursor1 + 3 - DJ] = a_pix;
                        row1[cursor1 + 3 + DJ] += (a_pix += k);
                        row1[cursor1 + 3] += (a_pix += k);
                        row0[cursor0 + 3 + DJ] += (a_pix + k);

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

            for (int i = 0; i < qPixels.Length; ++i)
                qPixels[i] = ditherable.DitherColorIndex(palette, nMaxColors, pixels[i]);            

            return qPixels;
        }
		
        public static Bitmap ProcessImagePixels(Bitmap dest, int[] qPixels, bool hasSemiTransparency, int transparentPixelIndex)
        {
            var bpp = Image.GetPixelFormatSize(dest.PixelFormat);
            if (bpp < 16)
                return dest;

            var w = dest.Width;
            var h = dest.Height;

            if (hasSemiTransparency && dest.PixelFormat < PixelFormat.Format32bppArgb)
                dest = new Bitmap(w, h, PixelFormat.Format32bppArgb);
            else if (transparentPixelIndex >= 0 && dest.PixelFormat < PixelFormat.Format16bppArgb1555)
                dest = new Bitmap(w, h, PixelFormat.Format16bppArgb1555);
            else if (dest.PixelFormat < PixelFormat.Format16bppRgb565)
                dest = new Bitmap(w, h, PixelFormat.Format16bppRgb565);

            bpp = Image.GetPixelFormatSize(dest.PixelFormat);
            var targetData = dest.LockBits(new Rectangle(0, 0, w, h), ImageLockMode.WriteOnly, dest.PixelFormat);

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

                if (bpp == 32)
                {
                    for (int y = 0; y < h; ++y)
                    {   // For each row...
                        for (int x = 0; x < w * 4;)
                        {
                            var c = Color.FromArgb(qPixels[pixelIndex++]);
                            pRowDest[x++] = c.B;
                            pRowDest[x++] = c.G;
                            pRowDest[x++] = c.R;
                            pRowDest[x++] = c.A;
                        }
                        pRowDest += strideDest;
                    }
                }
                else if (bpp == 16)
                {
                    for (int y = 0; y < h; ++y)
                    {   // For each row...
                        for (int x = 0; x < w * 2;)
                        {
                            var argb = qPixels[pixelIndex++];
                            pRowDest[x++] = (byte)(argb & 0xFF);
                            pRowDest[x++] = (byte)(argb >> 8);
                        }
                        pRowDest += strideDest;
                    }
                }
                else
                {
                    for (int y = 0; y < h; ++y)
                    {   // For each row...
                        for (int x = 0; x < w; ++x)
                            pRowDest[x] = (byte)qPixels[pixelIndex++];
                        pRowDest += strideDest;
                    }
                }
            }

            dest.UnlockBits(targetData);
            return dest;
        }
		
        public static Bitmap ProcessImagePixels(Bitmap dest, Color[] palettes, int[] qPixels, bool hasTransparent)
        {
            if (hasTransparent)
            {
                var propertyItem = (PropertyItem)FormatterServices.GetUninitializedObject(typeof(PropertyItem));
                propertyItem.Id = PropertyTagIndexTransparent;
                propertyItem.Len = 1;
                propertyItem.Type = PropertyTagTypeByte;
                propertyItem.Value = new byte[] { 0 };

                dest.SetPropertyItem(propertyItem);
            }

            var palette = dest.Palette;
            for (int i = 0; i < palettes.Length; ++i)
                palette.Entries[i] = palettes[i];
            dest.Palette = palette;

            var bpp = Image.GetPixelFormatSize(dest.PixelFormat);
            var w = dest.Width;
            var h = dest.Height;

            var targetData = dest.LockBits(new Rectangle(0, 0, w, h), ImageLockMode.WriteOnly, dest.PixelFormat);

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
                for (int y = 0; y < h; ++y)
                {	// For each row...
                    for (int x = 0; x < w; ++x)
                    {	// ...for each pixel...
                        byte nibbles = 0;
                        var index = (byte)qPixels[pixelIndex++];

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
            return dest;
        }        
		
        public static bool GrabPixels(Bitmap source, int[] pixels, ref bool hasSemiTransparency, ref Color transparentColor, ref int transparentPixelIndex, byte alphaThreshold, int nMaxColors = 2)
        {
            var bitmapWidth = source.Width;
            var bitmapHeight = source.Height;

            hasSemiTransparency = false;
            transparentPixelIndex = -1;

            int transparentIndex = -1;
            var palettes = source.Palette.Entries;
            if (nMaxColors > 2)
            {
                foreach (var pPropertyItem in source.PropertyItems)
                {
                    if (pPropertyItem.Id == PropertyTagIndexTransparent)
                    {
                        transparentIndex = pPropertyItem.Value[0];
                        var c = palettes[transparentIndex];
                        transparentColor = Color.FromArgb(0, c.R, c.G, c.B);
                        break;
                    }
                }
            }

            int pixelIndex = 0, strideSource;
            var data = source.LockBits(new Rectangle(0, 0, bitmapWidth, bitmapHeight), ImageLockMode.ReadOnly, PixelFormat.Format32bppArgb);

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
                for (int y = 0; y < bitmapHeight; ++y)
                {
                    var pPixelSource = pRowSource + (y * strideSource);
                    // For each row...
                    for (int x = 0; x < bitmapWidth; ++x)
                    {
                        var pixelBlue = *pPixelSource++;
                        var pixelGreen = *pPixelSource++;
                        var pixelRed = *pPixelSource++;
                        var pixelAlpha = *pPixelSource++;

                        var argb = Color.FromArgb(pixelAlpha, pixelRed, pixelGreen, pixelBlue);
                        var argb1 = Color.FromArgb(0, pixelRed, pixelGreen, pixelBlue);
                        if (transparentIndex > -1 && transparentColor.ToArgb() == argb1.ToArgb())
                        {
                            pixelAlpha = 0;
                            argb = argb1;
                        }

                        if (pixelAlpha < 0xE0)
                        {
                            if (pixelAlpha == 0)
                            {
                                transparentPixelIndex = pixelIndex;
                                if (nMaxColors > 2)
                                    transparentColor = argb;
                                else
                                    argb = transparentColor;
                            }
                            else if(pixelAlpha > alphaThreshold)
                                hasSemiTransparency = true;
                        }
                        pixels[pixelIndex++] = argb.ToArgb();
                    }
                }
                source.UnlockBits(data);
            }
            return true;
        }
    }
}
