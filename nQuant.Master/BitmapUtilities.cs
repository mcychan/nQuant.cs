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

        public static bool GrabPixels(Bitmap source, int[] pixels, ref bool hasSemiTransparency, ref Color transparentColor, ref int transparentPixelIndex)
        {
            int bitmapWidth = source.Width;
            int bitmapHeight = source.Height;

            hasSemiTransparency = false;
            transparentPixelIndex = -1;

            int transparentIndex = -1;
            var palettes = source.Palette.Entries;
            foreach (var pPropertyItem in source.PropertyItems)
            {
                if (pPropertyItem.Id == BitmapUtilities.PropertyTagIndexTransparent)
                {
                    transparentIndex = pPropertyItem.Value[0];
                    Color c = palettes[transparentIndex];
                    transparentColor = Color.FromArgb(0, c.R, c.G, c.B);
                }
            }

            int pixelIndex = 0;
            var data = source.LockBits(new Rectangle(0, 0, bitmapWidth, bitmapHeight), ImageLockMode.ReadOnly, PixelFormat.Format32bppArgb);

            // Declare an array to hold the bytes of the bitmap.
            int bytesLength = Math.Abs(data.Stride) * bitmapHeight;
            var rgbValues = new byte[bytesLength];

            // Copy the RGB values into the array.
            System.Runtime.InteropServices.Marshal.Copy(data.Scan0, rgbValues, 0, bytesLength);

            for (int i = 0; i < rgbValues.Length; i += 4)
            {
                byte pixelBlue = rgbValues[i];
                byte pixelGreen = rgbValues[i + 1];
                byte pixelRed = rgbValues[i + 2];
                byte pixelAlpha = rgbValues[i + 3];

                var argb = Color.FromArgb(pixelAlpha, pixelRed, pixelGreen, pixelBlue);
                var argb1 = Color.FromArgb(0, pixelRed, pixelGreen, pixelBlue);
                if (transparentIndex > -1 && transparentColor.ToArgb() == argb1.ToArgb())
                {
                    pixelAlpha = 0;
                    argb = argb1;
                }

                if (pixelAlpha < Byte.MaxValue)
                {

                    if (pixelAlpha == 0)
                    {
                        transparentPixelIndex = pixelIndex;
                        transparentColor = argb;
                    }
                    else
                        hasSemiTransparency = true;
                }
                pixels[pixelIndex++] = argb.ToArgb();
            }

            source.UnlockBits(data);
            return true;
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

            int bpp = Image.GetPixelFormatSize(dest.PixelFormat);
            int w = dest.Width;
            int h = dest.Height;

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
        public static Bitmap ProcessImagePixels(Bitmap dest, int[] qPixels, bool hasSemiTransparency, int transparentPixelIndex)
        {
            int bpp = Image.GetPixelFormatSize(dest.PixelFormat);
            if (bpp < 16)
                return dest;

            int w = dest.Width;
            int h = dest.Height;

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
    }
}
