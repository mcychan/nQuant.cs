using System;
using System.Collections.Generic;
using System.Drawing;

/* The Hilbert curve is a space filling curve that visits every point in a square grid with a size of any other power of 2.
Copyright (c) 2021 Miller Cy Chan
* It was first described by David Hilbert in 1892. Applications of the Hilbert curve are in image processing: especially image compression and dithering. */

using static nQuant.Master.HilbertCurve.Direction;

namespace nQuant.Master
{
    internal delegate ushort DitherFn(Color[] palette, int nMaxColors, int pixel);

    internal delegate int GetColorIndexFn(int pixel);

    class HilbertCurve
    {
        internal enum Direction { LEFT, RIGHT, DOWN, UP };

        internal sealed class ErrorBox
        {
            private readonly float[] p;
            internal ErrorBox()
            {
                p = new float[4];
            }

            internal ErrorBox(Color c)
            {
                p = new float[] {
                    c.R,
                    c.G,
                    c.B,
                    c.A
                };
            }

            internal float this[int i]
            {
                get { return p[i]; }
                set { p[i] = value; }
            }

            internal int Length
            {
                get { return p.Length; }
            }
        }

        private int x, y;
        private readonly int width;
        private readonly int height;
        private readonly int[] pixels;
	    private readonly Color[] palette;
	    private readonly int[] qPixels;
        private readonly DitherFn ditherFn;
        private readonly GetColorIndexFn getColorIndexFn;
        private readonly List<ErrorBox> errorq;
        private readonly float[] weights;
        private readonly int[] lookup;

        private const byte DITHER_MAX = 16;
        private const float BLOCK_SIZE = 256f;

        private HilbertCurve(int width, int height, int[] image, Color[] palette, int[] qPixels, DitherFn ditherFn, GetColorIndexFn getColorIndexFn)
        {
            x = 0;
            y = 0;
            this.width = width;
            this.height = height;
            this.pixels = image;
            this.palette = palette;
            this.qPixels = qPixels;
            this.ditherFn = ditherFn;
            this.getColorIndexFn = getColorIndexFn;
            errorq = new List<ErrorBox>();
            weights = new float[DITHER_MAX];
            lookup = new int[65536];
        }

        private void DitherCurrentPixel()
        {
            if (x >= 0 && y >= 0 && x < width && y < height)
            {
                Color pixel= Color.FromArgb(pixels[x + y * width]);
                ErrorBox error = new ErrorBox(pixel);
                for (int c = 0; c < DITHER_MAX; ++c)
                {
                    ErrorBox eb = errorq[c];
                    for (int j = 0; j < eb.Length; ++j)
                        error[j] += eb[j] * weights[c];
                }

                int r_pix = (int)Math.Min(BLOCK_SIZE - 1, Math.Max(error[0], 0.0));
                int g_pix = (int)Math.Min(BLOCK_SIZE - 1, Math.Max(error[1], 0.0));
                int b_pix = (int)Math.Min(BLOCK_SIZE - 1, Math.Max(error[2], 0.0));
                int a_pix = (int)Math.Min(BLOCK_SIZE - 1, Math.Max(error[3], 0.0));

                Color c2 = Color.FromArgb(a_pix, r_pix, g_pix, b_pix);
                if (palette.Length < 64)
                {
                    int offset = getColorIndexFn(c2.ToArgb());
                    if (lookup[offset] == 0)
                        lookup[offset] = (pixel.A == 0) ? 1 : ditherFn(palette, palette.Length, c2.ToArgb()) + 1;
                    qPixels[x + y * width] = lookup[offset] - 1;
                }
                else
                    qPixels[x + y * width] = ditherFn(palette, palette.Length, c2.ToArgb());

                errorq.RemoveAt(0);
                c2 = palette[qPixels[x + y * width]];
                error[0] = r_pix - c2.R;
                error[1] = g_pix - c2.G;
                error[2] = b_pix - c2.B;
                error[3] = a_pix - c2.A;

                for (int j = 0; j < error.Length; ++j)
                {
                    if (Math.Abs(error[j]) > DITHER_MAX)
                        error[j] = error[j] < 0 ? -DITHER_MAX : DITHER_MAX;
                }
                errorq.Add(error);
            }
        }

        private void Run()
        {
            /* Dithers all pixels of the image in sequence using
             * the Hilbert path, and distributes the error in
             * a sequence of 16 pixels.
             */
            x = y = 0;
            float weightRatio = (float)Math.Pow(BLOCK_SIZE + 1f, 1f / (DITHER_MAX - 1f));
            float weight = 1f, sumweight = 0f;
            for (int c = 0; c < DITHER_MAX; ++c)
            {
                errorq.Add(new ErrorBox());
                sumweight += (weights[DITHER_MAX - c - 1] = 1.0f / weight);
                weight *= weightRatio;
            }

            weight = 0f; /* Normalize */
            for (int c = 0; c < DITHER_MAX; ++c)
                weight += (weights[c] /= sumweight);
            weights[0] += 1f - weight;
            /* Walk the path. */
            int i = Math.Max(width, height), depth = 0;
            while (i > 0)
            {
                ++depth;
                i >>= 1;
            }

            Iter(depth, UP);
            DitherCurrentPixel();
        }

        private void NavTo(Direction dir)
        {
            DitherCurrentPixel();
            switch (dir)
            {
                case LEFT:
                    --x;
                    break;
                case RIGHT:
                    ++x;
                    break;
                case UP:
                    --y;
                    break;
                case DOWN:
                    ++y;
                    break;
            }
        }

        private void Curve(int level, Direction a, Direction b, Direction c, Direction d, Direction e, Direction f, Direction g)
        {
            Iter(level - 1, a);
            NavTo(e);
            Iter(level - 1, b);
            NavTo(f);
            Iter(level - 1, c);
            NavTo(g);
            Iter(level - 1, d);
        }

        private void Iter(int level, Direction dir)
        {
            if (level <= 0)
                return;

            switch (dir)
            {
                case LEFT:
                    Curve(level, UP, LEFT, LEFT, DOWN, RIGHT, DOWN, LEFT);
                    break;
                case RIGHT:
                    Curve(level, DOWN, RIGHT, RIGHT, UP, LEFT, UP, RIGHT);
                    break;
                case UP:
                    Curve(level, LEFT, UP, UP, RIGHT, DOWN, RIGHT, UP);
                    break;
                case DOWN:
                    Curve(level, RIGHT, DOWN, DOWN, LEFT, UP, LEFT, DOWN);
                    break;
            }
        }

        public static int[] Dither(int width, int height, int[] pixels, Color[] palette, DitherFn ditherFn, GetColorIndexFn getColorIndexFn)
        {
            var qPixels = new int[pixels.Length];
            new HilbertCurve(width, height, pixels, palette, qPixels, ditherFn, getColorIndexFn).Run();
            return qPixels;
        }
    }
}
