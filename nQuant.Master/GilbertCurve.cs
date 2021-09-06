using System;
using System.Collections.Generic;
using System.Drawing;

/* Generalized Hilbert ("gilbert") space-filling curve for rectangular domains of arbitrary (non-power of two) sizes.
Copyright (c) 2021 Miller Cy Chan
* A general rectangle with a known orientation is split into three regions ("up", "right", "down"), for which the function calls itself recursively, until a trivial path can be produced. */

namespace nQuant.Master
{
    class GilbertCurve
    {
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

        private const byte DITHER_MAX = 9;
        private const float BLOCK_SIZE = 243f;

        private GilbertCurve(int width, int height, int[] image, Color[] palette, int[] qPixels, DitherFn ditherFn, GetColorIndexFn getColorIndexFn)
        {
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
		
		private static int Sign(int x) {
    	    if(x < 0)
    		    return -1;
    	    return (x > 0) ? 1 : 0;
        }

        private void DitherPixel(int x, int y)
        {
            int bidx = x + y * width;
            Color pixel= Color.FromArgb(pixels[bidx]);
            ErrorBox error = new ErrorBox(pixel);
            for (int c = 0; c < DITHER_MAX; ++c)
            {
                ErrorBox eb = errorq[c];
                for (int j = 0; j < eb.Length; ++j)
                    error[j] += eb[j] * weights[c];
            }

            int r_pix = (int)Math.Min(Byte.MaxValue, Math.Max(error[0], 0.0));
            int g_pix = (int)Math.Min(Byte.MaxValue, Math.Max(error[1], 0.0));
            int b_pix = (int)Math.Min(Byte.MaxValue, Math.Max(error[2], 0.0));
            int a_pix = (int)Math.Min(Byte.MaxValue, Math.Max(error[3], 0.0));

            Color c2 = Color.FromArgb(a_pix, r_pix, g_pix, b_pix);
            if (palette.Length < 64)
            {
                int offset = getColorIndexFn(c2.ToArgb());
                if (lookup[offset] == 0)
                    lookup[offset] = (pixel.A == 0) ? 1 : ditherFn(palette, palette.Length, c2.ToArgb()) + 1;
                qPixels[bidx] = lookup[offset] - 1;
            }
            else
                qPixels[bidx] = ditherFn(palette, palette.Length, c2.ToArgb());

            errorq.RemoveAt(0);
            c2 = palette[qPixels[bidx]];
            if (palette.Length > 256)
                qPixels[bidx] = (short) getColorIndexFn(c2.ToArgb());

            error[0] = r_pix - c2.R;
            error[1] = g_pix - c2.G;
            error[2] = b_pix - c2.B;
            error[3] = a_pix - c2.A;

            for (int j = 0; j < error.Length; ++j)
            {
                if (Math.Abs(error[j]) < DITHER_MAX)
                    continue;

                error[j] /= 3.0f;			    
            }
            errorq.Add(error);
        }
		
	private void Generate2d(int x, int y, int ax, int ay, int bx, int by) {    	
		int w = Math.Abs(ax + ay);
		int h = Math.Abs(bx + by);
		int dax = Sign(ax);
		int day = Sign(ay);
		int dbx = Sign(bx);
		int dby = Sign(by);

		if (h == 1) {
			for (int i = 0; i < w; ++i){
				DitherPixel(x, y);
				x += dax;
				y += day;
			}
			return;
		}

		if (w == 1) {
			for (int i = 0; i < h; ++i){
				DitherPixel(x, y);
				x += dbx;
				y += dby;
			}
			return;
		}

		int ax2 = ax / 2;
		int ay2 = ay / 2;
		int bx2 = bx / 2;
		int by2 = by / 2;

		int w2 = Math.Abs(ax2 + ay2);
		int h2 = Math.Abs(bx2 + by2);

		if (2 * w > 3 * h) {
			if ((w2 % 2) != 0 && w > 2) {
				ax2 += dax;
				ay2 += day;
			}    		
			Generate2d(x, y, ax2, ay2, bx, by);
			Generate2d(x + ax2, y + ay2, ax - ax2, ay - ay2, bx, by);
			return;
		}

		if ((h2 % 2) != 0 && h > 2) {
			bx2 += dbx;
			by2 += dby;
		}

		Generate2d(x, y, bx2, by2, ax2, ay2);
		Generate2d(x + bx2, y + by2, ax, ay, bx - bx2, by - by2);
		Generate2d(x + (ax - dax) + (bx2 - dbx), y + (ay - day) + (by2 - dby), -bx2, -by2, -(ax - ax2), -(ay - ay2));    		
	}

        private void Run()
        {
            /* Dithers all pixels of the image in sequence using
             * the Gilbert path, and distributes the error in
             * a sequence of 9 pixels.
             */
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
			
            if (width >= height)
                Generate2d(0, 0, width, 0, 0, height);
            else
                Generate2d(0, 0, 0, height, width, 0);
        }

        public static int[] Dither(int width, int height, int[] pixels, Color[] palette, DitherFn ditherFn, GetColorIndexFn getColorIndexFn)
        {
            var qPixels = new int[pixels.Length];
            new GilbertCurve(width, height, pixels, palette, qPixels, ditherFn, getColorIndexFn).Run();
            return qPixels;
        }
    }
}
