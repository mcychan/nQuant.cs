using System.Drawing;

namespace nQuant.Master
{
    interface Ditherable
    {
		public int GetColorIndex(int pixel);

		public ushort DitherColorIndex(Color[] palette, int pixel, int pos);

	}
}
