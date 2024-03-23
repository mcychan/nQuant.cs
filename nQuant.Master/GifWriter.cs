using System;
using System.Collections.Generic;
using System.Drawing;
using System.Drawing.Imaging;
using System.IO;
using System.Runtime.Serialization;

namespace nQuant
{
	public class GifWriter {
		private const int PropertyTagFrameDelay = 0x5100;
		private const int PropertyTagLoopCount = 0x5101;
		private const short PropertyTagTypeLong = 4;
		private const short PropertyTagTypeShort = 3;

		private const int UintBytes = 4;

		private readonly bool _loop;
		private readonly string _destPath;
		private readonly uint _delay;

		public GifWriter(String destPath, uint delay = 850, bool loop = true)
		{
			_destPath = destPath;
			_loop = loop;
			_delay = delay;
		}

		private ImageCodecInfo GetEncoder(ImageFormat format)
		{
			var codecs = ImageCodecInfo.GetImageDecoders();
			foreach (ImageCodecInfo codec in codecs) {
				if (codec.FormatID == format.Guid)
					return codec;
			}
			return null;
		}

		private static void SetFrameDelay(Bitmap firstBitmap, uint delay, int count)
		{
			// PropertyItem for the frame delay (apparently, no other way to create a fresh instance).
			var frameDelay = (PropertyItem) FormatterServices.GetUninitializedObject(typeof(PropertyItem));
			frameDelay.Id = PropertyTagFrameDelay;
			frameDelay.Type = PropertyTagTypeLong;
			// Length of the value in bytes.
			frameDelay.Len = count * UintBytes;
			// The value is an array of 4-byte entries: one per frame.
			// Every entry is the frame delay in 1/100-s of a second, in little endian.
			frameDelay.Value = new byte[count * UintBytes];
			// E.g., here, we're setting the delay of every frame to 1 second.
			var frameDelayBytes = BitConverter.GetBytes(delay / 10);
			for (int j = 0; j < count; ++j)
				Array.Copy(frameDelayBytes, 0, frameDelay.Value, j * UintBytes, UintBytes);
			firstBitmap.SetPropertyItem(frameDelay);
		}

		private static void SetLoop(Bitmap firstBitmap, bool loop)
		{
			if (!loop)
				return;

			// PropertyItem for the number of animation loops.
			var loopPropertyItem = (PropertyItem)FormatterServices.GetUninitializedObject(typeof(PropertyItem));
			loopPropertyItem.Id = PropertyTagLoopCount;
			loopPropertyItem.Type = PropertyTagTypeShort;
			loopPropertyItem.Len = 2;
			// 0 means to animate forever.
			loopPropertyItem.Value = BitConverter.GetBytes((ushort)0);
			firstBitmap.SetPropertyItem(loopPropertyItem);
		}

		public void AddImages(List<Bitmap> bitmaps)
		{
			using var fs = new FileStream(_destPath, FileMode.Create);
			var firstBitmap = bitmaps[0];
			var gifEncoder = GetEncoder(ImageFormat.Gif);

			// Params of the first frame.
			var encoderParams = new EncoderParameters(1);
			encoderParams.Param[0] = new EncoderParameter(Encoder.SaveFlag, (long)EncoderValue.MultiFrame);
			SetFrameDelay(firstBitmap, _delay, bitmaps.Count);
			SetLoop(firstBitmap, _loop);
			firstBitmap.Save(fs, gifEncoder, encoderParams);

			// Params of other frames.
			encoderParams.Param[0] = new EncoderParameter(Encoder.SaveFlag, (long)EncoderValue.FrameDimensionTime);
			for (int i = 1; i < bitmaps.Count; ++i)
				firstBitmap.SaveAdd(bitmaps[i], encoderParams);

			// Params for the finalizing call.
			encoderParams.Param[0] = new EncoderParameter(Encoder.SaveFlag, (long)EncoderValue.Flush);
			firstBitmap.SaveAdd(encoderParams);
		}

	}
}