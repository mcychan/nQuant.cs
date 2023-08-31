using System;
using System.Collections.Generic;
using System.Drawing;
using System.Drawing.Imaging;
using System.IO;
using System.Linq;
using System.Runtime.Serialization;

namespace nQuant
{
	public class GifWriter	{	
		private const int PropertyTagFrameDelay = 0x5100;
		private const int PropertyTagLoopCount = 0x5101;
		private const short PropertyTagTypeLong = 4;
		private const short PropertyTagTypeShort = 3;

		private const int UintBytes = 4;

		private bool _hasAlpha, _loop;
		private string _destPath;
		private PropertyItem _frameDelay, _loopPropertyItem;

		public GifWriter(String destPath, bool hasAlpha, int count = 1, uint delay = 850, bool loop = true)
		{
			_destPath = destPath;
            _hasAlpha = hasAlpha;
            _loop = loop;            

            // PropertyItem for the frame delay (apparently, no other way to create a fresh instance).
            _frameDelay = (PropertyItem)FormatterServices.GetUninitializedObject(typeof(PropertyItem));
            _frameDelay.Id = PropertyTagFrameDelay;
            _frameDelay.Type = PropertyTagTypeLong;
            // Length of the value in bytes.
            _frameDelay.Len = count * UintBytes;
            // The value is an array of 4-byte entries: one per frame.
            // Every entry is the frame delay in 1/100-s of a second, in little endian.
            _frameDelay.Value = new byte[count * UintBytes];
            // E.g., here, we're setting the delay of every frame to 1 second.
            var frameDelayBytes = BitConverter.GetBytes(delay / 10);
            for (int j = 0; j < count; ++j)
                Array.Copy(frameDelayBytes, 0, _frameDelay.Value, j * UintBytes, UintBytes);

            // PropertyItem for the number of animation loops.
            _loopPropertyItem = (PropertyItem)FormatterServices.GetUninitializedObject(typeof(PropertyItem));
            _loopPropertyItem.Id = PropertyTagLoopCount;
            _loopPropertyItem.Type = PropertyTagTypeShort;
            _loopPropertyItem.Len = 2;
            // 0 means to animate forever.
            _loopPropertyItem.Value = BitConverter.GetBytes((ushort)0);
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

		public void AddImages(List<Bitmap> bitmaps)
		{
            using var fs = new FileStream(_destPath, FileMode.Create);
            var firstBitmap = bitmaps[0];            
            var gifEncoder = GetEncoder(ImageFormat.Gif);

            // Params of the first frame.
            var encoderParams1 = new EncoderParameters(1);
            encoderParams1.Param[0] = new EncoderParameter(Encoder.SaveFlag, (long)EncoderValue.MultiFrame);
            firstBitmap.SetPropertyItem(_frameDelay);
            if (_loop)
                firstBitmap.SetPropertyItem(_loopPropertyItem);
            firstBitmap.Save(fs, gifEncoder, encoderParams1);

            for (int i = 1; i < bitmaps.Count; ++i)
            {
                var bitmap = bitmaps[i];
                // Params of other frames.
                var encoderParamsN = new EncoderParameters(1);
                encoderParamsN.Param[0] = new EncoderParameter(Encoder.SaveFlag, (long)EncoderValue.FrameDimensionTime);
                firstBitmap.SaveAdd(bitmap, encoderParamsN);
            }

            // Params for the finalizing call.
            var encoderParamsFlush = new EncoderParameters(1);
            encoderParamsFlush.Param[0] = new EncoderParameter(Encoder.SaveFlag, (long)EncoderValue.Flush);
            firstBitmap.SaveAdd(encoderParamsFlush);            
        }

	}
}