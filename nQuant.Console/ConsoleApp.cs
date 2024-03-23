using System;
using System.Diagnostics;
using System.Drawing;
using System.Drawing.Imaging;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Text;

using nQuant.Master.Ga;
using System.Collections.Generic;
using PnnQuant;

namespace nQuant
{
	public class ConsoleApp
	{
		private static int maxColors = 256;
		private static bool dither = true;
		private static string targetPath = string.Empty;

		static ConsoleApp()
		{
			AppContext.SetSwitch("System.Drawing.EnableUnixSupport", true);
		}

		private static string ProcessArgs(string[] args)
		{
			string strAlgor = null;
			for (var index = 1; index < args.Length; ++index)
			{
				var currentArg = args[index].ToUpper();
				if (currentArg.Length > 1
				  && (currentArg.StartsWith("-", StringComparison.Ordinal)
				  || currentArg.StartsWith("–", StringComparison.Ordinal)
				  || currentArg.StartsWith("/", StringComparison.Ordinal)))
				{
					if (index >= args.Length - 1)
					{
						PrintUsage();
						Environment.Exit(1);
						return null;
					}

					currentArg = currentArg.Substring(1);
					switch (currentArg)
					{
						case "M":
							if (!Int32.TryParse(args[index + 1], out maxColors))
							{
								PrintUsage();
								Environment.Exit(1);
							}
							break;

						case "D":
							string strDither = args[index + 1].ToUpper();
							if (!(strDither == "Y" || strDither == "N"))
							{
								PrintUsage();
								Environment.Exit(1);
								break;
							}
							dither = strDither == "Y";
							break;
						case "O":
							targetPath = args[index + 1];
							break;

						case "A":
							strAlgor = args[index + 1].ToUpper();
							if (strAlgor != "OTSU" && !strAlgor.StartsWith("PNN"))
							{
								PrintUsage();
								Environment.Exit(1);
								break;
							}

							if (strAlgor == "OTSU")
							{
								maxColors = 2;
								return strAlgor;
							}
							break;

						default:
							PrintUsage();
							Environment.Exit(1);
							break;
					}
				}
			}
			return strAlgor;
		}

		private static void PrintUsage()
		{
			System.Console.WriteLine();
			System.Console.WriteLine("usage: nQuant <input image path> [options]");
			System.Console.WriteLine();
			System.Console.WriteLine("Valid options:");
			System.Console.WriteLine("  /m : Max Colors (pixel-depth) - Maximum number of colors for the output format to support. The default is 256 (8-bit).");
			System.Console.WriteLine("  /d : Dithering or not? y or n.");
			System.Console.WriteLine("  /a : Algorithm used - Choose one of them, otherwise give you the default pnn from [pnn, pnnlab, pnnlab+, otsu]");
			System.Console.WriteLine("  /o : Output image file path. The default is <source image path directory>\\<source image file name without extension>-quant<Max colors>.png");
		}


		private static void DoProcess(Bitmap source, string algorithm)
		{
			var copyright = Assembly.GetExecutingAssembly().GetCustomAttributes(typeof(AssemblyCopyrightAttribute), false)[0] as AssemblyCopyrightAttribute;
			var quantizer = algorithm == "PNN" ? new PnnQuant.PnnQuantizer() : new PnnQuant.PnnLABQuantizer();
			switch (algorithm)
			{
				case "PNNLAB":
					System.Console.WriteLine("nQuant Version {0} C# Color Quantizer. An adaptation of fast pairwise nearest neighbor based algorithm.", Assembly.GetExecutingAssembly().GetName().Version);
					System.Console.WriteLine(copyright.Copyright);
					using (var dest = quantizer.QuantizeImage(source, PixelFormat.Undefined, maxColors, dither))
					{
						dest.Save(targetPath, ImageFormat.Png);
						System.Console.WriteLine("Converted image: " + Path.GetFullPath(targetPath));
					}
				break;
				case "PNNLAB+":
					System.Console.WriteLine("nQuant Version {0} C# Color Quantizer. An adaptation of fast pairwise nearest neighbor based parallel genetic algorithm.", Assembly.GetExecutingAssembly().GetName().Version);
					System.Console.WriteLine(copyright.Copyright);
					var alg = new APNsgaIII<PnnQuant.PnnLABGAQuantizer>(new PnnQuant.PnnLABGAQuantizer(new PnnQuant.PnnLABQuantizer(), new List<Bitmap> { source }, maxColors));
					alg.Run(999, -Double.Epsilon);
					using (var pGAq = alg.Result) {
						System.Console.WriteLine("\n" + pGAq.Result);
						using (var dest = pGAq.QuantizeImage(dither)[0])
						{
							dest.Save(targetPath, ImageFormat.Png);
							System.Console.WriteLine("Converted image: " + Path.GetFullPath(targetPath));
						}
					}
				break;
				case "OTSU":
					System.Console.WriteLine("nQuant Version {0} C# Color Quantizer. An adaptation of Otsu's Image Segmentation Method.", Assembly.GetExecutingAssembly().GetName().Version);
					System.Console.WriteLine(copyright.Copyright);
					using (var dest = new OtsuThreshold.Otsu().ConvertGrayScaleToBinary(source))
					{
						dest.Save(targetPath, ImageFormat.Png);
						System.Console.WriteLine("Converted black and white image: " + Path.GetFullPath(targetPath));
					}
				break;
				default:
					System.Console.WriteLine("nQuant Version {0} C# Color Quantizer. An adaptation of fast pairwise nearest neighbor based algorithm.", Assembly.GetExecutingAssembly().GetName().Version);
					System.Console.WriteLine(copyright.Copyright);
					using (var dest = quantizer.QuantizeImage(source, PixelFormat.Undefined, maxColors, dither))
					{
						dest.Save(targetPath, ImageFormat.Png);
						System.Console.WriteLine("Converted image: " + Path.GetFullPath(targetPath));
					}
				break;
			}
		}

		private static void DoProcess(string sourcePath, string algorithm)
		{
			if (!File.Exists(sourcePath))
			{
				System.Console.WriteLine("The source file {0} you specified does not exist.", sourcePath);
				return;
			}

			if (string.IsNullOrEmpty(targetPath))
				targetPath = Path.GetFileNameWithoutExtension(sourcePath) + "-" + algorithm + "quant" + maxColors + ".png";

			Console.OutputEncoding = Encoding.UTF8;
			Stopwatch stopwatch = Stopwatch.StartNew();

			using (var bitmap = new Bitmap(sourcePath))
			{
				try
				{
					DoProcess(bitmap, algorithm);
				}
				catch (Exception q)
				{
#if (DEBUG)
					System.Console.WriteLine(q.StackTrace);
#else
					System.Console.WriteLine(q.Message);
					System.Console.WriteLine("Incorrect pixel format: {0} for {1} colors.", bitmap.PixelFormat.ToString(), maxColors);
#endif
				}
			}
			System.Console.WriteLine(@"Completed in {0:s\.fff} secs with peak memory usage of {1}.", stopwatch.Elapsed, Process.GetCurrentProcess().PeakWorkingSet64.ToString("#,#"));
		}


		private static void OutputImages(string sourceDir)
		{
			Console.OutputEncoding = Encoding.UTF8;
			var copyright = Assembly.GetExecutingAssembly().GetCustomAttributes(typeof(AssemblyCopyrightAttribute), false)[0] as AssemblyCopyrightAttribute;
			Stopwatch stopwatch = Stopwatch.StartNew();

			try {
				System.Console.WriteLine("nQuant Version {0} C# Color Quantizer. An adaptation of fast pairwise nearest neighbor based parallel genetic algorithm.", Assembly.GetExecutingAssembly().GetName().Version);
				System.Console.WriteLine(copyright.Copyright);

				if (string.IsNullOrEmpty(targetPath))
				{
					targetPath = sourceDir;
				}
				else if (!Directory.Exists(targetPath))
					Directory.CreateDirectory(targetPath);

				var paths = Directory.GetFiles(sourceDir).OrderBy(p => p).ToList();
				var bitmaps = paths.Select(p => {
					try
					{
						var path = Path.GetFullPath(p);
						return new Bitmap(path);
					}
					catch {
						return null;
					}
				}).Where(bitmap => bitmap != null).ToList();

				var alg = new APNsgaIII<PnnLABGAQuantizer>(new PnnLABGAQuantizer(new PnnLABQuantizer(), bitmaps, maxColors));
				alg.Run(999, -Double.Epsilon);
                using var pGAq = alg.Result;
                System.Console.WriteLine("\n" + pGAq.Result);
                var destPath = string.Empty;
                var imgs = pGAq.QuantizeImage(dither);
                if (maxColors > 256)
                {
                    for (int i = 0; i < imgs.Count; ++i)
                    {
                        var fname = Path.GetFileNameWithoutExtension(paths[i]);
                        destPath = Path.Combine(targetPath, fname) + " - PNNLAB+quant" + maxColors + ".png";
                        imgs[i].Save(destPath, ImageFormat.Png);
                        System.Console.WriteLine("Converted image: " + Path.GetFullPath(destPath));
                    }
                }
                else
                {
                    var fname = Path.GetFileNameWithoutExtension(paths[0]);
                    destPath = Path.Combine(targetPath, fname) + " - PNNLAB+quant" + maxColors + ".gif";
                    var gifWriter = new GifWriter(destPath, 850);
                    gifWriter.AddImages(imgs);
                }

                System.Console.WriteLine("Converted image: " + Path.GetFullPath(destPath));
            }
			catch (Exception q)
			{
				System.Console.WriteLine(q.StackTrace);
			}
			System.Console.WriteLine(@"Completed in {0:s\.fff} secs with peak memory usage of {1}.", stopwatch.Elapsed, Process.GetCurrentProcess().PeakWorkingSet64.ToString("#,#"));
		}
		
		public static void Main(string[] args)
		{
			string algorithm = "PNNLAB";
#if DEBUG
			var sourcePath = @"samples\SE5x9.jpg";
			maxColors = 256;
#else
			if (args.Length < 1)
			{
				PrintUsage();
				Environment.Exit(1);
			}
			var sourcePath = args[0];
			algorithm = ProcessArgs(args);
#endif
			if (Directory.Exists(sourcePath))
				OutputImages(sourcePath);
			else
				DoProcess(sourcePath, algorithm);
		}

	}
}