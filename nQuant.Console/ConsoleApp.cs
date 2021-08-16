using System;
using System.Diagnostics;
using System.Drawing;
using System.Drawing.Imaging;
using System.IO;
using System.Reflection;

namespace nQuant
{
    public class ConsoleApp
    {
        private static int maxColors = 256;
        private static bool dither = true;
        private static string targetPath = string.Empty;

        public static void Main(string[] args)
        {
            System.Console.WriteLine("nQuant Version {0} C# Color Quantizer. An adaptation of fast pairwise nearest neighbor based algorithm.", Assembly.GetExecutingAssembly().GetName().Version);
            System.Console.WriteLine("Copyright (C) 2018 - 2021 Miller Cy Chan.");

#if DEBUG
            var sourcePath = @"samples\SE5x9.jpg";
            maxColors = 128;
#else
            if (args.Length < 1)
            {
               PrintUsage();
               Environment.Exit(1);
            }
            var sourcePath = args[0];
            ProcessArgs(args);
#endif
            if (!File.Exists(sourcePath))
            {
                System.Console.WriteLine("The source file you specified does not exist.");
                Environment.Exit(1);
            }
            if (string.IsNullOrEmpty(targetPath))
            {
                var lastDot = sourcePath.LastIndexOf('.');
                if (lastDot == -1)
                    lastDot = sourcePath.Length;
                targetPath = sourcePath.Substring(0, lastDot) + "-quant" + maxColors + ".png";
            }

            Stopwatch stopwatch = Stopwatch.StartNew();
            PnnQuant.PnnQuantizer quantizer = new PnnQuant.PnnLABQuantizer();
            using (var bitmap = new Bitmap(sourcePath))
            {
                try
                {
                    using (var dest = quantizer.QuantizeImage(bitmap, PixelFormat.Undefined, maxColors, dither))
                    {
                        dest.Save(targetPath, ImageFormat.Png);
                        System.Console.WriteLine("Converted image: " + targetPath);
                    }
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

        private static void ProcessArgs(string[] args)
        {
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
                        return;
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

                        default:
                            PrintUsage();
                            Environment.Exit(1);
                            break;
                    }
                }
            }
        }

        private static void PrintUsage()
        {
            System.Console.WriteLine();
            System.Console.WriteLine("usage: nQuant <input image path> [options]");
            System.Console.WriteLine();
            System.Console.WriteLine("Valid options:");
            System.Console.WriteLine("  /m : Max Colors (pixel-depth) - Maximum number of colors for the output format to support. The default is 256 (8-bit).");
            System.Console.WriteLine("  /d : Dithering or not? y or n.");
            System.Console.WriteLine("  /o : Output image file path. The default is <source image path directory>\\<source image file name without extension>-quant<Max colors>.png");
        }

    }
}
