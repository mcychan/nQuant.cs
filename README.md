# nQuant.cs Color Quantizer
Fast pairwise nearest neighbor based algorithm with C# console

nQuant.cs is a C# color quantizer producing high quality 256 color 8 bit PNG images using an algorithm optimized for the highest quality possible.

Another advantage of nQuant is that it is a .net library that you can integrate nicely with your own C# code while many of the popular quantizers only provide command line implementations. nQuant.cs also provides a command line wrapper in case you want to use it from the command line. Less artifacts by using advanced dithering techniques such as Generalized Hilbert ("gilbert") space-filling curve and partial Blue noise distribution to diffuse the minimized errors. More importantly, a parallel genetic algorithm called PNNLAB+ is proposed for converting a sequence of similar images under the same palette.
If you are using C#, you would call nQuant as follows:

```cs
    bool dither = true;
    var quantizer = new PnnQuant.PnnQuantizer();
    using(var bitmap = new Bitmap(sourcePath))
    {
        try
        {                    
            using (var dest = quantizer.QuantizeImage(bitmap, pixelFormat, maxColors, dither))
            {
                dest.Save(targetPath, ImageFormat.Png);
                System.Console.WriteLine("Converted image: " + targetPath);
            }
        }
        catch (Exception q)
        {
            System.Console.WriteLine(q.StackTrace);
        }
    }
```

OTSU method (OTSU) is a global adaptive binarization threshold image segmentation algorithm. This algorithm takes the maximum inter class variance between the background and the target image as the threshold selection rule.
```cs
    var quantizer = new OtsuThreshold.Otsu();
    using(var bitmap = new Bitmap(sourcePath))
    {
        try
        {                    
            using (var dest = quantizer.ConvertGrayScaleToBinary(bitmap))
            {
                dest.Save(targetPath, ImageFormat.Png);
                System.Console.WriteLine("Converted black and white image: " + targetPath);
            }
        }
        catch (Exception q)
        {
            System.Console.WriteLine(q.StackTrace);
        }
    }
```
<p>Example image:<br /><img src="https://user-images.githubusercontent.com/26831069/142559831-f8f6f2ce-487e-4353-8aa1-7845706e7833.png" /></p>
<p>Resulted image:<br /><pre><img src="https://user-images.githubusercontent.com/26831069/142559920-88143e07-2787-46a2-a07c-cccf5a39065a.png" /></pre></p>

If you are using the command line. Assuming you are in the same directory as nQuant.exe and nQuant.Master.dll, you would enter:
`nQuant yourImage.jpg /o yourNewImage.png`

To switch algorithms, `/a otsu` can perform the above black and white conversion.

nQuant will quantize yourImage.jpg and create yourNewImage.png in the same directory.

There are a few configuration arguments you can optionaly use to try and influence how the image gets quantized. These are explained in the console application.
