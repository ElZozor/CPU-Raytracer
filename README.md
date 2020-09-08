# CPU-Raytracer
A little faculty project that implements a Raytracer which runs on the CPU.

The code is in C/C++.

## Usage

### Build

Just launch a terminal in the root directory and run the `make` command.  
Once done you can use the program with `./mrt [options]`.

### Needed

Two parameters are needed to run the program.  
First you need to provide an output file in which the render will be saved.  
This can be done with the `-filename` option.

Then, you need to provide which scene you want to render.  
This can be done with the `-scenenumber` parameter. The default value is 0.

### Options

Some options are available.  
 - width: Integer -> It sets the resulting iamge width.
 - height: Integer -> It sets the resulting image height.
 - aa: Integer -> The amount of Anti-aliasing for the render.
 - aamode: String ('none' | 'luma' | 'sobel') -> Which optimisation you want to apply to the anti-aliasing (see bellow).
 - focaldist: Float -> The focal distance for the render.
 - focalrange: Float -> The range within the render will be sharp.

#### AA optimizations  
This parameter can be sets to one of the following values : none, luma, sobel.  
 - `none` means that the anti-aliasing is just a basic superslamping.  
 - `luma` means that the objects edges are detected by their luminosity. The antialising will only be applied on them.
 - `sobel` means that the render will be parsed with the sobel operator to detect edges and a superslamping will be done in these areas.
 
## Some renders

<img src="https://raw.githubusercontent.com/ElZozor/CPU-Raytracer/master/previews/paysage.png"/>
<img src="https://raw.githubusercontent.com/ElZozor/CPU-Raytracer/master/previews/cube_minecraft.png"/>
