# A Simple ispc Example
Here is a walkthrough of a simple example of using ispc to compute an image of the Mandelbrot set. The full source code for this example is in the examples/mandelbrot directory of the ispc distribution; below is a lightly modified version of that example. (For example, we have elided some of the details related to benchmarking and writing the final output image that are in the actual mandelbrot.cpp file there.)First we have some code in the main program file, mandelbrot.cpp, which is implemented in C++. We mostly have a regular C++ source file that sets a number of variables that drive the overall computation and allocates memory to store the result of the Mandelbrot set computation. The call to mandelbrot_ispc is a regular subroutine call, though in this case it's calling out to code that was written in ispc rather than in C/C++. Because the ispc compiler generates regular object files using standard C/C++ calling conventions, this kind of close interoperation between the two languages is straightforward.int main() {
    unsigned int width = 768, height = 512;
    float x0 = -2.f, x1 = 1.f;
    float y0 = -1.f, y1 = 1.f;
    int maxIterations = 256;
    int *buf = new int[width*height];

    mandelbrot_ispc(x0, y0, x1, y1, width, height, maxIterations, buf);

    // write output...
}
Here is the start of the definition of the mandelbrot_ispc() function that the ispc compiler compiles. (This function is defined in the mandelbrot.ispc file.) There are a four main things to notice here.The export qualifier indicates that this function should be made available to be called by C/C++ code, not just from other ispc functions.This function directly takes parameters passed by the C/C++ code as regular function parameters. No API calls or a runtime library are necessary to set parameter values to pass to this function.The parameters are all tagged with the uniform qualifier. This qualifier indicates that a variable has the same value for all of the concurrently-executing program instances; it is discussed in more detail in the ispc User's Manual. Any values passed from the application to the ispc program must be uniform.The pointer to the output buffer is passed as an unsized array. Just like in C, arrays in function calls are converted to pointers to the first element of the array. (This parameter could also be declared to just be a pointer.) Like the other function parameters, this buffer is passed from the application to the ispc code without any copying or reformatting, as would be required with a driver model based implementation.export void mandelbrot_ispc(uniform float x0, uniform float y0,
                            uniform float x1, uniform float y1,
                            uniform int width, uniform int height,
                            uniform int maxIterations,
                            uniform int output[]) {
The function implementation starts by computing the step between pixels in the x and y directions.    float dx = (x1 - x0) / width;
    float dy = (y1 - y0) / height;
Next, the function loops over the pixels of the image with a doubly-nested loop. The outer loop is a regular loop over scanlines, while the inner loop introduces a new looping construct provided by ispc, foreach. In general, foreach specifies parallel iteration over an n-dimensional range of integer values. In this case, it specifies parallel iteration over a 1-D range of pixels in a scanline. Each time through the loop, the iteration variable i takes on values corresponding to a small range of the domain.For example, if ispc was generating code to run in 8-wide parallel fashion for an 8-wide AVX SIMD unit, i would have the values (0,1,2,3,4,5,6,7) across the group of executing program instances the first time the loop executed, (8,9,10,11,12,13,14,15) the second time the loop executed, and so forth.    for (uniform int j = 0; j < height; j++) {
        foreach (i = 0 ... width) {
Inside these loops, we first find the floating-point x and y values that correspond to the pixels at which to compute the Mandelbrot iteration counts.            float x = x0 + i * dx;
            float y = y0 + j * dy;
Once it has figured out the coordinates at which to compute the number of iterations, the program finds the appropriate index into the output array for them and calls a subroutine to do the actual computation.            int index = j * width + i;
            output[index] = mandel(x, y, maxIterations);
        }
    }
}
Here is the implementation of the mandel() routine. It takes real and imaginary coordinates and a maximum iteration count and applies the iterative computation that defines the Mandelbrot set. Although this appears to be a serial function, recall that ispc's execution model is SPMD; the program compiles down to run in parallel across the elements of the SIMD unit.static inline int mandel(float c_re, float c_im, int count) {
    float z_re = c_re, z_im = c_im;
    int i;
    for (i = 0; i < count; ++i) {
        if (z_re * z_re + z_im * z_im > 4.)
            break;

        float new_re = z_re * z_re - z_im * z_im;
        float new_im = 2.f * z_re * z_im;
        z_re = c_re + new_re;
        z_im = c_im + new_im;
    }
    return i;
}
To compile this example, install the ispc compiler in a directory that is in your PATH, and run make on Mac OS X or Linux, or build using the mandelbrot.vcxproj build file with MSVC 2010.% make
ispc -O2 --target=avx mandelbrot.ispc -o objs/mandelbrot_ispc.o -h objs/mandelbrot_ispc.h
g++ mandelbrot.cpp -Iobjs/ -O3 -Wall -c -o objs/mandelbrot.o
g++ mandelbrot_serial.cpp -Iobjs/ -O3 -Wall -c -o objs/mandelbrot_serial.o
g++ -Iobjs/ -O3 -Wall -o mandelbrot objs/mandelbrot.o objs/mandelbrot_ispc.o objs/mandelbrot_serial.o -lm
%
After the executable is built, running it benchmarks a serial implementation and the ispc version. On a Second Generation Intel® Core i7 CPU, doing this computation in parallel across the SIMD lanes of a single core gives a substantial speedup:% ./mandelbrot
[mandelbrot ispc]:         [58.570] million cycles
[mandelbrot serial]:     [335.005] million cycles
(5.72x speedup from ISPC)
%
Here is a downscaled version of the image that the program generates:Note that all of the exposition up until now has been related to harvesting the available parallelism from the SIMD unit in a single CPU core. To fully use the system's computational resources also requires parallelism across the CPU cores.Thanks to the fact that ispc cleanly inter-operates with application C/C++ code through a regular function call, one option would be to rewrite the application code to be multi-threaded and to then have parallel tasks that were in turn partially or fully implemented in ispc.The second option for parallelizing across cores is to use the built-in functionality for launching concurrent tasks that ispc offers. The ispc language provides constructs for launching asynchronous tasks that can execute concurrently with the rest of the program; see the directory examples/mandelbrot_tasks in the ispc distribution to see how to parallelize across both cores and the SIMD units.