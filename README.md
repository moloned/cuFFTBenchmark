# FFTBenchmark
Snippet of code to measure performance of cuFFT for images of different sizes

*Build*

%> make bench FROMTO=\<fft_type\>

<fft_type> should be one of Z2Z, Z2D, D2Z, C2C, C2R, R2C to specify the precision and any assumptions of Hermiticity.


*Run*
- Modify in.txt to specify sizes

     Each line should contain 2 integers to specify the image size
  
    The last line of in.txt should be "0 0"

- %> make test FROMTO=\<fft type\>


