function f = fft3(x, mrows, ncols, ppages)
%FFT3 Three-dimensional discrete Fourier Transform.
%   FFT3(X) returns the three-dimensional Fourier transform of array X.
%   If X is a vector, the result will have the same orientation.
%
%   FFT3(X,MROWS,NCOLS,PPAGES) pads array X with zeros to size
%   MROWS-by-NCOLS-by-PPAGES before transforming.
%
%   Class support for input X: 
%      float: double, single
%
%   See also FFT, FFT2, FFTN, FFTSHIFT, FFTW, IFFT, IFFT2, IFFT3, IFFTN.

s = size(x);
if length(s)==3
    if nargin==1
        f = fftn(x);
    else
        f = fftn(x,[mrows ncols, ppages]);
    end
else
    if nargin==1
        f = fft(fft(fft(x,[],3),[],2),[],1);
    else
        f = fft(fft(fft(x,ppages,3),ncols,2),mrows,1);
    end
end   
