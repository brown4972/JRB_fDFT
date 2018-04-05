function x = ifft3(varargin)
%IFFT3 Three-dimensional inverse discrete Fourier transform.
%   IFFT3(F) returns the three-dimensional inverse Fourier transform of
%   array F.  If F is a vector, the result will have the same orientation.
%
%   IFFT3(F,MROWS,NCOLS,PPAGES) pads matrix F with zeros to size 
%   MROWS-by-NCOLS-by-PPAGES before transforming.
%
%   IFFT3(..., 'symmetric') causes IFFT3 to treat F as conjugate symmetric
%   in three dimensions so that the output is purely real.  This option is
%   useful when F is not exactly conjugate symmetric merely because of
%   round-off error.  See the reference page for the specific mathematical
%   definition of this symmetry.
%
%   IFFT3(..., 'nonsymmetric') causes IFFT3 to make no assumptions about the
%   symmetry of F.
%
%   Class support for input F:
%      float: double, single
%
%   See also FFT, FFT2, FFT3, FFTN, FFTSHIFT, FFTW, IFFT, IFFT2, IFFTN.


narginchk(1,5)

f = varargin{1};
s = size(f);
m_in = s(1);
n_in = s(2);
p_in = s(3);
num_inputs = numel(varargin);
symmetry = 'nonsymmetric';
if ischar(varargin{end})
    symmetry = varargin{end};
    num_inputs = num_inputs - 1;
end

if num_inputs == 1
    m_out = m_in;
    n_out = n_in;
    p_out = p_in;

elseif num_inputs == 2 || num_inputs == 3
    error('If you specify MROWS, you also have to specify NCOLS and PPAGES.')
    
else
    m_out = double(varargin{2});
    n_out = double(varargin{3});
    p_out = double(varargin{4});
end

if ~isa(f, 'float')
    f = double(f);
end

if (m_out ~= m_in) || (n_out ~= n_in) || (p_out ~= p_in)
    out_size = size(f);
    out_size(1) = m_out;
    out_size(2) = n_out;
    out_size(3) = p_out;
    f2 = zeros(out_size, class(f));
    mm = min(m_out, m_in);
    nn = min(n_out, n_in);
    pp = min(p_out, p_in);
    f2(1:mm, 1:nn, 1:pp, :) = f(1:mm, 1:nn, 1:pp, :);
    f = f2;
end

if length(s) == 3
    x = ifftn(f, symmetry);
else
    x = ifft(ifft(ifft(f, [], 3), [], 2), [], 1, symmetry);
end   


