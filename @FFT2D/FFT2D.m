function  obj = FFT2D(imSize,mdSize, phase,mode)

%obj = p2DFT(imSize, mdSize, [ ,phase,mode])
%
%
%	Implementation of partial Fourier operator.
%	
%	input:
%			imSize - the image size (1x2)
%           mdSize - The sampling domain size (1x2)
%			phase - Phase of the image for phase correction
%			mode - 1- real, 2-cmplx
%
%	Output:
%			The operator
%
%	(c) Michael Lustig 2007

if nargin <3
	phase = 1;
end
if nargin <4
	mode = 2; % 0 - positive, 1- real, 2-cmplx
end

obj.adjoint = 0;
obj.imSize = imSize;
obj.mdSize = mdSize;
obj.ph = phase;
obj.mode = mode;
obj = class(obj,'FFT2D');

