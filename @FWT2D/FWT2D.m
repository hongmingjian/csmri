function res = FWT2D(imSize, sdSize, filterType, filterSize, wavScale)
%
% implements a wavelet operator
%
% (c) Michael Lustig 2007

res.adjoint = 0;
res.qmf = MakeONFilter(filterType, filterSize);
res.wavScale = wavScale;
res.imSize = imSize;
res.sdSize = sdSize;
res = class(res,'FWT2D');
