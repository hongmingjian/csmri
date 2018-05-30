function  res = TVOP(imSize)

%res = TVOP()
%
% Implements a spatial finite-differencing operator.
%
% (c) Michael Lustig 2007

res.adjoint = 0;
res.imSize=imSize;
res = class(res,'TVOP');

