function res = mtimes(a,b)


if a.adjoint
	res = adjD(b);

else
    b = reshape(b,a.imSize);
	res = D(b);

end



    
