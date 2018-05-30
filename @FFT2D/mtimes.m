function res = mtimes(this,b)

if this.adjoint
    b=reshape(b, this.mdSize);
    
    res = zpad(b,this.imSize(1),this.imSize(2));
    res = ifft2c(res);
    res = res.*conj(this.ph);
    switch this.mode
    case 0
		res = real(res);
   	case 1
		res = real(res);
    end
else
    b=reshape(b, this.imSize);

    switch this.mode
    case 0
		b = real(b);
   	case 1
		b = real(b);
    end
    
    b = b.*this.ph; % phase correct
    res = fft2c(b);
    res = crop(res,this.mdSize(1),this.mdSize(2));
end



    
