function res = mtimes(this,b)

if this.adjoint
    b=reshape(b, this.sdSize);
    res = IWT2_PO(real(b),this.wavScale,this.qmf) + 1i*IWT2_PO(imag(b),this.wavScale,this.qmf);
else
    b=reshape(b, this.imSize);
    res = FWT2_PO(real(b),this.wavScale,this.qmf) + 1i* FWT2_PO(imag(b),this.wavScale,this.qmf);
end


