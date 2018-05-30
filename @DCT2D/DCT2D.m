classdef DCT2D
    
    properties (SetAccess = private, GetAccess = public)
        adjoint;
        dim;
    end % properties
    
    methods
        function this = DCT2D(dim)
            this.adjoint = 0;
            this.dim=dim;
        end
        
        function this = ctranspose(this)
            this.adjoint = xor(this.adjoint,1);
        end
                
        function res = mtimes(this,b)            
            if this.adjoint
                b=reshape(b, this.dim);
		res=idct2(b);
            else
                b=reshape(b, this.dim);
		res=dct2(b);
            end
        end        
    end  % methods
end

