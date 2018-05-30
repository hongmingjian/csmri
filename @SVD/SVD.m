classdef SVD
    
    properties (SetAccess = private, GetAccess = public)
        adjoint;
        X;
    end % properties

    properties (SetAccess = private, GetAccess = public)
        U;
        V;
    end % properties
    
    methods
        function obj = SVD(X)
            obj.adjoint = 0;
            obj.X=X;
            [obj.U, S, obj.V]=svd(obj.X);
        end
        
        function this = ctranspose(this)
            this.adjoint = xor(this.adjoint,1);
        end
                
        function res = mtimes(this,b)            
            if this.adjoint
                b=reshape(b, size(this.X));
                res=this.U*b*this.V';
            else
                b=reshape(b, size(this.X));
                res=this.U'*b*this.V;
            end
        end
        
    end  % methods
end

