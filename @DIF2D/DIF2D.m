classdef DIF2D
    
    properties (SetAccess = private, GetAccess = public)
        adjoint;
        dim;
        Xfm;
        XfmInv;
    end % properties
    
    methods
        function this = DIF2D(dim)
            this.adjoint = 0;
            this.dim=dim;
            this.Xfm=eye(this.dim);
            this.Xfm(1, 1)=0;
            this.Xfm=[this.Xfm; this.Xfm(1,:)];
            this.Xfm(1, :)=[];
            this.Xfm=eye(this.dim)-this.Xfm;
            this.XfmInv=inv(this.Xfm);
        end
        
        function this = ctranspose(this)
            this.adjoint = xor(this.adjoint,1);
        end
                
        function res = mtimes(this,b)            
            if this.adjoint
                b=reshape(b, this.dim);
%                res=this.XfmInv*b; % H
                res=b*this.XfmInv;  % V
            else
                b=reshape(b, this.dim);
%                res=this.Xfm*b;    % H
                res=b*this.Xfm;     % V
            end
        end        
    end  % methods
end

