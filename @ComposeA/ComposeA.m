classdef ComposeA

   properties (SetAccess = public, GetAccess = public)
       Mx;
       Mxt;
       phi;
       psi;
   end % properties
   
   properties (SetAccess = private, GetAccess = public)
       adjoint;
   end % properties
   
   methods
      function obj = ComposeA(Mx, Mxt, phi, psi) 
         obj.adjoint = 0;
         obj.Mx=Mx;
         obj.Mxt=Mxt;
         obj.phi=phi;
         obj.psi=psi;
      end 
      
    function res = ctranspose(this)
        this.adjoint = xor(this.adjoint,1);
        res = this;      
    end
    
    function res = mtimes(this,b)
        if this.adjoint
            res=this.Mxt(b);
            res=this.phi'*res;
            res=this.psi*res;
            res=res(:);
        else
            res=this.psi'*b;
            res=this.phi*res;
            res=this.Mx(res);
            res=res(:);
        end
    end
    
   end  % methods
end
