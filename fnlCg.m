function x = fnlCg(x0,problem, params)
%-----------------------------------------------------------------------
%
% res = fnlCg(x0,problem, params)
%
% implementation of a L1 penalized non linear conjugate gradient reconstruction
%
% The function solves the following problem:
%
% given k-space measurments y, and a fourier operator F the function 
% finds the image x that minimizes:
%
% Phi(x) = ||F* W' *x - y||^2 + lambda1*|x|_1 + lambda2*TV(W'*x) 
%
%
% the optimization method used is non linear conjugate gradient with fast&cheap backtracking
% line-search.
% 
% (c) Michael Lustig 2007
%-------------------------------------------------------------------------
x = x0;


% line search parameters
maxlsiter = params.lineSearchItnlim ;
gradToll = params.gradToll ;
alpha = params.lineSearchAlpha;     beta = params.lineSearchBeta;
t0 = params.lineSearchT0;
k = 0;
t = 1;

% copmute g0  = grad(Phi(x))

g0 = wGradient(x,problem, params);

dx = -g0;


% iterations
while(1)

% backtracking line-search

	% pre-calculate values, such that it would be cheap to compute the objective
	% many times for efficient line-search
	[FTXFMtx, FTXFMtdx, DXFMtx, DXFMtdx] = preobjective(x, dx, problem);
	f0 = objective(FTXFMtx, FTXFMtdx, DXFMtx, DXFMtdx,x,dx, 0, problem, params);
	t = t0;
        [f1, ERRobj, RMSerr]  =  objective(FTXFMtx, FTXFMtdx, DXFMtx, DXFMtdx,x,dx, t, problem, params);
	
	lsiter = 0;

	while (f1 > f0 - alpha*t*abs(g0(:)'*dx(:)))^2 & (lsiter<maxlsiter)
		lsiter = lsiter + 1;
		t = t * beta;
		[f1, ERRobj, RMSerr]  =  objective(FTXFMtx, FTXFMtdx, DXFMtx, DXFMtdx,x,dx, t, problem, params);
	end

	if lsiter == maxlsiter
		disp('Reached max line search,.... not so good... might have a bug in operators. exiting... ');
		return;
	end

	% control the number of line searches by adapting the initial step search
	if lsiter > 2
		t0 = t0 * beta;
	end 
	
	if lsiter<1
		t0 = t0 / beta;
	end

	x = (x + t*dx);

	%--------- uncomment for debug purposes ------------------------	
	disp(sprintf('%d   , obj: %f, RMS: %f, L-S: %d', k,f1,RMSerr,lsiter));

	%---------------------------------------------------------------
	
    %conjugate gradient calculation
    
	g1 = wGradient(x,problem, params);
	bk = g1(:)'*g1(:)/(g0(:)'*g0(:)+eps);
	g0 = g1;
	dx =  - g1 + bk* dx;
	k = k + 1;
	
	%TODO: need to "think" of a "better" stopping criteria ;-)
	if (k > params.Itnlim) || (norm(dx(:)) < gradToll) 
		break;
	end

end


return;


function [FTXFMtx, FTXFMtdx, DXFMtx, DXFMtdx] = preobjective(x, dx, problem)

% precalculates transforms to make line search cheap

FTXFMtx = problem.A*x;
FTXFMtdx = problem.A*dx;

if problem.TVWeight
    DXFMtx = problem.TV*(problem.A.psi'*x);
    DXFMtdx = problem.TV*(problem.A.psi'*dx);
else
    DXFMtx = 0;
    DXFMtdx = 0;
end





function [res, obj, RMS] = objective(FTXFMtx, FTXFMtdx, DXFMtx, DXFMtdx, x,dx,t, problem, params)
%calculated the objective function

p = params.pNorm;

obj = FTXFMtx + t*FTXFMtdx - problem.y;
obj = obj(:)'*obj(:);

if problem.TVWeight
    w = DXFMtx(:) + t*DXFMtdx(:);
    TV = (w.*conj(w)+params.l1Smooth).^(p/2); 
else
    TV = 0;
end

if problem.xfmWeight
   w = x(:) + t*dx(:); 
   XFM = (w.*conj(w)+params.l1Smooth).^(p/2);
else
    XFM=0;
end



TV = sum(TV.*problem.TVWeight(:));
XFM = sum(XFM.*problem.xfmWeight(:));
RMS = sqrt(obj/sum(abs(problem.y(:))>0));

%fprintf('------- Obj=%5.2f XFM=%5.2f TV=%5.2f -------\n',obj,XFM,TV);

res = obj + (TV) + (XFM) ;

function grad = wGradient(x,problem,params)

gradXFM = 0;
gradTV = 0;

gradObj = gOBJ(x,problem);
if problem.xfmWeight
gradXFM = gXFM(x,params);
end
if problem.TVWeight
gradTV = gTV(x,problem, params);
end


grad = (gradObj +  problem.xfmWeight.*gradXFM + problem.TVWeight.*gradTV);



function grad = gOBJ(x,problem)
% computes the gradient of the data consistency

grad = problem.A'*(problem.A*x - problem.y);

grad = 2*grad ;

function grad = gXFM(x,params)
% compute gradient of the L1 transform operator

p = params.pNorm;

grad = p*x.*(x.*conj(x)+params.l1Smooth).^(p/2-1);


function grad = gTV(x,problem,params)
% compute gradient of TV operator

p = params.pNorm;

Dx = problem.TV*(problem.A.psi'*x);
G = p*Dx.*(Dx.*conj(Dx) + params.l1Smooth).^(p/2-1);
grad = problem.A.psi*(problem.TV'*G);
grad=grad(:);






