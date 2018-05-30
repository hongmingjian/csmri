function problem=csopt(problem, algo, niter)

A =@(z) problem.A*z;
At=@(z) problem.A'*z;

switch upper(algo)
    case 'GPSR.BB'
        xhat=problem.x0;
        for i=1:niter
            [xhat,x_debias,objective,times,debias_start,mses,taus]= ...
                GPSR_BB(problem.y, A, problem.xfmWeight,...
                'AT', At, 'INITIALIZATION', xhat, 'Verbose', 0);
        end
    
    case 'L1_LS' % works by changing 'MAX_NT_ITER' from the default(400) to 2.
        xhat = problem.x0;
        for i=1:1%niter
            [xhat,status,history] = l1_ls(problem.A,problem.A',problem.M,problem.N,...
                problem.y, problem.xfmWeight, 1e-3, false, 1e-3, 20);%, xhat);
        end
        
    case 'NLCG'
        % Parameters being fed into the Non-Linear CG
        params.Itnlim = 8;
        params.gradToll = 1e-30;
        params.l1Smooth = 1e-15;
        params.pNorm = 1;

        % and line search parameters
        params.lineSearchItnlim = 150;
        params.lineSearchAlpha = 0.01;
        params.lineSearchBeta = 0.6;
        params.lineSearchT0 = 1 ; 

        % Kick it off
        xhat = problem.x0;
        for n=1:niter
%            if mod(n, 2)==0,
%            tmp=abs(problem.A.psi'*xhat);
%            tmp=tmp/max(tmp(:));
%            problem.A.psi=SVD(tmp);  
%            xhat=problem.x0;
%            end

            xhat=fnlCg(xhat, problem, params);
            
%            info=sprintf('||xhat||_1/||tr(xhat)||_1=%.4f/%.4f=%.4f', norm(xhat, 1), norm(diag(reshape(xhat, 512, 512)), 1), norm(xhat, 1)/norm(diag(reshape(xhat, 512, 512)), 1));
%            if isfield(problem, 'xtrue')        
%                err=problem.xtrue-xhat;                      
%                decibels = 20*log10(sqrt(numel(err))/norm(err));
%                info=[info sprintf(', +%5.2f dB', decibels)];
%            end           
%            info=[info '\n'];            
%            fprintf(1, info);
        end
        
   case 'L1QC_LOGBARRIER'
        xhat = problem.x0;
        for i=1:niter
            xhat = l1qc_logbarrier(xhat, A, At, problem.y, 0.001);
        end
        
    case 'SPARSELAB.STOMP'
        xhat=SolveStOMP(@A_for_sparselab, problem.y, problem.N); 
        
    case 'SPARSELAB.LASSO'
        xhat=SolveLasso(@A_for_sparselab, problem.y, problem.N, 'lasso', 20); 
        
    case 'SPARSELAB.OMP'
        xhat=SolveOMP(@A_for_sparselab, problem.y, problem.N, 20);
        
    case 'SPARSELAB.BP'
        xhat=SolveBP(@A_for_sparselab, problem.y, problem.N, 20, problem.xfmWeight);
        
    otherwise
        error('Bad algorithm: %s.', algo);
end

problem.xhat=xhat;

    % operator for callback of SparseLab
    %     y = OperatorName(mode, m, n, x, I, dim)
    %   This function gets as input a vector x and an index set I, and returns
    %   y = A(:,I)*x if mode = 1, or y = A(:,I)'*x if mode = 2. 
    %   A is the m by dim implicit matrix implemented by the function. I is a
    %   subset of the columns of A, i.e. a subset of 1:dim of length n. x is a
    %   vector of length n if mode = 1, or a vector of length m if mode = 2.
    function y=A_for_sparselab(mode, m, n, x, I, dim)
    switch mode
        case 1
            xx=zeros(dim, 1);
            xx(I)=x;
            y=A(xx);
        case 2
            yy=At(x);
            y=yy(I);
    end
    end
end
