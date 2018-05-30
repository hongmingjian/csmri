function problem=csexp(dataset, phi, psi, algo, niter, mask, pdf1)

%% Formulate the problem
switch(upper(dataset))
    case 'WSHA'
        % Prepare the data
        K=7;      %  ???(?FFT?????)
        N=256;    %  ????
        M=64;     %  ???(M>=K*log(N/K),??40,???????)
        f1=50;    %  ????1
        f2=100;   %  ????2
        f3=200;   %  ????3
        f4=400;   %  ????4
        fs=800;   %  ????
        ts=1/fs;  %  ????
        Ts=1:N;   %  ????
        
        x=0.3*cos(2*pi*f1*Ts*ts)+...
            0.6*cos(2*pi*f2*Ts*ts)+...
            0.1*cos(2*pi*f3*Ts*ts)+...
            0.9*cos(2*pi*f4*Ts*ts);  %  ????
        
        x=x.'; % to a column vector
        
        % formulate the problem
        problem.M=M;
        problem.size=[1 N];
        
        S=randn(problem.M, prod(problem.size));
        Mx=@(z) S*z;
        Mxt=@(z) S'*z;
        
        phi=1;
        psi=FFT1D(prod(problem.size), prod(problem.size), 1, 0);
        
        problem.A=ComposeA(Mx, Mxt, phi, psi);
        problem.y = Mx(x);
        problem.x0= zeros(prod(problem.size), 1);       % Initial guess
        problem.xtrue=psi*(phi'*x);
        problem.TVWeight = 0;          % TV penalty
        problem.xfmWeight = 59;        % L1 penalty
        
    case 'LUSTIG'
        % Prepare the data
        tmp=load('/home/users/uqmhong1/work/cs/sparseMRI_v0.2/brain512.mat');
        mask=tmp.mask; pdf1=tmp.pdf; data=tmp.data;
        
        idx=find(mask==1);
        Mx=@(z) z(idx);
        Mxt=@(z) subsasgn(zeros(size(mask)), substruct('()', {idx}), z);
        
        if ischar(phi)
            switch(upper(phi))
                case {'FFT'}
                    phi=FFT2D(size(data), size(mask));
            end
        end
        
        tmp=data./pdf1;
        im_zf=phi'*tmp;
        data=data/max(abs(im_zf(:)));
        im_zf=im_zf/max(abs(im_zf(:)));
        
        if ischar(psi)
            switch(upper(psi))
                case {'DWT'}
                    psi=FWT2D(size(data), size(mask), 'Daubechies', 4, 4);
                case {'SVD'}
                    psi=SVD(im_zf);
            end
        end
        
        % formulate the problem
        problem.M=length(idx);
        problem.size=size(data);
        problem.A=ComposeA(Mx, Mxt, phi, psi);
        problem.y = problem.A.Mx(data);
        problem.x0 = psi*im_zf;  % Initial guess
        problem.x0 = problem.x0(:);
        problem.TV = TVOP(size(data));
        problem.TVWeight  = 0.002;   % TV penalty
        problem.xfmWeight = 0.005;   % L1 penalty
        
    case 'PHANTOM'
        full=phantom(512);
        full=full/max(full(:));
        
        if ischar(phi)
            switch(upper(phi))
                case {'FFT'}
                    phi=FFT2D(size(full), size(mask));
            end
        end
        
        full=phi*full;       % to k-spaced
        
        idx=find(mask==1);
        Mx=@(z) z(idx);
        Mxt=@(z) subsasgn(zeros(size(mask)), substruct('()', {idx}), z);
        
        data=Mxt(Mx(full));  % undersampling
        
        %        im_zf=phi'*(data./pdf1);
        %        data=data/max(abs(im_zf(:)));
        %        im_zf=abs(im_zf)/max(abs(im_zf(:)));
        im_zf=abs(phi'*data);
        im_zf=im_zf/max(im_zf(:));
        
        if ischar(psi)
            switch(upper(psi))
                case {'DWT', 'FWT'}
                    psi=FWT2D(size(data), size(mask), 'Daubechies', 4, 4);
                case 'SVD'
                    psi=SVD(im_zf);
            end
        end
        
        im_full=abs(phi'*full);
        im_full=im_full/max(im_full(:));
        
        % formulate the problem
        problem.M=length(idx);
        problem.size=size(data);
        problem.A=ComposeA(Mx, Mxt, phi, psi);
        problem.y=Mx(data);
        problem.x0=psi*im_zf;  % Initial guess
        problem.x0=problem.x0(:);
        problem.xtrue=psi*im_full;
        problem.xtrue=problem.xtrue(:);
        problem.TV=TVOP(size(data));
        problem.TVWeight=0.002;    % TV penalty
        problem.xfmWeight=0.005;   % L1 penalty
        
    otherwise
        tmp=load(dataset);
        full=tmp.data;
        full=full/max(full(:));
        
        if ischar(phi)
            switch(upper(phi))
                case {'FFT'}
                    phi=FFT2D(size(full), size(mask));
            end
        end
        
        full=phi*full;       % to k-spaced
        
        idx=find(mask==1);
        Mx=@(z) z(idx);
        Mxt=@(z) subsasgn(zeros(size(mask)), substruct('()', {idx}), z);
        
        data=Mxt(Mx(full));  % undersampling
        
        %        im_zf=phi'*(data./pdf1);
        %        data=data/max(abs(im_zf(:)));
        %        im_zf=abs(im_zf)/max(abs(im_zf(:)));
        im_zf=abs(phi'*data);
        im_zf=im_zf/max(im_zf(:));
        
        if ischar(psi)
            switch(upper(psi))
                case {'DWT', 'FWT'}
                    psi=FWT2D(size(data), size(mask), 'Daubechies', 4, 4);
                case 'SVD'
                    psi=SVD(im_zf);
            end
        end
        
        im_full=abs(phi'*full);
        im_full=im_full/max(im_full(:));
        
        % formulate the problem
        problem.M=length(idx);
        problem.size=size(data);
        problem.A=ComposeA(Mx, Mxt, phi, psi);
        problem.y=Mx(data);
        problem.x0=psi*im_zf;  % Initial guess
        problem.x0=problem.x0(:);
        problem.xtrue=psi*im_full;
        problem.xtrue=problem.xtrue(:);
        problem.TV=TVOP(size(data));
        problem.TVWeight=0.002;    % TV penalty
        problem.xfmWeight=0.005;   % L1 penalty
        
end

%% CS Optimization
tStart=tic;
problem=csopt(problem, algo, niter);
tElapsed=toc(tStart);

%% back to user domain
uhat=problem.A.psi'*problem.xhat;

%% Announce the result
switch(upper(dataset))
    case 'WSHA'
        figure(1);
        subplot(211);
        hold on;
        plot(x,'r'),plot(uhat,'k.-')
        legend('Original', 'Recovery');
        
        subplot(212)
        hold on
        plot(abs(problem.A.psi*x),'r'),plot(abs(problem.A.psi*uhat),'k.-')
        title('Sparse domain');
        
    otherwise
        u0=problem.A.psi'*problem.x0;
        u0=reshape(u0, problem.size);
        u0=abs(u0); u0=u0/max(u0(:));
        uhat=reshape(uhat, problem.size);
        uhat=abs(uhat); uhat=uhat/max(uhat(:));
        
        info=[];
        info=[info sprintf('%.2f%%', 100*length(find(mask==1))/numel(mask))];
        
        if isobject(problem.A.phi)
            info=[info ', \phi=' class(problem.A.phi)];
        end
        if isobject(problem.A.psi)
            info=[info ', \psi=' class(problem.A.psi)];
        end
        
        info=[info sprintf(', ||x||_1=%f, ||Ax-y||_2=%f', ...
            norm(problem.xhat, 1), ...
            norm(problem.A*problem.xhat-problem.y))];
        info=[info sprintf(', %s(%.4fs)', algo, tElapsed)];
        
        res=cat(2,u0,uhat);
        if isfield(problem, 'xtrue')
            utrue=problem.A.psi'*problem.xtrue;
            utrue=reshape(utrue, problem.size);
            utrue=abs(utrue); utrue=utrue/max(utrue(:));
            
            udiff=utrue-uhat;
            udiff=abs(udiff); udiff=udiff./max(udiff(:));
            
            res=cat(1, res, cat(2, udiff, utrue));
        end
        
        warning off Images:initSize:adjustingMag
        figure, imshow(res, 'InitialMagnification', 50);
        
        if isfield(problem, 'xtrue')
            err=problem.xtrue-problem.x0;
            decibels = 20*log10(sqrt(numel(err))/norm(err));
            text(0, 0, sprintf('Z/F(+%5.2f dB)', decibels), 'Color', 'w', 'VerticalAlignment', 'top');
            
            err=problem.xtrue-problem.xhat;
            decibels = 20*log10(sqrt(numel(err))/norm(err));
            info = [info sprintf(', %5.2f dB', decibels)];
            text(size(data, 2), 0, sprintf('CS(+%5.2f dB)', decibels), 'Color', 'w', 'VerticalAlignment', 'top');
            
            text(0, size(data, 1), '|Full-CS|', 'Color', 'w', 'VerticalAlignment', 'top');
            text(size(data, 2), size(data, 1), 'Full', 'Color', 'w', 'VerticalAlignment', 'top');
        else
            text(0, 0, 'Z/F', 'Color', 'w', 'VerticalAlignment', 'top');
            text(size(data, 2), 0, 'CS', 'Color', 'w', 'VerticalAlignment', 'top');
        end
        
        title(info);
        fprintf(1, '%s\n', info);
end

end
