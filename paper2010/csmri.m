function [problem,result] =csmri(dataset, lines, xform,border)
tmp=load(dataset);
eval(sprintf('mask=tmp.mask%d;', lines)); 
if nargin <4
    data_full= tmp.data;
else
    edata = conv2(tmp.data,border,'same');
    data_full = edata;
end
mask=mask'; %?????
clear tmp;

%%
phi=FFT2D(size(data_full),size(mask));
data=phi*data_full;%kspace
idx=find(mask==1);
Mx=@(z) z(idx);
Mxt=@(z) subsasgn(zeros(size(mask)), substruct('()', {idx}), z);
data=Mxt(Mx(data));%undersample
im_zf=abs(phi'*data);
im_zf=im_zf/max(im_zf(:));
im0 = im_zf;
tStart=tic;
for n=1:8
    switch(upper(xform))
    case 'SVD'
        psi=SVD(im0);
    case 'DWT'
        psi=FWT2D(size(data), size(mask), 'Daubechies', 4, 4);
    case 'DCT'
        psi=DCT2D(size(data));
    case 'IDT'
        psi=1;
    end
    
    problem.A=ComposeA(Mx,Mxt,phi,psi);
    problem.M=length(idx);
    problem.N=numel(mask);
    problem.size=size(data);
    problem.y = problem.A.Mx(data);
    problem.x0 = problem.A.psi*im0;  % Initial guess
    problem.x0 = problem.x0(:);
    problem.xtrue=problem.A.psi*data_full;
    problem.xtrue=problem.xtrue(:);
    problem.TV = TVOP(problem.size);
    
    problem.TVWeight  = 0.002;   % TV penalty
    problem.xfmWeight = 0.02;   % L1 penalty    
    
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
    
    xhat=fnlCg(xhat, problem, params);

    info=[];
    if isfield(problem, 'xtrue')
        err=problem.xtrue-xhat;
        decibels = 20*log10(sqrt(numel(err))/norm(err));
        info=sprintf('+%5.2f dB', decibels);
    end
    
    if strcmp(class(problem.A.psi), 'SVD')
        info=[info ', ' sprintf('||xhat||_1/||tr(xhat)||_1=%.4f/%.4f=%.4f', ...
            norm(xhat, 1), norm(diag(reshape(xhat, problem.size)), 1), ...
            norm(xhat, 1)/norm(diag(reshape(xhat, problem.size)), 1))];
    end
    %fprintf(1, '%s\n', info);
    im0=problem.A.psi'*xhat;
    
tElapsed=toc(tStart);
problem.time = tElapsed;
problem.xhat=xhat;
problem.img = problem.A.psi'* xhat;



%%
info=[];
info=[info sprintf('%.2f%% sampled', 100*length(find(mask==1))/numel(mask))];

if isobject(problem.A.phi)
    info=[info ', \phi=' class(problem.A.phi)];
end
if isobject(problem.A.psi)
    info=[info ', \psi=' class(problem.A.psi)];
end

info=[info sprintf(', ||x||_1=%f, ||Ax-y||_2=%f', ...
    norm(problem.xhat, 1), ...
    norm(problem.A*problem.xhat-problem.y))];
info=[info sprintf(', %.4f seconds', problem.time)];

if isfield(problem, 'xtrue')
    err=problem.xtrue-problem.xhat; 
    err = problem.A.psi'*err; err=err(:);
    decibels = 20*log10(sqrt(numel(err))/norm(err));
    info=[info sprintf(', +%5.2f dB', decibels)];
end

fprintf(1, '%s\n', info);
end



%%
uhat=problem.A.psi'*problem.xhat;
%uhat=uhat/max(abs(uhat(:)));
uhat=reshape(uhat, problem.size);

result.resim = uhat;
res=cat(2,im_zf,uhat);
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
title(info);

if isfield(problem, 'xtrue')
    im_zf= problem.A.psi*im_zf; im_zf = im_zf(:);
    err=problem.xtrue-im_zf;
    decibels = 20*log10(sqrt(numel(err))/norm(err));
    text(0, 0, sprintf('Z/F(+%5.2f dB)', decibels), 'Color', 'w', 'VerticalAlignment', 'top');
    
    err=problem.xtrue-problem.xhat;
    decibels = 20*log10(sqrt(numel(err))/norm(err));
    text(size(data, 2), 0, sprintf('CS(+%5.2f dB)', decibels), 'Color', 'w', 'VerticalAlignment', 'top');
    
    text(0, size(data, 1), '|Full-CS|', 'Color', 'w', 'VerticalAlignment', 'top');
    text(size(data, 2), size(data, 1), 'Full', 'Color', 'w', 'VerticalAlignment', 'top');
else
    text(0, 0, 'Z/F', 'Color', 'w', 'VerticalAlignment', 'top');
    text(size(data, 2), 0, 'CS', 'Color', 'w', 'VerticalAlignment', 'top');
end

end
