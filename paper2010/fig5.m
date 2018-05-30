
clear

tmp=load('f5.mat');
dataset=imresize(tmp.data, [32 32], 'bicubic');
radius=0.1;

figure;
nrow=2;
ncol=4;
i=1;
for it=1:nrow
    pdf1=genPDF([32 32], 5, 1/6, 2, radius, 0);
    if it == 1
        mask=genSampling(pdf1, 10, 60);
    else
        tmp=load('32mask1D.mat');
        mask=tmp.mask;
    end
    
    subplot(nrow,ncol,i)
    imshow(mask);
    i=i+1;

    pulse=zeros(size(mask));
    pulse(8,8)=1;
    result{i}=tpsf(pulse, 'fft', 'svd', mask, dataset);
    subplot(nrow,ncol,i);
    mesh(abs(result{i}), 'FaceColor',[1 1 1],'EdgeColor',[0.7 0.7 0.7]);
    set(gca, 'XTick', []);
    set(gca, 'YTick', []);
    xlim([1 size(result{i}, 1)]);
    ylim([1 size(result{i}, 1)]);
%    zlim([0 0.25]);
%    axis square;
    i=i+1;
    
    pulse=zeros(size(mask));
    pulse(8,16)=1;
    result{i}=tpsf(pulse, 'fft', 'svd', mask, dataset);
    subplot(nrow,ncol,i);
    mesh(abs(result{i}), 'FaceColor',[1 1 1],'EdgeColor',[0.7 0.7 0.7]);
    set(gca, 'XTick', []);
    set(gca, 'YTick', []);
    xlim([1 size(result{i}, 1)]);
    ylim([1 size(result{i}, 1)]);
%    zlim([0 0.25]);
%    axis square;
    i=i+1;
    
    pulse=zeros(size(mask));
    pulse(16,16)=1;
    result{i}=tpsf(pulse, 'fft', 'svd', mask, dataset);
    subplot(nrow,ncol,i);
    mesh(abs(result{i}), 'FaceColor',[1 1 1],'EdgeColor',[0.7 0.7 0.7]);
    set(gca, 'XTick', []);
    set(gca, 'YTick', []);
    xlim([1 size(result{i}, 1)]);
    ylim([1 size(result{i}, 1)]);
%    zlim([0 0.25]);    
%    axis square;
    i=i+1;
end
