
load('f5.mat')
A=imresize(data, [256 256], 'bicubic');
[U, S, V]=svd(A);
figure;
plot(log2(1:size(A, 1)), diag(S), '-s')
ylabel('\sigma_i')
set(gca,'XTickLabel',{'1';'2';'4';'8';'16';'32';'64';'128';'256'})
xlabel('i')

figure;
sum=zeros(size(A));
for k=1:size(A, 1)
    sum=sum+S(k, k)*U(:, k)*V(:, k)';
    mse(k)=norm(A(:)-sum(:))/sqrt(numel(A));    
end
plot(log2(1:size(A, 1)), mse, '-s')
ylabel('MSE')
set(gca,'XTickLabel',{'1';'2';'4';'8';'16';'32';'64';'128';'256'})
xlabel('k')

figure;
pc=[1 40 80];
S=diag(S);
res=A;
for i=1:length(pc)
    k=pc(i);
    Sk=zeros(size(S));
    Sk(1:k)=S(1:k);
    res=cat(2, res, U*diag(Sk)*V');
end
imshow(res);
