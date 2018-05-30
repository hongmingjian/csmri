%%
lines=[140 200 260];

%% Figure 6, 7, 8, 11, 12, 13
for i=1:length(lines)
    brain_dct{i}=csmri('brain.mat', lines(i), 'dct');
    brain_dwt{i}=csmri('brain.mat', lines(i), 'dwt');
    brain_svd{i}=csmri('brain.mat', lines(i), 'svd');
    angio_idt{i}=csmri('angio.mat', lines(i), 'idt');
    angio_dwt{i}=csmri('angio.mat', lines(i), 'dwt');
    angio_svd{i}=csmri('angio.mat', lines(i), 'svd');
end

%% Figure 10
figure; 
semilogx(abs(sort(abs(brain_dct{1}.xtrue), 1, 'descend')), 'b--','MarkerSize',5,'LineWidth',2)
hold on
semilogx(abs(sort(abs(brain_dwt{1}.xtrue), 1, 'descend')), 'g-.','MarkerSize',5,'LineWidth',2)
semilogx(abs(sort(abs(brain_svd{1}.xtrue), 1, 'descend')), 'r-','MarkerSize',5,'LineWidth',2)
legend('DCT', 'DWT', 'SVD');
xlabel('Index of pixels')
ylabel('Coefficients')

%% Figure 14
figure; 
semilogx(abs(sort(abs(angio_idt{1}.xtrue), 1, 'descend')), 'b--','MarkerSize',5,'LineWidth',2)
hold on
semilogx(abs(sort(abs(angio_dwt{1}.xtrue), 1, 'descend')), 'g-.','MarkerSize',5,'LineWidth',2)
semilogx(abs(sort(abs(angio_svd{1}.xtrue), 1, 'descend')), 'r-','MarkerSize',5,'LineWidth',2)
legend('IDT', 'DWT', 'SVD');
xlabel('Index of pixels')
ylabel('Coefficients')
