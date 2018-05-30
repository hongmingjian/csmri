function mask=genmask(imSize, pdf1, AF)

k=imSize(2)*AF;

indicators=zeros(1, imSize(2));
indicators(pdf1>=1)=1;

if k < sum(indicators)
    error('Bad parameter: pdf');
end

k=k-sum(indicators);

while k>0
    i=randi(imSize(2));
    if indicators(i) == 0
        xi=rand(1);
        if xi <= pdf1(i)
            indicators(i)=1;
            k=k-1;
        end
    end
end

mask=repmat(indicators, imSize(1), 1);