[U, S, V]=svd(im_zf);

res=cat(2, im_zf, 3*abs(U), 3*abs(S));

S1=U'*im_true*V;

res=cat(1, res, cat(2, im_true, 3*abs(V), 3*abs(S1)));

imshow(res);

