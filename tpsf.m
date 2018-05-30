function result=tpsf(delta, phi, psi, mask, dataset)

if ischar(phi)
    switch(upper(phi))
        case {'FFT'}
            phi=FFT2D(size(delta), size(delta));
    end
end

idx=find(mask==1);
Mx=@(z) z(idx);
Mxt=@(z) subsasgn(zeros(size(mask)), substruct('()', {idx}), z);

if ischar(psi)
    switch(upper(psi))
        case {'DWT'}
            psi=FWT2D(size(delta), size(delta), 'Daubechies', 4, 4);
        case {'SVD'}
            if ~exist('dataset', 'var')
                error('bad args');
            end
            if ischar(dataset)
                switch(upper(dataset))
                    case 'PHANTOM'
                        full=phantom('Modified Shepp-Logan', size(delta, 1));
                    otherwise
                        tmp=load(dataset);
                        full=tmp.data;
                end
            else
                full=dataset;
            end
            full=full/max(full(:));
            full=phi*full;
            data=Mxt(Mx(full));
            im_zf=abs(phi'*data);
            im_zf=im_zf/max(im_zf(:));
            
            psi=SVD(im_zf);
    end
end

result=psi*(phi'*(Mxt(Mx(phi*(psi'*delta)))));

end
