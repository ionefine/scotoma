function out = Gauss(pRF,x,y,dims)
%out = Gauss(pRF,x,y,dims)

if ~isfield(pRF,'ang')
    pRF.ang = 0;
end

rampx =  cos(pRF.ang)*(x-pRF.center(1)) + sin(pRF.ang)*(y-pRF.center(2));
rampy = -sin(pRF.ang)*(x-pRF.center(1)) + cos(pRF.ang)*(y-pRF.center(2));

out = exp(-(rampx.^2/(2*(pRF.ar*pRF.sig).^2)+rampy.^2/(2*pRF.sig.^2)));


if exist('dims','var')
    out= out(:);
end

