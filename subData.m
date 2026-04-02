function [bold,prfs] = subData(fullBold,fullPrfs,id)
% bold = subBold(fullBold,id)
bold.vertex = fullBold.vertex(id);
bold.varea = fullBold.varea(id);
bold.data = fullBold.data(:,id);
bold.hemisphere = fullBold.data(id);


prfs.vertex = fullPrfs.vertex(id);
prfs.varea = fullPrfs.varea(id);
prfs.x0 = fullPrfs.x0(id);
prfs.y0 = fullPrfs.y0(id);
prfs.sigma = fullPrfs.sigma(id);
prfs.vexpl = fullPrfs.vexpl(id);
prfs.hemisphere = fullPrfs.hemisphere(id);

