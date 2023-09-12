
% Material and geometric parameters: ae, ce, fe, he (he is element size and all
% parameters are assumed constant in the element)
% p is the number of nodes in element e. Order of element is p-1 (i.e. linear for p=2, quadratic for p=3)


function [feM, keM] = EM_321(ae, ce, fe, he, p)

if (p==2)
    keM = ae/he*[1, -1; -1, 1] + ce*he/6*[2, 1; 1, 2];
    feM = fe*he/2*[1, 1]';
end

if (p==3)
    keM = ae/(3*he)*[7, -8, 1; -8, 16, -8; 1, -8, 7] + ce*he/30*[4, 2, -1; 2, 16, 2; -1, 2, 4];
    feM = fe*he/6*[1, 4, 1]';
end

end
