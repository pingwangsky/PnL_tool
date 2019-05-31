% get alpha coordinates
function alp = getAlpha(BP,P)

    BP = [BP; 1 1 1 1];
    P = [P; ones(1,size(P,2))];
    alp = inv(BP) * P;	
