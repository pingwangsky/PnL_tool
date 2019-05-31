% size 12x1 or 3x4
function [BPN lm] = normBP(BP)
    
    flag = 0;
    if size(BP,1) == 12
        BP = reshape(BP,3,4);
        flag = 12;
    end
    
    % scale
    m = BP;
    D = [m(:,1)-m(:,4), m(:,2)-m(:,4), m(:,3)-m(:,4)]; 
    D2 = sum(D.^2);
    lm = sqrt(3/sum(D2));
    if BP(3,1) < 0
        lm = -lm;
    end
    BPN = BP * lm;    
    
    % vectorize
    if flag == 12
        BPN = BPN(:);
    end	
	
