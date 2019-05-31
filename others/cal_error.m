function [rotation_error,position_error]=cal_error(R_cw,T_cw,R,T)

min_error=inf;
for i=1:size(R_cw,3)
    temp=norm(R_cw(:,:,i)-R);
    if temp<min_error
        min_error=temp;
        optR=R_cw(:,:,i);
        optT=T_cw(:,i);
    end
end
rotation_error=norm(optR-R,2);
position_error=norm(optT-T,2);
end