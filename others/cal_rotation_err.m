function y = cal_rotation_err(R_estm, R_true)

X1= R_estm(:,1); X2= R_true(:,1);
Y1= R_estm(:,2); Y2= R_true(:,2);
Z1= R_estm(:,3); Z2= R_true(:,3);

X1= X1/norm(X1);
Y1= Y1/norm(Y1);
Z1= Z1/norm(Z1);
X2= X2/norm(X2);
Y2= Y2/norm(Y2);
Z2= Z2/norm(Z2);

exyz= [X1'*X2 Y1'*Y2 Z1'*Z2];
exyz(exyz>1)= 1;
exyz(exyz<-1)= -1;

y= max(abs(acos(exyz)))*180/pi;

return
