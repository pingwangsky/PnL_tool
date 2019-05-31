%******************************************************************************
% This programe is implemented via MATLAB 2013.                              *
% Author :  Ping Wang                                                        *
% Contact:  pingwangsky@gmail.com                                            *
% License:  Copyright (c) 2019 Ping Wang, All rights reserved.               *
% Address:  College of Electrical and Information Engineering,               *
%           Lanzhou University of Technology                                 *
% My site:  https://sites.google.com/view/ping-wang-homepage                 *
%*****************************************************************************/ 

function [ R_wc,t_wc ] = SRPnL( p1,p2,W1,W2 )

nLine = length(p1);
p1 = [p1; ones(1,nLine)];
p2 = [p2; ones(1,nLine)];

n_c=xnorm(cross(p1,p2));
%calculate the matrix Q in ahead;
Q=zeros(2*nLine,10);
N=zeros(2*nLine,3);
nx=n_c(1,:)'; ny=n_c(2,:)'; nz=n_c(3,:)';
Px=W1(1,:)'; Py=W1(2,:)'; Pz=W1(3,:)';
Q(1:2:end,:)=[Px.*nx+Py.*ny+Pz.*nz,...
          2*Py.*nz-2*Pz.*ny,...
          2*Pz.*nx-2*Px.*nz,...
          2*Px.*ny-2*Py.*nx,...
          Px.*nx-Py.*ny-Pz.*nz,...
          2.*Px.*ny+2*Py.*nx,...
          2.*Px.*nz+2*Pz.*nx,...
          Py.*ny-Px.*nx-Pz.*nz,...
          2*Py.*nz+2*Pz.*ny,...
          Pz.*nz-Px.*nx-Py.*ny];
Px=W2(1,:)'; Py=W2(2,:)'; Pz=W2(3,:)';
Q(2:2:end,:)=[Px.*nx+Py.*ny+Pz.*nz,...
          2*Py.*nz-2*Pz.*ny,...
          2*Pz.*nx-2*Px.*nz,...
          2*Px.*ny-2*Py.*nx,...
          Px.*nx-Py.*ny-Pz.*nz,...
          2.*Px.*ny+2*Py.*nx,...
          2.*Px.*nz+2*Pz.*nx,...
          Py.*ny-Px.*nx-Pz.*nz,...
          2*Py.*nz+2*Pz.*ny,...
          Pz.*nz-Px.*nx-Py.*ny];
N(1:2:end,:)=-n_c.';
N(2:2:end,:)=-n_c.';
%计算用于Gauss-Newton优化的参数;
CC=((N.'*N)\N.')*Q;
EE=Q-N*CC;
GG=EE.'*EE;

%计算直线方向和直线上任意一点;
Vw=W2-W1;
Vw = xnorm(Vw);
Pw = (W1 + W2) * 0.5;

xs=p1;
xe=p2;

n=nLine;
%choose the line with longest length in the image plane;
lineLenVec = sqrt(sum((xs-xe).^2));
[longestLineLen, LineID] = max(lineLenVec);

temp=n_c(:,1);   n_c(:,1)=n_c(:,LineID);  n_c(:,LineID)=temp;
temp = Vw(:,1);  Vw(:,1) = Vw(:, LineID); Vw(:,LineID) = temp;
temp = Pw(:,1);  Pw(:,1) = Pw(:, LineID); Pw(:,LineID) = temp;
lineLenVec(1) = 0;
l1 = xs(:,1) - xe(:,1); 
l1 = l1/norm(l1);

%first line is fixed. Find the second line
for i=2:n
    lineLenVec(LineID) = 0;
    [longestLineLen, LineID] = max(lineLenVec);%the current lonest line
    l2 = xs(:,LineID) - xe(:,LineID);
    l2 = l2/norm(l2);
    cosAngle = abs(l1'*l2);
    if cosAngle <  1.1 % 0<angle<180, 15<angle<165,or 30<angle<150
        break;
    end
end
% temp = xs(:,2);  xs(:,2) = xs(:, LineID); xs(:,LineID) = temp;
% temp = xe(:,2);  xe(:,2) = xe(:, LineID); xe(:,LineID) = temp;
temp=n_c(:,2);   n_c(:,2)=n_c(:,LineID);  n_c(:,LineID)=temp;
temp = Vw(:,2);  Vw(:,2) = Vw(:, LineID); Vw(:,LineID) = temp;
temp = Pw(:,2);  Pw(:,2) = Pw(:, LineID); Pw(:,LineID) = temp;
 
% Get the canonical configuration
% nc1 = cross(xs(:,1),xe(:,1));  nc1 = nc1/norm(nc1);
nc1=n_c(:,1);
Vw1 = Vw(:,1);    

    
Xm = cross(nc1,Vw1); Xm = Xm/norm(Xm); %the X axis of Model frame
Ym = nc1; %the Y axis of Model frame
Zm = cross(Xm,Ym);  Zm= Zm/norm(Zm);%the Z axis of Model frame;


Rot = [Xm, Ym, Zm]'; % Rot * [Xm, Ym, Zm] = I.

%rotate all the vector by Rot.    
% n_c = cross(xs,xe);
% n_c = xnorm(n_c);

nc_bar = Rot * n_c;
Vw_bar = Rot * Vw;
% Pw_bar = Rot * Pw;
    
% determine the angle psi, it is the angle between z axis and Vw_bar(:,1).
% the rotation matrix Rx(alpha) rotates Vw_bar(:,1) to z axis
cos_alpha= [0,0,1]*Vw_bar(:,1); %the angle between z axis and Vw_bar(:,1).

sin_alpha= sqrt(1 - cos_alpha*cos_alpha);
Rx= [1 0 0; 0 cos_alpha -sin_alpha; 0 sin_alpha cos_alpha];
Zaxis = Rx * Vw_bar(:,1); % should be the Z axis, i.e. [0, 0, 1]';
if 1 - abs(Zaxis(3)) > 1e-5
   Rx = Rx';
end
    
%fourth step, estimate the rotation angle phi by least square residual.
%i.e the rotation matrix Rz(phi)
Vm2 = Rx * Vw_bar(:,2);
A2= Vm2(1);      B2= Vm2(2);      C2= Vm2(3);
x2= nc_bar(1,2); y2= nc_bar(2,2); z2= nc_bar(3,2);

coef   = zeros(9,1); %coefficients of equation (7)
polyDF = zeros(16,1); %dF = ployDF(1) * t^15 + ployDF(2) * t^14 + ... + ployDF(15) * t + ployDF(16);
% construct the  polynomial F'
for i=3:n %The first two lines are included in every triplet.
    Vm3 = Rx*Vw_bar(:,i);
    A3= Vm3(1);      B3= Vm3(2);      C3= Vm3(3);
    x3= nc_bar(1,i); y3= nc_bar(2,i); z3= nc_bar(3,i);
    u11 = -z2*A2*y3*B3 + y2*B2*z3*A3;
    u12 = -y2*A2*z3*B3 + z2*B2*y3*A3;
    u13 = -y2*B2*z3*B3 + z2*B2*y3*B3 + y2*A2*z3*A3 - z2*A2*y3*A3;
    u14 = -y2*B2*x3*C3 + x2*C2*y3*B3;
    u15 =  x2*C2*y3*A3 - y2*A2*x3*C3;
    u21 = -x2*A2*y3*B3 + y2*B2*x3*A3;
    u22 = -y2*A2*x3*B3 + x2*B2*y3*A3;
    u23 =  x2*B2*y3*B3 - y2*B2*x3*B3 - x2*A2*y3*A3 + y2*A2*x3*A3;
    u24 =  y2*B2*z3*C3 - z2*C2*y3*B3;
    u25 =  y2*A2*z3*C3 - z2*C2*y3*A3;
    u31 = -x2*A2*z3*A3 + z2*A2*x3*A3;
    u32 = -x2*B2*z3*B3 + z2*B2*x3*B3;
    u33 =  x2*A2*z3*B3 - z2*A2*x3*B3 + x2*B2*z3*A3 - z2*B2*x3*A3;
    u34 =  z2*A2*z3*C3 + x2*A2*x3*C3 - z2*C2*z3*A3 - x2*C2*x3*A3;
    u35 = -z2*B2*z3*C3 - x2*B2*x3*C3 + z2*C2*z3*B3 + x2*C2*x3*B3;
    u36 = -x2*C2*z3*C3 + z2*C2*x3*C3;
        
    a4 =   u11*u11 + u12*u12 - u13*u13 - 2*u11*u12 +   u21*u21 + u22*u22 - u23*u23...
            -2*u21*u22 - u31*u31 - u32*u32 +   u33*u33 + 2*u31*u32;
    a3 =2*(u11*u14 - u13*u15 - u12*u14 +   u21*u24 -   u23*u25...
            - u22*u24 - u31*u34 + u33*u35 +   u32*u34);
    a2 =-2*u12*u12 + u13*u13 + u14*u14 -   u15*u15 + 2*u11*u12 - 2*u22*u22 + u23*u23...
            + u24*u24 - u25*u25 +2*u21*u22+ 2*u32*u32 -   u33*u33...
            - u34*u34 + u35*u35 -2*u31*u32- 2*u31*u36 + 2*u32*u36;
    a1 =2*(u12*u14 + u13*u15 +  u22*u24 +  u23*u25 -   u32*u34 - u33*u35 - u34*u36);
    a0 =   u12*u12 + u15*u15+   u22*u22 +  u25*u25 -   u32*u32 - u35*u35 - u36*u36 - 2*u32*u36;
    b3 =2*(u11*u13 - u12*u13 +  u21*u23 -  u22*u23 -   u31*u33 + u32*u33);
    b2 =2*(u11*u15 - u12*u15 +  u13*u14 +  u21*u25 -   u22*u25 + u23*u24 - u31*u35 + u32*u35 - u33*u34);
    b1 =2*(u12*u13 + u14*u15 +  u22*u23 +  u24*u25 -   u32*u33 - u34*u35 - u33*u36);
    b0 =2*(u12*u15 + u22*u25 -  u32*u35 -  u35*u36);
        
    d0 =    a0*a0 -   b0*b0;
    d1 = 2*(a0*a1 -   b0*b1);
    d2 =    a1*a1 + 2*a0*a2 +   b0*b0 - b1*b1 - 2*b0*b2;
    d3 = 2*(a0*a3 +   a1*a2 +   b0*b1 - b1*b2 -   b0*b3);
    d4 =    a2*a2 + 2*a0*a4 + 2*a1*a3 + b1*b1 + 2*b0*b2 - b2*b2 - 2*b1*b3;
    d5 = 2*(a1*a4 +   a2*a3 +   b1*b2 + b0*b3 -   b2*b3);
    d6 =    a3*a3 + 2*a2*a4 +   b2*b2 - b3*b3 + 2*b1*b3;
    d7 = 2*(a3*a4 +   b2*b3);
    d8 =    a4*a4 +   b3*b3;
        
    coef = coef + [a4, a3, a2, a1, a0, b3, b2, b1, b0]';
        
    polyDF(1) = polyDF(1) +                                8*d8*d8;
    polyDF(2) = polyDF(2) + 15* d7*d8;
    polyDF(3) = polyDF(3) + 14* d6*d8 +                    7*d7*d7;
    polyDF(4) = polyDF(4) + 13*(d5*d8 +  d6*d7);
    polyDF(5) = polyDF(5) + 12*(d4*d8 +  d5*d7)+           6*d6*d6;
    polyDF(6) = polyDF(6) + 11*(d3*d8 +  d4*d7 +  d5*d6);
    polyDF(7) = polyDF(7) + 10*(d2*d8 +  d3*d7 +  d4*d6) + 5*d5*d5;
    polyDF(8) = polyDF(8) + 9 *(d1*d8 +  d2*d7 +  d3*d6  +   d4*d5);
    polyDF(9) = polyDF(9) + 8 *(d1*d7 +  d2*d6 +  d3*d5) + 4*d4*d4 + 8*d0*d8;
    polyDF(10)= polyDF(10)+ 7 *(d1*d6 +  d2*d5 +  d3*d4) +           7*d0*d7;
    polyDF(11)= polyDF(11)+ 6 *(d1*d5 +  d2*d4)+           3*d3*d3 + 6*d0*d6;
    polyDF(12)= polyDF(12)+ 5 *(d1*d4 +  d2*d3)+                     5*d0*d5;
    polyDF(13)= polyDF(13)+ 4 * d1*d3 +                    2*d2*d2 + 4*d0*d4;
    polyDF(14)= polyDF(14)+ 3 * d1*d2 +                              3*d0*d3;
    polyDF(15)= polyDF(15)+                                  d1*d1 + 2*d0*d2;
    polyDF(16)= polyDF(16)+                                            d0*d1;
end

%solve polyDF
rs= roots(polyDF);
% retriving the local minima of the cost function.
maxreal= max(abs(real(rs)));
rs(abs(imag(rs))/maxreal > 0.01)= [];
minRoots = real(rs);

% poly    = (15:-1:1).*polyDF(1:15)';
% PolyVal = polyval(poly, minRoots);
% minRoots(PolyVal <= 0)= [];

count=1;
% minErr=inf;
numOfRoots = length(minRoots);
for rootId = 1 : numOfRoots
    cosbeta= minRoots(rootId);
    sign1 = sign(coef(1) * cosbeta^4 + coef(2) * cosbeta^3 + coef(3) * cosbeta^2 + coef(4) * cosbeta + coef(5));
    sign2 = sign(coef(6) * cosbeta^3 + coef(7) * cosbeta^2 + coef(8) * cosbeta   + coef(9));
    sinbeta= -sign1*sign2*sqrt(abs(1-cosbeta*cosbeta));
    Rz = [cosbeta -sinbeta 0; sinbeta cosbeta 0; 0 0 1];
    RzRxRot = Rz*Rx*Rot;
    
    nxi=nc_bar(1,:)'; nyi=nc_bar(2,:)'; nzi=nc_bar(3,:)';
    Vm=RzRxRot*Vw;
    Vxi=Vm(1,:)'; Vyi=Vm(2,:)'; Vzi=Vm(3,:)';
    A=[nxi .* Vxi + nzi .* Vzi,nxi .* Vzi - nzi .* Vxi,nyi .* Vyi];
    
    G=A.'*A;
    g11=G(1,1);
    g12=G(1,2);
    g13=G(1,3);
    g22=G(2,2);
    g23=G(2,3);
    F4=4*g12^2+g22^2+g11^2-2*g11*g22;
    F3=4*g12*g23+2*g11*g13-2*g13*g22;
    F2=g23^2+2*g11*g22+g13^2-4*g12^2-g11^2-g22^2;
    F1=2*g13*g22-2*g11*g13-2*g12*g23;
    F0=g12^2-g13^2;
    
    c= roots([F4,F3,F2,F1,F0]);
    maxreal1= max(abs(real(c)));
    c(abs(imag(c))/maxreal1 > 0.01)= [];
    
    c= real(c);
    c=c.';
    s=(2*g12*c.^2+g23*c-g12)./((g11-g22)*c+g13);
    
    for j=1:length(c)
        s1=c(j);
        s2=s(j);
        Ry = [s1, 0, s2; 0, 1, 0; -s2, 0, s1];
        R1 = (Rot') * (Ry * Rz * Rx) * Rot;
        %t1=Rot.'*V*[s1,s2,1].';
        
        solution = RefineGaussNewton(Cayley(R1),EE,GG);
        s1=solution(1);
        s2=solution(2);
        s3=solution(3);
        sr=[1,s1,s2,s3,s1^2,s1*s2,s1*s3,s2^2,s2*s3,s3^2].';
        factor=1/(1+s1^2+s2^2+s3^2); 
        Rr=factor*[1+s1^2-s2^2-s3^2,2*s1*s2-2*s3,2*s2+2*s1*s3;
                   2*s3+2*s1*s2,1-s1^2+s2^2-s3^2,2*s2*s3-2*s1;
                   2*s1*s3-2*s2,2*s1+2*s2*s3,1-s1^2-s2^2+s3^2];
        Rt=factor*CC*sr;
%         curErr=sr.'*GG*sr;
%         if curErr<minErr
%             optR=Rr;
%             optT=Rt;
%             minErr=curErr;
%         end        
        R_wc(:,:,count)=Rr;
        t_wc(:,count)=Rt;
        
        count=count+1;
    end
    
end
% R_wc=optR;
% t_wc=optT;

end

