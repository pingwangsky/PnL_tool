function [R_cw, T_cw] = Ansar(p1, p2, P1_w, P2_w)

	[V_w P_w] = getVP(P1_w, P2_w);
	
    [R_wc, T_wc] = NLL([p1; p2].', V_w.', P_w.');
    
    R_cw = R_wc.';
    T_cw = -R_cw * T_wc;
    
    

function [Rot_cw, Pos_cw] = NLL(lineEndPoints, V, P)
%This function follows the n line linear(NLL) algorithm proposed by Ansar
%and Daniilidis: Linear pose estimation from points or lines.
%input: lineEndPoints(i,:) = [startpointx, startpointy, endpointx endpointy];
%       V(i,:) = the direction of ith line in the world frame
%       P(i,:) = a point of ith line in the world frame
%output: Rot_cw = Rot_wc' which is parametrized as
%              [r1,r4,r7; r2,r5,r8; r3,r6,r9];  (V_w = Rot_cw V_c)
%        Pos_cw = the position of camera in world frame;
%              P_w = R_cw * P_c + Pos_cw;

%first, compute the image line direction alpha and orthogonal vector gama
n = size(lineEndPoints,1);
if(n~=size(V,1) || n~=size(P,1)) 
    error('Input data lineEndPoints, V and P are inconsistent'); 
end
if(n<4)
    error('The input data is too few to estimate the camera pose');
end
    
alpha = zeros(n,3);%the image line direction
gama  = zeros(n,3);%the vector which is orthogonal to alpha
d     = zeros(n,3);%the non-normalized vector  which is orthogonal to alpha
for i=1:n
 temp = [lineEndPoints(i,3) - lineEndPoints(i,1), lineEndPoints(i,4) - lineEndPoints(i,2), 0]; %[ex-sx, ey-sy, 0]
 alpha(i,:) = temp/norm(temp);
 c_i = [lineEndPoints(i,1), lineEndPoints(i,2), 1 ]; %(sx,sy,1);
 d_i = c_i - (  c_i * alpha(i,:)' ) * alpha(i,:);
 gama(i,:) = d_i/norm(d_i);
 d(i,:)    = d_i;
end

%second, build the system matrix M, we have M * xbar = 0 where 
% xbar = (r1r1, r1r2, r1r3, r1r4, ... r1r9, r2r2, r2r3,...r2r9, r3r3,r3r4,...r3r9, ... , r8r8, r8r9, r9r9)' ;
M = zeros(n*(2*n-1)+12, 46);
k = 1;
for i = 1:n
    Vi      = V(i,:);
    alpha_i = alpha(i,:);
    gama_i  = gama(i,:);
    ag11i = alpha_i(1) * alpha_i(1) + gama_i(1) * gama_i(1);
    ag12i = alpha_i(1) * alpha_i(2) + gama_i(1) * gama_i(2);     ag21i = ag12i;
    ag13i = alpha_i(1) * alpha_i(3) + gama_i(1) * gama_i(3);     ag31i = ag13i;
    ag22i = alpha_i(2) * alpha_i(2) + gama_i(2) * gama_i(2);    
    ag23i = alpha_i(2) * alpha_i(3) + gama_i(2) * gama_i(3);     ag32i = ag23i;
    ag33i = alpha_i(3) * alpha_i(3) + gama_i(3) * gama_i(3);    
    for j = i:n
        Vj      = V(j,:);
        alpha_j = alpha(j,:);
        gama_j  = gama(j,:);
        ag11j = alpha_j(1) * alpha_j(1) + gama_j(1) * gama_j(1);
        ag12j = alpha_j(1) * alpha_j(2) + gama_j(1) * gama_j(2);     ag21j = ag12j;
        ag13j = alpha_j(1) * alpha_j(3) + gama_j(1) * gama_j(3);     ag31j = ag13j;
        ag22j = alpha_j(2) * alpha_j(2) + gama_j(2) * gama_j(2);
        ag23j = alpha_j(2) * alpha_j(3) + gama_j(2) * gama_j(3);     ag32j = ag23j;
        ag33j = alpha_j(3) * alpha_j(3) + gama_j(3) * gama_j(3);
        %according to equation (9)
        M(k,1)  = Vi(1)*Vj(1)*ag11i*ag11j + Vi(1)*Vj(1)*ag12i*ag21j + Vi(1)*Vj(1)*ag13i*ag31j; %r1r1
        M(k,2)  = Vi(1)*Vj(1)*ag11i*ag12j + Vi(1)*Vj(1)*ag21i*ag11j + Vi(1)*Vj(1)*ag12i*ag22j + Vi(1)*Vj(1)*ag22i*ag21j + Vi(1)*Vj(1)*ag13i*ag32j + Vi(1)*Vj(1)*ag23i*ag31j; %r1r2
        M(k,3)  = Vi(1)*Vj(1)*ag11i*ag13j + Vi(1)*Vj(1)*ag31i*ag11j + Vi(1)*Vj(1)*ag12i*ag23j + Vi(1)*Vj(1)*ag32i*ag21j + Vi(1)*Vj(1)*ag13i*ag33j + Vi(1)*Vj(1)*ag33i*ag31j; %r1r3
        M(k,4)  = Vi(1)*Vj(2)*ag11i*ag11j + Vi(2)*Vj(1)*ag11i*ag11j + Vi(1)*Vj(2)*ag12i*ag21j + Vi(2)*Vj(1)*ag12i*ag21j + Vi(1)*Vj(2)*ag13i*ag31j + Vi(2)*Vj(1)*ag13i*ag31j; %r1r4
        M(k,5)  = Vi(1)*Vj(2)*ag11i*ag12j + Vi(2)*Vj(1)*ag21i*ag11j + Vi(1)*Vj(2)*ag12i*ag22j + Vi(2)*Vj(1)*ag22i*ag21j + Vi(1)*Vj(2)*ag13i*ag32j + Vi(2)*Vj(1)*ag23i*ag31j; %r1r5
        M(k,6)  = Vi(1)*Vj(2)*ag11i*ag13j + Vi(2)*Vj(1)*ag31i*ag11j + Vi(1)*Vj(2)*ag12i*ag23j + Vi(2)*Vj(1)*ag32i*ag21j + Vi(1)*Vj(2)*ag13i*ag33j + Vi(2)*Vj(1)*ag33i*ag31j; %r1r6
        M(k,7)  = Vi(1)*Vj(3)*ag11i*ag11j + Vi(3)*Vj(1)*ag11i*ag11j + Vi(1)*Vj(3)*ag12i*ag21j + Vi(3)*Vj(1)*ag12i*ag21j + Vi(1)*Vj(3)*ag13i*ag31j + Vi(3)*Vj(1)*ag13i*ag31j; %r1r7
        M(k,8)  = Vi(1)*Vj(3)*ag11i*ag12j + Vi(3)*Vj(1)*ag21i*ag11j + Vi(1)*Vj(3)*ag12i*ag22j + Vi(3)*Vj(1)*ag22i*ag21j + Vi(1)*Vj(3)*ag13i*ag32j + Vi(3)*Vj(1)*ag23i*ag31j; %r1r8
        M(k,9)  = Vi(1)*Vj(3)*ag11i*ag13j + Vi(3)*Vj(1)*ag31i*ag11j + Vi(1)*Vj(3)*ag12i*ag23j + Vi(3)*Vj(1)*ag32i*ag21j + Vi(1)*Vj(3)*ag13i*ag33j + Vi(3)*Vj(1)*ag33i*ag31j; %r1r9
        M(k,10) = Vi(1)*Vj(1)*ag21i*ag12j + Vi(1)*Vj(1)*ag22i*ag22j + Vi(1)*Vj(1)*ag23i*ag32j; %r2r2
        M(k,11) = Vi(1)*Vj(1)*ag21i*ag13j + Vi(1)*Vj(1)*ag31i*ag12j + Vi(1)*Vj(1)*ag22i*ag23j + Vi(1)*Vj(1)*ag32i*ag22j + Vi(1)*Vj(1)*ag23i*ag33j + Vi(1)*Vj(1)*ag33i*ag32j; %r2r3
        M(k,12) = Vi(1)*Vj(2)*ag21i*ag11j + Vi(2)*Vj(1)*ag11i*ag12j + Vi(1)*Vj(2)*ag22i*ag21j + Vi(2)*Vj(1)*ag12i*ag22j + Vi(1)*Vj(2)*ag23i*ag31j + Vi(2)*Vj(1)*ag13i*ag32j; %r2r4
        M(k,13) = Vi(1)*Vj(2)*ag21i*ag12j + Vi(2)*Vj(1)*ag21i*ag12j + Vi(1)*Vj(2)*ag22i*ag22j + Vi(2)*Vj(1)*ag22i*ag22j + Vi(1)*Vj(2)*ag23i*ag32j + Vi(2)*Vj(1)*ag23i*ag32j; %r2r5
        M(k,14) = Vi(1)*Vj(2)*ag21i*ag13j + Vi(2)*Vj(1)*ag31i*ag12j + Vi(1)*Vj(2)*ag22i*ag23j + Vi(2)*Vj(1)*ag32i*ag22j + Vi(1)*Vj(2)*ag23i*ag33j + Vi(2)*Vj(1)*ag33i*ag32j; %r2r6
        M(k,15) = Vi(1)*Vj(3)*ag21i*ag11j + Vi(3)*Vj(1)*ag11i*ag12j + Vi(1)*Vj(3)*ag22i*ag21j + Vi(3)*Vj(1)*ag12i*ag22j + Vi(1)*Vj(3)*ag23i*ag31j + Vi(3)*Vj(1)*ag13i*ag32j; %r2r7
        M(k,16) = Vi(1)*Vj(3)*ag21i*ag12j + Vi(3)*Vj(1)*ag21i*ag12j + Vi(1)*Vj(3)*ag22i*ag22j + Vi(3)*Vj(1)*ag22i*ag22j + Vi(1)*Vj(3)*ag23i*ag32j + Vi(3)*Vj(1)*ag23i*ag32j; %r2r8
        M(k,17) = Vi(1)*Vj(3)*ag21i*ag13j + Vi(3)*Vj(1)*ag31i*ag12j + Vi(1)*Vj(3)*ag22i*ag23j + Vi(3)*Vj(1)*ag32i*ag22j + Vi(1)*Vj(3)*ag23i*ag33j + Vi(3)*Vj(1)*ag33i*ag32j; %r2r9
        M(k,18) = Vi(1)*Vj(1)*ag31i*ag13j + Vi(1)*Vj(1)*ag32i*ag23j + Vi(1)*Vj(1)*ag33i*ag33j; %r3r3
        M(k,19) = Vi(1)*Vj(2)*ag31i*ag11j + Vi(2)*Vj(1)*ag11i*ag13j + Vi(1)*Vj(2)*ag32i*ag21j + Vi(2)*Vj(1)*ag12i*ag23j + Vi(1)*Vj(2)*ag33i*ag31j + Vi(2)*Vj(1)*ag13i*ag33j; %r3r4
        M(k,20) = Vi(1)*Vj(2)*ag31i*ag12j + Vi(2)*Vj(1)*ag21i*ag13j + Vi(1)*Vj(2)*ag32i*ag22j + Vi(2)*Vj(1)*ag22i*ag23j + Vi(1)*Vj(2)*ag33i*ag32j + Vi(2)*Vj(1)*ag23i*ag33j; %r3r5
        M(k,21) = Vi(1)*Vj(2)*ag31i*ag13j + Vi(2)*Vj(1)*ag31i*ag13j + Vi(1)*Vj(2)*ag32i*ag23j + Vi(2)*Vj(1)*ag32i*ag23j + Vi(1)*Vj(2)*ag33i*ag33j + Vi(2)*Vj(1)*ag33i*ag33j; %r3r6
        M(k,22) = Vi(1)*Vj(3)*ag31i*ag11j + Vi(3)*Vj(1)*ag11i*ag13j + Vi(1)*Vj(3)*ag32i*ag21j + Vi(3)*Vj(1)*ag12i*ag23j + Vi(1)*Vj(3)*ag33i*ag31j + Vi(3)*Vj(1)*ag13i*ag33j; %r3r7
        M(k,23) = Vi(1)*Vj(3)*ag31i*ag12j + Vi(3)*Vj(1)*ag21i*ag13j + Vi(1)*Vj(3)*ag32i*ag22j + Vi(3)*Vj(1)*ag22i*ag23j + Vi(1)*Vj(3)*ag33i*ag32j + Vi(3)*Vj(1)*ag23i*ag33j; %r3r8
        M(k,24) = Vi(1)*Vj(3)*ag31i*ag13j + Vi(3)*Vj(1)*ag31i*ag13j + Vi(1)*Vj(3)*ag32i*ag23j + Vi(3)*Vj(1)*ag32i*ag23j + Vi(1)*Vj(3)*ag33i*ag33j + Vi(3)*Vj(1)*ag33i*ag33j; %r3r9
        M(k,25) = Vi(2)*Vj(2)*ag11i*ag11j + Vi(2)*Vj(2)*ag12i*ag21j + Vi(2)*Vj(2)*ag13i*ag31j; %r4r4
        M(k,26) = Vi(2)*Vj(2)*ag11i*ag12j + Vi(2)*Vj(2)*ag21i*ag11j + Vi(2)*Vj(2)*ag12i*ag22j + Vi(2)*Vj(2)*ag22i*ag21j + Vi(2)*Vj(2)*ag13i*ag32j + Vi(2)*Vj(2)*ag23i*ag31j; %r4r5
        M(k,27) = Vi(2)*Vj(2)*ag11i*ag13j + Vi(2)*Vj(2)*ag31i*ag11j + Vi(2)*Vj(2)*ag12i*ag23j + Vi(2)*Vj(2)*ag32i*ag21j + Vi(2)*Vj(2)*ag13i*ag33j + Vi(2)*Vj(2)*ag33i*ag31j; %r4r6
        M(k,28) = Vi(2)*Vj(3)*ag11i*ag11j + Vi(3)*Vj(2)*ag11i*ag11j + Vi(2)*Vj(3)*ag12i*ag21j + Vi(3)*Vj(2)*ag12i*ag21j + Vi(2)*Vj(3)*ag13i*ag31j + Vi(3)*Vj(2)*ag13i*ag31j; %r4r7
        M(k,29) = Vi(2)*Vj(3)*ag11i*ag12j + Vi(3)*Vj(2)*ag21i*ag11j + Vi(2)*Vj(3)*ag12i*ag22j + Vi(3)*Vj(2)*ag22i*ag21j + Vi(2)*Vj(3)*ag13i*ag32j + Vi(3)*Vj(2)*ag23i*ag31j; %r4r8
        M(k,30) = Vi(2)*Vj(3)*ag11i*ag13j + Vi(3)*Vj(2)*ag31i*ag11j + Vi(2)*Vj(3)*ag12i*ag23j + Vi(3)*Vj(2)*ag32i*ag21j + Vi(2)*Vj(3)*ag13i*ag33j + Vi(3)*Vj(2)*ag33i*ag31j; %r4r9
        M(k,31) = Vi(2)*Vj(2)*ag21i*ag12j + Vi(2)*Vj(2)*ag22i*ag22j + Vi(2)*Vj(2)*ag23i*ag32j; %r5r5
        M(k,32) = Vi(2)*Vj(2)*ag21i*ag13j + Vi(2)*Vj(2)*ag31i*ag12j + Vi(2)*Vj(2)*ag22i*ag23j + Vi(2)*Vj(2)*ag32i*ag22j + Vi(2)*Vj(2)*ag23i*ag33j + Vi(2)*Vj(2)*ag33i*ag32j; %r5r6
        M(k,33) = Vi(2)*Vj(3)*ag21i*ag11j + Vi(3)*Vj(2)*ag11i*ag12j + Vi(2)*Vj(3)*ag22i*ag21j + Vi(3)*Vj(2)*ag12i*ag22j + Vi(2)*Vj(3)*ag23i*ag31j + Vi(3)*Vj(2)*ag13i*ag32j; %r5r7
        M(k,34) = Vi(2)*Vj(3)*ag21i*ag12j + Vi(3)*Vj(2)*ag21i*ag12j + Vi(2)*Vj(3)*ag22i*ag22j + Vi(3)*Vj(2)*ag22i*ag22j + Vi(2)*Vj(3)*ag23i*ag32j + Vi(3)*Vj(2)*ag23i*ag32j; %r5r8
        M(k,35) = Vi(2)*Vj(3)*ag21i*ag13j + Vi(3)*Vj(2)*ag31i*ag12j + Vi(2)*Vj(3)*ag22i*ag23j + Vi(3)*Vj(2)*ag32i*ag22j + Vi(2)*Vj(3)*ag23i*ag33j + Vi(3)*Vj(2)*ag33i*ag32j; %r5r9     
        M(k,36) = Vi(2)*Vj(2)*ag31i*ag13j + Vi(2)*Vj(2)*ag32i*ag23j + Vi(2)*Vj(2)*ag33i*ag33j; %r6r6
        M(k,37) = Vi(2)*Vj(3)*ag31i*ag11j + Vi(3)*Vj(2)*ag11i*ag13j + Vi(2)*Vj(3)*ag32i*ag21j + Vi(3)*Vj(2)*ag12i*ag23j + Vi(2)*Vj(3)*ag33i*ag31j + Vi(3)*Vj(2)*ag13i*ag33j; %r6r7
        M(k,38) = Vi(2)*Vj(3)*ag31i*ag12j + Vi(3)*Vj(2)*ag21i*ag13j + Vi(2)*Vj(3)*ag32i*ag22j + Vi(3)*Vj(2)*ag22i*ag23j + Vi(2)*Vj(3)*ag33i*ag32j + Vi(3)*Vj(2)*ag23i*ag33j; %r6r8
        M(k,39) = Vi(2)*Vj(3)*ag31i*ag13j + Vi(3)*Vj(2)*ag31i*ag13j + Vi(2)*Vj(3)*ag32i*ag23j + Vi(3)*Vj(2)*ag32i*ag23j + Vi(2)*Vj(3)*ag33i*ag33j + Vi(3)*Vj(2)*ag33i*ag33j; %r6r9  
        M(k,40) = Vi(3)*Vj(3)*ag11i*ag11j + Vi(3)*Vj(3)*ag12i*ag21j + Vi(3)*Vj(3)*ag13i*ag31j; %r7r7
        M(k,41) = Vi(3)*Vj(3)*ag11i*ag12j + Vi(3)*Vj(3)*ag21i*ag11j + Vi(3)*Vj(3)*ag12i*ag22j + Vi(3)*Vj(3)*ag22i*ag21j + Vi(3)*Vj(3)*ag13i*ag32j + Vi(3)*Vj(3)*ag23i*ag31j; %r7r8
        M(k,42) = Vi(3)*Vj(3)*ag11i*ag13j + Vi(3)*Vj(3)*ag31i*ag11j + Vi(3)*Vj(3)*ag12i*ag23j + Vi(3)*Vj(3)*ag32i*ag21j + Vi(3)*Vj(3)*ag13i*ag33j + Vi(3)*Vj(3)*ag33i*ag31j; %r7r9  
        M(k,43) = Vi(3)*Vj(3)*ag21i*ag12j + Vi(3)*Vj(3)*ag22i*ag22j + Vi(3)*Vj(3)*ag23i*ag32j; %r8r8
        M(k,44) = Vi(3)*Vj(3)*ag21i*ag13j + Vi(3)*Vj(3)*ag31i*ag12j + Vi(3)*Vj(3)*ag22i*ag23j + Vi(3)*Vj(3)*ag32i*ag22j + Vi(3)*Vj(3)*ag23i*ag33j + Vi(3)*Vj(3)*ag33i*ag32j; %r8r9  
        M(k,45) = Vi(3)*Vj(3)*ag31i*ag13j + Vi(3)*Vj(3)*ag32i*ag23j + Vi(3)*Vj(3)*ag33i*ag33j; %r9r9
        M(k,46) = -Vi * Vj';%rho
        k = k+1;
        %according to equation (10)
        if j>i
            %first element of Rvi x Rvj
            M(k,1)  = Vi(1)*Vj(1)*ag12i*ag13j - Vi(1)*Vj(1)*ag13i*ag12j;%r1r1
            M(k,2)  = Vi(1)*Vj(1)*ag12i*ag23j - Vi(1)*Vj(1)*ag23i*ag12j + Vi(1)*Vj(1)*ag22i*ag13j - Vi(1)*Vj(1)*ag13i*ag22j;%r1r2
            M(k,3)  = Vi(1)*Vj(1)*ag12i*ag33j - Vi(1)*Vj(1)*ag33i*ag12j + Vi(1)*Vj(1)*ag32i*ag13j - Vi(1)*Vj(1)*ag13i*ag32j;%r1r3
            M(k,4)  = Vi(1)*Vj(2)*ag12i*ag13j - Vi(2)*Vj(1)*ag13i*ag12j + Vi(2)*Vj(1)*ag12i*ag13j - Vi(1)*Vj(2)*ag13i*ag12j;%r1r4
            M(k,5)  = Vi(1)*Vj(2)*ag12i*ag23j - Vi(2)*Vj(1)*ag23i*ag12j + Vi(2)*Vj(1)*ag22i*ag13j - Vi(1)*Vj(2)*ag13i*ag22j;%r1r5
            M(k,6)  = Vi(1)*Vj(2)*ag12i*ag33j - Vi(2)*Vj(1)*ag33i*ag12j + Vi(2)*Vj(1)*ag32i*ag13j - Vi(1)*Vj(2)*ag13i*ag32j;%r1r6
            M(k,7)  = Vi(1)*Vj(3)*ag12i*ag13j - Vi(3)*Vj(1)*ag13i*ag12j + Vi(3)*Vj(1)*ag12i*ag13j - Vi(1)*Vj(3)*ag13i*ag12j;%r1r7
            M(k,8)  = Vi(1)*Vj(3)*ag12i*ag23j - Vi(3)*Vj(1)*ag23i*ag12j + Vi(3)*Vj(1)*ag22i*ag13j - Vi(1)*Vj(3)*ag13i*ag22j;%r1r8
            M(k,9)  = Vi(1)*Vj(3)*ag12i*ag33j - Vi(3)*Vj(1)*ag33i*ag12j + Vi(3)*Vj(1)*ag32i*ag13j - Vi(1)*Vj(3)*ag13i*ag32j;%r1r9
            M(k,10) = Vi(1)*Vj(1)*ag22i*ag23j - Vi(1)*Vj(1)*ag23i*ag22j;%r2r2
            M(k,11) = Vi(1)*Vj(1)*ag22i*ag33j - Vi(1)*Vj(1)*ag33i*ag22j + Vi(1)*Vj(1)*ag32i*ag23j - Vi(1)*Vj(1)*ag23i*ag32j;%r2r3
            M(k,12) = Vi(1)*Vj(2)*ag22i*ag13j - Vi(2)*Vj(1)*ag13i*ag22j + Vi(2)*Vj(1)*ag12i*ag23j - Vi(1)*Vj(2)*ag23i*ag12j;%r2r4
            M(k,13) = Vi(1)*Vj(2)*ag22i*ag23j - Vi(2)*Vj(1)*ag23i*ag22j + Vi(2)*Vj(1)*ag22i*ag23j - Vi(1)*Vj(2)*ag23i*ag22j;%r2r5
            M(k,14) = Vi(1)*Vj(2)*ag22i*ag33j - Vi(2)*Vj(1)*ag33i*ag22j + Vi(2)*Vj(1)*ag32i*ag23j - Vi(1)*Vj(2)*ag23i*ag32j;%r2r6
            M(k,15) = Vi(1)*Vj(3)*ag22i*ag13j - Vi(3)*Vj(1)*ag13i*ag22j + Vi(3)*Vj(1)*ag12i*ag23j - Vi(1)*Vj(3)*ag23i*ag12j;%r2r7
            M(k,16) = Vi(1)*Vj(3)*ag22i*ag23j - Vi(3)*Vj(1)*ag23i*ag22j + Vi(3)*Vj(1)*ag22i*ag23j - Vi(1)*Vj(3)*ag23i*ag22j;%r2r8
            M(k,17) = Vi(1)*Vj(3)*ag22i*ag33j - Vi(3)*Vj(1)*ag33i*ag22j + Vi(3)*Vj(1)*ag32i*ag23j - Vi(1)*Vj(3)*ag23i*ag32j;%r2r9
            M(k,18) = Vi(1)*Vj(1)*ag32i*ag33j - Vi(1)*Vj(1)*ag33i*ag32j;%r3r3
            M(k,19) = Vi(1)*Vj(2)*ag32i*ag13j - Vi(2)*Vj(1)*ag13i*ag32j + Vi(2)*Vj(1)*ag12i*ag33j - Vi(1)*Vj(2)*ag33i*ag12j;%r3r4
            M(k,20) = Vi(1)*Vj(2)*ag32i*ag23j - Vi(2)*Vj(1)*ag23i*ag32j + Vi(2)*Vj(1)*ag22i*ag33j - Vi(1)*Vj(2)*ag33i*ag22j;%r3r5
            M(k,21) = Vi(1)*Vj(2)*ag32i*ag33j - Vi(2)*Vj(1)*ag33i*ag32j + Vi(2)*Vj(1)*ag32i*ag33j - Vi(1)*Vj(2)*ag33i*ag32j;%r3r6
            M(k,22) = Vi(1)*Vj(3)*ag32i*ag13j - Vi(3)*Vj(1)*ag13i*ag32j + Vi(3)*Vj(1)*ag12i*ag33j - Vi(1)*Vj(3)*ag33i*ag12j;%r3r7
            M(k,23) = Vi(1)*Vj(3)*ag32i*ag23j - Vi(3)*Vj(1)*ag23i*ag32j + Vi(3)*Vj(1)*ag22i*ag33j - Vi(1)*Vj(3)*ag33i*ag22j;%r3r8
            M(k,24) = Vi(1)*Vj(3)*ag32i*ag33j - Vi(3)*Vj(1)*ag33i*ag32j + Vi(3)*Vj(1)*ag32i*ag33j - Vi(1)*Vj(3)*ag33i*ag32j;%r3r9
            M(k,25) = Vi(2)*Vj(2)*ag12i*ag13j - Vi(2)*Vj(2)*ag13i*ag12j;%r4r4
            M(k,26) = Vi(2)*Vj(2)*ag12i*ag23j - Vi(2)*Vj(2)*ag23i*ag12j + Vi(2)*Vj(2)*ag22i*ag13j - Vi(2)*Vj(2)*ag13i*ag22j;%r4r5
            M(k,27) = Vi(2)*Vj(2)*ag12i*ag33j - Vi(2)*Vj(2)*ag33i*ag12j + Vi(2)*Vj(2)*ag32i*ag13j - Vi(2)*Vj(2)*ag13i*ag32j;%r4r6
            M(k,28) = Vi(2)*Vj(3)*ag12i*ag13j - Vi(3)*Vj(2)*ag13i*ag12j + Vi(3)*Vj(2)*ag12i*ag13j - Vi(2)*Vj(3)*ag13i*ag12j;%r4r7
            M(k,29) = Vi(2)*Vj(3)*ag12i*ag23j - Vi(3)*Vj(2)*ag23i*ag12j + Vi(3)*Vj(2)*ag22i*ag13j - Vi(2)*Vj(3)*ag13i*ag22j;%r4r8
            M(k,30) = Vi(2)*Vj(3)*ag12i*ag33j - Vi(3)*Vj(2)*ag33i*ag12j + Vi(3)*Vj(2)*ag32i*ag13j - Vi(2)*Vj(3)*ag13i*ag32j;%r4r9
            M(k,31) = Vi(2)*Vj(2)*ag22i*ag23j - Vi(2)*Vj(2)*ag23i*ag22j;%r5r5
            M(k,32) = Vi(2)*Vj(2)*ag22i*ag33j - Vi(2)*Vj(2)*ag33i*ag22j + Vi(2)*Vj(2)*ag32i*ag23j - Vi(2)*Vj(2)*ag23i*ag32j;%r5r6
            M(k,33) = Vi(2)*Vj(3)*ag22i*ag13j - Vi(3)*Vj(2)*ag13i*ag22j + Vi(3)*Vj(2)*ag12i*ag23j - Vi(2)*Vj(3)*ag23i*ag12j;%r5r7
            M(k,34) = Vi(2)*Vj(3)*ag22i*ag23j - Vi(3)*Vj(2)*ag23i*ag22j + Vi(3)*Vj(2)*ag22i*ag23j - Vi(2)*Vj(3)*ag23i*ag22j;%r5r8
            M(k,35) = Vi(2)*Vj(3)*ag22i*ag33j - Vi(3)*Vj(2)*ag33i*ag22j + Vi(3)*Vj(2)*ag32i*ag23j - Vi(2)*Vj(3)*ag23i*ag32j;%r5r9
            M(k,36) = Vi(2)*Vj(2)*ag32i*ag33j - Vi(2)*Vj(2)*ag33i*ag32j;%r6r6
            M(k,37) = Vi(2)*Vj(3)*ag32i*ag13j - Vi(3)*Vj(2)*ag13i*ag32j + Vi(3)*Vj(2)*ag12i*ag33j - Vi(2)*Vj(3)*ag33i*ag12j;%r6r7
            M(k,38) = Vi(2)*Vj(3)*ag32i*ag23j - Vi(3)*Vj(2)*ag23i*ag32j + Vi(3)*Vj(2)*ag22i*ag33j - Vi(2)*Vj(3)*ag33i*ag22j;%r6r8
            M(k,39) = Vi(2)*Vj(3)*ag32i*ag33j - Vi(3)*Vj(2)*ag33i*ag32j + Vi(3)*Vj(2)*ag32i*ag33j - Vi(2)*Vj(3)*ag33i*ag32j;%r6r9
            M(k,40) = Vi(3)*Vj(3)*ag12i*ag13j - Vi(3)*Vj(3)*ag13i*ag12j;%r7r7
            M(k,41) = Vi(3)*Vj(3)*ag12i*ag23j - Vi(3)*Vj(3)*ag23i*ag12j + Vi(3)*Vj(3)*ag22i*ag13j - Vi(3)*Vj(3)*ag13i*ag22j;%r7r8
            M(k,42) = Vi(3)*Vj(3)*ag12i*ag33j - Vi(3)*Vj(3)*ag33i*ag12j + Vi(3)*Vj(3)*ag32i*ag13j - Vi(3)*Vj(3)*ag13i*ag32j;%r7r9            
            M(k,43) = Vi(3)*Vj(3)*ag22i*ag23j - Vi(3)*Vj(3)*ag23i*ag22j;%r8r8
            M(k,44) = Vi(3)*Vj(3)*ag22i*ag33j - Vi(3)*Vj(3)*ag33i*ag22j + Vi(3)*Vj(3)*ag32i*ag23j - Vi(3)*Vj(3)*ag23i*ag32j;%r8r9  
            M(k,45) = Vi(3)*Vj(3)*ag32i*ag33j - Vi(3)*Vj(3)*ag33i*ag32j;%r9r9            
            M(k,20) = M(k,20)+ Vi(1)*Vj(2) - Vi(2)*Vj(1);%r3r5
            M(k,23) = M(k,23)+ Vi(1)*Vj(3) - Vi(3)*Vj(1);%r3r8
            M(k,14) = M(k,14)+ Vi(2)*Vj(1) - Vi(1)*Vj(2);%r2r6
            M(k,38) = M(k,38)+ Vi(2)*Vj(3) - Vi(3)*Vj(2);%r6r8
            M(k,17) = M(k,17)+ Vi(3)*Vj(1) - Vi(1)*Vj(3);%r2r9
            M(k,35) = M(k,35)+ Vi(3)*Vj(2) - Vi(2)*Vj(3);%r5r9
            k = k+1;
            
            %second element of Rvi x Rvj
            M(k,1)  = Vi(1)*Vj(1)*ag13i*ag11j - Vi(1)*Vj(1)*ag11i*ag13j;%r1r1
            M(k,2)  = Vi(1)*Vj(1)*ag13i*ag21j - Vi(1)*Vj(1)*ag21i*ag13j + Vi(1)*Vj(1)*ag23i*ag11j - Vi(1)*Vj(1)*ag11i*ag23j;%r1r2
            M(k,3)  = Vi(1)*Vj(1)*ag13i*ag31j - Vi(1)*Vj(1)*ag31i*ag13j + Vi(1)*Vj(1)*ag33i*ag11j - Vi(1)*Vj(1)*ag11i*ag33j;%r1r3
            M(k,4)  = Vi(1)*Vj(2)*ag13i*ag11j - Vi(2)*Vj(1)*ag11i*ag13j + Vi(2)*Vj(1)*ag13i*ag11j - Vi(1)*Vj(2)*ag11i*ag13j;%r1r4
            M(k,5)  = Vi(1)*Vj(2)*ag13i*ag21j - Vi(2)*Vj(1)*ag21i*ag13j + Vi(2)*Vj(1)*ag23i*ag11j - Vi(1)*Vj(2)*ag11i*ag23j;%r1r5
            M(k,6)  = Vi(1)*Vj(2)*ag13i*ag31j - Vi(2)*Vj(1)*ag31i*ag13j + Vi(2)*Vj(1)*ag33i*ag11j - Vi(1)*Vj(2)*ag11i*ag33j;%r1r6
            M(k,7)  = Vi(1)*Vj(3)*ag13i*ag11j - Vi(3)*Vj(1)*ag11i*ag13j + Vi(3)*Vj(1)*ag13i*ag11j - Vi(1)*Vj(3)*ag11i*ag13j;%r1r7
            M(k,8)  = Vi(1)*Vj(3)*ag13i*ag21j - Vi(3)*Vj(1)*ag21i*ag13j + Vi(3)*Vj(1)*ag23i*ag11j - Vi(1)*Vj(3)*ag11i*ag23j;%r1r8
            M(k,9)  = Vi(1)*Vj(3)*ag13i*ag31j - Vi(3)*Vj(1)*ag31i*ag13j + Vi(3)*Vj(1)*ag33i*ag11j - Vi(1)*Vj(3)*ag11i*ag33j;%r1r9
            M(k,10) = Vi(1)*Vj(1)*ag23i*ag21j - Vi(1)*Vj(1)*ag21i*ag23j;%r2r2
            M(k,11) = Vi(1)*Vj(1)*ag23i*ag31j - Vi(1)*Vj(1)*ag31i*ag23j + Vi(1)*Vj(1)*ag33i*ag21j - Vi(1)*Vj(1)*ag21i*ag33j;%r2r3
            M(k,12) = Vi(1)*Vj(2)*ag23i*ag11j - Vi(2)*Vj(1)*ag11i*ag23j + Vi(2)*Vj(1)*ag13i*ag21j - Vi(1)*Vj(2)*ag21i*ag13j;%r2r4
            M(k,13) = Vi(1)*Vj(2)*ag23i*ag21j - Vi(2)*Vj(1)*ag21i*ag23j + Vi(2)*Vj(1)*ag23i*ag21j - Vi(1)*Vj(2)*ag21i*ag23j;%r2r5
            M(k,14) = Vi(1)*Vj(2)*ag23i*ag31j - Vi(2)*Vj(1)*ag31i*ag23j + Vi(2)*Vj(1)*ag33i*ag21j - Vi(1)*Vj(2)*ag21i*ag33j;%r2r6
            M(k,15) = Vi(1)*Vj(3)*ag23i*ag11j - Vi(3)*Vj(1)*ag11i*ag23j + Vi(3)*Vj(1)*ag13i*ag21j - Vi(1)*Vj(3)*ag21i*ag13j;%r2r7
            M(k,16) = Vi(1)*Vj(3)*ag23i*ag21j - Vi(3)*Vj(1)*ag21i*ag23j + Vi(3)*Vj(1)*ag23i*ag21j - Vi(1)*Vj(3)*ag21i*ag23j;%r2r8
            M(k,17) = Vi(1)*Vj(3)*ag23i*ag31j - Vi(3)*Vj(1)*ag31i*ag23j + Vi(3)*Vj(1)*ag33i*ag21j - Vi(1)*Vj(3)*ag21i*ag33j;%r2r9
            M(k,18) = Vi(1)*Vj(1)*ag33i*ag31j - Vi(1)*Vj(1)*ag31i*ag33j;%r3r3
            M(k,19) = Vi(1)*Vj(2)*ag33i*ag11j - Vi(2)*Vj(1)*ag11i*ag33j + Vi(2)*Vj(1)*ag13i*ag31j - Vi(1)*Vj(2)*ag31i*ag13j;%r3r4
            M(k,20) = Vi(1)*Vj(2)*ag33i*ag21j - Vi(2)*Vj(1)*ag21i*ag33j + Vi(2)*Vj(1)*ag23i*ag31j - Vi(1)*Vj(2)*ag31i*ag23j;%r3r5
            M(k,21) = Vi(1)*Vj(2)*ag33i*ag31j - Vi(2)*Vj(1)*ag31i*ag33j + Vi(2)*Vj(1)*ag33i*ag31j - Vi(1)*Vj(2)*ag31i*ag33j;%r3r6
            M(k,22) = Vi(1)*Vj(3)*ag33i*ag11j - Vi(3)*Vj(1)*ag11i*ag33j + Vi(3)*Vj(1)*ag13i*ag31j - Vi(1)*Vj(3)*ag31i*ag13j;%r3r7
            M(k,23) = Vi(1)*Vj(3)*ag33i*ag21j - Vi(3)*Vj(1)*ag21i*ag33j + Vi(3)*Vj(1)*ag23i*ag31j - Vi(1)*Vj(3)*ag31i*ag23j;%r3r8
            M(k,24) = Vi(1)*Vj(3)*ag33i*ag31j - Vi(3)*Vj(1)*ag31i*ag33j + Vi(3)*Vj(1)*ag33i*ag31j - Vi(1)*Vj(3)*ag31i*ag33j;%r3r9
            M(k,25) = Vi(2)*Vj(2)*ag13i*ag11j - Vi(2)*Vj(2)*ag11i*ag13j;%r4r4
            M(k,26) = Vi(2)*Vj(2)*ag13i*ag21j - Vi(2)*Vj(2)*ag21i*ag13j + Vi(2)*Vj(2)*ag23i*ag11j - Vi(2)*Vj(2)*ag11i*ag23j;%r4r5
            M(k,27) = Vi(2)*Vj(2)*ag13i*ag31j - Vi(2)*Vj(2)*ag31i*ag13j + Vi(2)*Vj(2)*ag33i*ag11j - Vi(2)*Vj(2)*ag11i*ag33j;%r4r6
            M(k,28) = Vi(2)*Vj(3)*ag13i*ag11j - Vi(3)*Vj(2)*ag11i*ag13j + Vi(3)*Vj(2)*ag13i*ag11j - Vi(2)*Vj(3)*ag11i*ag13j;%r4r7
            M(k,29) = Vi(2)*Vj(3)*ag13i*ag21j - Vi(3)*Vj(2)*ag21i*ag13j + Vi(3)*Vj(2)*ag23i*ag11j - Vi(2)*Vj(3)*ag11i*ag23j;%r4r8
            M(k,30) = Vi(2)*Vj(3)*ag13i*ag31j - Vi(3)*Vj(2)*ag31i*ag13j + Vi(3)*Vj(2)*ag33i*ag11j - Vi(2)*Vj(3)*ag11i*ag33j;%r4r9
            M(k,31) = Vi(2)*Vj(2)*ag23i*ag21j - Vi(2)*Vj(2)*ag21i*ag23j;%r5r5
            M(k,32) = Vi(2)*Vj(2)*ag23i*ag31j - Vi(2)*Vj(2)*ag31i*ag23j + Vi(2)*Vj(2)*ag33i*ag21j - Vi(2)*Vj(2)*ag21i*ag33j;%r5r6
            M(k,33) = Vi(2)*Vj(3)*ag23i*ag11j - Vi(3)*Vj(2)*ag11i*ag23j + Vi(3)*Vj(2)*ag13i*ag21j - Vi(2)*Vj(3)*ag21i*ag13j;%r5r7
            M(k,34) = Vi(2)*Vj(3)*ag23i*ag21j - Vi(3)*Vj(2)*ag21i*ag23j + Vi(3)*Vj(2)*ag23i*ag21j - Vi(2)*Vj(3)*ag21i*ag23j;%r5r8
            M(k,35) = Vi(2)*Vj(3)*ag23i*ag31j - Vi(3)*Vj(2)*ag31i*ag23j + Vi(3)*Vj(2)*ag33i*ag21j - Vi(2)*Vj(3)*ag21i*ag33j;%r5r9
            M(k,36) = Vi(2)*Vj(2)*ag33i*ag31j - Vi(2)*Vj(2)*ag31i*ag33j;%r6r6
            M(k,37) = Vi(2)*Vj(3)*ag33i*ag11j - Vi(3)*Vj(2)*ag11i*ag33j + Vi(3)*Vj(2)*ag13i*ag31j - Vi(2)*Vj(3)*ag31i*ag13j;%r6r7
            M(k,38) = Vi(2)*Vj(3)*ag33i*ag21j - Vi(3)*Vj(2)*ag21i*ag33j + Vi(3)*Vj(2)*ag23i*ag31j - Vi(2)*Vj(3)*ag31i*ag23j;%r6r8
            M(k,39) = Vi(2)*Vj(3)*ag33i*ag31j - Vi(3)*Vj(2)*ag31i*ag33j + Vi(3)*Vj(2)*ag33i*ag31j - Vi(2)*Vj(3)*ag31i*ag33j;%r6r9
            M(k,40) = Vi(3)*Vj(3)*ag13i*ag11j - Vi(3)*Vj(3)*ag11i*ag13j;%r7r7
            M(k,41) = Vi(3)*Vj(3)*ag13i*ag21j - Vi(3)*Vj(3)*ag21i*ag13j + Vi(3)*Vj(3)*ag23i*ag11j - Vi(3)*Vj(3)*ag11i*ag23j;%r7r8
            M(k,42) = Vi(3)*Vj(3)*ag13i*ag31j - Vi(3)*Vj(3)*ag31i*ag13j + Vi(3)*Vj(3)*ag33i*ag11j - Vi(3)*Vj(3)*ag11i*ag33j;%r7r9            
            M(k,43) = Vi(3)*Vj(3)*ag23i*ag21j - Vi(3)*Vj(3)*ag21i*ag23j;%r8r8
            M(k,44) = Vi(3)*Vj(3)*ag23i*ag31j - Vi(3)*Vj(3)*ag31i*ag23j + Vi(3)*Vj(3)*ag33i*ag21j - Vi(3)*Vj(3)*ag21i*ag33j;%r8r9  
            M(k,45) = Vi(3)*Vj(3)*ag33i*ag31j - Vi(3)*Vj(3)*ag31i*ag33j;%r9r9 
            M(k,6)  = M(k,6) + Vi(1)*Vj(2) - Vi(2)*Vj(1);%r1r6
            M(k,9)  = M(k,9) + Vi(1)*Vj(3) - Vi(3)*Vj(1);%r1r9
            M(k,19) = M(k,19)+ Vi(2)*Vj(1) - Vi(1)*Vj(2);%r3r4
            M(k,30) = M(k,30)+ Vi(2)*Vj(3) - Vi(3)*Vj(2);%r4r9
            M(k,22) = M(k,22)+ Vi(3)*Vj(1) - Vi(1)*Vj(3);%r3r7
            M(k,37) = M(k,37)+ Vi(3)*Vj(2) - Vi(2)*Vj(3);%r6r7
            k = k+1;
            
            %third element of Rvi x Rvj
            M(k,1)  = Vi(1)*Vj(1)*ag11i*ag12j - Vi(1)*Vj(1)*ag12i*ag11j;%r1r1
            M(k,2)  = Vi(1)*Vj(1)*ag11i*ag22j - Vi(1)*Vj(1)*ag22i*ag11j + Vi(1)*Vj(1)*ag21i*ag12j - Vi(1)*Vj(1)*ag12i*ag21j;%r1r2
            M(k,3)  = Vi(1)*Vj(1)*ag11i*ag32j - Vi(1)*Vj(1)*ag32i*ag11j + Vi(1)*Vj(1)*ag31i*ag12j - Vi(1)*Vj(1)*ag12i*ag31j;%r1r3
            M(k,4)  = Vi(1)*Vj(2)*ag11i*ag12j - Vi(2)*Vj(1)*ag12i*ag11j + Vi(2)*Vj(1)*ag11i*ag12j - Vi(1)*Vj(2)*ag12i*ag11j;%r1r4
            M(k,5)  = Vi(1)*Vj(2)*ag11i*ag22j - Vi(2)*Vj(1)*ag22i*ag11j + Vi(2)*Vj(1)*ag21i*ag12j - Vi(1)*Vj(2)*ag12i*ag21j;%r1r5
            M(k,6)  = Vi(1)*Vj(2)*ag11i*ag32j - Vi(2)*Vj(1)*ag32i*ag11j + Vi(2)*Vj(1)*ag31i*ag12j - Vi(1)*Vj(2)*ag12i*ag31j;%r1r6
            M(k,7)  = Vi(1)*Vj(3)*ag11i*ag12j - Vi(3)*Vj(1)*ag12i*ag11j + Vi(3)*Vj(1)*ag11i*ag12j - Vi(1)*Vj(3)*ag12i*ag11j;%r1r7
            M(k,8)  = Vi(1)*Vj(3)*ag11i*ag22j - Vi(3)*Vj(1)*ag22i*ag11j + Vi(3)*Vj(1)*ag21i*ag12j - Vi(1)*Vj(3)*ag12i*ag21j;%r1r8
            M(k,9)  = Vi(1)*Vj(3)*ag11i*ag32j - Vi(3)*Vj(1)*ag32i*ag11j + Vi(3)*Vj(1)*ag31i*ag12j - Vi(1)*Vj(3)*ag12i*ag31j;%r1r9
            M(k,10) = Vi(1)*Vj(1)*ag21i*ag22j - Vi(1)*Vj(1)*ag22i*ag21j;%r2r2
            M(k,11) = Vi(1)*Vj(1)*ag21i*ag32j - Vi(1)*Vj(1)*ag32i*ag21j + Vi(1)*Vj(1)*ag31i*ag22j - Vi(1)*Vj(1)*ag22i*ag31j;%r2r3
            M(k,12) = Vi(1)*Vj(2)*ag21i*ag12j - Vi(2)*Vj(1)*ag12i*ag21j + Vi(2)*Vj(1)*ag11i*ag22j - Vi(1)*Vj(2)*ag22i*ag11j;%r2r4
            M(k,13) = Vi(1)*Vj(2)*ag21i*ag22j - Vi(2)*Vj(1)*ag22i*ag21j + Vi(2)*Vj(1)*ag21i*ag22j - Vi(1)*Vj(2)*ag22i*ag21j;%r2r5
            M(k,14) = Vi(1)*Vj(2)*ag21i*ag32j - Vi(2)*Vj(1)*ag32i*ag21j + Vi(2)*Vj(1)*ag31i*ag22j - Vi(1)*Vj(2)*ag22i*ag31j;%r2r6
            M(k,15) = Vi(1)*Vj(3)*ag21i*ag12j - Vi(3)*Vj(1)*ag12i*ag21j + Vi(3)*Vj(1)*ag11i*ag22j - Vi(1)*Vj(3)*ag22i*ag11j;%r2r7
            M(k,16) = Vi(1)*Vj(3)*ag21i*ag22j - Vi(3)*Vj(1)*ag22i*ag21j + Vi(3)*Vj(1)*ag21i*ag22j - Vi(1)*Vj(3)*ag22i*ag21j;%r2r8
            M(k,17) = Vi(1)*Vj(3)*ag21i*ag32j - Vi(3)*Vj(1)*ag32i*ag21j + Vi(3)*Vj(1)*ag31i*ag22j - Vi(1)*Vj(3)*ag22i*ag31j;%r2r9
            M(k,18) = Vi(1)*Vj(1)*ag31i*ag32j - Vi(1)*Vj(1)*ag32i*ag31j;%r3r3
            M(k,19) = Vi(1)*Vj(2)*ag31i*ag12j - Vi(2)*Vj(1)*ag12i*ag31j + Vi(2)*Vj(1)*ag11i*ag32j - Vi(1)*Vj(2)*ag32i*ag11j;%r3r4
            M(k,20) = Vi(1)*Vj(2)*ag31i*ag22j - Vi(2)*Vj(1)*ag22i*ag31j + Vi(2)*Vj(1)*ag21i*ag32j - Vi(1)*Vj(2)*ag32i*ag21j;%r3r5
            M(k,21) = Vi(1)*Vj(2)*ag31i*ag32j - Vi(2)*Vj(1)*ag32i*ag31j + Vi(2)*Vj(1)*ag31i*ag32j - Vi(1)*Vj(2)*ag32i*ag31j;%r3r6
            M(k,22) = Vi(1)*Vj(3)*ag31i*ag12j - Vi(3)*Vj(1)*ag12i*ag31j + Vi(3)*Vj(1)*ag11i*ag32j - Vi(1)*Vj(3)*ag32i*ag11j;%r3r7
            M(k,23) = Vi(1)*Vj(3)*ag31i*ag22j - Vi(3)*Vj(1)*ag22i*ag31j + Vi(3)*Vj(1)*ag21i*ag32j - Vi(1)*Vj(3)*ag32i*ag21j;%r3r8
            M(k,24) = Vi(1)*Vj(3)*ag31i*ag32j - Vi(3)*Vj(1)*ag32i*ag31j + Vi(3)*Vj(1)*ag31i*ag32j - Vi(1)*Vj(3)*ag32i*ag31j;%r3r9
            M(k,25) = Vi(2)*Vj(2)*ag11i*ag12j - Vi(2)*Vj(2)*ag12i*ag11j;%r4r4
            M(k,26) = Vi(2)*Vj(2)*ag11i*ag22j - Vi(2)*Vj(2)*ag22i*ag11j + Vi(2)*Vj(2)*ag21i*ag12j - Vi(2)*Vj(2)*ag12i*ag21j;%r4r5
            M(k,27) = Vi(2)*Vj(2)*ag11i*ag32j - Vi(2)*Vj(2)*ag32i*ag11j + Vi(2)*Vj(2)*ag31i*ag12j - Vi(2)*Vj(2)*ag12i*ag31j;%r4r6
            M(k,28) = Vi(2)*Vj(3)*ag11i*ag12j - Vi(3)*Vj(2)*ag12i*ag11j + Vi(3)*Vj(2)*ag11i*ag12j - Vi(2)*Vj(3)*ag12i*ag11j;%r4r7
            M(k,29) = Vi(2)*Vj(3)*ag11i*ag22j - Vi(3)*Vj(2)*ag22i*ag11j + Vi(3)*Vj(2)*ag21i*ag12j - Vi(2)*Vj(3)*ag12i*ag21j;%r4r8
            M(k,30) = Vi(2)*Vj(3)*ag11i*ag32j - Vi(3)*Vj(2)*ag32i*ag11j + Vi(3)*Vj(2)*ag31i*ag12j - Vi(2)*Vj(3)*ag12i*ag31j;%r4r9
            M(k,31) = Vi(2)*Vj(2)*ag21i*ag22j - Vi(2)*Vj(2)*ag22i*ag21j;%r5r5
            M(k,32) = Vi(2)*Vj(2)*ag21i*ag32j - Vi(2)*Vj(2)*ag32i*ag21j + Vi(2)*Vj(2)*ag31i*ag22j - Vi(2)*Vj(2)*ag22i*ag31j;%r5r6
            M(k,33) = Vi(2)*Vj(3)*ag21i*ag12j - Vi(3)*Vj(2)*ag12i*ag21j + Vi(3)*Vj(2)*ag11i*ag22j - Vi(2)*Vj(3)*ag22i*ag11j;%r5r7
            M(k,34) = Vi(2)*Vj(3)*ag21i*ag22j - Vi(3)*Vj(2)*ag22i*ag21j + Vi(3)*Vj(2)*ag21i*ag22j - Vi(2)*Vj(3)*ag22i*ag21j;%r5r8
            M(k,35) = Vi(2)*Vj(3)*ag21i*ag32j - Vi(3)*Vj(2)*ag32i*ag21j + Vi(3)*Vj(2)*ag31i*ag22j - Vi(2)*Vj(3)*ag22i*ag31j;%r5r9
            M(k,36) = Vi(2)*Vj(2)*ag31i*ag32j - Vi(2)*Vj(2)*ag32i*ag31j;%r6r6
            M(k,37) = Vi(2)*Vj(3)*ag31i*ag12j - Vi(3)*Vj(2)*ag12i*ag31j + Vi(3)*Vj(2)*ag11i*ag32j - Vi(2)*Vj(3)*ag32i*ag11j;%r6r7
            M(k,38) = Vi(2)*Vj(3)*ag31i*ag22j - Vi(3)*Vj(2)*ag22i*ag31j + Vi(3)*Vj(2)*ag21i*ag32j - Vi(2)*Vj(3)*ag32i*ag21j;%r6r8
            M(k,39) = Vi(2)*Vj(3)*ag31i*ag32j - Vi(3)*Vj(2)*ag32i*ag31j + Vi(3)*Vj(2)*ag31i*ag32j - Vi(2)*Vj(3)*ag32i*ag31j;%r6r9
            M(k,40) = Vi(3)*Vj(3)*ag11i*ag12j - Vi(3)*Vj(3)*ag12i*ag11j;%r7r7
            M(k,41) = Vi(3)*Vj(3)*ag11i*ag22j - Vi(3)*Vj(3)*ag22i*ag11j + Vi(3)*Vj(3)*ag21i*ag12j - Vi(3)*Vj(3)*ag12i*ag21j;%r7r8
            M(k,42) = Vi(3)*Vj(3)*ag11i*ag32j - Vi(3)*Vj(3)*ag32i*ag11j + Vi(3)*Vj(3)*ag31i*ag12j - Vi(3)*Vj(3)*ag12i*ag31j;%r7r9            
            M(k,43) = Vi(3)*Vj(3)*ag21i*ag22j - Vi(3)*Vj(3)*ag22i*ag21j;%r8r8
            M(k,44) = Vi(3)*Vj(3)*ag21i*ag32j - Vi(3)*Vj(3)*ag32i*ag21j + Vi(3)*Vj(3)*ag31i*ag22j - Vi(3)*Vj(3)*ag22i*ag31j;%r8r9  
            M(k,45) = Vi(3)*Vj(3)*ag31i*ag32j - Vi(3)*Vj(3)*ag32i*ag31j;%r9r9          
            M(k,12) = M(k,12)+ Vi(1)*Vj(2) - Vi(2)*Vj(1);%r2r4
            M(k,15) = M(k,15)+ Vi(1)*Vj(3) - Vi(3)*Vj(1);%r2r7
            M(k,5)  = M(k,5) + Vi(2)*Vj(1) - Vi(1)*Vj(2);%r1r5
            M(k,33) = M(k,33)+ Vi(2)*Vj(3) - Vi(3)*Vj(2);%r5r7
            M(k,8)  = M(k,8) + Vi(3)*Vj(1) - Vi(1)*Vj(3);%r1r8
            M(k,29) = M(k,29)+ Vi(3)*Vj(2) - Vi(2)*Vj(3);%r4r8
            k = k+1;
        end
    end
end
%according to equation (11), R^T*R = I
M(k,1)  = 1;   M(k,10) = 1;    M(k,18) = 1;  M(k,46) = -1;     k = k+1; %r1*r1+r2*r2+r3*r3 = 1;
M(k,25) = 1;   M(k,31) = 1;    M(k,36) = 1;  M(k,46) = -1;     k = k+1; %r4*r4+r5*r5+r6*r6 = 1;
M(k,40) = 1;   M(k,43) = 1;    M(k,45) = 1;  M(k,46) = -1;     k = k+1; %r7*r7+r8*r8+r9*r9 = 1;
M(k,4)  = 1;   M(k,13) = 1;    M(k,21) = 1;                    k = k+1; %r1*r4+r2*r5+r3*r6 = 0;
M(k,7)  = 1;   M(k,16) = 1;    M(k,24) = 1;                    k = k+1; %r1*r7+r2*r8+r3*r9 = 0;
M(k,28) = 1;   M(k,34) = 1;    M(k,39) = 1;                    k = k+1; %r4*r7+r5*r8+r6*r9 = 0;
%according to equation (11), R*R^T = I
M(k,1)  = 1;   M(k,25) = 1;    M(k,40) = 1;  M(k,46) = -1;     k = k+1; %r1*r1+r4*r4+r7*r7 = 1;
M(k,10) = 1;   M(k,31) = 1;    M(k,43) = 1;  M(k,46) = -1;     k = k+1; %r2*r2+r5*r5+r8*r8 = 1;
M(k,18) = 1;   M(k,36) = 1;    M(k,45) = 1;  M(k,46) = -1;     k = k+1; %r3*r3+r6*r6+r9*r9 = 1;
M(k,2)  = 1;   M(k,26) = 1;    M(k,41) = 1;                    k = k+1; %r1*r2+r4*r5+r7*r8 = 0;
M(k,3)  = 1;   M(k,27) = 1;    M(k,42) = 1;                    k = k+1; %r1*r3+r4*r6+r7*r9 = 0;
M(k,11) = 1;   M(k,32) = 1;    M(k,44) = 1;                    k = k+1; %r2*r3+r5*r6+r8*r9 = 0;

%solve the linear system M * xbar = 0  using SVD
[U, S, W] = svd(M);
if n>=5 %for 5 or more lines, ker(M) is one-dimensional. 
    xbar = W(:,46); % the last column of W;
    xbar = xbar/xbar(46); % the condition rho = 1; rho is the last element of xbar.    
else %for four lines,  ker(M) is 8-dimensional (46 - 38).
    % xbar = lamda1 * w1 + ... + lamda6 * w6;
    rankM = rank(M);
    if rankM < 38        
        error('rankM<38, Input data is degenerated');
    end
    w = zeros(46,8);
    w(:, 1:8) =  W(:,39:46);
    % according to equation (5), build the K matrix
    % we have K * lamda = 0 where lamda = (lamda1*lamda1, lamda1*lamda2, lamda1*lamda3, ..., lamda1*lamda8,
    %   lamda2*lamda2, lamda2*lamda3, ..., lamda2*lamda8, ..., lamda7*lamda7, lamda7*lamda8, lamda8*lamda8);
    K = zeros(126*2, 36);
    row = 1;
    for i = 1:9
        for j = i+1:9
            for k = j+1:9
                for l = k+1:9
                    index_ij = 9*(i-1) - i*(i-1)/2 + j;% the index of element rirj in the vector xbar
                    index_kl = 9*(k-1) - k*(k-1)/2 + l;% the index of element rkrl in the vector xbar
                    index_ik = 9*(i-1) - i*(i-1)/2 + k;% the index of element rirj in the vector xbar
                    index_jl = 9*(j-1) - j*(j-1)/2 + l;% the index of element rjrl in the vector xbar
                    index_il = 9*(i-1) - i*(i-1)/2 + l;% the index of element rirl in the vector xbar
                    index_jk = 9*(j-1) - j*(j-1)/2 + k;% the index of element rjrk in the vector xbar  
                    %rirj * rkrl = rirk * rjrl
                    K(row,1)   = w(index_ij, 1) * w(index_kl, 1) - w(index_ik, 1) * w(index_jl, 1); %lamda1*lamda1; 
                    K(row,2)   = w(index_ij, 1) * w(index_kl, 2) - w(index_ik, 1) * w(index_jl, 2) + w(index_ij, 2) * w(index_kl, 1) - w(index_ik, 2) * w(index_jl, 1); %lamda1*lamda2; 
                    K(row,3)   = w(index_ij, 1) * w(index_kl, 3) - w(index_ik, 1) * w(index_jl, 3) + w(index_ij, 3) * w(index_kl, 1) - w(index_ik, 3) * w(index_jl, 1); %lamda1*lamda3; 
                    K(row,4)   = w(index_ij, 1) * w(index_kl, 4) - w(index_ik, 1) * w(index_jl, 4) + w(index_ij, 4) * w(index_kl, 1) - w(index_ik, 4) * w(index_jl, 1); %lamda1*lamda4; 
                    K(row,5)   = w(index_ij, 1) * w(index_kl, 5) - w(index_ik, 1) * w(index_jl, 5) + w(index_ij, 5) * w(index_kl, 1) - w(index_ik, 5) * w(index_jl, 1); %lamda1*lamda5; 
                    K(row,6)   = w(index_ij, 1) * w(index_kl, 6) - w(index_ik, 1) * w(index_jl, 6) + w(index_ij, 6) * w(index_kl, 1) - w(index_ik, 6) * w(index_jl, 1); %lamda1*lamda6; 
                    K(row,7)   = w(index_ij, 1) * w(index_kl, 7) - w(index_ik, 1) * w(index_jl, 7) + w(index_ij, 7) * w(index_kl, 1) - w(index_ik, 7) * w(index_jl, 1); %lamda1*lamda7; 
                    K(row,8)   = w(index_ij, 1) * w(index_kl, 8) - w(index_ik, 1) * w(index_jl, 8) + w(index_ij, 8) * w(index_kl, 1) - w(index_ik, 8) * w(index_jl, 1); %lamda1*lamda8; 
                    K(row,9)   = w(index_ij, 2) * w(index_kl, 2) - w(index_ik, 2) * w(index_jl, 2); %lamda2*lamda2; 
                    K(row,10)  = w(index_ij, 2) * w(index_kl, 3) - w(index_ik, 2) * w(index_jl, 3) + w(index_ij, 3) * w(index_kl, 2) - w(index_ik, 3) * w(index_jl, 2); %lamda2*lamda3; 
                    K(row,11)  = w(index_ij, 2) * w(index_kl, 4) - w(index_ik, 2) * w(index_jl, 4) + w(index_ij, 4) * w(index_kl, 2) - w(index_ik, 4) * w(index_jl, 2); %lamda2*lamda4; 
                    K(row,12)  = w(index_ij, 2) * w(index_kl, 5) - w(index_ik, 2) * w(index_jl, 5) + w(index_ij, 5) * w(index_kl, 2) - w(index_ik, 5) * w(index_jl, 2); %lamda2*lamda5; 
                    K(row,13)  = w(index_ij, 2) * w(index_kl, 6) - w(index_ik, 2) * w(index_jl, 6) + w(index_ij, 6) * w(index_kl, 2) - w(index_ik, 6) * w(index_jl, 2); %lamda2*lamda6; 
                    K(row,14)  = w(index_ij, 2) * w(index_kl, 7) - w(index_ik, 2) * w(index_jl, 7) + w(index_ij, 7) * w(index_kl, 2) - w(index_ik, 7) * w(index_jl, 2); %lamda2*lamda7; 
                    K(row,15)  = w(index_ij, 2) * w(index_kl, 8) - w(index_ik, 2) * w(index_jl, 8) + w(index_ij, 8) * w(index_kl, 2) - w(index_ik, 8) * w(index_jl, 2); %lamda2*lamda8; 
                    K(row,16)  = w(index_ij, 3) * w(index_kl, 3) - w(index_ik, 3) * w(index_jl, 3); %lamda3*lamda3;                     
                    K(row,17)  = w(index_ij, 3) * w(index_kl, 4) - w(index_ik, 3) * w(index_jl, 4) + w(index_ij, 4) * w(index_kl, 3) - w(index_ik, 4) * w(index_jl, 3); %lamda3*lamda4; 
                    K(row,18)  = w(index_ij, 3) * w(index_kl, 5) - w(index_ik, 3) * w(index_jl, 5) + w(index_ij, 5) * w(index_kl, 3) - w(index_ik, 5) * w(index_jl, 3); %lamda3*lamda5; 
                    K(row,19)  = w(index_ij, 3) * w(index_kl, 6) - w(index_ik, 3) * w(index_jl, 6) + w(index_ij, 6) * w(index_kl, 3) - w(index_ik, 6) * w(index_jl, 3); %lamda3*lamda6; 
                    K(row,20)  = w(index_ij, 3) * w(index_kl, 7) - w(index_ik, 3) * w(index_jl, 7) + w(index_ij, 7) * w(index_kl, 3) - w(index_ik, 7) * w(index_jl, 3); %lamda3*lamda7; 
                    K(row,21)  = w(index_ij, 3) * w(index_kl, 8) - w(index_ik, 3) * w(index_jl, 8) + w(index_ij, 8) * w(index_kl, 3) - w(index_ik, 8) * w(index_jl, 3); %lamda3*lamda8; 
                    K(row,22)  = w(index_ij, 4) * w(index_kl, 4) - w(index_ik, 4) * w(index_jl, 4); %lamda4*lamda4; 
                    K(row,23)  = w(index_ij, 4) * w(index_kl, 5) - w(index_ik, 4) * w(index_jl, 5) + w(index_ij, 5) * w(index_kl, 4) - w(index_ik, 5) * w(index_jl, 4); %lamda4*lamda5; 
                    K(row,24)  = w(index_ij, 4) * w(index_kl, 6) - w(index_ik, 4) * w(index_jl, 6) + w(index_ij, 6) * w(index_kl, 4) - w(index_ik, 6) * w(index_jl, 4); %lamda4*lamda6; 
                    K(row,25)  = w(index_ij, 4) * w(index_kl, 7) - w(index_ik, 4) * w(index_jl, 7) + w(index_ij, 7) * w(index_kl, 4) - w(index_ik, 7) * w(index_jl, 4); %lamda4*lamda7; 
                    K(row,26)  = w(index_ij, 4) * w(index_kl, 8) - w(index_ik, 4) * w(index_jl, 8) + w(index_ij, 8) * w(index_kl, 4) - w(index_ik, 8) * w(index_jl, 4); %lamda4*lamda8; 
                    K(row,27)  = w(index_ij, 5) * w(index_kl, 5) - w(index_ik, 5) * w(index_jl, 5); %lamda5*lamda5; 
                    K(row,28)  = w(index_ij, 5) * w(index_kl, 6) - w(index_ik, 5) * w(index_jl, 6) + w(index_ij, 6) * w(index_kl, 5) - w(index_ik, 6) * w(index_jl, 5); %lamda5*lamda6; 
                    K(row,29)  = w(index_ij, 5) * w(index_kl, 7) - w(index_ik, 5) * w(index_jl, 7) + w(index_ij, 7) * w(index_kl, 5) - w(index_ik, 7) * w(index_jl, 5); %lamda5*lamda7; 
                    K(row,30)  = w(index_ij, 5) * w(index_kl, 8) - w(index_ik, 5) * w(index_jl, 8) + w(index_ij, 8) * w(index_kl, 5) - w(index_ik, 8) * w(index_jl, 5); %lamda5*lamda8; 
                    K(row,31)  = w(index_ij, 6) * w(index_kl, 6) - w(index_ik, 6) * w(index_jl, 6); %lamda6*lamda6; 
                    K(row,32)  = w(index_ij, 6) * w(index_kl, 7) - w(index_ik, 6) * w(index_jl, 7) + w(index_ij, 7) * w(index_kl, 6) - w(index_ik, 7) * w(index_jl, 6); %lamda6*lamda7; 
                    K(row,33)  = w(index_ij, 6) * w(index_kl, 8) - w(index_ik, 6) * w(index_jl, 8) + w(index_ij, 8) * w(index_kl, 6) - w(index_ik, 8) * w(index_jl, 6); %lamda6*lamda8; 
                    K(row,34)  = w(index_ij, 7) * w(index_kl, 7) - w(index_ik, 7) * w(index_jl, 7); %lamda7*lamda7; 
                    K(row,35)  = w(index_ij, 7) * w(index_kl, 8) - w(index_ik, 7) * w(index_jl, 8) + w(index_ij, 8) * w(index_kl, 7) - w(index_ik, 8) * w(index_jl, 7); %lamda7*lamda8;
                    K(row,36)  = w(index_ij, 8) * w(index_kl, 8) - w(index_ik, 8) * w(index_jl, 8); %lamda8*lamda8; 
                    row = row + 1;
                    %rirj * rkrl = rirl * rjrk
                    K(row,1)   = w(index_ij, 1) * w(index_kl, 1) - w(index_il, 1) * w(index_jk, 1); %lamda1*lamda1; 
                    K(row,2)   = w(index_ij, 1) * w(index_kl, 2) - w(index_il, 1) * w(index_jk, 2) + w(index_ij, 2) * w(index_kl, 1) - w(index_il, 2) * w(index_jk, 1); %lamda1*lamda2; 
                    K(row,3)   = w(index_ij, 1) * w(index_kl, 3) - w(index_il, 1) * w(index_jk, 3) + w(index_ij, 3) * w(index_kl, 1) - w(index_il, 3) * w(index_jk, 1); %lamda1*lamda3; 
                    K(row,4)   = w(index_ij, 1) * w(index_kl, 4) - w(index_il, 1) * w(index_jk, 4) + w(index_ij, 4) * w(index_kl, 1) - w(index_il, 4) * w(index_jk, 1); %lamda1*lamda4; 
                    K(row,5)   = w(index_ij, 1) * w(index_kl, 5) - w(index_il, 1) * w(index_jk, 5) + w(index_ij, 5) * w(index_kl, 1) - w(index_il, 5) * w(index_jk, 1); %lamda1*lamda5; 
                    K(row,6)   = w(index_ij, 1) * w(index_kl, 6) - w(index_il, 1) * w(index_jk, 6) + w(index_ij, 6) * w(index_kl, 1) - w(index_il, 6) * w(index_jk, 1); %lamda1*lamda6; 
                    K(row,7)   = w(index_ij, 1) * w(index_kl, 7) - w(index_il, 1) * w(index_jk, 7) + w(index_ij, 7) * w(index_kl, 1) - w(index_il, 7) * w(index_jk, 1); %lamda1*lamda7; 
                    K(row,8)   = w(index_ij, 1) * w(index_kl, 8) - w(index_il, 1) * w(index_jk, 8) + w(index_ij, 8) * w(index_kl, 1) - w(index_il, 8) * w(index_jk, 1); %lamda1*lamda8; 
                    K(row,9)   = w(index_ij, 2) * w(index_kl, 2) - w(index_il, 2) * w(index_jk, 2); %lamda2*lamda2; 
                    K(row,10)  = w(index_ij, 2) * w(index_kl, 3) - w(index_il, 2) * w(index_jk, 3) + w(index_ij, 3) * w(index_kl, 2) - w(index_il, 3) * w(index_jk, 2); %lamda2*lamda3; 
                    K(row,11)  = w(index_ij, 2) * w(index_kl, 4) - w(index_il, 2) * w(index_jk, 4) + w(index_ij, 4) * w(index_kl, 2) - w(index_il, 4) * w(index_jk, 2); %lamda2*lamda4; 
                    K(row,12)  = w(index_ij, 2) * w(index_kl, 5) - w(index_il, 2) * w(index_jk, 5) + w(index_ij, 5) * w(index_kl, 2) - w(index_il, 5) * w(index_jk, 2); %lamda2*lamda5; 
                    K(row,13)  = w(index_ij, 2) * w(index_kl, 6) - w(index_il, 2) * w(index_jk, 6) + w(index_ij, 6) * w(index_kl, 2) - w(index_il, 6) * w(index_jk, 2); %lamda2*lamda6; 
                    K(row,14)  = w(index_ij, 2) * w(index_kl, 7) - w(index_il, 2) * w(index_jk, 7) + w(index_ij, 7) * w(index_kl, 2) - w(index_il, 7) * w(index_jk, 2); %lamda2*lamda7; 
                    K(row,15)  = w(index_ij, 2) * w(index_kl, 8) - w(index_il, 2) * w(index_jk, 8) + w(index_ij, 8) * w(index_kl, 2) - w(index_il, 8) * w(index_jk, 2); %lamda2*lamda8; 
                    K(row,16)  = w(index_ij, 3) * w(index_kl, 3) - w(index_il, 3) * w(index_jk, 3); %lamda3*lamda3;                     
                    K(row,17)  = w(index_ij, 3) * w(index_kl, 4) - w(index_il, 3) * w(index_jk, 4) + w(index_ij, 4) * w(index_kl, 3) - w(index_il, 4) * w(index_jk, 3); %lamda3*lamda4; 
                    K(row,18)  = w(index_ij, 3) * w(index_kl, 5) - w(index_il, 3) * w(index_jk, 5) + w(index_ij, 5) * w(index_kl, 3) - w(index_il, 5) * w(index_jk, 3); %lamda3*lamda5; 
                    K(row,19)  = w(index_ij, 3) * w(index_kl, 6) - w(index_il, 3) * w(index_jk, 6) + w(index_ij, 6) * w(index_kl, 3) - w(index_il, 6) * w(index_jk, 3); %lamda3*lamda6; 
                    K(row,20)  = w(index_ij, 3) * w(index_kl, 7) - w(index_il, 3) * w(index_jk, 7) + w(index_ij, 7) * w(index_kl, 3) - w(index_il, 7) * w(index_jk, 3); %lamda3*lamda7; 
                    K(row,21)  = w(index_ij, 3) * w(index_kl, 8) - w(index_il, 3) * w(index_jk, 8) + w(index_ij, 8) * w(index_kl, 3) - w(index_il, 8) * w(index_jk, 3); %lamda3*lamda8; 
                    K(row,22)  = w(index_ij, 4) * w(index_kl, 4) - w(index_il, 4) * w(index_jk, 4); %lamda4*lamda4; 
                    K(row,23)  = w(index_ij, 4) * w(index_kl, 5) - w(index_il, 4) * w(index_jk, 5) + w(index_ij, 5) * w(index_kl, 4) - w(index_il, 5) * w(index_jk, 4); %lamda4*lamda5; 
                    K(row,24)  = w(index_ij, 4) * w(index_kl, 6) - w(index_il, 4) * w(index_jk, 6) + w(index_ij, 6) * w(index_kl, 4) - w(index_il, 6) * w(index_jk, 4); %lamda4*lamda6; 
                    K(row,25)  = w(index_ij, 4) * w(index_kl, 7) - w(index_il, 4) * w(index_jk, 7) + w(index_ij, 7) * w(index_kl, 4) - w(index_il, 7) * w(index_jk, 4); %lamda4*lamda7; 
                    K(row,26)  = w(index_ij, 4) * w(index_kl, 8) - w(index_il, 4) * w(index_jk, 8) + w(index_ij, 8) * w(index_kl, 4) - w(index_il, 8) * w(index_jk, 4); %lamda4*lamda8; 
                    K(row,27)  = w(index_ij, 5) * w(index_kl, 5) - w(index_il, 5) * w(index_jk, 5); %lamda5*lamda5; 
                    K(row,28)  = w(index_ij, 5) * w(index_kl, 6) - w(index_il, 5) * w(index_jk, 6) + w(index_ij, 6) * w(index_kl, 5) - w(index_il, 6) * w(index_jk, 5); %lamda5*lamda6; 
                    K(row,29)  = w(index_ij, 5) * w(index_kl, 7) - w(index_il, 5) * w(index_jk, 7) + w(index_ij, 7) * w(index_kl, 5) - w(index_il, 7) * w(index_jk, 5); %lamda5*lamda7; 
                    K(row,30)  = w(index_ij, 5) * w(index_kl, 8) - w(index_il, 5) * w(index_jk, 8) + w(index_ij, 8) * w(index_kl, 5) - w(index_il, 8) * w(index_jk, 5); %lamda5*lamda8; 
                    K(row,31)  = w(index_ij, 6) * w(index_kl, 6) - w(index_il, 6) * w(index_jk, 6); %lamda6*lamda6; 
                    K(row,32)  = w(index_ij, 6) * w(index_kl, 7) - w(index_il, 6) * w(index_jk, 7) + w(index_ij, 7) * w(index_kl, 6) - w(index_il, 7) * w(index_jk, 6); %lamda6*lamda7; 
                    K(row,33)  = w(index_ij, 6) * w(index_kl, 8) - w(index_il, 6) * w(index_jk, 8) + w(index_ij, 8) * w(index_kl, 6) - w(index_il, 8) * w(index_jk, 6); %lamda6*lamda8; 
                    K(row,34)  = w(index_ij, 7) * w(index_kl, 7) - w(index_il, 7) * w(index_jk, 7); %lamda7*lamda7; 
                    K(row,35)  = w(index_ij, 7) * w(index_kl, 8) - w(index_il, 7) * w(index_jk, 8) + w(index_ij, 8) * w(index_kl, 7) - w(index_il, 8) * w(index_jk, 7); %lamda7*lamda8;
                    K(row,36)  = w(index_ij, 8) * w(index_kl, 8) - w(index_il, 8) * w(index_jk, 8); %lamda8*lamda8;                   
                    row = row + 1;
                end
            end
        end
    end    
    %solve the linear system K * lamda = 0  using SVD
    [uk,sk,wk] = svd(K);
    lamda = wk(:,36); %the last comlun of wk.
    if lamda(1) < 0
        lamda = -lamda;
    end
    lamda1 = sqrt(lamda(1)); %lamda1*lamda1 = lamda(1);
    lamda2 = sign(lamda1) * sign(lamda(2)) * sqrt(abs(lamda(9)));  % lamda1*lamda2 = lamda(2); lamda2*lamda2 =  lamda(9); 
    lamda3 = sign(lamda1) * sign(lamda(3)) * sqrt(abs(lamda(16))); % lamda1*lamda3 = lamda(3); lamda3*lamda3 =  lamda(16); 
    lamda4 = sign(lamda1) * sign(lamda(4)) * sqrt(abs(lamda(22))); % lamda1*lamda4 = lamda(4); lamda4*lamda4 =  lamda(22); 
    lamda5 = sign(lamda1) * sign(lamda(5)) * sqrt(abs(lamda(27))); % lamda1*lamda5 = lamda(5); lamda5*lamda5 =  lamda(27); 
    lamda6 = sign(lamda1) * sign(lamda(6)) * sqrt(abs(lamda(31))); % lamda1*lamda6 = lamda(6); lamda6*lamda6 =  lamda(31); 
    lamda7 = sign(lamda1) * sign(lamda(7)) * sqrt(abs(lamda(34))); % lamda1*lamda7 = lamda(7); lamda7*lamda7 =  lamda(34); 
    lamda8 = sign(lamda1) * sign(lamda(8)) * sqrt(abs(lamda(36))); % lamda1*lamda8 = lamda(8); lamda8*lamda8 =  lamda(36); 
    scale = lamda1 *  w(46, 1) + lamda2 *  w(46, 2) + lamda3 *  w(46, 3) + lamda4 *  w(46, 4) ...
          + lamda5 *  w(46, 5) + lamda6 *  w(46, 6) + lamda7 *  w(46, 7) + lamda8 *  w(46, 8);
    lamda1 = lamda1/scale; lamda2 = lamda2/scale; lamda3 = lamda3/scale; lamda4 = lamda4/scale;...
    lamda5 = lamda5/scale; lamda6 = lamda6/scale; lamda7 = lamda7/scale; lamda8 = lamda8/scale;
    xbar = lamda1 * w(:, 1) + lamda2 * w(:, 2) + lamda3 * w(:, 3) + lamda4 * w(:, 4)...
         + lamda5 * w(:, 5) + lamda6 * w(:, 6) + lamda7 * w(:, 7) + lamda8 * w(:, 8);
end

if xbar(1) <0
    xbar = -1 * xbar;
end
r1 = sqrt(xbar(1));%r1*r1 = xbar(1);
r2 = sign(r1) * sign(xbar(2)) * sqrt(abs(xbar(10))); % r1*r2 = xbar(2); r2*r2 =  xbar(10);
r3 = sign(r1) * sign(xbar(3)) * sqrt(abs(xbar(18))); % r1*r3 = xbar(3); r3*r3 =  xbar(18);
r4 = sign(r1) * sign(xbar(4)) * sqrt(abs(xbar(25))); % r1*r4 = xbar(4); r4*r4 =  xbar(25);
r5 = sign(r1) * sign(xbar(5)) * sqrt(abs(xbar(31))); % r1*r5 = xbar(5); r5*r5 =  xbar(31);
r6 = sign(r1) * sign(xbar(6)) * sqrt(abs(xbar(36))); % r1*r6 = xbar(6); r6*r6 =  xbar(36);
r7 = sign(r1) * sign(xbar(7)) * sqrt(abs(xbar(40))); % r1*r7 = xbar(7); r7*r7 =  xbar(40);
r8 = sign(r1) * sign(xbar(8)) * sqrt(abs(xbar(43))); % r1*r8 = xbar(8); r8*r8 =  xbar(43);
r9 = sign(r1) * sign(xbar(9)) * sqrt(abs(xbar(45))); % r1*r9 = xbar(9); r9*r9 =  xbar(45);
rot = [r1, r4, r7; r2,r5,r8; r3,r6,r9];
rot = rot/det(rot); % make sure det(rot) = 1;

%now, estimate the translation vector; Build the matrix A;
A = zeros(3*n, 4);
row = 1;
for i =1:n
    rp      = rot * P(i,:)';
    gama_i  = gama(i,:);
    d_i     = d(i,:);
    A(row, 1) = gama_i(1) * gama_i(1); 
    A(row, 2) = gama_i(2) * gama_i(1);
    A(row, 3) = gama_i(3) * gama_i(1) - d_i(1);
    A(row, 4) = rp(1) * gama_i(1) * gama_i(1) + rp(2) * gama_i(2) * gama_i(1) + rp(3) * gama_i(3) * gama_i(1) - rp(3) * d_i(1);
    row = row + 1;
    A(row, 1) = gama_i(1) * gama_i(2); 
    A(row, 2) = gama_i(2) * gama_i(2);
    A(row, 3) = gama_i(3) * gama_i(2) - d_i(2);
    A(row, 4) = rp(1) * gama_i(1) * gama_i(2) + rp(2) * gama_i(2) * gama_i(2) + rp(3) * gama_i(3) * gama_i(2) - rp(3) * d_i(2);
    row = row + 1;
    A(row, 1) = gama_i(1) * gama_i(3); 
    A(row, 2) = gama_i(2) * gama_i(3);
    A(row, 3) = gama_i(3) * gama_i(3) - d_i(3);
    A(row, 4) = rp(1) * gama_i(1) * gama_i(3) + rp(2) * gama_i(2) * gama_i(3) + rp(3) * gama_i(3) * gama_i(3) - rp(3) * d_i(3);
    row = row + 1;
end

[ua, sa, va] = svd(A);
pos1 = va(1,4)/va(4,4);
pos2 = va(2,4)/va(4,4);
pos3 = va(3,4)/va(4,4);
Pos_wc = [pos1, pos2, pos3]';%Point_c = R_wc * Point_w + Pos_wc

Rot_cw = rot';% rot is estimated rotation matrix Rot_wc, so Rot_cw = Rot_wc'.
Pos_cw = -Rot_cw*Pos_wc;%in this case: Point_w = R_cw * Point_c + Pos_cw, so Pos_cw = -Rot_wc' * Pos_wc
return
