function [R_cw, T_cw] = Mirzaei(p1, p2, P1_w, P2_w)

	[V_w P_w] = getVP(P1_w, P2_w);

	nLine = length(p1);
	p1 = [p1; ones(1,nLine)];
	p2 = [p2; ones(1,nLine)];
	
    [R_wc, T_wc] = AlgLS(p1, p2, V_w, P_w);
    
    R_cw = R_wc.';
    T_cw = -R_cw * T_wc;


function [rot_cw, pos_cw] = AlgLS(xs,xe,Vw,Pw)
%This function follows the AlgLS algorithm proposed by Faraz and Stergios: 
% Globally optimal pose estimation from line corrsepondences.
%input: xs(:, i) = the start point of the ith image line [startpointx, startpointy, 1];
%       xe(:, i) = the end point of the ith image line   [endpointx,   endpointy, 1];
%       Vw(:, i) = the direction of ith line in the world frame
%       Pw(:, i) = a point of ith line in the world frame
%output: rot_cw = the orientation of camera in rotation matrix parametrization
%                 (W_w = rot_cw * V_c)
%        pos_cw = the position of camera in global frame;
%                 (P_w = rot_cw * P_c + pos_cw;

addpath('AlgLSUtils','AlgLSUtils/Cexp','AlgLSUtils/robotics3D');%add the necessary file path

n = size(xs,2);
if(n~=size(xe,2) || n~=size(Pw,2) || n~=size(Vw,2)) 
    error('Input data xs, xe, P1w and P2w are inconsistent'); 
end

if(n<4)
    error('The input data is too few to determine a unique camera pose');
end

%compute the normal of the interpretation plane passing through to camera center and the line 
nc = zeros(3,n);
for i=1:n
    temp = cross(xs(:,i),xe(:,i));  
    nc(:,i) = temp/norm(temp); 
end
%call Faraz's EstimateVanishingPoints() function, which returns the
%rotation matrix rot_cw in quaternion form.
I_moment_n = nc;
G_line = Vw;
[q_sols, residuals, info] = EstimateVanishingPoints(I_moment_n, G_line, 'relaxed');
rot_cw = quat2rot(q_sols);

%estimate translation vector.
Mat = zeros(2*n, 4);
for i = 1:n
    nw = rot_cw * nc(:,i);
    Uw = cross(Pw(:,i),Vw(:,i));
    Vxi = Vw(1,i);     Vyi = Vw(2,i);    Vzi = Vw(3,i);
    nxi = nw(1);       nyi = nw(2);      nzi = nw(3);
    Uxi = Uw(1);       Uyi = Uw(2);      Uzi = Uw(3);
    Mat(2*i-1, 1) =  nxi*Vzi;             Mat(2*i, 1) =  nyi*Vyi + nzi*Vzi;
    Mat(2*i-1, 2) =  nyi*Vzi;             Mat(2*i, 2) = -nyi*Vxi;
    Mat(2*i-1, 3) = -nxi*Vxi - nyi*Vyi;   Mat(2*i, 3) = -nzi*Vxi;
    Mat(2*i-1, 4) =  nxi*Uyi - nyi*Uxi;   Mat(2*i, 4) =  nzi*Uyi - nyi*Uzi;
end
%solve the linear system Mat * [tx, ty, tz, 1]' = 0  using SVD,
[UMat, SMat, VMat] = svd(Mat);
vec = VMat(:,4);% the last column of Vmat;
vec = vec/vec(4); %the condition that the last element of vec should be 1.
%now we get the translation vector pos_cw
pos_cw = vec(1:3);

return
