function [R_cw, T_cw] = LPnL_Bar_LS(p1, p2, P1_w, P2_w)
% the line is expressed by the start and end points
% inputs:
%	 p1: 2d projection of the start point
%	 p2: 2d projection of the end point
%	 P1_w: 3d coordinates of the start point in world frame
%	 P2_w: 3d coordinates of the end point in world frame
% outputs:
%	 R: estimated rotation
%	 T: estimated translation

	n = length(p1);	
	
	% get base points
	ct_w = mean([P1_w P2_w],2);
	BP_w = [eye(3) zeros(3,1)] + kron([1 1 1 1],ct_w);
	alp1_w = getAlpha(BP_w,P1_w);
	alp2_w = getAlpha(BP_w,P2_w);

	% normal of porjected lines
	nl = getProjNorm(p1,p2);

	% matrix
	M1 = kron([1 1 1 1], [nl nl]');
	M2 = kron([alp1_w alp2_w]', [1 1 1]);
	M = M1 .* M2;
	
	MTM = M' * M;
	[XV XD] = xeig(MTM);
	X = normBP(XV(:,1));
	
    BP_c = reshape(X,3,4);
	[R_cw T_cw] = getRT(BP_w.',BP_c.');
	
