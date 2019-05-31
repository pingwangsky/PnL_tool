function [R T] = DLT_Lines(p1, p2, P1_w, P2_w)
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
	
	nl = getProjNorm(p1,p2);

	M1 = kron([1 1 1], [nl nl]');
	M2 = kron([P1_w P2_w]', [1 1 1]);
	M = [M1 .* M2 [nl nl]'];
	D = ones(2*n,1);

	MTM = M' * diag(D) * M;
	[XV XD] = xeig(MTM);
	
	X = XV(:,1);	
	[R T] = getRT9(X);

