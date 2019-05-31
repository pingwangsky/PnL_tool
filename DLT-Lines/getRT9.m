function [R T] = getRT9(X)

	X = X * sign(X(12));
	
	R = reshape(X(1:9),3,3);
	sc = mean(sqrt(sum(R.^2))) * sign(det(R));
	lm = 1 / sc;
	R = R * lm;
	
%	R = rodrigues(rodrigues(R));
	
	T = X(10:12) * lm;
