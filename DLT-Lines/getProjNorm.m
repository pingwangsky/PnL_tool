% norm of a projection plane of a 2D line
function nl = getProjNorm(p1,p2)

	n = size(p1,2);
	d1 = [p1; ones(1,n)];
	d2 = [p2; ones(1,n)];
	nl = cross(d1,d2,1);
	nl = xnorm(nl);
