% eig -- the eigenvalue is sorted
function [V,D] = xeig(A)

	[V,D] = eig(A);
	if ~issorted(diag(D))
		[V,D] = eig(A);
		[D,I] = sort(diag(D));
		V = V(:, I);
	end
