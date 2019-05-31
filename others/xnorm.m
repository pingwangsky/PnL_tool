% normalize each column
function V1 = xnorm(V0)

	V1 = V0 .* kron(ones(3,1), 1./sqrt(sum(V0.^2)));
