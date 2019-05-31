function [V P] = getVP(P1,P2)

	V = P2 - P1;
	V = xnorm(V);
	P = (P1 + P2) * 0.5;
