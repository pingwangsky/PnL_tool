function L = Left(x) 

L = [ x(4)*eye(3) - skewsymm(x(1:3)) x(1:3) ; -x(1:3)' x(4)];