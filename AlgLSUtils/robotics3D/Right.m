function R = Right(x) 

R = [ x(4)*eye(3) + skewsymm(x(1:3)) x(1:3) ; -x(1:3)' x(4)];