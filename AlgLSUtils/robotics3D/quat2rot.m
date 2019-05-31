function [A]=quat2rot(q)

% quaternions_to_rot_matrix(q) takes a unit quaternion and returns the corresponding rotational matrix

A(1,1) = q(1,1)^2 - q(2,1)^2 - q(3,1)^2 + q(4,1)^2  ;
A(1,2) = 2*( q(1,1)*q(2,1) + q(3,1)*q(4,1) ); 
A(1,3) = 2*( q(1,1)*q(3,1) - q(2,1)*q(4,1) );

A(2,1) = 2*( q(1,1)*q(2,1) - q(3,1)*q(4,1) );
A(2,2) = -q(1,1)^2 + q(2,1)^2 - q(3,1)^2 + q(4,1)^2;
A(2,3) = 2*( q(2,1)*q(3,1) + q(1,1)*q(4,1) );

A(3,1) = 2*( q(1,1)*q(3,1) + q(2,1)*q(4,1) );
A(3,2) = 2*( q(2,1)*q(3,1) - q(1,1)*q(4,1) );
A(3,3) = -q(1,1)^2 - q(2,1)^2 + q(3,1)^2 + q(4,1)^2;