function dth = quat_diff(q1,q2)

dq = quat_mul(quat_inv(q1),q2);

dth = dq(1:3,:)*2;