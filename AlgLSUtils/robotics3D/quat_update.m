function qUpdated = quat_update(a, b)

if length(a) == 3 && length(b) == 4
    q = [a/2; 1];
    q = q/norm(q);
    qUpdated = quat_mul(q,b);
    
elseif length(a) == 4 && length(b) == 3
    q = [b/2; 1];
    q = q/norm(q);
    qUpdated = quat_mul(a,q);
    
else
    error('one input should be a small angle 3D vector and the other should be a quaternion!')
end