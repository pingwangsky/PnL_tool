function err = projectionErr(xs, xe, P1w, P2w, R0)

%first compute the weight of each line and the normal of the interpretation plane passing through to camera center and the line 
n = length(xs);
w  = zeros(n,1);
nc = zeros(3,n);
for i=1:n
    w(i) = 1/norm(xs(:,i)-xe(:,i)); % the weight of a line is the inverse of its image length
    temp = xcross(xs(:,i),xe(:,i));  
    nc(:,i) = temp/norm(temp);   
end

s = CGR(R0);

A= zeros(n*2,12);
for i= 1:n
    nx= nc(1,i); ny= nc(2,i); nz= nc(3,i);
    px= P1w(1,i); py= P1w(2,i); pz= P1w(3,i);
    A(i*2-1,:)= [nx*px ny*px nz*px nx*py ny*py nz*py nx*pz ny*pz nz*pz nx ny nz];
    px= P2w(1,i); py= P2w(2,i); pz= P2w(3,i);
    A(i*2,:)= [nx*px ny*px nz*px nx*py ny*py nz*py nx*pz ny*pz nz*pz nx ny nz];
end

A1= A(:,1:9);
A2= A(:,10:12);

A4= -inv(A2.'*A2)*A2.'*A1;
A3= A1+A2*A4;
M= A3.'*A3;

R= CGR2(s);
x = R(:);
t= A4*x;

err = x' * M * x;

return

function c = xcross(a,b)

c = [a(2)*b(3)-a(3)*b(2);
     a(3)*b(1)-a(1)*b(3);
     a(1)*b(2)-a(2)*b(1)];
 
return
 
function s= CGR(R)

A= R.';

q4= sqrt(1+A(1,1)+A(2,2)+A(3,3))/2;

if q4 > 0.01
	q1= (A(3,2)-A(2,3))/q4/4;
	q2= (A(1,3)-A(3,1))/q4/4;
	q3= (A(2,1)-A(1,2))/q4/4;
else
	q1= sqrt(1+A(1,1)-A(2,2)-A(3,3))/2;
	q2= (A(1,2)+A(2,1))/q1/4;
	q3= (A(1,3)+A(3,1))/q1/4;
	q4= (A(3,2)-A(2,3))/q1/4;
end

s= [q1; q2; q3]/q4;

return

function R= CGR2(s)

s1= s(1);
s2= s(2);
s3= s(3);

R= [ s1^2 - s2^2 - s3^2 + 1,   2*s3 + 2*s1*s2,           2*s1*s3 - 2*s2;
     2*s1*s2 - 2*s3, - s1^2 + s2^2 - s3^2 + 1,           2*s1 + 2*s2*s3;
     2*s2 + 2*s1*s3,           2*s2*s3 - 2*s1, - s1^2 - s2^2 + s3^2 + 1];
 
R= R/(1+s1*s1+s2*s2+s3*s3);

return


