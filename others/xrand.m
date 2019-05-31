function y= xrand(m,n,r)
%function description:
%产生范围在r(1)到r(2)之间的m行，n列的随机数;
y= r(1)+rand(m,n)*(r(2)-r(1));

return