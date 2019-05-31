function plotframe(q,p)

% get a quaternion and plot a frame rotated by it.

or = [0;0;0];

plotarrow(or,[1;0;0],'r');
hold on
plotarrow(or,[0;1;0],'b');
plotarrow(or,[0;0;1],'g');

R = quat2rot(q);

if nargin < 2
    p = or;
end

plotarrow(p,R'*[1;0;0],'r--');
plotarrow(p,R'*[0;1;0],'b--');
plotarrow(p,R'*[0;0;1],'g--');
 
hold off

% R = quat2rot(q);
% 
% or = [0;0;0];
% if nargin < 2
%     p = or;
% end
% 
% e = eye(3);
% 
% X = [or(1)*ones(3,1) ; p(1)*ones(3,1)];
% Y = [or(2)*ones(3,1) ; p(2)*ones(3,1)];
% Z = [or(3)*ones(3,1) ; p(3)*ones(3,1)];
% 
% U = [e(:,1); R'*e(:,1)];
% V = [e(:,2); R'*e(:,2)];
% W = [e(:,3); R'*e(:,3)];

%quiver3(X,Y,Z,U,V,W,'LineWidth',2,...
%    'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',12);