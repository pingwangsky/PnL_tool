%ROTZ Rotation about Z axis
%
%	TR = ROTZ(theta)
%
% Returns a homogeneous transformation representing a rotation of theta 
% about the Z axis.
%
% See also: ROTX, ROTY, ROTVEC.

% $Log: rotz.m,v $
% Revision 1.1.1.1  2007/06/28 03:45:21  faraz
% this is a first try to set bundle adjustment on CVS (version 6) 
%
% Revision 1.2  2002/04/01 11:47:16  pic
% General cleanup of code: help comments, see also, copyright, remnant dh/dyn
% references, clarification of functions.
%
% $Revision: 145 $
% Copyright (C) 1993-2002, by Peter I. Corke

function r = rotz(t)
	ct = cos(t);
	st = sin(t);
	r =    [ct	-st	0
		st	ct	0
		0	0	1];
