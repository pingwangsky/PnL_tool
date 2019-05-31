%ROTY Rotation about Y axis
%
%	TR = ROTY(theta)
%
% Returns a homogeneous transformation representing a rotation of theta 
% about the Y axis.
%
% See also: ROTX, ROTZ, ROTVEC.

% $Log: roty.m,v $
% Revision 1.1.1.1  2007/06/28 03:45:21  faraz
% this is a first try to set bundle adjustment on CVS (version 6) 
%
% Revision 1.2  2002/04/01 11:47:16  pic
% General cleanup of code: help comments, see also, copyright, remnant dh/dyn
% references, clarification of functions.
%
% $Revision: 493 $
% Copyright (C) 1993-2002, by Peter I. Corke

function r = roty(t)
	ct = cos(t);
	st = sin(t);
	r =    [ct	0	st
		0	1	0
		-st	0	ct];
