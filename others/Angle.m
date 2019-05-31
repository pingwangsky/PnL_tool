function Angles = Angle( rotation )
%calculation the  three components of the rotation and translation;

% Angles(1,1)=atan(rotation(1,2)/rotation(1,1))/pi*180; %yaw;
% Angles(2,1)=-asin(rotation(1,3))/pi*180; %pitch;
% Angles(3,1)=atan(rotation(2,3)/rotation(3,3))/pi*180; %roll;

Angles(1,1)=atan(rotation(1,2)/rotation(1,1)); %yaw;
Angles(2,1)=-asin(rotation(1,3)); %pitch;
Angles(3,1)=atan(rotation(2,3)/rotation(3,3)); %roll;

end

