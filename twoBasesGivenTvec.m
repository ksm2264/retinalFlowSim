function [basis1,basis2] = twoBasesGivenTvec(gravity_angle,about_y_angle,tvec,eyeHeight)
%% function gives two eye coordinate basis vectors given an eye height, and eye orientation in space
% INPUTS:
% gravity_angle: angle of gaze in world relative to gravity (straight
% down), which is [0, -1, 0]
% about_y_angle: angle of gaze in world about Y axis, with 0 degrees being
% straight ahead, which is [0, 0, 1]. Positive rotates towards +X ([1, 0,
% 0])
% tvec: translation vector of eye through space while fixating
% eyeHeight: height of the eye above the ground plane



%% first function computes the basis vectors for the first eye pose, rows are vectors
% first row = right direction (we constrain this to always be in the XZ
% plane, although this isn't accurate since the eye can inherit head tilts,
% like when you put your ear towards your shooulder)
% second row = up
% third row = eye direction
R = [cos(about_y_angle) 0 -sin(about_y_angle);
    cos(gravity_angle)*sin(about_y_angle) sin(gravity_angle) cos(gravity_angle)*cos(about_y_angle);
    sin(gravity_angle)*sin(about_y_angle) -cos(gravity_angle) sin(gravity_angle)*cos(about_y_angle)];

%% next we compute the gaze on ground intersection (a), 
% then shift it by tvec, capturing the movement of the gaze on ground location in the eye coordinate frame

% a is gaze on ground intersection, which we compute by dividing eye height
% by the Y component of the eye direction (third row of R), then scaling
% the eye direction unit vector by the result. This result is a scalar
% which tells you how many eye direction unit vectors you need to repeat to
% reach the ground plane
a = eyeHeight/cos(gravity_angle)*R(3,:);

% b is where the ground intersection location translates to in eye relative
% coordinates, so just negate the eye in world translation vector and apply
% that to a
b = a - tvec;

% normalize both so now they are unit direction vectors
a = normr(a);
b = normr(b);

% w is the axis about which the eye will be rotating, here we're applying
% Listing's law (https://en.wikipedia.org/wiki/Listing%27s_law). This
% states that the eye physically rotates around a single axis of rotation,
% which is achieved through specific synergistic activation of the
% different oculmotor muscles. Thus we only need this one axis as well as
% an angle to rotate about this axis to specify a rotation. The cross
% product tells us what this axis is because it gives us a direction
% perpendicular to both a and b
w = normr(cross(a,b));

% phi is the actual angle between a and b, in other words it is the angle
% that the eye needs to rotate around the axis w, in order to fixate on a
% as it moves due to translation tvec
phi = 2*atan2(norm(a-b),norm(a+b));

% here we can specify a gain factor, epsilon. This will let us scale how
% well the eye is fixating the initial location. For example epsilon=0.5
% would have the eye only rotate halfway between wherever a translates to
% due to tvec
epsilon = 1;

% K is the 'symmetric skew' matrix, which is used in the formula for a
% rotation matrix given an axis angle rotation (see:
% https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle)
K = [0 -w(3) w(2);
     w(3) 0 -w(1);
     -w(2) w(1) 0];

% here we're just plugging everything into that formula. The resulting
% rotation matrix R2 can be used to left multiply column vectors, or the
% transpose can be used to right multiply row vectors, in order to rotate
% them by the rotation needed to rotate the eye from the intial ground
% fixated location, the the new location of that initially fixated
% location after applying translation tvec
R2 = cos(epsilon*phi)*eye(3) + sin(epsilon*phi)*K + (1-cos(epsilon*phi))*(w'*w);

% we then apply R2 to the initial eye pose, to yield the second eye pose,
% which is the eye pose that fixates a as it translates to b. We'll call
% the two poses basis1 and basis2.
basis2 = R*R2';
basis1 = R;

end






