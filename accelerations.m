close all
clearvars

%% params
walking_speed = 1500; %walking speed in mm / s
calc_displacement = 1; % how small you want to discretize and calculate velocities over (I think you can just leave this as is)
gravity_angle = 45; % gaze angle relative to straight downwards (degrees)
about_y_angle = 15; % gaze angle about the vertical axis, with 0 being straight ahead (in plane with [0,0,1]) (degrees)
eyeHeight = 1800; % height of the eye above the ground plane (mm)

%% compute / convert
tvec = [0 0 calc_displacement]; % translation vector due to walking
gravity_angle = deg2rad(gravity_angle); % convert to radians
about_y_angle = deg2rad(about_y_angle);

% calculate sequential eye basis vectors
[basis1,basis2] = twoBasesGivenTvec(gravity_angle,about_y_angle,tvec,eyeHeight);

% calculate a scale factor depending on desired walking speed and
% discretization parameter
scale_factor = walking_speed/calc_displacement;

% calculate retinal flow based on eye basis vectors, translation, and scale
% factor
[flow,~] = cam2camFlow(eyeHeight,tvec,basis1,basis2,scale_factor);

%% for acceleration
gravity_angle = asin(basis2(2,2)); % new gravity angle is just the eye orientation 2 gravity angle from previous section
 % this will allow calculation of the next flow field so you can approximate acceleration

about_y_angle = acos(basis2(1,1)); % new about y angle also just uses the value from eye orientation 2

% same as before, but now for the next set of eye orientations
[basis1,basis2] = twoBasesGivenTvec(gravity_angle,about_y_angle,tvec,eyeHeight);
[flow2,rhoGrid,thetaGrid] = cam2camFlow(eyeHeight,tvec,basis1,basis2,scale_factor);

%% visualization
figure(1)
clf
subplot(1,2,1)
imagesc(flow.Magnitude)
title('Frame 1 flow magnitude (deg/s)')
colorbar
subplot(1,2,2)
plot(flow,'decimationfactor',[3 3]);
axis ij
xlim([1 size(flow.Magnitude,1)]);
ylim([1 size(flow.Magnitude,1)]);
title('Frame 1 flow directions');

figure(2)
clf
subplot(1,2,1)
imagesc(flow2.Magnitude)
title('Frame 2 flow magnitude (deg/s)')
colorbar
subplot(1,2,2)
plot(flow2,'decimationfactor',[3 3]);
axis ij
xlim([1 size(flow.Magnitude,1)]);
ylim([1 size(flow.Magnitude,1)]);
title('Frame 2 flow directions');

figure(3)
clf
imagesc((flow2.Magnitude-flow.Magnitude)*scale_factor);
title('Acceleration from frame1 to frame2 in deg/s/s')
colorbar

% retinal coordinates for reference
figure(4)
clf
subplot(1,2,1)
imagesc(rad2deg(rhoGrid));
colorbar
colormap(gca,parula);
title('Eccentricity (deg)');
subplot(1,2,2);
imagesc(rad2deg(thetaGrid));
colorbar
colormap(gca,hsv)
title('Polar angle (deg), 0 = right');
