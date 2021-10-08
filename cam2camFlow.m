function [flow,rhoGrid,thetaGrid] = cam2camFlow(eyeHeight,tvec,cam1_basis,cam2_basis,scale_factor)

%% function takes two eye orientation matrices / basis vectors, eye height, tvec, and the scale factor and computes the resulting retinal flow
%INPUTS:
% eyeHeight: height of the eye above the ground planee
% tvec: translation vector of the eye through space
% cam1_basis: the orientation of the eye before tvec
% cam2_basis: the orientation of the eye after tvec
% scale_factor: how much to scale the speeds by, since we have to
% discretize, this is to properly scale things depending on walking speed

%% params
grid_res = 50; % resolution of the retinal flow array in pixels
max_ecc = 45; % maximum eccentricity in degrees of visual angle

%% create retinal coordinate grid to compute vectors
% we use meshgrid to create a [grid_res x grid_res] array
[xx, yy] = meshgrid(1:grid_res,1:grid_res);
yy = flipud(yy); % flip the y so that +y is upwards, so that in the resulting array upwards is up in world coordinates

% assign corresponding eccentricities and polar angles of each point
% in the retinal array (the center will be 0 or close depending on rounding
% shenanigans). Rho (eccentricity) is assigned by taking the distance from
% the center of the [xx,yy] grid in pixels, dividing it by grid_res/2, and
% then multiplying by the maximum eccentricity. This maps the pixel units
% linearly to degrees of eccentricity, when measuring distance from the
% center. So the middle is 0, and then all the way to the right is max_ecc.
% Theta is just the polar angle about the center.
rhoGrid = reshape(vecnorm([xx(:) yy(:)]-[grid_res/2 grid_res/2],2,2)/(grid_res/2)*(deg2rad(max_ecc)),[grid_res grid_res]);
thetaGrid = reshape(atan2(yy(:)-grid_res/2,xx(:)-grid_res/2),[grid_res grid_res]);


%% convert retinal coordinate grids into 3D vectors

% here we convert (theta,rho) into their corresponding 3D vectors, using
% spherical coordinates. So the Z-component is just cos(rho) (at 0
% eccentricity, cos(0) = 1, so the vector is straight ahead ([0,0,1]). The
% x and y components are proportional to cos(theta) and sin(theta)
% respectively, but then scaled by sin(rho), so the more eccentric they
% are, the more influenced by theta.
shootVecs = [cos(thetaGrid(:)).*sin(rhoGrid(:)) sin(thetaGrid(:)).*sin(rhoGrid(:)) cos(rhoGrid(:))];


%% intersect the 3D vectors with the ground

% set origin of these vectors at eyeHeight above the ground
cam1_pos = [0 eyeHeight 0];

% rotate these vectors by the eye orientation, so they are now in world
% coordinates (they were in eye coordinates before with [0, 0, 1] being eye
% direction, now eye direction is row 3 of cam1_basis
shootVecs_cam1 = shootVecs*cam1_basis;

% calculate indeces of intersections that we will effectively set at
% infinity (above the horizon). This will eventually result in motion in
% the visual field that is only influenced by rotation, not by translation
% (points at infinity don't move when you move)
nan_dex = shootVecs_cam1(:,2)>-1e-12;

% calculate ground intersection of these vectors. This is done by taking
% the Y component of each vector, and dividing the eye height by it. Then
% the resulting scale factor is used to scale up each of the vectors so
% that it intersects with the ground plane.
gp_intersect_cam1 = -cam1_pos(:,2)./shootVecs_cam1(:,2).*shootVecs_cam1;

% then translate each of these ground intersection points by negated eye
% translation vector (points move in the opposite direction of the eye in
% eye coordinates)
gp_intersect_cam2 = gp_intersect_cam1 - tvec;

%% compute new positions of ground points in eye orientation #2 relative coordinates

% normalize the new ground point intersection points so that they are unit
% vectors
receiveVecs_cam2 = normr(gp_intersect_cam2);

% rotate them into the reference frame of eye orientation #2 (this
% incorporates the rotation caused by fixating the ground)
receiveVecs_cam2 = receiveVecs_cam2*cam2_basis';

%% calculate movement of points in visual space by referencing corresponding points in eye orientation #1 and #2

% eye #1 vectors
points1 = shootVecs;
% eye #2 vectors
points2 = receiveVecs_cam2;

% substitute in the points at infinity calculate earlier, by instead of
% intersecting with the ground and translating, just use the direction
% vectors from before and rotate them by eye#1 to eye#2 rotation
points2(nan_dex,:) = normr(shootVecs_cam1(nan_dex,:)*cam2_basis');

% calculate the displacement speed in degrees of visual angle. This is just
% the angle between each corresponding direction vector for eye#1 and eye#2
displacementMag = 2*atan2(vecnorm(normr(points1)-normr(points2),2,2),...
    vecnorm(normr(points1)+normr(points2),2,2));

% calculate the displacement direction, this is done in the XY plane, so
% we're ignoring the eye direction dimension for defining direction
displacementXY = normr(points2)-normr(points1);
displacementXY = normr(displacementXY(:,1:2)); % normalize so we have unit vector

% then we scale each displacementXY by the magnitude, so that the magnitude
% of the resulting 2D vector is the speed in degrees/second
displacementXY = displacementXY.*rad2deg(displacementMag)*scale_factor;

% reshape the arrays such that we have two 2D arrays, one for velocity X
% and one for velocity Y. These will match the original polar eye
% coordinate grids (rhoGrid, thetaGrid)
retVx = reshape(displacementXY(:,1),[grid_res grid_res]);
retVy = -reshape(displacementXY(:,2),[grid_res grid_res]);

% create matlab optical flow object (nice for plotting, and automatically
% stores the direction array, as Orientation, and speed, as Magnitude, as
% well as Vx and Vy).
flow = opticalFlow(retVx,retVy);

end


