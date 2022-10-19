% Generate a vector 
x = rand(1, 10^7);
y = rand(1, 10^7);
z = rand(1, 10^7);
% x, y and z should be mapped (with care)
% if not multiplied by the same number, they are not really random?
x = pi * x;
y = 2 * pi * y;
z = 2 * ( z - 0.5 );
minAllowableDistance = 0.0005; % Particles diameter
numberOfParticles    = 100000;  % Number of particles
numberOfPoints = 10^7;
% Initialize first point.
keeperX = x(1);
keeperY = y(1);
keeperZ = z(1);
% Try dropping down more points.
counter = 2;
for k = 2 : numberOfPoints
	% Get a trial point.
	thisX = x(k);
	thisY = y(k);
    thisZ = z(k);
	% See how far it is away from existing keeper points.
	distances = sqrt((thisX-keeperX).^2 + (thisY - keeperY).^2 + (thisZ - keeperZ).^2);
	minDistance = min(distances);
    
	if minDistance >= minAllowableDistance
		keeperX(counter) = thisX;
		keeperY(counter) = thisY;
        keeperZ(counter) = thisZ;
		counter = counter + 1;
    end
    
    if counter > numberOfParticles 
        break
    end
end
plot3(keeperX, keeperY, keeperZ, 'b*');
grid on;

fileID = fopen('RandomParticles_100000.txt','w');
fprintf(fileID,'%u\r\n',numberOfParticles);
for ii = 1:numberOfParticles
    fprintf(fileID, '%E %E %E\r\n', keeperX(ii), keeperY(ii), keeperZ(ii));
end
fclose(fileID);
