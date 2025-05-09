function theta = getVectorAngle(x,y)
%GETVECTORANGLE returns the angle between two input vectors x and y in
%degrees
%   Detailed explanation goes here
nx = norm(x);
ny = norm(y);

theta = acosd(dot(x,y)/(nx*ny));
end

