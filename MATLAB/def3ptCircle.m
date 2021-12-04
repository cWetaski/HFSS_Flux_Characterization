function [center,radius] = def3ptCircle(a,b,c)
%DEF3POINTCIRCLE obtain center and radius of circle given 3 reference points (a,b,c)
% Each point must have form p = [x, y]
% Method from J. Riggs, https://www.mathworks.com/matlabcentral/answers/uploaded_files/207922/Circle%20Calculations.pdf

% Generate some variables for ease of coding
dy_ab = b(2) - a(2); % y distance from a to b
dx_ab = b(1) - a(1); % x distance from a to b
dy_bc = c(2) - b(2); % y distance from b to c
dx_bc = c(1) - b(1); % x distance from b to c
dy_ac = c(2) - a(2); % y distance from a to c
dx_ac = c(1) - a(1);  % x distance from a to c

% Get slope of line ab
if dx_ab == 0 % a and b have same x coordinate, slope is infinite
    m_ab = inf;
else
    m_ab = dy_ab/dx_ab; % slope = rise over run
end

% Get slope of line bc
if dx_bc == 0 
    m_bc = inf;
else
    m_bc = dy_bc/dx_bc;
end

% Get slope of line ac
if dx_ac == 0
    m_ac = inf;
else
    m_ac = dy_ac/dx_ac;
end

% check for colinearity
if m_ab == m_bc
    error('points are colinear, circle cannot be defined from three colinear points')
end

% choose two points to define a line that passes thru center of circle. By
% default, points a and b are selected. However, if slope of ab is 0, then
% then points a and c are used instead.

p1 = a;

if m_ab ~= 0 % slope of ab is not 0, so we can use points a and b
    p2 = b;
    m0 = m_ab;
    p3 = c;
else % we must use points a and c instead
    p2 = c;
    m0 = m_ac;
    p3 = b;
end

m = -1/m0;

% compute midpoint of line from p1 to p2;

mid = [(p1(1) + p2(1))/2, (p1(2) + p2(2))/2]; % midpoint is average of x and y coordinates of p1 and p2, respectively

% compute y intercept of line with slope m which passes thru midpoint;
yint = mid(2) - m*mid(1);

% find center point, d - the point which is equidistant to points a, b, and c
% and is on line y = mx + yint
d = [0 0];

d(1) = (p3(1)^2 + p3(2)^2 - p1(1)^2 - p1(2)^2 + 2*yint*(p1(2) - p3(2)))/(2*(p3(1) - p1(1) + m*(p3(2)-p1(2))));
d(2) = m*d(1) + yint;

center = d;

% check radii are equal
d_ad = sqrt((d(2) - a(2))^2 + (d(1) - a(1))^2);
d_bd = sqrt((d(2) - b(2))^2 + (d(1) - b(1))^2);
d_cd = sqrt((d(2) - c(2))^2 + (d(1) - c(1))^2);

if abs(d_ad - d_bd) < 0.001 && abs(d_ad - d_cd) < 0.001
    radius = d_ad;
else
    error('Something went wrong, centerpoint is not equidistant from points');
end

end
