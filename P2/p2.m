function [sp, sq, slope_segment, slope, p] = p2(xSample, ySample, xp, xq, xt)
% Problem 2
% INPUT
%  xSample: column vector with the x values of the table
%  ySample: column vector with the y value of the table
%    diffy: approximate values of the derivatives at each point x of the
%           sample       
%       xp: abscissa of the point (being not in the sample) at which we 
%           approximation by cubic splines in part (a) is computed
%       xq: abscissa of the point (being not in the sample) at which hints 
%           in parts (a) and (c) are given
%       xt: point of at which at the approximation for the slope is
%           computed.
% OUTPUT
%            sp: spline approximation at point xp
%            sq: spline approximation at point xq
% slope_segment: slope of the segments joining the points (x(4),y(4)) and
%                (x(8),y(8)) of the sample.
%         slope: approximation of the slope at xt computed from the points
%                of the sample 

%Part (a)
%Cubic splines
sp = spline(xSample,ySample,xp);
sq = spline(xSample,ySample,xq);

%Part(b)
id1 = find(abs(xt-xSample) < 2*eps);
if id1 == 1 || id1 == length(xSample)
     error('Cannot approximate slope neither for the 1st nor for the last point of the sample\n')
end
slope_segment = (ySample(8)-ySample(4))/(xSample(8)-xSample(4));
slope = (ySample(id1+1)-ySample(id1-1))/(xSample(id1+1)-xSample(id1-1));

%Part (c)
%Hermite interpolation
id1 = find(xp > xSample, 1, 'last');  % xSample(1) < ... < xSample(end)
id2 = find(xp < xSample, 1, 'first'); % the point xp is placed at interval [xSample(id1), xSample(id2)]
if id1 == 1 || id2 == length(xSample)
    error('Cannot approximate slope neither for the 1st nor for the last point of the sample\n')
end
diffyId1 = (ySample(id2)- ySample(id1-1))/(xSample(id2)- xSample(id1-1));
diffyId2 = (ySample(id2+1) - ySample(id1))/(xSample(id2+1)- xSample(id1));
A = [xSample(id1)^3 xSample(id1)^2 xSample(id1) 1;
     xSample(id2)^3 xSample(id2)^2 xSample(id2) 1;
     3*xSample(id1)^2 2*xSample(id1) 1 0;
     3*xSample(id2)^2 2*xSample(id2) 1 0];
%dy(3) = 1/4/log(2);
%dy(4) = 1/8/log(2);
b = [ySample(id1); ySample(id2); diffyId1; diffyId2];
p =A\b;

end

