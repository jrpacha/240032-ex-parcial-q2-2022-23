%Problem 2

clearvars
close all
clc

%Data
perm = 1; xp = 7.13; xq = 4.92; xt = 4.0;

x = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024];
y = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10];

%Compute the slopes
dy1 = 1/0.69315; dyn = dy1/1024;
n = length(x);
dy = [dy1, (y(3:n)-y(1:n-2))./(x(3:n)-x(1:n-2)), dyn];

fprintf('\tProblem 2 (%d)\n',perm);

[sp, sq, slope_segment, slope, p] = p2(x, y, xp, xq, xt);

%Part (a)
%Cubic splines
fprintf('Part (a)\n');
fprintf('Spline approximation for log2(%.2f),\n', xp(perm));
fprintf('       log2(%.2f) \x2248 %.4e\n', xp(perm), sp);
fprintf('Hint: the same approximation for x = log2(%.2f) gives,\n', xq(perm));
fprintf('       log2(%.2f) \x2248 %.4e\n\n', xq(perm), sq);

%Part (b)
%Compute the slopes
fprintf('Part (b)\n');
fprintf('%10s%6s%12s\n','X', 'Y', 'approx Y''');
for i = 1:n
    if i == find(abs(xt(perm)-x) < 2*eps)
        mark = '->';
    else
        mark = '';
    end
    fprintf('%3s%7.1f%6.1f%12.4e\n', mark, x(i), y(i), dy(i));
end
fprintf('\nValue estimated for the slope at x = %.1f: %.4e\n',...
    xt(perm), slope);
fprintf('Hint: the slope of the segment from (x3,y3) to\n');
fprintf('        (x7,y7) is m = %.4e\n', slope_segment);

%Part (c)
yp = polyval(p, xp);
yq = polyval(p, xq);
fprintf('\nPart (c)\n');
fprintf('The value of log2(%.2f) for this approximation is\n', xp(perm));
fprintf('       log2(%.2f) \x2248 %.4e\n', xp(perm), yp);
fprintf('Hint. Using this estimation method\n');
fprintf('       log2(%.2f) \x2248 %.4e\n\n', xq(perm), yq);

%Plot spline vs Hermite interpolation 
%for the interval 4 < x < 8
xx = linspace(x(3), x(4), 4001);
%Function f(x) = log2(x)
fxx = log2(xx);
%Cubic Spline approximation
sp = spline(x, y, xx);
% Hermite interpolation
yy = polyval(p, xx);
figure()
hold on
plot(xx, fxx, '-', 'Color', 'black', 'LineWidth', 1)
plot(xx, sp, '-', 'Color', 'red', 'Linewidth', 1)
plot(xx, yy, '-', 'Color', 'blue', 'LineWidth', 1)
xlabel('$x$', 'Interpreter', 'LaTex', 'FontSize', 14)
ylabel('$y\qquad$', 'Interpreter', 'latex', 'FontSize', 14, 'Rotation', 360)
legend('$y =\log_{2}(x)$', 'Cubic spline', 'Hermite interp.',...
    'Location', 'SouthEast', 'interpreter', 'LaTeX', 'fontSize', 10)
title('Cubic spline v.s.~deg.~3 Hermite interpolation for $y = \log_{2}(x)$',...
    'FontSize', 12, 'Interpreter', 'LaTeX')
hold off
