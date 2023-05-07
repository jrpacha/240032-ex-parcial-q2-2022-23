clearvars
close all

% BVP
% (a1(x)u')' + a0(x) u = f(x), x \in (a,b)
% u(a) = 1,
% u(b) = alpha
%
% with,
% a = 0, b = pi/4,
% a1(x) = 1+(tan(x))^2,
% a0(x) = a0,     constant value passed to the function as parameter
% f(x) = x,
% alpha = 2.0,    for parts (a) and (b)     

% Data
a0 = 6;
n = 100;
xp = 0.388;
a = 0; b = pi/4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n1 = n+1;     %number of nodes
h = (b-a)/n; %length of the elements

nodes = linspace(a,b,n1)'; %position of the nodes
elem = [(1:n)',(2:n1)'];   %connectivity matrix

numNodes = size(nodes,1);
numElem = size(elem,1);

K = zeros(n1);
F = zeros(n1,1);
Q = zeros(n1,1);

for e=1:numElem
    rows = [elem(e,1), elem(e,2)];
    cols = rows;
    x1 = nodes(rows(1,1),1); x2 = nodes(rows(1,2),1);
    Ke = (tan(x2)-tan(x1))*[1, -1;-1, 1]/h^2 + a0*h*[2,1;1,2]/6;
    Fe = h*[2*x1+x2; x1+2*x2]/6;
    K(rows,cols) = K(rows,cols) + Ke;
    F(rows)=F(rows)+Fe;
end

%BC
fixedNodes = [1,numNodes];
freeNodes = setdiff(1:numNodes,fixedNodes);

%Natural BC
Q(freeNodes) = 0; %Not necessary, since Q was initialised to 0

%Essential BC
u = zeros(numNodes,1);
u(1) = 1;
u(numNodes) = 2;

%Reduced system
Fm = F(freeNodes) + Q(freeNodes) - K(freeNodes, fixedNodes)*u(fixedNodes);
Km = K(freeNodes,freeNodes);

%Solve the reduced system:
um = Km\Fm;
u(freeNodes) = um;

%Part (a)
interpU = interp1(nodes,u,xp);
u15 = u(15);
clc
fprintf('Problem 3\n')
fprintf('Part (a) Interpolated value of U: %.4e\n',interpU)
fprintf('Check: straightforward calculation %.4e\n',...
    u(50)+(u(51)-u(50))*(xp-nodes(50))/h)
fprintf('Hint: u(15) = %.4e\n',u15)

%Part (b)
avU = sum(u)/numNodes;
fprintf('Part (b) Averaged value of u, <u> = %.4e\n',avU)

%Part (c)
fixedNodes = 1;
freeNodes = setdiff(1:numNodes, fixedNodes);
Fm = F(freeNodes) + Q(freeNodes) - K(freeNodes, fixedNodes)*u(fixedNodes);
Km = K(freeNodes, freeNodes);

Km(end,:) = 1;
Fm(end) = numNodes - u(1);

um = Km\Fm;
u(freeNodes) = um;
alpha = u(end);

fprintf('Part (c) alpha = %.4e\n', alpha)
fprintf('Check: value of the new average, <u> = %e\n', ...
    sum(u)/numNodes)
