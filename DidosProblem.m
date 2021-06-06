%{ 
Dido's optimization problem: Using Matlab's constrained optimization function, fmincon, 
implement an approximation of Dido's Problem to find an approximate 
solution using finite dimensional optimization/nonlinear programming. 
%}

%% clear workspace
clear all 
clc

%% specify problem parameters
N = 50; % number of perturbations (segments of string length)
L = 10; % total fixed length of string (initial given, fixed for all trials)
theta = linspace(0,90,N); % create a vector of segmented angles between 0-90 (degrees)

xi = [];
yi = [];
R = []; % empty vector for R at each angle (theta)

maxRn = L; % initialize variable for max value of R for each iteration

for i=1:N
    minRn = L/4;
    rn = (maxRn-minRn).*rand(1,1) + minRn;
    xn = rn*cosd(theta(i));
    yn = rn*sind(theta(i));
    xi = [xi; xn];
    yi = [yi; yn];
    R = [R; rn];
    maxRn = sqrt(xn^2 + yn^2);
end

len = size(R);
len = len(1);

TestArea = DidoArea(R); % calculate area for guess values
TestPerimeter = DidoPerimeter(R);

%% optimize the vector R to maximize the area for a given length
A = [];
b = [];
Aeq = [];
beq = [];
lb = zeros([1 len]); % set lower bound of r to be 0
ub = [];
for i=1:(len-1)
    ub = [ub;L];
end
ub = [ub;L]; % set upper bound of r at full string length L

nonlcon = @constraints;
options = optimoptions('fmincon','Algorithm','interior-point','Display','iter');
X = fmincon(@DidoArea,R,A,b,Aeq,beq,lb,ub,nonlcon,options) % optimal solution (X) for R

Xi = [];
Yi = [];
for i=1:N
    % collect X and Y cartesian coordinates from optimal R for plotting 
    Xn = X(i)*cosd(theta(i));
    Yn = X(i)*sind(theta(i));
    Xi = [Xi; Xn];
    Yi = [Yi; Yn];
end

Area = DidoArea(X)
Perimeter = DidoPerimeter(X)

%% plot the initial guess
f1 = figure(1);
clf
subplot(1,2,1)
grid on 
hold on
plot(xi,yi)

for i=1:(N)
    hold on
    scatter(xi(i),yi(i),25,'filled')
end

xlabel('X')
ylabel('Y')
str = sprintf('Dido Problem Initial Guess: n = %d', N);
title(str)

H1=area(xi,yi,'FaceColor',[1 0.85 0]);
hold on

Areastr = sprintf('Area = %d', TestArea);
Lenstr = sprintf('Length = %d', L);
Perimstr = sprintf('Perimeter = %d', TestPerimeter);
text(1,1,Areastr,'HorizontalAlignment','left')
text(1,0.75,Perimstr,'HorizontalAlignment','left')
text(1,0.5,Lenstr,'HorizontalAlignment','left')


%% plot the optimized solution
%f2 = figure(2);
%clf
subplot(1,2,2)
grid on 
hold on
plot(Xi,Yi)
for i=1:(N)
    hold on
    scatter(Xi(i),Yi(i),25,'filled')
end

xlabel('X')
ylabel('Y')
str = sprintf('Dido Problem Optimized Solution: n = %d', N);
title(str)

H1=area(Xi,Yi,'FaceColor',[1 0.85 0]);
hold on

Areastr = sprintf('Area = %d', Area);
Lenstr = sprintf('Length = %d', L);
Perimstr = sprintf('Perimeter = %d', Perimeter);

text(1,1,Areastr,'HorizontalAlignment','left')
text(1,0.75,Perimstr,'HorizontalAlignment','left')
text(1,0.5,Lenstr,'HorizontalAlignment','left')

%% save figures and workspace variables
path = "PhD/OptimalControl/"; % specify directory to save files in
params = "L" + L + "N" + N;

pic1 = path + params + 'TestDidoV3.png';
pic2 = path + params + 'DidoSolutionV3.png';

%saveas(f1,pic1);
saveas(f1,pic2);

%fig1 = path + params + 'TestDidoV3.fig';
fig2 = path + params + 'DidoSolutionV3.fig';

%saveas(f1,fig1);
saveas(f1,fig2);

Rfile = path + params + 'Rguess.mat';
Xfile = path + params + 'Ropt.mat';

save(Rfile,'R');
save(Xfile,'X');

%% functions
function Area = DidoArea(R)
% given a vector R (radius length for every angle theta),
% compute the Area (sum of all X-Y discretized trapezoidal area segments)
len = size(R);
len = len(1);
theta = linspace(0,90,len);
xi = [];
yi = [];
for i=1:len
    rn = R(i);
    th = theta(i);
    xn = rn*cosd(th);
    yn = rn*sind(th);
    xi = [xi; xn];
    yi = [yi; yn];
end
Area = trapz(xi,yi);
end

function Perimeter = DidoPerimeter(R)
% given a vector R (radius length for every angle theta),
% compute the perimeter (sum of all X-Y discretized line segments)
len = size(R);
len = len(1);
theta = linspace(0,90,len);
xi = [];
yi = [];
for i=1:len
    xn = R(i)*cosd(theta(i));
    yn = R(i)*sind(theta(i));
    xi = [xi; xn];
    yi = [yi; yn];
end
d = diff([xi(:) yi(:)]);
Perimeter = sum(sqrt(sum(d.*d,2)));
end

function [c1,ceq] = constraints(R)
% specify equality constraint P - L = 0, where L = 10
Perimeter = DidoPerimeter(R);
c1 = [];
ceq = Perimeter - 10;
end
