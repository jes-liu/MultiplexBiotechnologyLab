clear all; close all; clc;

% testing small sample of spillover & getting the equation
dir = 'C:\Users\jesse\OneDrive\Documents\MATLAB\spillover.xlsx';

data = readmatrix(dir);

data = sortrows(data, 1);

x = data(:,1);
y = data(:,2);

y = (y*1366)+1084.703;

y(54)=1366;
x(54)=0;

p = polyfit(x,y,2);

new_x = -8:.001:8;
new_y = polyval(p,new_x);


xfix=[0];
yfix=[1366];
n=4;
xder=[-8 8];
yder=[0 0];

% polyfix(x,y,n,xfix,yfix,xder,dydx)
P = polyfix(new_x,new_y,n,xfix,yfix,xder,yder)
fixed_y = polyval(P,new_x);

figure()
hold on
plot(x, y, 'o')
plot(new_x, fixed_y)
eq = sprintf('y = %fx^4 %fx^3 %fx^2 +%fx +%.0f',P(1), P(2), P(3), P(4), P(5));
text(-6.5, 1400, eq)
xlabel('Distance from bright bead')
ylabel('Intensity')
hold off

%% 
clear all; close all; clc;

folder = 'C:\Users\jesse\OneDrive\Documents\MATLAB\spillover';

filename = dir(folder);
filename = filename(3:12);

bright_beads=[1350, 1450, 1550, 1650, 1750, 1850, 1950, 2250, 2750, 3500];
n=4;
eq_coef=zeros(length(bright_beads), n+1);

for i=1:length(filename)
    dir = fullfile(folder, filename(i).name);
    data = readmatrix(dir);

    % exceuting polyfit
    x = data(:,1);
    y = data(:,2);

    p = polyfit(x,y,4);
    new_x = -8:.001:8;
    new_y = polyval(p,new_x);

    %xfix=[0];
    %yfix=[bright_beads(i)];
    xfix=[];
    yfix=[];
    xder=[-8 8];
    yder=[0 0];
    % polyfix(x,y,n,xfix,yfix,xder,dydx)
    P = polyfix(new_x,new_y,n,xfix,yfix,xder,yder);
    fixed_y = polyval(P,new_x);
    
    figure()
    hold on
    plot(x, y, 'o')
    plot(new_x, fixed_y, 'k-', 'Linewidth', 3)
    ylim([min(fixed_y)-100, 1400])
    eq = sprintf('y = %fx^4 %fx^3 %fx^2 +%fx +%.0f',P(1), P(2), P(3), P(4), P(5));
    text(-6.5, 1350, eq)
    xlabel('Distance from bright bead')
    ylabel('Intensity')
    hold off
    
    eq_coef(i,:)=P;
    
end

% export the equations
writematrix(eq_coef,'parabolic.xlsx')
