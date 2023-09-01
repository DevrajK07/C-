%mypanelgen function to find x and z arrays of points on an aerofoil based
%off of 3 inputs from the user; the 4-digit NACA aerofoil code, the angle
%of attack, and the number of panels. In general, it is recommended to use
%150+ panels to generate a smooth airfoil and that level of discretisation
%is especially important if this function is to be used in accordance with
%the PanelStrength() function as well as the main script. 
function [xp, zp] = mypanelgen(NACAcode,AoA,n)
%Extract features from the NACAcode whilst simultaneiously converting it
%to 'double' type
m = str2double(NACAcode(1))/100;
p = str2double(NACAcode(2))/10;
t = str2double(NACAcode(3:4))/100;

%Number of panels is converted to even if odd input given
if mod(n,2) == 1
    n = n + 1;
end
%Number of points is the number of panels plus 1 
ngp = n+1;
x = ones(1,ceil(ngp/2));
%how finely to discretise distance from leading edge to trailing edge and
%the x coordinate points without taking into consideration orthogonality
%with the camber line
for i = 1:ceil(ngp/2)
    x(i) = 1 - 0.5*(1-cos(2*pi*((i-1)/n)));
end
l = length(x);
%Camber function
yc = ones(1,l);
%Camber function to generate yc which is dependednt on x location before or
%after p
%Equation used is from handout
for i=1:l
    if (x(i) < p)
        yc(i) = (m/p^2)*(2*p*x(i) - (x(i))^2);

    else 
        yc(i) = (m/((1-p)^2))*((1-2*p)+2*p*x(i) - (x(i))^2);

    end
end

%Thickness function
c1 = 0.2969;
c2 = -0.1260;
c3 = -0.3516;
c4 = 0.2843;
%c5 = -0.1015; This value of c5 doesn't allow a perfect connection at the
%end point and after research, using c5=-0.1036 is a much nicer solution
c5 = -0.1036;
yt = ones(1,l);
%Generating thickness points from formula in handout
for j=1:l    
    yt(j) = 5*t*(c1*sqrt(x(j)) + c2*x(j) + c3*(x(j))^2 + c4*(x(j))^3 + c5*((x(j)^4)));
end

%Coordinates
%Pre-allocating derivative and theta arrays for speed
dycdx = ones(1,l);
theta = ones(1,l);
%For loop to generate theta values to be used to obtain final x and z
%arrays based off of orthogonality to camber line 
for k=1:l
    if x(k) < p
        dycdx(k) = (2*m/(p^2))*(p-x(k));

    else
        dycdx(k) = (2*m/((1-p)^2))*(p-(x(k)));
    
    end
    theta(k) = atan(dycdx(k));
        
end
%Generates upper and lower panel points based off of orhtogonality to
%camber line
xu = ones(1,l);
xl = ones(1,l);
zu = ones(1,l);
zl = ones(1,l);
for a=1:l
    xu(a) = x(a) - yt(a)*sin((theta(a)));
    xl(a) = x(a) + yt(a)*sin((theta(a)));
    zu(a) = yc(a) + yt(a)*cos((theta(a)));
    zl(a) = yc(a) - yt(a)*cos((theta(a)));
end
xl = xl(1:end-1);
zl = zl(1:end-1);

%Incorporating the Angle of Attack and appending wake panel onto end of
%aerofoil
mag = 999999;
xend2 = 999999;
zend2 = mag*tand(AoA)*-1;
xp = [xl(1:floor(ngp/2)) flip(xu(1:ceil(ngp/2))) 1 xend2];
zp = [zl(1:floor(ngp/2)) flip(zu(1:ceil(ngp/2))) 0 zend2];
% xp = [xl(1,1:n/2),xu(1,n/2+1:end),xend2];
% zp = [zl(1,1:n/2),zu(1,n/2+1:end),zend2];

%Code for plotting if using this function on its own, otherwise commented
%out to prevent an extra figure being generated.
% f1 = figure(1);
% hold on
% grid on
% axis([-0.2 1.2 -0.7 0.7])
% plot(xp(1:end-2),zp(1:end-2),'r-','LineWidth',2);
% hold off
xp = xp;
zp= zp;