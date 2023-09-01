%PanelStrength function to compute and display panel strengths and a crude
%estimation of the coefficient of lift generated from four inputs;
%The 4 digit NACA aerofoil, freestream velocity, angle of attack of aerofoil
%and also number of panels. In general a good guideline for fairly accurate
%estimations of the CL value is to use a panel number of 150+. Within this
%function, mypanelgen() has also been called to generate x and z arrays of
%points on the aerofoil for calcualtions of midpoints and the Beta
%matrix.
function [mu, C_L, xp, zp] = PanelStrength(NACAcode,Uinf,AoA,n)
tic
%Uses mypanelgen() function to output the x and z arrays
[xp,zp]  = mypanelgen(NACAcode,AoA,n);
%Appends large magnitude wake panel to end of aerofoil and generates final
%x and z array
mag = 999999;
xend2 = mag;
zend2 = mag*tand(AoA)*-1;
xp = [xp(1:end-2) xend2];
zp = [zp(1:end-2) zend2];
A = zeros(n+1);
B = zeros(n+1,1);
%Using matrix operators to compute the beta matrix and in turn the B matrix
Beta = atan((zp(2:end) - zp(1:end-1))./((xp(2:end) - xp(1:end-1))));
B = -Uinf*sin(AoA*pi/180 - Beta);

for i=1:n
    %Loop used to generate all the midpoints of the panels as this is where
    %the velocities are evaluated
    xmid(i) = 0.5*(xp(i)+xp(i+1));
    zmid(i) = 0.5*(zp(i)+zp(i+1));
    for j = 1:n+1
        %Using cdoubletfin() function (MATLAB had an issue with me saving
        %it as just cdoublet(), to generate velocities and there the matrix
        %A to be used for the final mu calculation
        [uip,vip] = cdoubletfin([xmid(i),zmid(i)],[xp(j),zp(j)],[xp(j+1),zp(j+1)]);
        A(i,j) = vip*cos(Beta(i)) - uip*sin(Beta(i));
    end
end
% A(n+1,:) = zeros(n+1,n+1);
%Appends kutta condition onto A and B matrix to generate correct mu values
A(n+1,1) = 1;
A(n+1,n) = -1;
A(n+1,n+1) = 1;
B(n+1) = 0;
%mu calculated using matrix division
mu = A\B';
%Estimation of lift coefficient from mu values and given freeestream
%velocity
C_L = -2*(mu(end)/(Uinf));


mu = mu;
C_L = C_L;
end