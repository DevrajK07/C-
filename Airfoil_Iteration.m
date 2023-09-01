% Code to iterate through all 9999 NACA airfoils to determine which one has
% the highest CL

%MATLAB Aerofoil Script
clear
clc
%Asks user to input NACA code, then the if statement sets initialise to 1 
%only if the 2412 aerofoil is submitted (which then carries out additional
%operations)
NACAcode = input('What NACA 4-series aerofoil do you want? ','s');
code = str2double(NACAcode);
if code == 2412
    initialise = 1;
else
    initialise = 0;
end
%Generates AoA array to match the ones given in X-foil data to plot against
AoAarray = linspace(-17.5,19.25,145);
AoAarray2 = linspace(0,10,50);
%Generates array of number of panels
nArray = [50,100,200];
%Uses switch-case statement to carry out aditional operations for the 2412
%case
switch initialise
    case 1
        %Generates different colour lines for different Panel lines on
        %CL-AoA plot
        colour = 'rgm';
        fig1 = figure;
        %For loop for number of panels wanted to be evaluated
        for i=1:length(nArray)
            muarray = [];
            C_Larray = [];
            %Nested j loop to use PanelStrength() function for each angle
            %of attack with each length of panel requested
            for j=1:length(AoAarray)
                [mu, C_L] = PanelStrength(NACAcode,15,AoAarray(j),nArray(i));
                mu = mu;
                C_Larray(j) = C_L;
            end
                %Generates plots for AoA against CL for different numbers
                %of panels while holding on to plot them against each other
                plot(AoAarray,C_Larray,colour(i),'LineWidth',1.5)
                hold on
        end
        %Now generating all Cl values for a range of values between 0 and
        %10 degrees angle of attack
        CLvals = ones(length(AoAarray2),1);
        for val=1:length(AoAarray2)
           [mu CLvals(val)] = PanelStrength(NACAcode,15,AoAarray2(val),200);
        end
        disp(CLvals)
        table = readtable("xf-naca2412-il-1000000.txt");
        tbl = table(:, {'alpha','CL'});
        %Adds X-Foil data plot for comparison against N panel values
        plot(tbl.alpha,tbl.CL,'k','LineWidth',1.5);
        title('CL v.s. AoA plot')
        xlabel('AoA')
        ylabel('CL')
        legend('50 points','100 points','200 points','X-FOIL Data','Location','best')
        axis([0 10 0 5])
        hold off
        saveas(fig1,'CLvAoAgraph.png');
        %Now continues with the regular process of generating the flow
        %field around the aerofoil
        AoA = 10;
        n = 200;
        Uinf = 15;
        %Brings x and z points from panelgen() into this script
        [xp,zp]  = mypanelgen(NACAcode,AoA,n);
        %For wake panel large magnitude and appending those x and z points
        %onto this array of x and z points
        mag = 999999;
        xend2 = mag;
        zend2 = mag*tand(AoA)*-1;
        xp = [xp(1:end-2) xend2];
        zp = [zp(1:end-2) zend2];
        A = zeros(n+1);
        B = zeros(n+1,1);
        %Computing Beta and B with array operators
        Beta = atan((zp(2:end) - zp(1:end-1))./((xp(2:end) - xp(1:end-1))));
        B = -Uinf*sin(AoA*pi/180 - Beta);
        
        %Defines midpoints of panels as this is where velocities are
        %evaluated
        for i=1:n
            xmid(i) = 0.5*(xp(i)+xp(i+1));
            zmid(i) = 0.5*(zp(i)+zp(i+1));
            for j = 1:n+1
                %Now using cdoubletfin() function to determine velocities
                %at those points
                [uip,vip] = cdoubletfin([xmid(i),zmid(i)],[xp(j),zp(j)],[xp(j+1),zp(j+1)]);
                A(i,j) = vip*cos(Beta(i)) - uip*sin(Beta(i));
            end
        end
        %Implements kutta condition into the matrix A
        A(n+1,1) = 1;
        A(n+1,n) = -1;
        A(n+1,n+1) = 1;
        B(n+1) = 0;
        mu = A\B';
        C_L = -2*(mu(end)/(Uinf));
        
        %resolution determines mesh size and how accurate my streamlines will be as
        %it evaluates at more points along each of the streamlines
        resolution = 150;
        xgrid = linspace(-0.2,1.2,resolution);
        ygrid = linspace(-0.7,0.7,resolution);
        U(1:resolution,1:resolution) = 0;
        V(1:resolution,1:resolution) = 0;  
        %Nested for loop to generate streamlines to follow curve of aerofoil
        for i=1:resolution
            for j=1:resolution
                %Adding the constants for the final matrix beforehand for
                %better efficiency
                U(i,j) = Uinf*cosd(AoA);
                V(i,j) = Uinf*sind(AoA);
                for k = 1:n+1
                    %Generating streamlines at each panel
                    [u,v] = cdoubletfin([xgrid(i), ygrid(j)],[xp(k),zp(k)],[xp(k+1),zp(k+1)]);
                    U(i,j) = U(i,j) + mu(k)*u;
                    V(i,j) = V(i,j) + mu(k)*v;
                end
            end
        end
        