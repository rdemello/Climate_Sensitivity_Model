clc;

% Read in data from a NetCDF file (name.nc)
  forcingFile='inputForcing.nc';
  %display the variables and metadata stored in the NetCDF file
  ncdisp(forcingFile); 
  %read in the variable called Aerosol_cloud_lifetime
  aclf = ncread(forcingFile,'Aerosol_cloud_lifetime');
  Year = ncread(forcingFile,'Year');
% Note that the Skeie et al. Radiative Forcing is with respect to 1750, but we start in 1850.

%% Declare All Variables

F_total = load('RCP8.5_Forcing.txt');

ECS = 1;            %climate sensitivity 
a = 3.74/ECS;       %defined alpha for climate feedback parameter

density = 1027;     %"p" density of water in kg/m3
c_p = 4218;         %specific heat capacity of water J/kg/K

k = 1;              %vertical diffusivity cm^2/s
h_u = 100;          %upper height m
h_d = 900;          %lower height m

C_u = density*c_p*h_u; %thermal interia upper J/(m^2 K^1 s^1/2)
C_d = density*c_p*h_d; %thermal interia deep J/(m^2 K^1 s^1/2)

g = (2*k*c_p*density)/(h_u+h_d); %heat diffusion

num_years = 161;
timestep = 1;
F_diff = 1;         %variable from code snippet

T_d = zeros(num_years,1); %empty array for deep temp
T_u = zeros(num_years,1); %empty array for upper temp

%% Run Loop

for i = 1:num_years-1;
    T_u(i+1,1) = (F_total(i,2) - (a*(T_u(i,1))) - (g*(T_u(i,1) - T_d(i,1))))/C_u; 
    T_d(i+1,1) = (g*(T_u(i,1) - T_d(i,1)))/C_d;
end
  
%% Plot Graph

  figure(1);
  plot(Year,T_u,'r','LineWidth',3);
  title('Global Mean Temperatures','FontWeight','bold','FontSize',14);
  ylabel('Temp Variation from Total Forcing on Upper Ocean','FontSize',12);
  xlabel('Year','FontWeight','bold','FontSize',12);
  hold all
  plot(Year,T_d,'r','LineWidth',3);
  title('Global Mean Temperatures','FontWeight','bold','FontSize',14);
  ylabel('Temp Variation from Total Forcing on Lower Ocean','FontSize',12);
  xlabel('Year','FontWeight','bold','FontSize',12);

%% Old Equations

% F_atmos = F_total - a*T_u;
% F_mixed = F_atmos - F_diff;
% dT_u = (F_mixed*time)/c_p; 
% F_diff = k * dT/dz;

%% Example Code

% % An example of forward timestepping
%   Td = zeros(num_years,1); %create a new empty array 
% % Iterate the model
%   for i = 1:num_years-1;
%    % determine annual energy into each box (Joules)
%    deep_energy = timestep * diffusive_flux;
%    % determine new temperature (old temp + energy as temperature change)
%    Td(i+1,1)=Td(i,1)+deep_energy/(dens*cp*hd);
%   end

% % Compare this to the HadCRUT4 global (NH+SH)/2 data
% % Downloaded as txt file from http://www.metoffice.gov.uk/hadobs/hadcrut4/data/current/download.html
% % NOTE: expressed as anomalies w.r.t 1961-1990
%   HadCRUT4_file=fopen('HadCRUT.filename.txt'); 
%   HadCRUT4_all=textscan(HadCRUT4_file,'%*d %f %*f %*f %*f %*f %*f %*f %*f %*f %f %f');
%   T_median = HadCRUT4_all{1,1}(1:161,1);
%   T_95percent = HadCRUT4_all{1,2}(1:161,1);
%   T_5percent = HadCRUT4_all{1,3}(1:161,1);


% a plot of your results compared to data
%   figure(1);
%   plot(Year,T_u,'r','LineWidth',3);
%   title('Global Mean Temperatures','FontWeight','bold','FontSize',14);
%   ylabel('Temperature Anomaly (oC w.r.t. 1961-1990)','FontSize',12);
%   xlabel('Year','FontWeight','bold','FontSize',12);
%   hold all;
%   plot(Year,T_median,'g','LineWidth',3);
%   plot(Year,T_5percent,'g','LineWidth',1);
%   plot(Year,T_95percent,'g','LineWidth',1);
%   legend('Modelled','Observed','Location','NorthWest');
%   hold off;

