clc;
  
%% Define Variables

forcingFile='inputForcing.nc';
%ncdisp(forcingFile); %display all variables for the forcing file
  
% set variables for all forcings
F_aerosolCloudLife = ncread(forcingFile,'Aerosol_cloud_lifetime');
F_ghg = ncread(forcingFile,'Greenhouse_gases');
F_stratOzone = ncread(forcingFile,'Strato_ozone');
F_tropOzone = ncread(forcingFile,'Tropo_ozone');
F_stratoWater = ncread(forcingFile,'Strato_water');
F_aerosolDirect = ncread(forcingFile,'Aerosol_direct');
F_aerosolCloudAlbedo = ncread(forcingFile,'Aerosol_cloud_albedo');
F_landUse = ncread(forcingFile,'Land_use');
F_snowAlbedo = ncread(forcingFile,'Snow_albedo');
F_solar = ncread(forcingFile,'Solar');
F_volcanic = ncread(forcingFile,'Volcanic');
Year = ncread(forcingFile,'Year');

num_years = 161;
secInYear = 60*60*24*365;
timestep =  secInYear;

F_total = zeros(num_years,1);

F_RCPForcing = load('RCP8.5_Forcing.txt');

F_selected = F_RCPForcing;

ECS = 1.5;          %climate sensitivity 
a = 3.74/ECS;       %defined alpha for climate feedback parameter

density = 1027;     %"p" density of water in kg/m3
c_p = 4218;         %specific heat capacity of water J/kg/K

k = 0.0001;         %vertical diffusivity m^2/s
h_u = 100;         	%upper height m
h_d = 900;         	%lower height m

C_u = density*c_p*h_u; %thermal interia upper J/(m^2 K^1 s^1/2)
C_d = density*c_p*h_d; %thermal interia deep J/(m^2 K^1 s^1/2)

g = (2*k*c_p*density)/(h_u+h_d); %heat diffusion m2 * 1/s * J * 1/kg * 1/K * kg * 1/m3 * 1/m = J/m^2 * K * s

T_d = zeros(num_years,1); %empty array for deep temp
T_u = zeros(num_years,1); %empty array for upper temp

%% Run Loop

for i = 1:num_years-1;
    upper_energy = timestep * (F_selected(i,2) - (a*(T_u(i))) - (g*(T_u(i) - T_d(i))));  %solved from equation
    T_u(i+1) = T_u(i) + upper_energy/C_u; 
    
    deep_energy = timestep * g*(T_u(i) - T_d(i));
    T_d(i+1) = T_d(i) + deep_energy/C_d;
end

%% Plot Graph

figure(1);
plot(Year,T_u,'c','LineWidth',2);
title('Temperature Variation Based on RCP8.5 Forcing','FontWeight','bold','FontSize',14);
ylabel('Temperature Variation','FontSize',12);
xlabel('Year','FontWeight','bold','FontSize',12);
hold all
plot(Year,T_d,'b','LineWidth',2);
legend('Upper Ocean','Deep Ocean');

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

