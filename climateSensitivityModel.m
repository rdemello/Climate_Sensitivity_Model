clc;
  
%% Define Variables

arctic_depth = 1038;
med_depth = 1500;


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

HadCRUT4_file=fopen('HadCRUT.txt'); 
HadCRUT4_all=textscan(HadCRUT4_file,'%*d %f %*f %*f %*f %*f %*f %*f %*f %*f %f %f');
HadCRUT4s_file=fopen('HadCRUTs.txt'); 
HadCRUT4s_all=textscan(HadCRUT4s_file,'%*d %f %*f %*f %*f %*f %*f %*f %*f %*f %f %f');
T_median = HadCRUT4_all{1,1}(1:161,1);
T_medians = HadCRUT4s_all{1,1}(1:161,1);
T_adjustedMedian = zeros(161,1);

for i=1:161
   T_adjustedMedian(i) = (T_medians(i)-T_medians(1));
end

num_years = 161;
secInYear = 60*60*24*365;
timestep =  secInYear;

Individual_Forcing = [ F_aerosolCloudLife F_ghg F_stratOzone F_tropOzone F_stratoWater F_aerosolDirect F_aerosolCloudAlbedo F_landUse F_snowAlbedo F_solar F_volcanic];
F_total = sum(Individual_Forcing,2);
F_allForcing = [Individual_Forcing F_total];
F_RCPForcing = load('RCP8.5_Forcing.txt');

F_selected = F_allForcing;
h_d = 900;
T_d = zeros(num_years,5); %empty array for deep temp
T_u = zeros(num_years,5); %empty array for upper temp


%% Run Loop
%     for j = 1:12
%         for i = 1:num_years-1
            %seaTempCalc(i,j,T_u,T_d,timestep,F_selected)
            %upper_energy = timestep * (F_selected(i,j) - (a*(T_u(i,j))) - (g*(T_u(i,j) - T_d(i,j))));  %solved from equation
            %T_u(i+1,j) = T_u(i,j) + upper_energy/C_u; 

            %deep_energy = timestep * g*(T_u(i,j) - T_d(i,j));
            %T_d(i+1,j) = T_d(i,j) + deep_energy/C_d;
%         end
%     end

seaTempCalc(T_u,T_d,num_years,timestep,F_selected,h_d);

%% Plot Graph

figure(1);
plot(Year,T_u(:,5),'LineWidth',2);
title('Temperature Variation on Upper Ocean Based on Input Forcing','FontWeight','bold','FontSize',14);
ylabel('Temperature Variation','FontSize',12);
xlabel('Year','FontWeight','bold','FontSize',12);
%legend('AerosolCloudLife', 'ghg', 'stratOzone', 'tropOzone', 'stratoWater', 'aerosolDirect', 'aerosolCloudAlbedo', 'landUse', 'snowAlbedo', 'solar', 'volcanic');
hold all;
plot(Year,T_adjustedMedian,'g','LineWidth',1);
%plot(Year,T_medians,'g','LineWidth',1);
%plot(Year,T_u(:,1:11),'LineWidth',0.1);
%plot(Year,T_d(:,12),'LineWidth',2);

% % figure(2);
% plot(Year,F_allForcing,'LineWidth',2);

% figure(2);
% plot(Year,T_d,'LineWidth',2);
% title('Temperature Variation on Deep Ocean Based on Input Forcing','FontWeight','bold','FontSize',14);
% ylabel('Temperature Variation','FontSize',12);
% xlabel('Year','FontWeight','bold','FontSize',12);
%legend('AerosolCloudLife', 'ghg', 'stratOzone', 'tropOzone', 'stratoWater', 'aerosolDirect', 'aerosolCloudAlbedo', 'landUse', 'snowAlbedo', 'solar', 'volcanic');


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

