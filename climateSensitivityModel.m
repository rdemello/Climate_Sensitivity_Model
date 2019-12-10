clc;
  
%% Define Variables

med_depth = 1300;
arctic_depth = 838;
pacific_depth = 4080;
atlantic_depth = 3139;

forcingFile='inputForcing.nc';
Year = ncread(forcingFile,'Year');

F_selected = forcing();
F_noVolcano = noVolcano()

% HadCRUT4_file=fopen('HadCRUT.txt'); 
% HadCRUT4_all=textscan(HadCRUT4_file,'%*d %f %*f %*f %*f %*f %*f %*f %*f %*f %f %f');
% HadCRUT4s_file=fopen('HadCRUTs.txt'); 
% HadCRUT4s_all=textscan(HadCRUT4s_file,'%*d %f %*f %*f %*f %*f %*f %*f %*f %*f %f %f');
% T_median = HadCRUT4_all{1,1}(1:161,1);
% T_medians = HadCRUT4s_all{1,1}(1:161,1);
% T_adjustedMedian = zeros(161,1);
% 
% for i=1:161
%    T_adjustedMedian(i) = (T_medians(i)-T_medians(1));
% end

%% Run Loop

T_med = seaTempCalc(F_selected,med_depth);
T_pacific = seaTempCalc(F_selected,pacific_depth);
T_arctic = seaTempCalc(F_selected,arctic_depth);
T_atlantic = seaTempCalc(F_selected,atlantic_depth);

T_upper_all = [T_med(:,1),T_arctic(:,1),T_atlantic(:,1),T_pacific(:,1)];
T_deep_all = [T_med(:,2),T_arctic(:,2),T_atlantic(:,2),T_pacific(:,2)];

T_med_nv = seaTempCalc(F_noVolcano,med_depth);
T_pacific_nv = seaTempCalc(F_noVolcano,pacific_depth);
T_arctic_nv = seaTempCalc(F_noVolcano,arctic_depth);
T_atlantic_nv = seaTempCalc(F_noVolcano,atlantic_depth);

T_upper_all_nv = [T_med_nv(:,1),T_arctic_nv(:,1),T_atlantic_nv(:,1),T_pacific_nv(:,1)];
T_deep_all_nv = [T_med_nv(:,2),T_arctic_nv(:,2),T_atlantic_nv(:,2),T_pacific_nv(:,2)];

%% Plot Graph

figure(1);
plot(Year,T_upper_all,'LineWidth',2);
title('Temperature Variation on Upper Ocean Based at different Ocean depths','FontWeight','bold','FontSize',14);
ylabel('Temperature Variation','FontSize',12);
xlabel('Year','FontWeight','bold','FontSize',12);
%legend('AerosolCloudLife', 'ghg', 'stratOzone', 'tropOzone', 'stratoWater', 'aerosolDirect', 'aerosolCloudAlbedo', 'landUse', 'snowAlbedo', 'solar', 'volcanic');
legend('Med', 'Arctic', 'Atlantic', 'Pacific')
hold all;
%plot(Year,T_adjustedMedian,'g','LineWidth',1);
%plot(Year,T_medians,'g','LineWidth',1);
%plot(Year,T_u(:,1:11),'LineWidth',0.1);
%plot(Year,T_d(:,5),'LineWidth',2);

figure(2);
plot(Year,T_deep_all,'LineWidth',2);
title('Temperature Variation on Deep Ocean Based at different Ocean depths','FontWeight','bold','FontSize',14);
ylabel('Temperature Variation','FontSize',12);
xlabel('Year','FontWeight','bold','FontSize',12);
legend('Med', 'Arctic', 'Atlantic', 'Pacific')


figure(3);
plot(Year,F_selected,'LineWidth',2);
title('Input Forcings','FontWeight','bold','FontSize',14);
legend('AerosolCloudLife', 'ghg', 'stratOzone', 'tropOzone', 'stratoWater', 'aerosolDirect', 'aerosolCloudAlbedo', 'landUse', 'snowAlbedo', 'solar', 'volcanic');

figure(4);
plot(Year,T_upper_all_nv,'LineWidth',2);
title('Temperature Variation on Upper Ocean Based with different Ocean depths without Volcanic Forcing','FontWeight','bold','FontSize',14);
ylabel('Temperature Variation','FontSize',12);
xlabel('Year','FontWeight','bold','FontSize',12);
%legend('AerosolCloudLife', 'ghg', 'stratOzone', 'tropOzone', 'stratoWater', 'aerosolDirect', 'aerosolCloudAlbedo', 'landUse', 'snowAlbedo', 'solar', 'volcanic');
legend('Med', 'Arctic', 'Atlantic', 'Pacific')


figure(5);
plot(Year,T_deep_all_nv,'LineWidth',2);
title('Temperature Variation on Deep Ocean Based with different Ocean depths without Volcanic Forcing','FontWeight','bold','FontSize',14);
ylabel('Temperature Variation','FontSize',12);
xlabel('Year','FontWeight','bold','FontSize',12);
legend('Med', 'Arctic', 'Atlantic', 'Pacific')

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

