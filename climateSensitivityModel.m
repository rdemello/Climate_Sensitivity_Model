clc;

% Read in data from a NetCDF file (name.nc)
  forcingFile='inputForcing.nc';
  %display the variables and metadata stored in the NetCDF file
  ncdisp(forcingFile); 
  %read in the variable called Aerosol_cloud_lifetime
  aclf = ncread(forcingFile,'Aerosol_cloud_lifetime');    
% Note that the Skeie et al. Radiative Forcing is with respect to 1750,
% but we start in 1850.

F_total = load('RCP8.5_Forcing.txt');

ECS = 1; %climate sensitivity 
a = 3.74/ECS; %define alpha for climate feedback parameter

density = 1027; %p density of water in kg/m3
c_p = 4218; %specific heat capacity of water J/kg/K

k_d = 0.0001; %vertical diffusivity m^2/s
h_u = 100; %upper height m
h_d = 900; %lower height m

C_u = density*c_p*h_u; %thermal interia upper J/(m^2 K^1 s^1/2)
C_d = density*c_p*h_d; %thermal interia deep J/(m^2 K^1 s^1/2)

k = k_d*c_p*density; %heat diffusion

num_years = 161;
timestep = 1;

F_atmos = F_total - a*T_u;
F_mixed = F_atmos - F_diff;
dT_u = (F_mixed*time)/c_p; 
F_diff = k * dT/dz;

% An example of forward timestepping
  Td = zeros(num_years,1); %create a new empty array 
% Iterate the model
  for i = 1:num_years-1;
   % determine annual energy into each box (Joules)
   deep_energy = timestep * F_diff;
   % determine new temperature (old temp + energy as temperature change)
   Td(i+1,1)=Td(i,1)+deep_energy/(density*c_p*h_d);
  end