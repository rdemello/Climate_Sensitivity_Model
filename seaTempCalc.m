function [T_filled] = seaTempCalc(F_selected,h_d,length)

    num_years = 161;
    secInYear = 60*60*24*365;
    timestep =  secInYear;
    j=1;
    
    T_d = zeros(num_years,1); %empty array for deep temp
    T_u = zeros(num_years,1); %empty array for upper temp

    ECS = 3;            %climate sensitivity 
    a = 3.74/ECS;       %defined alpha for climate feedback parameter

    density = 1027;     %"p" density of water in kg/m3
    c_p = 4218;         %specific heat capacity of water J/kg/K

    k = 0.0001;         %vertical diffusivity m^2/s
    h_u = 200;         	%upper height m

    C_u = density*c_p*h_u; %thermal interia upper J/(m^2 K^1 s^1/2)
    C_d = density*c_p*h_d; %thermal interia deep J/(m^2 K^1 s^1/2)

    g = (2*k*c_p*density)/(h_u+h_d); %heat diffusion m2 * 1/s * J * 1/kg * 1/K * kg * 1/m3 * 1/m = J/m^2 * K * s
    
    for i = 1:num_years-1
        upper_energy = timestep * (F_selected(i,12) - (a*(T_u(i,j))) - (g*(T_u(i,j) - T_d(i,j))));  %solved from equation
        T_u(i+1,j) = T_u(i,j) + upper_energy/C_u;

        deep_energy = timestep * g*(T_u(i,j) - T_d(i,j));
        T_d(i+1,j) = T_d(i,j) + deep_energy/C_d;
    end
    
    T_filled=[T_u, T_d];
    
end

