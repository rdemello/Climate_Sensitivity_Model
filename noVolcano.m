function [forcingTotal] = noVolcano()
    forcingFile='inputForcing.nc';
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
    F_volcano = zeros(161,1);

    Individual_Forcing = [ F_aerosolCloudLife F_ghg F_stratOzone F_tropOzone F_stratoWater F_aerosolDirect F_aerosolCloudAlbedo F_landUse F_snowAlbedo F_solar F_volcano];
    F_total = sum(Individual_Forcing,2);
    F_allForcing = [Individual_Forcing F_total];

    forcingTotal = F_allForcing;

end


