classdef Water
    properties(Constant)
        % mu = 8.89*10^-4; %[Pa*s or (N/m^2)*s or (kg*m*s/s^2)/(m^2) or kg/(m*s))]
        % rho = 997; %[kg/m^3]
        % nu = B2KConstants.Water.mu ./ B2KConstants.Water.rho; %[N*m*s/kg or m^2/s]

        mu = DimVar(UC(8.89*10^-4,1*10^-6),'Pa-s'); %Dynamic viscosity
        rho = DimVar(UC(997,1),'kg/m^3'); %Density
        nu = B2KConstants.Water.mu ./ B2KConstants.Water.rho; %Kinematic viscosity
        epsilon_r = UC(78.39,0.05); %Dielectric constant or relative permittivity
    end
end