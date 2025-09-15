classdef Isopropanol
    properties(Constant)
        % mu = 2.07*10^-3; %[Pa*s or (N/m^2)*s or (kg*m*s/s^2)/(m^2) or kg/(m*s))]
        % rho = 785; %[kg/m^3]
        % nu = B2KConstants.Isopropanol.mu ./ B2KConstants.Isopropanol.rho; %[N*m*s/kg or m^2/s]

        mu = DimVar(UC(2.07*10^-3,1*10^-5),'Pa-s'); %Dynamic viscosity
        rho = DimVar(UC(785,1),'kg/m^3'); %Density
        nu = B2KConstants.Isopropanol.mu ./ B2KConstants.Isopropanol.rho; %Kinematic viscosity
        epsilon_r = UC(19.92,0.05); %Dielectric constant or relative permittivity
    end
end