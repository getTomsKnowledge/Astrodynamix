classdef  OrbitConstants
    properties ( Constant = true )
        
        % Earth:
        siderealDay_seconds = 86164.090524;
        oRateEarth_radsec = 7.292115*1e-5;
        R_earth_km = 6378.1; % km
        R_moon_km = 1738.1;
        r_earthSun_km = 1.495978707*1e8;
        L1_earthSun_km = 148.11*1e6;
        % from Earth along Earth-moon axis
        r_earthMoon_km = 384399;
        rCM_earthMoon_km = 4672.62;
        xi_L1_km = 321699.759; % barycentric position of L1
        r1_moonToL1_km = 58027;
        L1_earthMoon_km =  326372; 
        L2_earthMoon_km =  448923;
        L3_earthMoon_km = -381673;
        L4_Coords_km    =  [192200,  332899];
        L5_Coords_km    =  [192200, -332899];
        mu_moon_km = 4.9048695*1e3;
        mu_earth_km = 3.986004418e5; % earth grav parameter, km^3/s^2
        mu_sun_km = 1.32712440042e11; % sun grav param, km^3/s^2
        fineStruct = 0.007297352525693; % Fine-structure constant
        planck = 6.62607015e-34; % Planck's constant
        eFundamental = 1.602176634e-19; % Fundamental charge
        c_ms = 299792458; % lightspeed m/s
        G_Nmkg = 6.67430*1e-11; % Universal gravitational constant

        %{
        mu_0 = 2*fineStruct*planck / ((eFundamental^2)*c);
        epsilon_0 = 1/(mu_0*(c^2));
        %}
    end
end