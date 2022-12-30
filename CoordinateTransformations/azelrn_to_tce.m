function rhoTCE = azelrn_to_tce(az, el, rn, timeUTC, llaInput)

    rho_ENZ = get_enz_from_azelrn;
    lstInput_deg = siderealtime(timeUTC, llaInput);
    
    % Convert AzElRange to EastNorthZenith (sph. pol.):
    thetaLST = -(90 + lstInput_deg(1));
    phiInput = llaInput.latitude_deg - 90;
    cLST = cosd(thetaLST);
    sLST = sind(thetaLST);
    cPHi = cosd(phiInput);
    sPhi = sind(phiInput);
    
    ROTx = [1  0   0;
            0  cLST sLST;
            0  -sLST  cLST];
    
    ROTz = [cPHi sPhi 0;
            -sPhi  cPHi 0;
            0   0  1];
    
    rho_TCE = ROTz * ROTx * rho_ENZ;
    R_snag = [R_a_snag.position_km(1) R_a_snag.position_km(2) R_a_snag.position_km(3)]';
    
    % Apply vector addition to get the final ECI vector of the satellite:
    r_Sat = rho_TCE + R_snag;
    
    % Transform to spherical coordinates:
    rho = norm(r_Sat);
    dec_Sat = acosd(sqrt(r_Sat(1)^2 + r_Sat(2)^2) / rho);
    rAsc_Sat = acosd(r_Sat(1) / (rho*cosd(dec_Sat)) );

end