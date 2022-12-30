function rhoENZ = get_enz_from_azelrn(azIn_deg, elIn_deg, rho)

    cAz = cosd(azIn_deg);
    sAz = sind(azIn_deg);
    cEl = cosd(elIn_deg);
    sEl = sind(elIn_deg);

    rhoENZ = rho * [cEl*cAz;
                    cEl*sAz;
                    sEl];

end