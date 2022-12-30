%% Determines inclination range based on a set of limits

%{

% Part a:
phiWallops = 37.940194; % degrees
zetaMinA = 90;
zetaMaxA = 160;

% Part b:
phiVandenberg = 34.7483;
zetaMinB = 153;
zetaMaxB = 240;

% Part c:
phiKapustin = 48.4;
zetaMinC = 350;
zetaMaxC = 450;

% Part d:
phiSriharikota = 13.7;
zetaMinD = 100;
zetaMaxD = 290;


%}

function rangeOut = get_inclination_range(phiInput, azMin, azMax)

    cInput = cosd(phiInput);
    zetaRange = azMin:0.01:azMax;
    rangeOut.iMin = 10000;
    rangeOut.iMax = -1;
    rangeOut.zMin = 10000;
    rangeOut.zMax = -1;
    auxiliaryNow = 0;
    %flag = 0;
    
    for n = 1:length(zetaRange)
        if (n >= 360)
            zetaRange(n) = wrapTo360(zetaRange(n));
        end
        
        iNow = acosd(sind(zetaRange(n))*cInput);
        auxiliaryNow = acosd(cosd(zetaRange(n)) / sind(iNow));
    
        if ((iNow > rangeOut.iMax) && (auxiliaryNow >= 0))
            rangeOut.iMax = iNow;
            rangeOut.zMax = zetaRange(n);
        end
        if ((iNow < rangeOut.iMin) && (auxiliaryNow >= 0))
            rangeOut.iMin = iNow;
            rangeOut.zMin = zetaRange(n);
        end
    end

    % Problem 1:
    fprintf('  Inclination Ranges\n\n');
    fprintf('        iMin = %0.3f deg., azMin = %0.3f deg.\n', ...
        rangeOut.iMin, rangeOut.zMin);
    fprintf('        iMax = %0.3f deg., azMax = %0.3f deg.\n\n', ...
        rangeOut.iMax, rangeOut.zMax);

end


end

    
