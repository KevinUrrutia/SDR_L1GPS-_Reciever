clc; clear; close all; 

load("nav_sol.mat");

truella = [40.6899, -104.9561, 2.6776e+04]';

X_ECEF = [navSol.X; navSol.Y; navSol.Z];
%calculate displacement magnitude vector over time 
dispVecMag = sqrt(sum(diff(X_ECEF(:,3:end-1),[],2).^2));
figure; plot(dispVecMag)
title('Displacement Vector Magnitude')
xlabel('Solution count')
ylabel('meters')
grid on

% Calculate and plot error magnitude for each solution
trueECEF = LLA_to_ECEF(truella);
recalc_ECEF = zeros(3, size(X_ECEF, 2));
for ii = 1:size(recalc_ECEF, 2)
    recalc_ECEF(:,ii) = LLA_to_ECEF([-navSol.latitude(ii); navSol.longitude(ii); -navSol.height(ii)]);
end

errorECEF = sqrt(sum((recalc_ECEF(1:2,2:end-1) - trueECEF(1:2)).^2));
figure; plot(errorECEF)
title('ECEF Error')
xlabel('Solution count')
ylabel('meters')
grid on


% Plot clock estimate
figure; plot(navSol.bias_t(1,3:end-1));
title('Clock Estimate')
xlabel('Solution count')
ylabel('meters')
grid on


function r_ECEF = LLA_to_ECEF(r)
    
    lat = r(1);
    lon = r(2);
    alt = r(3);
   
    a = 6378137;     % semi-major axis of the earth in meters (WGS 84)
    f = 1/298.257223563; % flattening factor of wgs84 (oblate) spheriod of Earth 
    e_sq = (2*f - f^2); % eccentricity squared
    Rn = a./sqrt(1-e_sq*sind(lat).^2);
    h = alt;
    x = (Rn+h).*cosd(lat).*cosd(lon);
    y = (Rn+h).*cosd(lat).*sind(lon);
    z = ((1-e_sq)*Rn + h).*sind(lat);

    
    r_ECEF = [x;y;z];
   
end