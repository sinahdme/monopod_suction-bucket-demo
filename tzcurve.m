function [tz_in_final, tz_out_final] = tzcurve(P_i_out, fi, gamma, do_b, di_b, z, displacement, l_b, n_l, n, fi_crit,delta,D_r,e_max,e_min)
% In this function delta can be either 0 or a specific amount.
    % Validate inputs
    if any([gamma, do_b, l_b] <= 0) || any(size(z) ~= size(fi))
        error('Invalid input values');
    end

    %% fi_crit_model calculation
if  e_max ==0 && e_min ==0
    n_max = 0.4764;
    n_min = 0.2595;
    e_max = n_max/(1-n_max);
    e_min = n_min/(1-n_min);
else
n_max = e_max/(1+e_max);  %porosity ratio
n_min= e_min/(1+e_min);  %porosity ratio
end
 e_c  = (e_max-(e_max-e_min)*D_r);
 n_c = e_c/(1+e_c);

n_c_r = (n_max - n_c) / (n_max - n_min); % equation 15
n_c_model = 0.4764-n_c_r*(0.4764-0.2595);
fi_crit_model =asin( ( pi) / (6 * (1 - n_c_model)))*180/pi;

if fi_crit_model >45
    fi_crit_model  = 90 - (asin( ( pi) / (6 * (1 - n_c_model)))*180/pi); % deg
end
%% 
    % Fit a polynomial to soil properties
    soil_properties = [15, 0.48; 20, 67; 25, 81; 30, 96; 35, 115];
    p = polyfit(soil_properties(:, 1), soil_properties(:, 2), 1);
    f = polyval(p, atand(2/3 * tand(fi))); % Function of friction angle

    % Compute soil stiffness and max allowable settlement
    sigma_p = gamma .* z;
    DRR = D_r*100;
    ocr=max(1,(19.267*DRR.*sigma_p.^(-1.146) -1302.2.*sigma_p.^(-1.172)));
   fi_critical_ocr=atan(tand(fi_crit_model) ./ocr).*180/pi; %deg
KO_OCR=((3/4)^0.5-sind(fi_critical_ocr))./(2*cosd(fi_critical_ocr)-1);
lambda_ocr=0.95;
K0=KO_OCR.*ocr-KO_OCR.*(ocr-1)*lambda_ocr;
%     K0 =1 - sind(fi);
%     if delta==0
%     t_max_in = K0 .* sigma_p .* (2/3 * tand(fi));
%     t_max_out = P_i_out .* (2/3 * tand(fi));
%     else
     t_max_in = K0 .* sigma_p .* (tand(delta));
    t_max_out = P_i_out .* (tand(delta));
%     end
    % Preallocate arrays for efficiency
    t_out = zeros(size(displacement));
    t_in = zeros(size(displacement));

    % Compute settlement values using vectorized operations
    mask = abs(displacement) <= 2.5;
    t_out(mask) = (displacement(mask) / 2.5) .* t_max_out(mask);
    t_out(~mask) = sign(displacement(~mask)) .* t_max_out(~mask);

    t_in(mask) = (displacement(mask) / 2.5) .* t_max_in(mask);
    t_in(~mask) = sign(displacement(~mask)) .* t_max_in(~mask);

    % Scale settlement values
    scale_factor = (l_b / n_l) * (pi / n);
    tz_in_final = t_in * scale_factor * di_b; % kN
    tz_out_final = t_out * scale_factor * do_b; % kN
% disp(size(tz_in_final));
% disp(size(tz_out_final));
end
