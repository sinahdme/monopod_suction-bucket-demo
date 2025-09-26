function transformedPoint = xtransformation(P, Q, eulerAngles)
    % Extract Euler angles
    phi = eulerAngles(1); % Roll angle
    theta = eulerAngles(2); % Pitch angle
    psi = eulerAngles(3); % Yaw angle
    
    % Calculate rotation matrices
    R_phi = [1, 0, 0;
             0, cos(phi), -sin(phi);
             0, sin(phi), cos(phi)];

    R_theta = [cos(theta), 0, sin(theta);
               0, 1, 0;
              -sin(theta), 0, cos(theta)];

    R_psi = [cos(psi), -sin(psi), 0;
             sin(psi), cos(psi), 0;
             0, 0, 1];

    % Combine rotation matrices (in reverse order)
    R = (R_phi * R_theta * R_psi);
    
    % Translate point P to the origin by subtracting Q
    P_translated = P - Q;
    
    % Apply the rotation matrix to the translated point
    transformedPoint = R * P_translated + Q;
end
