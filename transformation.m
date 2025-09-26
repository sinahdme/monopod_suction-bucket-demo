function transformedPoint = transformation(P, Q, eulerAngles)
% this function calculates the position of the point P in 3D in the
% cartesian frame with respect to center of rotation Q 
    % Extract Euler angles
    phi = eulerAngles(1); % Roll angle around X-axis
    theta = eulerAngles(2); % Pitch angle around Y-axis
    psi = eulerAngles(3); % Yaw angle around Z-axis
    
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

    % Combine rotation matrices
    R =R_psi' * R_theta'* R_phi';
    
    % Translate point P to the origin by subtracting Q
    P_translated = P - Q;
    
    % Apply the rotation matrix to the translated point
    transformedPoint = (R * P_translated) + Q;
end