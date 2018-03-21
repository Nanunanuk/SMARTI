function [ R,T,A ] = AR_RAT( komplex_n1,komplex_nAR,komplex_n2,wavelength,d_AR,theta)

% General form of Fresnels equation when no AR coating is involved
% (Altermatt lectures)
if komplex_n1==komplex_nAR||komplex_n2==komplex_nAR
    
    cos_theta2=1/komplex_n2*sqrt(komplex_n2^2-komplex_n1^2*sin(theta)^2);
    rs=(komplex_n1*cos(theta)-komplex_n2*cos_theta2)/(komplex_n1*cos(theta)+komplex_n2*cos_theta2);
    rp=(komplex_n2*cos(theta)-komplex_n1*cos_theta2)/(komplex_n2*cos(theta)+komplex_n1*cos_theta2);
    R=0.5*rs*conj(rs)+0.5*rp*conj(rp);
    
    T=1-R;
    A=0;
    
    
else
    %   Transfer Matrix Methode by Macleod
    
    phi1 = asin(komplex_n1 * sin(theta)/komplex_nAR);
    phi2 = asin(komplex_nAR * sin(phi1)/komplex_n2);
    
    
    k0 =  (2 * pi / wavelength);
    delta = k0 * komplex_nAR * d_AR * cos (phi1);
    
    y = 2.3544e-3; % [Si] "optical admittance" im vacuum
    
    % "optical admittance" im Medium mit index n - % S-Polarisation (TE)
    Y0_s = komplex_n1 * y * cos(theta);
    Y1_s = komplex_nAR * y * cos(phi1);
    Y2_s = komplex_n2 * y * cos(phi2);
    
    % "optical admittance" im Medium mit index n - % P-Polarisation (TM)
    Y0_p = komplex_n1 * y * (1/cos(theta));
    Y1_p = komplex_nAR * y * (1/cos(phi1));
    Y2_p = komplex_n2 * y * (1/cos(phi2));
    
    Matrix_s = [ cos(delta), (1i/Y1_s)*sin(delta)
        1i*Y1_s*sin(delta), cos(delta) ];
    
    Matrix_p = [ cos(delta), (1i/Y1_p)*sin(delta)
        1i*Y1_p*sin(delta), cos(delta) ];
    
    BC_s = Matrix_s*[1;Y2_s];
    B_s = BC_s(1);
    C_s = BC_s(2);
    
    BC_p = Matrix_p*[1;Y2_p];
    B_p = BC_p(1);
    C_p = BC_p(2);
    
    R_s = ((Y0_s * B_s - C_s)/(Y0_s * B_s + C_s)) * ((Y0_s * B_s - C_s)/(Y0_s * B_s + C_s))'; % TE - Reflexion
    T_s = (4 * Y0_s * real(Y2_s) ) / ((Y0_s * B_s + C_s) * (Y0_s * B_s + C_s)');              % TE - Transmission
    A_s = (4 * Y0_s * real(B_s*C_s' - Y2_s)) / ((Y0_s * B_s + C_s) * (Y0_s * B_s + C_s)');    % TE - Absorption
    
    R_p = ((Y0_p * B_p - C_p)/(Y0_p * B_p + C_p)) * ((Y0_p * B_p - C_p)/(Y0_p * B_p + C_p))'; % TM - Reflexion
    T_p = (4 * Y0_p * real(Y2_p) ) / ((Y0_p * B_p + C_p) * (Y0_p * B_p + C_p)');              % TM - Transmission
    A_p = (4 * Y0_p * real(B_p*C_p' - Y2_p)) / ((Y0_p * B_p + C_p) * (Y0_p * B_p + C_p)');    % TM - Absorption
    
    R = (R_s+R_p)/2;
    A = (A_s+A_p)/2;
    T = (T_s+T_p)/2;
    
end

end