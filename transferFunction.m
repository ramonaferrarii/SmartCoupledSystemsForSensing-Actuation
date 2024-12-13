function [tf, tau] = transferFunction(k, wi, R_L, Cpi, csi, w)
%% Formulas
% Variables
tau = -R_L * Cpi;
% Transfer function V/Fi
% tf = (k * tau *(w * k^2 * tau + 2 * tau * w))/((w * k^2 * tau + 2 * tau * w)^2 + 4) ...
%     + 1i * (k * tau)/((w * k^2 * tau + 2 * tau * w)^2 + 4);
tf = -(1i*tau*w*wi*k)/(wi^2-w^2-2*csi*wi*w^2*tau+1i*(2*csi*wi*w+wi^2*w*tau-tau*w^3+wi^2*w*k^2*tau));



