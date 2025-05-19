function area = Nasmyth_area(k_nd)
%
% Ryan Walter, 10/25/2024
%
% Function returns the integral of the Nasmyth spectrum from k = 0 to k =
%  k_nd divided by the complete integral from k = 0 to k = infinity. The
%  wavenumber k_nd must be non-dimensional where k_nd = k * L_K where k is
%  in units of cpm (cycles per metre) and L_k = (\nu^3/\epsilon)^{1/4} is
%  the Kolmogorov length.
% The area is a fraction that can be used to boost an estimate of \epsilon
% derived from an integral of a shear spectrum that uses a finite upper
% limit of spectral integration.
% 
%



% Equation 20 from Lueck 2022b (Statistics of Oceanic Turbulence Part II)   

a =  61.5;
b =  18.1;
c =  52.5;

x = k_nd.^(4/3);

area = tanh(a * x) - b * x .* exp(-c * x);  



% %Lueck Spectrum
% a = 65.5; % 65.5
% b = 19.0; %18.5
% c = 54.5; % 53.5
% 
% x = k_nd.^(4/3); % makes the syntax simpler
% area = tanh(a * x) - b * x .* exp(-c * x);  


