%_________________________________________________________________________%
% Adaptive Accelerated Gravitational Search Algorithm (AAGSA)             %
%                                                                         %
% Developed in MATLAB R2018b                                              %
%                                                                         %
% Inventor and programmer: Farshad Rezaei, PhD                            %
%                                                                         %
% e-Mail: farshad.rezaei@gmail.com                                        %
%         f.rezaei@alumni.iut.ac.ir                                       %
%                                                                         %
% Homepage: https://www.linkedin.com/in/farshad-rezaei-5a92559a/          %
%                                                                         %
% Main paper: Kamran, S.; Safavi, H.R.; Golmohammadi, M.H.; Rezaei, F.;   %
% Abd Elaziz, M.; Forestiero, A.; Lu, S. Maximizing Sustainability in     %
% Reservoir Operation under Climate Change Using a Novel Adaptive         %
% Accelerated Gravitational Search Algorithm. Water 2022,                 %
% 14, 905. https://doi.org/10.3390/w14060905                              %
%_________________________________________________________________________%

% This function is to initialize the position and velocity of the wolves to start the optimization process
function [pp,pv] = Initialization(np,nx,varmax,varmin,velmax,velmin)
pp=zeros(np,nx);
pv=zeros(np,nx);
for j=1:np
    pp(j,1:nx)=(varmax-varmin).*rand(1,nx)+varmin;
    pv(j,1:nx)=(velmax-velmin).*rand(1,nx)+velmin;
end