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

% The initial parameters that you need are:
%_________________________________________________________________________
% fobj=@YourCostFunction
% nx=number of your variables
% lb=the lower bound of variables which can generally be a fixed number or a vector
% ub=the upper bound of variables which can generally be a fixed number or a vector
% notice: if the lower and upper bounds are not fixed for all variables, 
% they appear in the forms of the vectors "varmin" and "varmax", as illustrated in following

% To run AAGSA: [z_iter,z_final,pos_final]=AAGSA(np,nx,maxit,varmax,varmin,velmax,velmin,thr_pow_max,thr_pow_min,c_random,final_per,g_initial,alpha,mut_max,mut_min,fobj);

%_________________________________________________________________________
% Set the required parameters to run the AAGSA algorithm

% This code is for solving the minimization problems. To maximize a desired 
% cost function,please implement this code upon inverting the sign of the cost function

clc
clear
close all
tic
run=30; % Maximum number of the algorithm runnings conducted
np=30; % Number of search agents
Function_name='F1'; % Name of the test function that can be from F1 to F13 
maxit=1000; % Maximum number of iterations
thr_pow_max=4; % Upper bound set for the power in Eq.(17) 
thr_pow_min=1; % Lower bound set for the power in Eq.(17) 
c_random=2; % A coefficient set in Eq.(21)
[lb,ub,nx,fobj]=Objective_Function(Function_name); % Load details of the selected benchmark function
varmax=ub*ones(1,nx); % Upper bound defined for the positions which can generally be a desired vector
varmin=lb*ones(1,nx); % Lower bound defined for the positions which can generally be a desired vector
limvel=0.1; % A ratio of the maximum distance in the search space to form the maximum velocity 
velmax=limvel*(varmax(1,1:nx)-varmin(1,1:nx)); % Upper bound defined for the velocities
velmin=-velmax; % Lower bound defined for the velocities
final_per=2; % Final percentage of the population that are included in the First Kbest elite agents
g_initial=100; % Initial gravitational constant
alpha=20; % A constant helping the gravitational constant be highly reduced
mut_max=0.9; % Maximum mutation coefficient used in the bound constraint handling technique
mut_min=0.1; % Minimum mutation coefficient used in the bound constraint handling technique
z_iter_main=zeros(run,maxit);
z_final_main=zeros(run);
pos_final_main=zeros(run,nx);
x1=zeros(maxit);
y1=zeros(maxit);

% Run the AAGSA algorithm for "run" times 
for nrun=1:run
    [z_iter,z_final,pos_final]=AAGSA(np,nx,maxit,varmax,varmin,velmax,velmin,thr_pow_max,thr_pow_min,c_random,final_per,g_initial,alpha,mut_max,mut_min,fobj);
     z_iter_main(nrun,1:maxit)=z_iter(1:maxit);
     z_final_main(nrun)=z_final;
     pos_final_main(nrun,1:nx)=pos_final(1:nx);
end

% Display the comprehensive results
disp(['The final statistical results calculated when implementing the AAGSA algorithm for ',num2str(run),' times are as follows:']);
disp(['The average of the final objective function values calculated over ',num2str(run),' times = ',num2str(mean(z_final_main(1:run)))]);
disp(['The median of the final objective function values calculated over ',num2str(run),' times = ',num2str(median(z_final_main(1:run)))]);
disp(['The best of the final objective function values calculated over ',num2str(run),' times = ',num2str(min(z_final_main(1:run)))]);
disp(['The standard deviation of the final objective function values calculated over ',num2str(run),' times = ',num2str(std(z_final_main(1:run)))]);

% Plot the convergence curve of the AAGSA over the course of iterations
for i=1:maxit
    x1(i)=i;sum1=0;
    for j=1:run
        sum1=sum1+z_iter_main(j,i);
    end
    y1(i)=sum1/run;
end
semilogy(x1,y1,'-r')
xlabel('Iteration');
ylabel('Average best-so-far');
legend('AAGSA');
hold on
time_aagsa = toc;
disp(['Elapsed time of running the AAGSA for ',num2str(run),' times = ',num2str(time_aagsa),' seconds']);