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
% AAGSA algorithm                                                                  
function [z_iter,z_final,pos_final] = AAGSA(np,nx,maxit,varmax,varmin,velmax,velmin,thr_pow_max,thr_pow_min,c_random,final_per,g_initial,alpha,mut_max,mut_min,fobj)
% disp(['Number of Iterations = ',num2str(it)]);
z_gbest=inf;
mass=zeros(np);
index=zeros(np);
optimal_pos=zeros(1,nx);
z=zeros(np);
pos_final=zeros(nx);
z_iter=zeros(maxit);

% Initialization process of the algorithm
[pp,pv]=Initialization(np,nx,varmax,varmin,velmax,velmin);

% Start the optimization process
it=1;

% Calculating the maximum distance between two agents in the search space
r_max=sqrt(sum((varmax(1,:)-varmin(1,:)).^2));

% Objective function evaluation
for j=1:np
    z(j)=fobj(pp(j,1:nx));
end

% Calculating Masses
g_constant=g_initial*exp(-alpha*it/maxit);
for j=1:np
    if j==1
        z_min=z(j);
        z_max=z(j);
    elseif z(j)<z_min 
        z_min=z(j);
    elseif z(j)>z_max
        z_max=z(j);
    end
end
for j=1:np
    mass(j)=(z_max-z(j))/(z_max-z_min);
    index(j)=j;
end
sum_mass=sum(mass(1:np));
mass(1:np)=mass(1:np)/sum_mass;

% Calculating Kbest nad Determining the First Kbest search Agents
kbest=final_per+(1-it/maxit)*(100-final_per);
kbest=round(np*kbest/100);
for j=1:np-1
    for jj=j+1:np
        if z(jj)<z(j)
            c1=z(j);
            z(j)=z(jj);
            z(jj)=c1;
            c2=index(j);
            index(j)=index(jj);
            index(jj)=c2;
        end
    end
end

% Calculating the Acceleration
thres_power=thr_pow_max-(thr_pow_max-thr_pow_min)*(it/maxit); % Eq.(19)
mu1=1-(it/maxit)^thres_power; % Eq.(17)
c_thres1=1-(1-it/maxit)/(1-mu1*(it/maxit)); % Eq.(15)
mu2=1-(it/maxit); % Eq.(18)
c_thres2=1-(1-it/maxit)/(1-mu2*(it/maxit)); % Eq.(16)
aa=zeros(np,nx);
for j=1:np
    for jj=1:kbest
        if index(jj)~=j
            p1(1:nx)=pp(index(jj),1:nx); 
            p2(1:nx)=pp(j,1:nx);
            rr=norm(p1-p2);
            random=c_thres1+rand(1,1)*(c_thres2-c_thres1); % Eq.(20)
            aa(j,1:nx)=aa(j,1:nx)+c_random*random*g_constant*...
                mass(index(jj))/(rr/r_max+eps)*(p1-p2);
        end
    end
end

% Determining the best-so-far objective value and the best-so-far agent
for j=1:np
    if z(j)<z_gbest
        z_gbest=z(j);
        pp_gbest(1:nx)=pp(j,1:nx);
    end
end
z_optimal(it)=z_gbest;
optimal_pos(it,:)=pp_gbest(:);

% Save the best-so-far objective value in the current run
z_iter(it)=z_optimal(it);

% The Main Loop
while it<maxit
    it=it+1;
    mut=mut_max-(mut_max-mut_min)*(it/maxit); % Eq.(27)
%     disp(['Number of Iterations= ',num2str(it)]);
    for j=1:np     
        pv(j,1:nx)=rand(1,nx).*pv(j,1:nx)+aa(j,1:nx); % Eq.(8)
        
        % Return back the velocity of the particles if going beyond the velocity boundaries
        flag4lbv=pv(j,:)<velmin(1,:);
        flag4ubv=pv(j,:)>velmax(1,:);
        pv(j,:)=pv(j,:).*(~(flag4lbv+flag4ubv))+velmin.*flag4lbv+velmax.*flag4ubv;
        pp(j,:)=pp(j,:)+pv(j,:); % Eq.(9)
        
        % Return back the position of the particles if going beyond the position boundaries
        flag4lbp=pp(j,:)<varmin(1,:);
        flag4ubp=pp(j,:)>varmax(1,:);
        varmin_new(1,:)=varmin(1,:)+rand(1,nx).*mut.*(varmax(1,:)-varmin(1,:)); % Eq.(25)
        varmax_new(1,:)=varmax(1,:)-rand(1,nx).*mut.*(varmax(1,:)-varmin(1,:)); % Eq.(26)
        pp(j,:)=pp(j,:).*(~(flag4lbp+flag4ubp))+varmin_new.*flag4lbp+varmax_new.*flag4ubp;
        
        % Objective function evaluations and determining of the personal best solutions and objectives
        z(j)=fobj(pp(j,:));
    end
    
    % Calculating Masses
    g_constant=g_initial*exp(-alpha*it/maxit);
    for j=1:np
        if j==1
            z_min=z(j);
            z_max=z(j);
        elseif z(j)<z_min
            z_min=z(j);
        elseif z(j)>z_max
            z_max=z(j);
        end
    end
    for j=1:np
        mass(j)=(z_max-z(j))/(z_max-z_min);
        index(j)=j;
    end
    sum_mass=sum(mass(1:np));
    mass(1:np)=mass(1:np)/sum_mass;

    % Calculating Kbest nad Determining the First Kbest search Agents
    kbest=final_per+(1-it/maxit)*(100-final_per);
    kbest=round(np*kbest/100);
    for j=1:np-1
        for jj=j+1:np
            if z(jj)<z(j)
                c1=z(j);
                z(j)=z(jj);
                z(jj)=c1;
                c2=index(j);
                index(j)=index(jj);
                index(jj)=c2;
            end
        end
    end
    
    % Calculating the Acceleration
    thres_power=thr_pow_max-(thr_pow_max-thr_pow_min)*(it/maxit); % Eq.(19)
    mu1=1-(it/maxit)^thres_power; % Eq.(17)
    c_thres1=1-(1-it/maxit)/(1-mu1*(it/maxit)); % Eq.(15)
    mu2=1-(it/maxit); % Eq.(18)
    c_thres2=1-(1-it/maxit)/(1-mu2*(it/maxit)); % Eq.(16)
    aa=zeros(np,nx);
    for j=1:np
        for jj=1:kbest
            if index(jj)~=j
                p1(1:nx)=pp(index(jj),1:nx); 
                p2(1:nx)=pp(j,1:nx);
                rr=norm(p1-p2);
                random=c_thres1+rand(1,1)*(c_thres2-c_thres1); % Eq.(20)
                aa(j,1:nx)=aa(j,1:nx)+c_random*random*g_constant*...
                    mass(index(jj))/(rr/r_max+eps)*(p1-p2);
            end
        end
    end
    
    % Determining the best-so-far objective value and the best-so-far agent
    for j=1:np
        if z(j)<z_gbest
            z_gbest=z(j);
            pp_gbest(1:nx)=pp(j,1:nx);
        end
    end
    z_optimal(it)=z_gbest;
    optimal_pos(it,:)=pp_gbest(:);

    % Save the best-so-far objective value in the current run
    z_iter(it)=z_optimal(it);
end

% Save the final best solution and objective revealed upon the end of the optimization process
z_final=z_optimal(maxit);
pos_final(1:nx)=optimal_pos(maxit,1:nx);
end