DR = 120; % Number of disorder realizations to average over
N_steps = 3*10^3; % Number of discrete time steps

sigma_Bz = .5; % MHz, Spread in zeeman frequencies

% tau = .2; %us, 2tau=time between pulses. 
tt = 3; 
tau_min = 0.05;
tau_max = 0.35;
tau = linspace(tau_min, tau_max, tt);

t_max = 10; % us

% NN = 3; % This is how many densities to simulate. Each one generates a figure of coherence vs epsilon. 
% n_max = -5; % 10^n_max, nm^-3, max density
% n_min=-7;
% n = logspace(n_min,n_max, NN)
nn = 10^-5;


nnn=0; % linear index for the 2nd sweep variable

%note at a density of 10^-3 nm^-3 the probability
% for a spin to have a nearest neighbor greater than 2.3 nm
% is ~95%. Use this to determine the appropriate time step
% s.t. most of the simulations are accurate. 

mm=15; %# of epsilon steps
eps_min = -pi/2;
eps_max = pi/2;
Sig_Sweep_eps = zeros(tt, mm+1); % row signifies density, column signifies epsilon
Std_Sweep_eps = zeros(tt, mm+1); % row signifies density, column signifies epsilon
                
                
% sweep some parameter (tau or density or other)
for time=tau
    nnn = nnn+1;
    tautau = ['tau = ', num2str(time), ' us'];
    disp(tautau)
    
    avg_r = 0.55/(nn^(1/3)); % nm, avg spacing between nearest neighbors in 3d
    sigma_r = 0.21/(nn^(1/3)); % nm, std dev of spacing between nearest neighbors in 3d
    
    %sweep epsilon for a given density
    for f=0:1:mm
        epseps = ['epsilon = ', num2str((180/pi)*(eps_min+f*(eps_max-eps_min)/mm)), ' degrees'];
        disp(epseps)
        Coherence = zeros(DR, 1); % One row holds the dynamics of Sx, each row is a different DR
        Coherence_Avgd = 0; % averaged decay
        Coherence_Stdev = 0; % Stdev 
        % average over disorder realizations for a given density and
        % epsilon
        for ii=1:1:DR
            % define random coordinates for grid of 4 spins. Put spin 1 at the origin. 
            r12 = avg_r+sigma_r*randn;
            if r12<0
                r12=avg_r;
            end
            x2 = 2*r12*(rand-1/2);
            y2 = 2*((r12^2-x2^2)^(1/2))*(rand-1/2);
            z2 = 2*(r12^2-x2^2-y2^2)^(1/2);
            theta12 = atan(((x2^2+y2^2)^(1/2))/z2);

            r13 = avg_r+sigma_r*randn;
            if r13<0
                r13=avg_r;
            end
            x3 = 2*r13*(rand-1/2);
            y3 = 2*((r13^2-x3^2)^(1/2))*(rand-1/2);
            z3 = 2*(r13^2-x3^2-y3^2)^(1/2);
            theta13 = atan(((x3^2+y3^2)^(1/2))/z3);
            
            r14 = avg_r+sigma_r*randn;
            if r14<0
                r14=avg_r;
            end
            x4 = 2*r14*(rand-1/2);
            y4 = 2*((r14^2-x4^2)^(1/2))*(rand-1/2);
            z4 = 2*(r14^2-x4^2-y4^2)^(1/2);
            theta14 = atan(((x4^2+y4^2)^(1/2))/z4);
                   
            x23 = x3-x2;
            y23 = y3-y2;
            z23 = z3-z2;
            r23 = (x23^2+y23^2+z23^2)^(1/2);
            theta23 = atan(((x23^2+y23^2)^(1/2))/z23);
            
            x24 = x4-x2;
            y24 = y4-y2;
            z24 = z4-z2;
            r24 = (x24^2+y24^2+z24^2)^(1/2);
            theta24 = atan(((x24^2+y24^2)^(1/2))/z24);
            
            x34 = x3-x4;
            y34 = y3-y4;
            z34 = z3-z4;
            r34 = (x34^2+y34^2+z34^2)^(1/2);
            theta34 = atan(((x34^2+y34^2)^(1/2))/z34);
            
            % define random fields such that their average is still zero
            B1 = sigma_Bz*randn;
            B2 = sigma_Bz*randn;
            B3 = sigma_Bz*randn;
            B4 = -(B1+B2+B3);
            
            
            Coherence(ii,1) = CPMG_Sense_Simulation_4Spins(time,r12,r13,r14,r23,r24,r34,theta12,...
                theta13,theta14,theta23,theta24,theta34,B1,B2,B3,B4,t_max,N_steps,eps_min+f*(eps_max-eps_min)/mm);
        end
        
        %average over DRs
        Coherence_Avgd = mean(Coherence);
        Coherence_Stdev = std(Coherence);       
        
        
        Sig_Sweep_eps(nnn, f+1) = Coherence_Avgd
        Std_Sweep_eps(nnn, f+1) = Coherence_Stdev/(DR^(1/2))
    end
end 


eps = linspace(eps_min*180/pi, eps_max*180/pi, mm+1);
box on;
hold on;
grid on;
for jj=1:1:tt
%     rdn = -(floor(log10(n(jj)*10^21))-1);
%     txt = [num2str(round(n(jj)*10^21, rdn)), 'cm^{-3}'];
    txt = [num2str(tau(jj)), 'us'];
    errorbar(eps, Sig_Sweep_eps(jj,:),Std_Sweep_eps(jj,:), '-s','MarkerSize',10,'DisplayName',txt);
end
xlabel('Epsilon (degrees)','FontSize',36,'FontName', 'Century Gothic') 
ylabel('\langle S_x(t_f)\rangle','FontSize',36,'FontName', 'Century Gothic')
legend show
