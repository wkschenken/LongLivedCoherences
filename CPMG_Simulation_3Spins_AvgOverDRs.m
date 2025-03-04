DR = 100; % Number of disorder realizations to average over
N_steps = 5*10^3; % Number of discrete time steps

sigma_Bz = 1; % MHz, Spread in zeeman frequencies

% tau = .2; %us, 2tau=time between pulses. 
t_max = 20; % us

NN = 4; % This is how many densities or taus to simulate. Each one generates a figure of coherence vs epsilon. 
% n_max = -3; % 10^n_max, nm^-3, max density
% n_min=-6;
% n = logspace(n_min,n_max, NN)
nnn=0; % linear index for the density

tau_min = 0.05; 
tau_max = 0.35; 
tau = linspace(tau_min, tau_max, NN);

density = 10^-4; 

%note at a density of 10^-3 nm^-3 the probability
% for a spin to have a nearest neighbor greater than 2.3 nm
% is ~95%. Use this to determine the appropriate time step
% s.t. most of the simulations are accurate. 

mm=10; %# of epsilon steps
eps_min = -pi/30;
eps_max = pi/30;
Sig_Sweep_eps = zeros(NN, mm+1); % row signifies density, column signifies epsilon
Std_Sweep_eps = zeros(NN, mm+1); % row signifies density, column signifies epsilon
                
                
% sweep density or tau 
for nn=tau
    tautau = ['tau = ', num2str(time), ' us'];
    disp(tautau)
    nnn=nnn+1;
    
    avg_r = 0.55/(density^(1/3)); % nm, avg spacing between nearest neighbors in 3d
    sigma_r = 0.21/(density^(1/3)); % nm, std dev of spacing between nearest neighbors in 3d
    
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
            % define random coordinates for grid of 3 spins. Put spin 1 at the origin. 
            r12 = avg_r+sigma_r*randn;
            if r12<0
                r12=avg_r;
            end
            x2 = r12*(rand-1/2);
            y2 = ((r12^2-x2^2)^(1/2))*(rand-1/2);
            z2 = (r12^2-x2^2-y2^2)^(1/2);
            theta12 = atan(((x2^2+y2^2)^(1/2))/z2);

            r13 = avg_r+sigma_r*randn;
            if r13<0
                r13=avg_r;
            end
            x3 = r13*(rand-1/2);
            y3 = ((r13^2-x3^2)^(1/2))*(rand-1/2);
            z3 = (r13^2-x3^2-y3^2)^(1/2);
            theta13 = atan(((x3^2+y3^2)^(1/2))/z3);

            x23 = x3-x2;
            y23 = y3-y2;
            z23 = z3-z2;
            r23 = (x23^2+y23^2+z23^2)^(1/2);
            theta23 = atan(((x23^2+y23^2)^(1/2))/z23);
            
            % define random fields such that their average is still zero
            B1 = sigma_Bz*randn;
            B2 = sigma_Bz*randn;
            B3 = -(B1+B2);
            
            
            Coherence(ii,1) = CPMG_Sense_Simulation_3Spins(nn, r12, r23, r13, theta12,...
                theta23, theta13,B1,B2,B3, t_max, N_steps,eps_min+f*(eps_max-eps_min)/mm);
        end
        
        %average over DRs
        Coherence_Avgd = mean(Coherence);
        Coherence_Stdev = std(Coherence);       
        
        
        Sig_Sweep_eps(nnn, f+1) = Coherence_Avgd/3
        Std_Sweep_eps(nnn, f+1) = Coherence_Stdev/(DR^(1/2))/3

    end
end 


eps = linspace(eps_min*180/pi, eps_max*180/pi, mm+1);
box on;
hold on;
grid on;
for jj=1:1:NN
%     rdn = -(floor(log10(n(jj)*10^21))-1);
%     txt = [num2str(round(n(jj)*10^21, rdn)), 'cm^{-3}'];
    txt = ['\tau = ', num2str(tau(jj)*(10^3)), ' ns'];
    errorbar(eps, Sig_Sweep_eps(jj,:),Std_Sweep_eps(jj,:), '-s','MarkerSize',10,'DisplayName',txt);
end
xlabel('Epsilon (degrees)','FontSize',36,'FontName', 'Century Gothic') 
ylabel('Coherence','FontSize',36,'FontName', 'Century Gothic')
legend show
