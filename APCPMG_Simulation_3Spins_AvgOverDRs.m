DR = 200; % Number of disorder realizations to average over
N_steps = 5*10^3; % Number of discrete time steps

sigma_Bz = 1; % MHz, Spread in zeeman frequencies
Bz_avg = 2; % detuning, effectively a global field

tau = .2; %us, 2tau=time between pulses. 
t_max = 20; % us

nPulses = round(t_max/(2*tau));


NN = 5; % This is how many densities or taus to simulate. Each one generates a figure of coherence vs epsilon. 
n_max = -2; % 10^n_max, nm^-3, max density
n_min=-6;
n = logspace(n_min,n_max, NN);
nnn=0; % linear index for the density

%note at a density of 10^-3 nm^-3 the probability
% for a spin to have a nearest neighbor greater than 2.3 nm
% is ~95%. Use this to determine the appropriate time step
% s.t. most of the simulations are accurate. 

eps = 0.6;
Sig_Sweep_PulseNo = zeros(NN, nPulses); % row signifies density, column signifies pulse number
Std_Sweep_PulseNo = zeros(NN, nPulses); % row signifies density, column signifies pulse number
                
                
% sweep density
for nn=n
    nn
    nnn=nnn+1;

    avg_r = 0.55/(nn^(1/3)); % nm, avg spacing between nearest neighbors in 3d
    sigma_r = 0.21/(nn^(1/3)); % nm, std dev of spacing between nearest neighbors in 3d

    % average over disorder realizations for a given density and
    % epsilon

    Coherence = zeros(DR, nPulses);

    for ii=1:1:DR
        % define random coordinates for grid of 3 spins. Put spin 1 at the origin. 
        r12 = avg_r+sigma_r*randn;
        if r12<0
            r12=avg_r;
        end
        x2 = 2*r12*(rand-1/2);
        y2 = 2*((r12^2-x2^2)^(1/2))*(rand-1/2);
        z2 = (r12^2-x2^2-y2^2)^(1/2);
        theta12 = atan(((x2^2+y2^2)^(1/2))/z2);

        r13 = avg_r+sigma_r*randn;
        if r13<0
            r13=avg_r;
        end
        x3 = 2*r13*(rand-1/2);
        y3 = 2*((r13^2-x3^2)^(1/2))*(rand-1/2);
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
        B3 = Bz_avg-(B1+B2);
        
        
        Coherence(ii,:) = APCPMG_Sense_Simulation_3Spins(tau, r12, r23, r13, theta12,...
            theta23, theta13,B1,B2,B3, t_max, N_steps,eps);

        100*((nnn-1)/NN + ii/(NN*DR))
        
        %average over DRs
    end

    Coherence_Avgd = mean(Coherence)/2;
    Coherence_Stdev = (std(Coherence)/(DR^(1/2)))/2;       
        
    Sig_Sweep_PulseNo(nnn, :) = Coherence_Avgd; % row signifies density, column signifies pulse number
    Std_Sweep_PulseNo(nnn, :) = Coherence_Stdev; % row signifies density, column signifies pulse number
end 


time = linspace(0, 2*nPulses*tau, nPulses);
box on;
hold on;
grid on;
for jj=1:1:NN
%     rdn = -(floor(log10(n(jj)*10^21))-1);
%     txt = [num2str(round(n(jj)*10^21, rdn)), 'cm^{-3}'];
    txt = ['density = ', num2str(n(jj)*10^21), ' /cm^3'];
    errorbar(time, Sig_Sweep_PulseNo(jj,:),Std_Sweep_PulseNo(jj,:), '-s','MarkerSize',10,'DisplayName',txt);
end
xlabel('Time (us)','FontSize',36,'FontName', 'Century Gothic') 
ylabel('APCPMG Coherence','FontSize',36,'FontName', 'Century Gothic')
legend show
