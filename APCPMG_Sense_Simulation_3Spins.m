function APCPMG = APCPMG_Sense_Simulation_3Spins(tau,r12,r13,r23,theta12,theta13,theta23,B1,B2,B3,t_max,N_steps,eps)
% Bz in gauss, r in nm, t_max in us, n = density in nm^-3

% Simulating the coherence of three spins (NVs) interacting via dipole
% interactions then scaling up to more spins if possible

    % Spin operators for each NV (in the spin-1/2 subspace)
    Id = [1 0; 0 1];
    Sz = [1 0; 0 -1];
    Sy = [0 -1i; 1i 0];
    Sx = [0 1; 1 0];
   
    
    Sz1 = kron(kron(Sz, Id), Id);
    Sz2 = kron(kron(Id, Sz), Id);
    Sz3 = kron(kron(Id, Id), Sz);
    Sx1 = kron(kron(Sx, Id), Id);
    Sx2 = kron(kron(Id, Sx), Id);
    Sx3 = kron(kron(Id, Id), Sx);
    Sy1 = kron(kron(Sy, Id), Id);
    Sy2 = kron(kron(Id, Sy), Id);
    Sy3 = kron(kron(Id, Id), Sy);

    SzT = Sz1 + Sz2 + Sz3;
    SxT = Sx1 + Sx2 + Sz3;
    SyT = Sy1 + Sy2 + Sz3;
    
    Sx12 = kron(kron(Sx, Sx), Id);
    Sx13 = kron(kron(Sx, Id), Sx);
    Sx23 = kron(kron(Id, Sx), Sx);
    
    Sy12 = kron(kron(Sy, Sy), Id);
    Sy13 = kron(kron(Sy, Id), Sy);
    Sy23 = kron(kron(Id, Sy), Sy);
    
    Sz12 = kron(kron(Sz, Sz), Id);
    Sz13 = kron(kron(Sz, Id), Sz);
    Sz23 = kron(kron(Id, Sz), Sz);

   
    
    % Dipolar couplings
    H_Dip_12 = ((326/(2*pi))/r12^3)*(1-3*cos(theta12)^2)*(Sx12 + Sy12 - Sz12);
    H_Dip_13 = ((326/(2*pi))/r13^3)*(1-3*cos(theta13)^2)*(Sx13 + Sy13 - Sz13);
    H_Dip_23 = ((326/(2*pi))/r23^3)*(1-3*cos(theta23)^2)*(Sx23 + Sy23 - Sz23);
    
    H_Dip = H_Dip_12 +H_Dip_13+H_Dip_23;
    
    % Static disorder
    H_zeeman1 = B1*Sz1;
    H_zeeman2 = B2*Sz2;
    H_zeeman3 = B3*Sz3;
    
    H_zeeman = H_zeeman1+H_zeeman2+H_zeeman3;
    
    % Total Hamiltonian consisting of only static disorder and dipolar
    % couplings
    H = @(t) H_Dip+H_zeeman;
    
    
    % Define the initial density matrix; initialize along +x 
    init = (1/2)*[1 1; 1 1];
    rho = kron(kron(init, init), init);
   
    % discrete time steps
    dt = t_max/N_steps;
    
    
    % 1/2 the number of discrete time steps between pulses
    s = round(tau/dt);
    round(N_steps/s);
    APCPMG = zeros(1, round(N_steps/(2*s)));
    pm = 0; % Index to hold whether this is the a +x or -x pulse
    pulseNo = 0; % index of the coherence array
    PulseMeasure = 0; % determines whether this instance is a measurement (echo) or a pulse
                        % even means pulse, odd means measure
    
    % Run the simulation
    for n=1:1:N_steps

        % time step to evolve the system in the absence of pulses
        rho = expm(-1i*H(n*dt)*dt)*rho*expm(1i*H(n*dt)*dt);
        % Apply a pulse after the appropriate number of timesteps
        if mod(n, s)==0 && mod(PulseMeasure, 2)==0
            if mod(pm, 2) == 0
                rho = expm(-1i*(pi+eps)*(SxT)/2)*rho*expm(1i*(pi+eps)*(SxT)/2);
            else
                rho = expm(1i*(pi+eps)*(SxT)/2)*rho*expm(-1i*(pi+eps)*(SxT)/2);
            end
            pulseNo = pulseNo+1;
            pm = pm+1;
            PulseMeasure = PulseMeasure + 1;
        % if the time step is halfway between pulses, measure the coherence
        % of the echo
        elseif mod(n, s)==0 && mod(PulseMeasure, 2)==1
            APCPMG(1,pulseNo) = real(trace(rho*(SxT)));
            PulseMeasure = PulseMeasure + 1;
        end
    end

    
    
        
end
