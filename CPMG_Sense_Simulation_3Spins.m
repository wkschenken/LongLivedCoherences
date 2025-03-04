function CPMG = CPMG_Sense_Simulation_3Spins(tau,r12,r13,r23,theta12,theta13,theta23,B1,B2,B3,t_max,N_steps,eps)
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
    
    % Matrix to hold the simulation data for these runs; the first row
    % holds the expectation value of Sx, the second holds the entanglement
    % entropy 
%     CPMG = zeros(1, N_steps);
    % But if we only need the last datapoint at the end of the evolution
    % (as we do when we sweep e.g. epsilon or signal frequency and monitor
    % coherence) then we only need to keep track of the density matrix
    % until the final point i.e. no need to extract <Sx> at every time. 
    CPMG=0;
    
    % Number of discrete time steps between pulses
    s = round(2*tau/dt);
    
    % Run the simulation
    for n=1:1:N_steps
        % Calculate expectation value of Sx as a function of time
%         CPMG(1,n) = real(trace(rho*(SxT)));
        % Calculate entanglement entropy
%         if n>1
%             rho1 = [(rho(1,1) + rho(2,2) + rho(3,3) + rho(4,4)) (rho(1,5)+rho(2,6)+rho(3,7)+rho(4,8)); (rho(5,1)+rho(6,2)+rho(7,3)+rho(8,4)) (rho(5,5)+rho(6,6)+rho(7,7)+rho(8,8))];
%             CPMG(2,n) = real(trace(rho1*logm(rho1)));
%         end
        % time step to evolve the system in the absence of pulses
        rho = expm(-1i*H(n*dt)*dt)*rho*expm(1i*H(n*dt)*dt);
        % Apply a pulse after the appropriate number of timesteps
        if mod(n, s)==0
            rho = expm(-1i*(pi+eps)*(SxT)/2)*rho*expm(1i*(pi+eps)*(SxT)/2);
        end
    end
    % Calculate expectation value of Sx at the final time
    CPMG = real(trace(rho*(SxT)));

    
    
        
end
