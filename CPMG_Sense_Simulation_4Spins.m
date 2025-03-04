function CPMG = CPMG_Sense_Simulation_4Spins(tau,r12,r13,r14,r23,r24,r34,theta12,theta13,theta14,theta23,theta24,...
    theta34,B1,B2,B3,B4,t_max,N_steps,eps)
% Bz in gauss, r in nm, t_max in us, n = density in nm^-3

% Simulating the coherence of 4 spins (NVs) interacting via dipole

    % Spin operators for each NV (in the spin-1/2 subspace)
    Id = [1 0; 0 1];
    Z = [1 0; 0 -1];
    Y = [0 -1i; 1i 0];
    X = [0 1; 1 0];
    
    Nspins=4;
    
    
    % Defining the desired matrices; first, in the order which to take the
    % kron product
    Z1 = {Z,Id,Id,Id};
    Z2 = {Id,Z,Id,Id};
    Z3 = {Id,Id,Z,Id};
    Z4 = {Id,Id,Id,Z};
    X1 = {X,Id,Id,Id};
    X2 = {Id,X,Id,Id};
    X3 = {Id,Id,X,Id};
    X4 = {Id,Id,Id,X};
    Y1 = {Y,Id,Id,Id};
    Y2 = {Id,Y,Id,Id};
    Y3 = {Id,Id,Y,Id};
    Y4 = {Id,Id,Id,Y};
    
    Z1Z2={Z,Z,Id,Id};
    Z1Z3={Z,Id,Z,Id};
    Z1Z4={Z,Id,Id,Z};
    Z2Z3={Id,Z,Z,Id};
    Z2Z4={Id,Z,Id,Z};
    Z3Z4={Id,Id,Z,Z};
    X1X2={X,X,Id,Id};
    X1X3={X,Id,X,Id};
    X1X4={X,Id,Id,X};
    X2X3={Id,X,X,Id};
    X2X4={Id,X,Id,X};
    X3X4={Id,Id,X,X};
    Y1Y2={Y,Y,Id,Id};
    Y1Y3={Y,Id,Y,Id};
    Y1Y4={Y,Id,Id,Y};
    Y2Y3={Id,Y,Y,Id};
    Y2Y4={Id,Y,Id,Y};
    Y3Y4={Id,Id,Y,Y};
    
    
    % first element of the kron product
    z1 = Z1{1};
    z2 = Z2{1};
    z3 = Z3{1};
    z4 = Z4{1};
    x1 = X1{1};
    x2 = X2{1};
    x3 = X3{1};
    x4 = X4{1};
    y1 = Y1{1};
    y2 = Y2{1};
    y3 = Y3{1};
    y4 = Y4{1};
    z1z2 = Z1Z2{1};
    z1z3 = Z1Z3{1};
    z1z4 = Z1Z4{1};
    z2z3 = Z2Z3{1};
    z2z4 = Z2Z4{1};
    z3z4 = Z3Z4{1};
    x1x2 = X1X2{1};
    x1x3 = X1X3{1};
    x1x4 = X1X4{1};
    x2x3 = X2X3{1};
    x2x4 = X2X4{1};
    x3x4 = X3X4{1};
    y1y2 = Y1Y2{1};
    y1y3 = Y1Y3{1};
    y1y4 = Y1Y4{1};
    y2y3 = Y2Y3{1};
    y2y4 = Y2Y4{1};
    y3y4 = Y3Y4{1};
    
    % multiply out the rest of the kron product
    for k = 2:Nspins
        z1 = kron(z1,Z1{k});
        z2 = kron(z2,Z2{k});
        z3 = kron(z3,Z3{k});
        z4 = kron(z4,Z4{k});
        x1 = kron(x1,X1{k});
        x2 = kron(x2,X2{k});
        x3 = kron(x3,X3{k});
        x4 = kron(x4,X4{k});
        y1 = kron(y1,Y1{k});
        y2 = kron(y2,Y2{k});
        y3 = kron(y3,Y3{k});
        y4 = kron(y4,Y4{k});
        z1z2 = kron(z1z2,Z1Z2{k});
        z1z3 = kron(z1z3,Z1Z3{k});
        z1z4 = kron(z1z4,Z1Z4{k});
        z2z3 = kron(z2z3,Z2Z3{k});
        z2z4 = kron(z2z4,Z2Z4{k});
        z3z4 = kron(z3z4,Z3Z4{k});
        x1x2 = kron(x1x2,X1X2{k});
        x1x3 = kron(x1x3,X1X3{k});
        x1x4 = kron(x1x4,X1X4{k});
        x2x3 = kron(x2x3,X2X3{k});
        x2x4 = kron(x2x4,X2X4{k});
        x3x4 = kron(x3x4,X3X4{k});
        y1y2 = kron(y1y2,Y1Y2{k});
        y1y3 = kron(y1y3,Y1Y3{k});
        y1y4 = kron(y1y4,Y1Y4{k});
        y2y3 = kron(y2y3,Y2Y3{k});
        y2y4 = kron(y2y4,Y2Y4{k});
        y3y4 = kron(y3y4,Y3Y4{k});
        
    end
    
    
    zT = z1+z2+z3+z4;
    xT = x1+x2+x3+x4;
    yT = y1+y2+y3+y4;
    

   
    
    % Dipolar couplings
    H_Dip_12 = ((326/(2*pi))/r12^3)*(1-3*cos(theta12)^2)*(x1x2 + y1y2 - z1z2);
    H_Dip_13 = ((326/(2*pi))/r13^3)*(1-3*cos(theta13)^2)*(x1x3 + y1y3 - z1z3);
    H_Dip_14 = ((326/(2*pi))/r14^3)*(1-3*cos(theta14)^2)*(x1x4 + y1y4 - z1z4);
    H_Dip_23 = ((326/(2*pi))/r23^3)*(1-3*cos(theta23)^2)*(x2x3 + y2y3 - z2z3);
    H_Dip_24 = ((326/(2*pi))/r24^3)*(1-3*cos(theta24)^2)*(x2x4 + y2y4 - z2z4);
    H_Dip_34 = ((326/(2*pi))/r34^3)*(1-3*cos(theta34)^2)*(x3x4 + y3y4 - z3z4);

    H_Dip = H_Dip_12 +H_Dip_13+H_Dip_14+H_Dip_23+H_Dip_24+H_Dip_34;
    
    % Static disorder
    H_zeeman1 = B1*z1;
    H_zeeman2 = B2*z2;
    H_zeeman3 = B3*z3;
    H_zeeman4 = B4*z4;
    
    H_zeeman = H_zeeman1+H_zeeman2+H_zeeman3+H_zeeman4;
    
    % Total Hamiltonian consisting of only static disorder and dipolar
    % couplings
    H = @(t) H_Dip+H_zeeman;
    
    
    % Define the initial density matrix; initialize along +x 
    init = (1/2)*[1 1; 1 1];
    rho = kron(kron(kron(init, init), init),init);
   
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
            rho = expm(-1i*(pi+eps)*(xT)/2)*rho*expm(1i*(pi+eps)*(xT)/2);
        end
    end
    % Calculate expectation value of Sx at the final time
    CPMG = real(trace(rho*(xT)));
    
    
        
end


