function [GPxhatOut] = SimPropPolyMGPs(N_GP,GPxhatIn,F,B,u,P,Q,origin,max_speed,debugFlag)


% F - the system matrix
% F - the system matrix
% F = [   1 0     delT    0
%     0 1     0       delT
%     0 0     1       0
%     0 0     0       1 ];
% B - Control matrix
% B = [   0       0
%     0       0
%     1       0
%     0       1  ];
%F = [1 0 1 0; 0 1 0 1; 0 0 1 0; 0 0 0 1];
% u - wind velocity
%u = [0.1 0.1]';
%Q - Model Uncertainty
% Qfactor = 1; % process noise mult factor
% sigmaQ = Qfactor*sqrt(0.1);
% Q = sigmaQ^2 * [ 1 0 0 0% process noise covariance matrix
%     0 1 0 0
%     0 0 1 0
%     0 0 0 1];
%configure scenario physical parameters
%center location of the phenomenon
originX = origin(1);
originY = origin(2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_state = 4;
% GPxhatOut ~ [ x y progress rate]
GPxhatOut = zeros(N_GP,N_state); % state prediction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1: N_GP
    %P = Pgp(i,:,:);
    [GPxhatOut(i,:) P]= predict_update(GPxhatIn(i,:)', P, F,B,Q, u);
    %Pgp(i,:,:) = P;
    %point #i
    if i==N_GP
        leftNeighbour = 1;
    else
        leftNeighbour = i+1;
    end
    if i==1
        rightNeighbour = N_GP;
    else
        rightNeighbour = i-1;
    end
    
    %TODO: consider slope of terain
    neighboursVel = (GPxhatIn(leftNeighbour,[3 4]) + GPxhatIn(rightNeighbour,[3 4]))/2;
    % Set the heuristic of spreading rate consider environment params
    updatedVel = 2/3 * GPxhatOut(i,[3 4])+ 1/3 * neighboursVel;
    GPxhatOut(i,3:4) = updatedVel;
    GPxhatOut(i,:)= CheckOutboundDirection(GPxhatOut(i,:)',GPxhatIn(i,:)',GPxhatIn(leftNeighbour,:),GPxhatIn(rightNeighbour,:),originY,originX);
    %check for maximum speed
    speed = norm(GPxhatOut(i,[3 4]));
    %max_speed = 0.1;
    %normalized vector is vector with a magnitude of 1. This ensures that we're not messing with the
    %direction of the vector, only its magnitude.
    if speed > max_speed
        GPxhatOut(i,[3 4]) = GPxhatOut(i,[3 4])*max_speed/speed;
    end
    %     if debugFlag
    %         figure(1)
    % %         xlim([-2 35]);
    % %         ylim([-15 35]);
    % %         plot(GPxhatOut(i,1), GPxhatOut(i,2), 's');
    %         hold on;
    %     end
end
% if debugFlag
%     figure(1)
%     plot( [ GPxhatOut(N_GP,1) ; GPxhatOut(:,1)],[GPxhatOut(N_GP,2); GPxhatOut(:,2)]);
% end
% if ~debugFlag
%     figure(1)
%     for i = 1:T
%         plot( [ GPxhatOut(N_GP,1,i) ; GPxhatOut(:,1,i)],[ GPxhatOut(N_GP,2,i); GPxhatOut(:,2,i)]);
%     end
% end


