function [xhat,Pout,Cov_a_b,sig2_Z,distCP2Obs] = UpdateFGEstimator(t,delT,origin,x,z,F,Pin,Q,H, R,detectStatus,MLE_flag,MAP_flag,QKF_flag,QKF_OptK_flag,avoidPredFlag,avoidEstFlag,passedBoundaryTimetag,limitCalcSpreadRateFlag,fixed_direction_enableFlag,decayFactor,corr_corr_flag,MaxSpreadRate,MinSpreadRate,MaxAffectedDist,MinAffectedDist,MaxCP_Heading,MinCP_Heading,MinLocError,DEBUG_FLAG,REC_CNT,TestedCFG)
%%% Matlab script to generate some "true" data for later assessment
%%% Generates:
%%% x: the state history which evolves according to
%%% x(k+1) = Fx(k) + w(k)
%%% w: the process noise history (randomly generated)
%%% z: a set of observations on the state corrupted by noise
%%% v: the noise on each observation (randomly generated)
%MaxSpreadRate = 0.2;

% Uncertainty with measurments
%set detected status
DetectedStatus=detectStatus;
%%define the measurements status
LABELED_IN = 1;
LABELED_OUT = 0;

curMeas = z;
%In case we have measurments in the same time step we should avoid
%prediction phase to prevent adding process noise in the same time step
%multiple times
% P_prior=P;
% xhat_prior = x;

if avoidPredFlag==1
    xpred = x;
    Ppred = Pin;
else
    [xpred, Ppred] = predict(x, Pin, F, Q);
end

P_prior=Ppred;
xhat_prior = xpred;


[T,psi] = Calc2DTransformationMatrix([xpred(1) xpred(2)]',curMeas);
%evaluate the distance for translation of the measurements along the new
%coordinate system
dist_OBS_CPpred = sqrt((curMeas(1)-xpred(1))^2 + (curMeas(2)-xpred(2))^2);
%debug
persistent cnt;

% Upon the first call we need to initialize the variables.
if isempty(cnt)
    
    cnt = 0;
end

% Compute the accumulated sum and the number of values so far.
cnt= cnt+1;


vel_a_b  =  T* xpred([3 4]);%TODO test

if DetectedStatus==LABELED_IN
    %TODO: change back to vel_a_b(2)
    z_a_b = [dist_OBS_CPpred 0 vel_a_b(1) 0]';% a,b coordinate system - 1D
else
    z_a_b = [dist_OBS_CPpred 0 vel_a_b(1) 0]';% a,b coordinate system - 1D
end
Cov_a_b = zeros(4,4);
BigT = zeros(4,4);
BigT([1 2],[1 2]) = T;
BigT([3 4],[3 4]) = T;
%transform error covariance
Cov_a_b =  BigT*Ppred*BigT';%TODO test
%%
%set location error to minimum in case of crossing the periphery
if passedBoundaryTimetag==t
    %%%heuristic solution
    %if corr_corr_flag
        corr_factor = sqrt(MinLocError^2/(Cov_a_b(1,1)));%ratio between the estimated and prior
        Cov_a_b(1,[2,3,4]) = Cov_a_b(1,[2,3,4])*corr_factor;% scale up all the associated entries
        Cov_a_b([2,3,4],1) = Cov_a_b(1,[2,3,4]);% set first column
    %end
    Cov_a_b(1,1) = MinLocError^2;%only along one axes  (Cov_aa)
    Ppred =  BigT'*Cov_a_b*BigT;
    xpred(1:2) = curMeas(1:2);
    z_a_b(1) = 0;%the distance between is zero at crossing time
    z_a_b(3) = 0;%the spread rate zeroed at crossing time

    P_prior=Ppred;
    xhat_prior = xpred;
end
%%
%SR
% ind =(find(abs(Cov_a_b)< 0.0001));
% if(~isnan(ind))
%     Cov_a_b(ind) =0;
% end
if 1%DEBUG_FLAG
    figure (98)
    hold on;
    if ~all(eig(Cov_a_b([1 2],[1 2]))>0)
        warning('The covariance matrix (Cov_a_b) must be positive definite (it has non-positive eigenvalues) (@time=%d)',t);
    else
        error_ellipse(Cov_a_b([1 2],[1 2]),z_a_b([1 2]));
    end
    plot(0,0,'>');
    plot(z_a_b(1),z_a_b(2),'sq');
    
    line([z_a_b(1),0],[z_a_b(2),0]);
    axis equal;
end

% Substitute the variance along the connecting line (1D)
sig2_Z = Cov_a_b(1,1);

cp_origVec = xpred([1 2]) - origin;

% Update the mean value of the measurment in a-b coordinate system
if (dist_OBS_CPpred >= MinAffectedDist  &&  dist_OBS_CPpred < MaxAffectedDist)
    % Set the higher limit of current interval
    
    %The first step for determining the conditional mean and covariance
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Conditional mean of the measurment vector
    expectedOutDir = atan2d(cp_origVec(2),cp_origVec(1));
    %angle of vector connect control point and current observation location
    %(x,y) coordinate system
    DP_CP_vec = xpred(1:2) - curMeas(1:2)';% vector connecting observation (DP) and CP = Vec_rel
    relSpatialAngel = AngleBitweenVectors(cp_origVec,DP_CP_vec);%angle between Vec_CP and Vec_rel
    behind_originFlag = abs(AngleBitweenVectors(cp_origVec,curMeas)) > MaxCP_Heading;%check if the angle between vectors is bigger than MaxCP_Heading

    firstCondFlag = abs(relSpatialAngel)>=MaxCP_Heading; % 0 for non conflict in case it detects IN
    secondCondFlag = wrapTo360(relSpatialAngel)< MaxCP_Heading; % 0 for non conflict in case it detects OUT
    firstConflictFlag = (DetectedStatus)&(firstCondFlag);% in case drone observe burned but it is not expected
    secondConflictFlag = (~DetectedStatus)&(secondCondFlag)&(~behind_originFlag);% in case drone observe unburned but it is not expected
    %     relSpatialAngel = atan2d(xpred(2)-curMeas(2),xpred(1)-curMeas(1));
    %Wrap angle in degrees to [-180 180]
    %     relTotalAngel =(expectedOutDir-relSpatialAngel);%wrapTo180(relSpatialAngel-spreadAngel);%abs(spreadAngel);
    %     firstCondFlag = abs(expectedOutDir)>=MaxCP_Heading; %the expected angle toward spread direction
    %     secondCondFlag = abs(relSpatialAngel)>=MaxCP_Heading;%the angle toward observation
    %     firstConflictFlag = (~DetectedStatus)&(firstCondFlag==secondCondFlag);% in case drone observe unburned but it is not expected
    %     secondConflictFlag = (DetectedStatus)&(((~firstCondFlag)&(secondCondFlag))|((firstCondFlag)&(~secondCondFlag)));% in case drone observe burned area but it is not expected
    conflictFlag = firstConflictFlag | secondConflictFlag;
    if  ~firstConflictFlag && ~secondConflictFlag% the periphery is ahead of the current measurment
        % Set the higher limit of current interval
        low_interval_limit=0;%z_a_b(1);
        high_interval_limit= 10000;% xpred(1);
        
        %%conflicted cases
        %could be heuristic that have a common sense
    elseif firstConflictFlag || secondConflictFlag% the periphery is behind of the current measurment
        %the predicted periphery grid point move along the connected line
        % but since the measurement status is unburned it implies that the burned part is much
        % more probable to be backward
        % Set the higher limit of current interval
        low_interval_limit=-10000;%the lowest limit of the scalar function is related to the variance along the line
        high_interval_limit= 0;%2*z_a_b(1);%hueristic method: keep it symetric to the first case (cutting tail)
        warning('Measurments are conflicted with expectation (@time=%d)',t);

    else
        warning('Measurments is probably behind the second side of the periphery (@time=%d)',t);
        
        xhat= xpred;
        Pout =  Ppred;
        distCP2Obs = sqrt((curMeas(1)-xhat(1))^2 + (curMeas(2)-xhat(2))^2);
        
        return;
    end
else
    warning('Measurments is not in the affected distance (@time=%d)',t);
    xhat= xpred;
    Pout =  Ppred;
    distCP2Obs = sqrt((curMeas(1)-xhat(1))^2 + (curMeas(2)-xhat(2))^2);
    
    return;
    
end


[EzA , Cov_zA] = CalcScalarExpectation(t,low_interval_limit,high_interval_limit,z_a_b(1),sqrt(sig2_Z),MLE_flag,MAP_flag,QKF_flag,DEBUG_FLAG,REC_CNT,TestedCFG);%the function get sigma and not sigma^2
%Conditional covariance of the measurment vector
%[Cov_zA] = CalcScalarQuantCov(low_interval_limit,high_interval_limit,z_a_b(1),sqrt(sigZ));%the function get sigma and not sigma^2
%TODO: We check if the expected and variance is not valid for next step of
%the estimator. If it is INVALID we will skip over the quantized update
if isnan(EzA) || isinf(EzA) || isnan (Cov_zA) || isinf(Cov_zA)
    warning('INVALID expected value and variance (@time=%d)',t);
    xhat= xpred;
    Pout =  Ppred;
    distCP2Obs = sqrt((curMeas(1)-xhat(1))^2 + (curMeas(2)-xhat(2))^2);
    
    return;
end
%%INVALID variance
if Cov_zA <= 0
    warning('INVALID variance (@time=%d)',t);
    Cov_zA = sig2_Z;
    xhat= xpred;
    Pout =  Ppred;
    distCP2Obs = sqrt((curMeas(1)-xhat(1))^2 + (curMeas(2)-xhat(2))^2);
    
    return;
end
%TODO: fix the problem with Cov_a_b(1,2) that get very big on close distance
% if Cov_zA <= 100
%     warning('limited variance (@time=%d)',t);
%     Cov_zA = 20;
% end


%Plot for Debug
%PlotPdf(EzA,Cov_zA,'b');

% Use the Kalman Filter to determine the conditional mean and covariance of
% the desired state vector
%%
%Evaluate the innovation of the
if QKF_flag
    
    [nu, S] = quant_innovation(EzA, Cov_a_b,z_a_b, H(1,1:4), R(1,1),Cov_zA,QKF_OptK_flag);
    [Ex_zA, Pz_a_b , K] = quant_innovation_update(z_a_b,Cov_a_b, nu, S, H(1,1:4));
    %The second step for determining the conditional mean and covariance
    Cov_x_zA = Pz_a_b + K*Cov_zA*K';
    %Back  Translation & Transformation
    xhat = xpred;
    %%TODO: decide what is the best way to translate back
    %xhat = xpred + BigT'*(Ex_zA - z_a_b) * decayFactor; %Pull mean value over to measurment;%TODO test
    xhat([1 2]) = xpred(1:2) + T'*(Ex_zA([1 2]) - z_a_b([1 2])) * decayFactor;
    xhat([3 4]) =  T'*Ex_zA([3 4]);%TODO test
    figure(99)
    hold on; grid on;
    plot(t,psi,'b.');
    %%%heuristic solution
    if corr_corr_flag
        corr_factor = sqrt(Cov_zA/(Cov_a_b(1,1)));%ratio between the estimated and prior
        Cov_x_zA(1,[2,3,4]) = Cov_x_zA(1,[2,3,4])*corr_factor;% scale up all the associated entries
        Cov_x_zA([2,3,4],1) = Cov_x_zA(1,[2,3,4]);% set first column
    end
    %%%
    Pout = BigT'*Cov_x_zA*BigT;%TODO test
    if ~all(eig(Pout)>0)
        warning('The covariance matrix must be positive definite (it has non-positive eigenvalues) (@time=%d)',t);
    end
    
    %%%%%%%
elseif MLE_flag %Fitted CDF
    xhat_ml = zeros(4,1);
    xhat_ml(1:2) = xpred(1:2) + T*([EzA 0 ]' - z_a_b(1:2)) * decayFactor ; % translation
    xhat_ml(3:4) = xpred(3:4);
    corr_factor = Cov_zA/sqrt(Cov_a_b(1,1));%ratio between the estimated and prior
    Cov_a_b(1,1) = Cov_zA;%the CDF fit returns sigma
    if corr_corr_flag
        Cov_a_b(1,[2,3,4]) = Cov_a_b(1,[2,3,4])*corr_factor;% scale up all the associated entries
        Cov_a_b([2,3,4],1) = Cov_a_b(1,[2,3,4]);% set first column
    end
    P_ml = BigT*Cov_a_b*BigT';%transformation
    %set the approximated maximum likelihood estimator as the one for next
    %iteration
    Pout= P_ml;
    xhat=xhat_ml;
    %DEBUG
    %     figure(95); hold on; grid on;
    %         plot(t,xhat_ml(1),'.r',t,xhat_prior(1),'.b','MarkerSize',12);
    %         figure(96); hold on; grid on;
    %         plot(t,P_ml(1,1),'.r',t,P_prior(1,1),'.b',t,P_ml(3,3),'.k','MarkerSize',12);
    
    if MAP_flag % in case we combine prior with likelihood to have the a posterior result
        %xhat_p prior P_p prior covariance
        [nu, S] = innovation(xhat_prior, P_prior,xhat_ml(1:2), H(1:2,1:4), P_ml(1:2,1:2));%xpred, Ppred, z, H, R
        [xhat_map, P_map , K] = innovation_update(xhat_prior,P_prior, nu, S, H(1:2,1:4));%xpred, Ppred, nu, S, H
        %                 xhat_map = (P_ml/(P_ml+P_prior))*xhat_prior + (P_prior/(P_ml+P_prior))*xhat_ml;
        %                 P_map = (P_ml*P_prior) / (P_ml+P_prior);
        %         P_map = (P_ml(1,1)*P_prior(1,1)) / (P_ml(1,1)+P_prior(1,1));
        %DEBUG
        %         figure(97); hold on; grid on;
        %         plot(t,xhat_map(1),'.r',t,xhat_prior(1),'.b',t,xhat_ml(1),'.k','MarkerSize',12);
        %         figure(98); hold on; grid on;
        %         plot(t,P_map(1,1),'.r',t,P_prior(1,1),'.b',t,P_ml(1,1),'.k','MarkerSize',12);
        %
        Pout= P_map;
        xhat=xhat_map;
        
    end
else %test optimization for MAP or MLE
    %from QKF
    [nu, S] = quant_innovation(EzA, Cov_a_b,z_a_b, H(1,1:4), R(1,1));
    [Ex_zA, P , K] = quant_innovation_update(z_a_b,Cov_a_b, nu, S, H(1,1:4));
    %The second step for determining the conditional mean and covariance
    Cov_x_zA = P + K*Cov_zA*K';
    
    xhat_opt = zeros(4,1);
    %from MLE And MAP
    %     xhat_opt (1:2) = xpred(1:2) + T*([EzA 0 ]' - z_a_b(1:2)) * decayFactor ; % translation
    %     xhat_opt (3:4) = xpred(3:4);
    %from QKF
    xhat_opt = xpred + BigT*(Ex_zA - z_a_b) * decayFactor; %Pull mean value over to measurment
    xhat_opt([3 4]) =  T*Ex_zA([3 4]);
    
    %     corr_factor = Cov_zA/sqrt(Cov_a_b(1,1));%ratio between the estimated and prior
    %     Cov_a_b(1,1) = Cov_zA^2;%the CDF fit returns sigma
    %     if corr_corr_flag
    %         Cov_a_b(1,[2,3,4]) = Cov_a_b(1,[2,3,4])*corr_factor;% scale up all the associated entries
    %         Cov_a_b([2,3,4],1) = Cov_a_b(1,[2,3,4]);% set first column
    %     end
    %transformation of
    %for MLE And MAP
    %P_opt  = BigT*Cov_a_b*BigT';
    %for QKF
    P_opt  = BigT'*Cov_x_zA*BigT;
    
    %set the approximated maximum likelihood estimator as the one for next
    %iteration
    Pout= P_opt ;
    xhat=xhat_opt ;
end

%direction between grid point velocity vector and tangent line
relVelAngleInDegrees = AngleBitweenVectors(xhat([3 4]),cp_origVec);
if fixed_direction_enableFlag
if abs(relVelAngleInDegrees)>= 40 %aspect angle of change is limited in between 
    %warning('Avoid update velocity (@time=%d)',t);
    curSpreadRate = norm(xpred([3 4]));
    estSpreadRate = norm(xhat([3 4]));
    proj_est = dot(xpred([3 4]),xhat([3 4]))/curSpreadRate^2 * xpred([3 4]); %projection of estimation along fix cp_origVec
    xhat([3 4]) = proj_est;%avoid motion
    %keep the same location as before
end
end
%jump allowed only when conflict between expected to observation comes up
%this condition suppose to prevent position jump in case of estimated
%velocity change in direction
% if ~conflictFlag && (abs(relVelAngleInDegrees)> 100)
%     warning('Avoid update state vector (@time=%d)',t);
%     %xhat = xpred;%avoid jump
%     %     Pout =  Ppred;
% end
% avoid updating position when it skips the origin of fire
relPosAngleInDegrees = AngleBitweenVectors(xhat([1 2]),xpred([1 2]));
if (abs(relPosAngleInDegrees)> 170)
    %warning('Avoid update position (@time=%d)',t);
    xhat = xpred;%avoid jump
end
curSpreadRate = norm(xhat([3 4]));
calcSpreadRate = sqrt(((x(1) - xhat(1))/delT).^2+((x(2) - xhat(2))/delT).^2);

%limit estimated spread rate
if curSpreadRate~=0 %avoid zero divisor
    if curSpreadRate > MaxSpreadRate
        %warning('spread rate reach the maximum (@time=%d)',t);
        xhat([3 4]) = xhat([3 4])*MaxSpreadRate/curSpreadRate;
    elseif curSpreadRate < MinSpreadRate
        xhat([3 4]) = xhat([3 4])*MinSpreadRate/curSpreadRate;
        warning('spread rate reach the minimum (@time=%d)',t);
    end
end

%limit calculated spread rate
if limitCalcSpreadRateFlag
  if calcSpreadRate > MaxSpreadRate
    %warning('Calculated spread rate is above limiter (@time=%d)',t);    
      %psi = atan2((xhat(2) - x(2)),(xhat(1) - x(1)));
      psi = atan2(xhat(4),xhat(3));
      xhat([1 2]) = x([1 2]) + [ MaxSpreadRate*cos(psi); MaxSpreadRate*sin( psi)]* delT;
  end
end

%No estimation for the last cycles after crossing the boundary
if avoidEstFlag
    warning('External overide to avoid update state and covariance (@time=%d)',t);
    xhat= xpred;
    Pout =  Ppred;
    
    Pout = BigT'*Cov_a_b*BigT;%TODO test
    
end
distCP2Obs = sqrt((curMeas(1)-xhat(1))^2 + (curMeas(2)-xhat(2))^2);


% if DEBUG_FLAG
%     figure(99)
%     subplot(3,1,1);
%     hold on
%     plot(cnt,xhat(1),'.-r');
%     plot(cnt,xhat(2),'.-b');
%     xlabel('iter');ylabel('pos(x_r,y_b)');
%     subplot(3,1,2);
%     hold on
%     plot(cnt,xhat(3),'.-r');
%     plot(cnt,xhat(4),'.-b');
%     plot(cnt,norm(xhat([3 4])),'.-g');
%     xlabel('iter');ylabel('vel(x_r,y_b,tot_g)');
%     subplot(3,1,3);
%     hold on
%     plot(cnt,dist,'.-b');
%     plot(cnt,sqrt(sigZ),'.-g');
%     if QKF
%     plot(cnt,Ex_zA(1)-EzA,'.-r');
%     end
%     xlabel('iter');ylabel('\sigma (g), dist(d_b),inovation(a_r)');
%
%     %plot(cnt,,'.-b');
% end







