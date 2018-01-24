function UT_UpdateQKF1x1TrajReshape(varargin)
if nargin==0
    clear all;
    close all;
    testCase = 1;
elseif nargin==1
    testCase = varargin{1};
else
    error('Incorrect number of input arguments.')
end
TestedCFG = 'FGE_1CP_mObs_test';%fitted Gaussian Estimator
TestedMainTitle = 'Fitted Gaussian Estimator (QKF)';
recFlag =1;
DEBUG_FLAG =0;
REC_CNT =0;%the record counter is used to initiating the time step to record the resulted calculation
plottools('on') ;plotbrowser('off'); set(0,'DefaultFigureWindowStyle','docked')

%To create a 95% confidence ellipse from the 1$\sigma$ error ellipse, we must enlarge it by a factor of 2.4477
ellipse_factor = 2.4477;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Setup and initialization %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 25;
delT = 1;
MaxAffectedDist = 2000;
MinAffectedDist = 100;
MinLocError = 20;
MaxSpreadRate = 6;% 22kmh ~ 6m/s
MinSpreadRate = 0;
NominalSpreadRate = 3;%10.8 kmh ~ 3m/s
droneNominalSpeed = 20;% 20m/sec
MaxCP_Heading = 90;
MinCP_Heading = 90;
sigma2Q = 0.01;
sigma2R = 0.1;
var_sys = sigma2Q;
var_obs = sigma2R;
% Setting the random seed, so the same example can be run several times
s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s);
%% Set up the movie.
if recFlag
    writerObj = VideoWriter(['.\figures\',TestedCFG,'Simulation',num2str(8),'.avi']); % Name it.
    writerObj.FrameRate = 2; % How many frames per second.
    open(writerObj);
end
%% Setup the process model (transition matrix, observation matrix, process noise and observation noise)
F = [1 0 delT 0
    0 1 0 delT
    0 0 1 0
    0 0 0 1 ];
H = [1 0 0 0;0 1 0 0;0 0 0 0;0 0 0 0]; % measurement matrix
Qfactor = 1; % process noise mult factor
Rfactor = 0; % measurement noise mult factor
sigmaQ = Qfactor*sqrt(0.1);
sigmaR = Rfactor*sqrt(0.1);
Q_NoProcessNoise = 0*sigmaQ^2 * [ 1 0% process noise covariance matrix
    0 1];
Q = sigmaQ^2 * [ 1 0 0 0% process noise covariance matrix
    0 1 0 0
    0 0 1 0
    0 0 0 1];

P = [];
P(:,:,1) = 10*Q;

%% Initialize states vectors and parameters
x = zeros(4,N);% [ x y  rate]
w = zeros(4,N);
Cov_a_b =[];
% z = zeros(2,N);% [ x y status] the status is quantization of the measurments -1 for inside the convex hull, +1 for outside of the convex hull
v = zeros(2,N);
%number of cycles to avoid estimation after crossing actual boundary
avoidEstAftCrossDelay = 1;
%flag to enable limiter on calculated spread rate
limitCalcSpreadRateFlag = 1;
%flag to enable fixed line motion of the control point
fixed_direction_enableFlag = 1;
% propotional gain on the innovation of expected value
decayFactor = 1;% the size of innovation we allowed to the estimator given the observation
MAP_flag=0;
MLE_flag=0;
corr_corr_flag = 0;%enable correlation correction based on calculated ratio
% flag to indicate CP updated due to drone crossing actual boundary
updateOnCrossingFlag=0;
% flag to enable CP update feature
enableUpdateOnCrossingFlag=0;
%% Setup the scenario (Fire origin, actual periphery, actual spread rate, fleet deployment)
origin = [0 ;0];
actualBoundary = [360;0];%after 30 minutes of fire start with nominal spread rate of 0.2 m/sec (0.2*30*60)
CP_Behind = actualBoundary(1)-180;% if the uncertianty of the spread rate is going down to 0.1 (0.1*30*60)
CP_Infront = actualBoundary(1)+180;% if the uncertianty of the spread rate is going up to 0.3 (0.3*30*60)
%process noise is set to be around the +- point of the initial location
Phi_s = 0.5;%error in spread rate
Q = Phi_s*[delT^3/3 0 delT^2/2 0; 0 delT^3/3 0 delT^2/2; delT^2/2 0 delT 0; 0 delT^2/2 0 delT];
%Q = Phi_s*eye(4);

%Observations noise

R = 20^2 * [ 1 0;0 1];

%start with high uncertainty of the prior
%1000 meter of standard deviation and 2 meter per second in spread rate
P(:,:,1) = [1000^2 0 0 0; 0 1000^2 0 0; 0 0 MaxSpreadRate^2 0; 0 0 0 MaxSpreadRate^2];
%common observation locations for first set of unit tests
drone_initial_location =[-1000+CP_Behind 0];
distObs2Bound(1) = sqrt((drone_initial_location(1)-actualBoundary(1,1))^2 + (drone_initial_location(2)-actualBoundary(2,1))^2);
time_to_go = distObs2Bound(1)/(droneNominalSpeed - NominalSpreadRate);
%% Setup a single observer (default)
z1(1,:)= [drone_initial_location:droneNominalSpeed*delT:time_to_go*droneNominalSpeed+500];%approaching from the left
N= length(z1);
z1(2,:) = zeros(1,N);
z1(3,:) = zeros(1,N);

switch testCase
    case 1
        %CASE I: drone on BURNED area approaching MOVING control point from the left
        MLE_flag = 0;
        MAP_flag = 0;
        QKF_flag = 1;
        QKF_OptK_flag = 0;
        TestedCFG = 'FGE_1CP_mObs_test';%fitted Gaussian Estimator
        TestedMainTitle = 'Fitted Gaussian Estimator (QKF)';
        set(0,'defaultTextInterpreter','latex'); %trying to set the default
        fig_ct = 1;
        %%%%Initialize observation configuration
        OnBurnedArea=1;%set detected status
        NumberOfObserver = 1; %set the number of observer for interlacing to one sequential measurement
        staticTraj1Flag = 0;%set first observer to a fixed location
        staticTraj2Flag = 0;%set first observer to a fixed location
        staticTraj3Flag = 1;%set first observer to a fixed location
        
        if staticTraj1Flag
            z1(1,:)= (-1000+CP_Behind)*ones(1,N);%static observer
        end
        z1(2,:)= 10*ones(1,N);
        %set compass direction of each drone (e.g. faced toward CP or with
        %compass angle)
        %the special case where the major axis of the ellipse is aligned with the X-axis
        theta1 = -30;
        theta2 = 30;
        theta3 = 10;
        
        T1 = [cosd(theta1) sind(theta1) ;-sind(theta1) cosd(theta1)];
        T2 = [cosd(theta2) sind(theta2) ;-sind(theta2) cosd(theta2)];
        T3 = [cosd(theta3) sind(theta3) ;-sind(theta3) cosd(theta3)];
        
        z1([1:2],:) = (T1*z1([1:2],:));%rotate 45 deg and test from different angle
        z1(3,:) = ones(1,N)*OnBurnedArea;
        
        if NumberOfObserver==2
            %set second observer
            %static observer - (TODO: comment out for dynamic)
            if staticTraj2Flag
                z2(1,:) = (1000+CP_Behind)*ones(1,N);
            else
                %dynamic observer case
                z2(1,:)= [(1000+CP_Behind):(-droneNominalSpeed*delT):-(time_to_go*droneNominalSpeed)-140];
            end
            z2(2,:)= 10*ones(1,N);
            z2([1:2],:) = (T2*z2([1:2],:));%rotate 45 deg and test from different angle
            z2(3,:) = (~OnBurnedArea)*ones(1,N);
            % interleaves two same sized matrices by row
            z1=z1';
            z2=z2';
            %         col_interleave = reshape([z1(:) z2(:)]',2*size(z1,1), []);
            % Note that the reshape requires that a and b be the same size.
            z = zeros(N,3,NumberOfObserver);
            z(:,:,1)= [z1];% col_interleave';
            z(:,:,2)= [z2];
        elseif NumberOfObserver==3
            %set second observer
            if staticTraj2Flag
                z2(1,:) = (1000+CP_Behind)*ones(1,N);
            else
                %dynamic observer case
                z2(1,:)= [(1000+CP_Behind):(-droneNominalSpeed*delT):-(time_to_go*droneNominalSpeed)-140];
            end
            z2(2,:)= 10*ones(1,N);
            %rotate trajectory with selected heading angle
            z2([1:2],:) = (T2*z2([1:2],:));%rotate theta [deg] and test from different angle
            z2(3,:) = (~OnBurnedArea)*ones(1,N);
            %set third observer
            if staticTraj3Flag
                z3(1,:) = (1000+CP_Behind)*ones(1,N);
            else
                %dynamic observer case
                z3(1,:)= [(1000+CP_Behind):(-droneNominalSpeed*delT):-(time_to_go*droneNominalSpeed)-140];
            end
            z3(2,:)= 10*ones(1,N);
            %rotate trajectory with selected heading angle
            z3([1:2],:) = (T3*z3([1:2],:));%rotate theta [deg] and test from different angle
            z3(3,:) = (~OnBurnedArea)*ones(1,N);
            % Note that the reshape requires that a and b be the same size.
            z = zeros(N,3,NumberOfObserver);
            z(:,:,1)= [z1'];% col_interleave';
            z(:,:,2)= [z2'];
            z(:,:,3)= [z3'];
        else
            z = zeros(N,3,1);%default 2 to avoid indexes problems
            z1=z1';
            z(:,:,1)= z1;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Uncertainty with measurments
        xhat = zeros(4,N); % state estimate
        xhat(1,1) = CP_Behind;%initial x position
        xhat(2,1) = 10;%initial y position
        xhat(3,1) = NominalSpreadRate;%initial vx velocity
        xhat(4,1) = 0;%initial vy velocity
        plotDir = 'r-->';
        heading = pi/2;
        actualBorderVel = [NominalSpreadRate;0];
        %rotate initial CP and boundary and actual spread rate
        if 0
            xhat([1:2],1) = (T1*xhat([1:2],1));%rotate 45 deg and test from different angle
            xhat([3:4],1) = (T1*xhat([3:4],1));%rotate 45 deg and test from different angle
            actualBoundary = T1*actualBoundary;
            actualBorderVel = T1*actualBorderVel;
        end
        
        avoidPredFlag = 0;
        avoidEstFlag = 0;
        avoidEstAftCrossFlag = 0;
        enableUpdateOnCrossingFlag =1;
        enableTrajReshape = 1;
        
        TestedMainSubTitle = ' - Single CP and Fixed Observers ( Predict-';
        if avoidPredFlag
            TestedMainSubTitle = [TestedMainSubTitle, 'OFF'];
        else
            TestedMainSubTitle = [TestedMainSubTitle, 'ON'];
        end
        if avoidEstFlag
            TestedMainSubTitle = [TestedMainSubTitle, ', Estimation-OFF'];
        else
            TestedMainSubTitle = [TestedMainSubTitle, ', Estimation-ON'];
        end
        if xhat(1,1) == CP_Behind
            TestedMainSubTitle = [TestedMainSubTitle, ', CP-Behind) '];
        else
            TestedMainSubTitle = [TestedMainSubTitle, ', CP-Ahead) '];
        end
        
        
end

%z(1,:)= [-10:1:N-11];

%z(1,:)= [0:1:N-1];


h = [];
avoidEstIndicator = avoidEstFlag*ones(NumberOfObserver,1);
avoidPredIndicator = avoidPredFlag*ones(NumberOfObserver,1);
avoidEstAftCrossIndicator = avoidEstAftCrossFlag*ones(NumberOfObserver,1);
%% Initialize recording vector after selected configuration (unit-test)
distCP2Obs = inf*ones(NumberOfObserver,N);%each observer in a seperate row
distCP2Bound = inf*ones(1,N);%single CP
distObs2Bound= inf*ones(NumberOfObserver,N);%each observer in a seperate row

distCP2Bound(1) = sqrt((actualBoundary(1,1)-xhat(1,1))^2 + (actualBoundary(2,1)-xhat(2,1))^2);
distCP2Obs(:,1) = sqrt((z(1,1,1:NumberOfObserver)-xhat(1,1)).^2 + (z(1,2,1:NumberOfObserver)-xhat(2,1)).^2);
distObs2Bound(:,1) = sqrt((z(1,1,1:NumberOfObserver)-actualBoundary(1,1)).^2 + (z(1,2,1:NumberOfObserver)-actualBoundary(2,1)).^2);

P([1 2],[1 2],1) =T1'*P([1 2],[1 2],1)*T1;
P = fix_nonPosCov(P);
Cov_a_b = P(:,:,1);
sigma2_z(1) = P(1,1,1);
%running configuration
once=1;
passedBoundaryFlag = zeros(1,NumberOfObserver);
passedBoundaryTimetag = inf*ones(1,NumberOfObserver);
%%Init measurement
GPxhat=[];
[GPxhat] = UpdateSimPropPolyMGPs(1,N,GPxhat,actualBoundary(1),NominalSpreadRate);
for i=2:N
    %Update complete actual periphery
    [GPxhat] = UpdateSimPropPolyMGPs(i,N,GPxhat,NominalSpreadRate);
    
    %Update boundary
    actualBoundary(:,i) = actualBoundary(:,i-1) + actualBorderVel*delT;
    %save prior and posterior between observations
    x_prior = xhat(:,i-1);
    x_est = xhat(:,i-1);
    P_prior = P(:,:,i-1);
    P_est = P(:,:,i-1);
    avoidPredInSeqFlag = avoidPredFlag;% start any sequence with predict active
    for j=1:NumberOfObserver
        x_prior = x_est;% the same value for the begining of each sequence of different observer, but  take the last estimated to improve with next observer
        P_prior = P_est;
        % check for trajectory reshape mode
        if enableTrajReshape
            [cov_angle x_tip]=eval_cov_angle(P_prior([1 2],[1 2]),x_prior,[z(i,1,1),z(i,2,1)],~DEBUG_FLAG);
            obs_heading_angle = atan2(z(i,2,1)-z(i-1,2,1),z(i,1,1)-z(i-1,1,1));
            pointA = [ z(i,1,1), z(i,2,1),   obs_heading_angle ];
            pointB = [ x_tip(1), x_tip(2), pi+cov_angle ]; %add 180 to head toward the eigen vector
            TurnRadius = 100;
            PathStep = droneNominalSpeed*delT;
            %tic
            updated_traj = dubins_curve(pointA,pointB, TurnRadius, PathStep,DEBUG_FLAG);%DEBUG_FLAG~quiet
            if i<N
                z(i+1,1:2,j) = updated_traj(2,1:2);%next point
            end
            figure(101)
            hold on;grid on;
            plot(i, cov_angle,'.');
            
            curMeas = z(i,1:3,j);
        else
            curMeas = z(i,1:3,j);%sequence of measurements from both side of the grid point
        end
        [in,on] = inpolygon(curMeas(1),curMeas(2),GPxhat(:,1,i),GPxhat(:,2,i));
        z(i,3,j) = in;
        OnBurnedArea = in;%curMeas(3);
        %Update tracking status
        distObs2Bound(j,i) = sqrt((curMeas(1)-actualBoundary(1,i))^2 + (curMeas(2)-actualBoundary(2,i))^2);
        boundTol = 20;
        if distObs2Bound(j,i) < boundTol && (~passedBoundaryFlag(j)) % passing by the  boundary from the left
            OnBurnedArea= ~OnBurnedArea;%set detected status
            z(i:N,3,j)= ones(1,N-i+1)*OnBurnedArea;% z is the sequential measuerement
            passedBoundaryFlag(j) = 1;
            passedBoundaryTimetag(j) = i;
        end
        
        if enableUpdateOnCrossingFlag && passedBoundaryFlag(j) &&~updateOnCrossingFlag
            %x_prior([1 2]) = curMeas([1 2]);
            %TODO set the spread rate...
            updateOnCrossingFlag =1;%indicate that CP location has been updated
        end
        %Avoid estimation after passing by boundary for several cycles TBD
        if passedBoundaryFlag(j) == 1 && i < (passedBoundaryTimetag(j) + avoidEstAftCrossDelay)
            avoidEstIndicator(j) = 1;
        end
        % the following condition is to prevent estimation for couple cycles
        % after drone croses the boundary and only in case we want to support
        % that feature in the current unit test
        if (i > (passedBoundaryTimetag(j) + avoidEstAftCrossDelay) ) && (~(avoidEstAftCrossIndicator(j)))
            avoidEstIndicator(j) = 0;
        end
        %initiate timestep to record the resulted computation
        %     if i==2
        %         REC_CNT=2;
        %     else
        %         REC_CNT=0;
        %     end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%MAIN CALL
        [x_est,P_est,Cov_a_b(:,:,i),sigma2_z(i),distCP2Obs(j,i)] = UpdateFGEstimator(i,delT,origin,x_prior,curMeas,F,P_prior,Q,H,R, OnBurnedArea,MLE_flag,MAP_flag,QKF_flag,QKF_OptK_flag,avoidPredInSeqFlag,avoidEstIndicator(j),passedBoundaryTimetag(j),limitCalcSpreadRateFlag,fixed_direction_enableFlag,decayFactor,corr_corr_flag,MaxSpreadRate,MinSpreadRate,MaxAffectedDist,MinAffectedDist,MaxCP_Heading,MinCP_Heading,MinLocError,DEBUG_FLAG,REC_CNT,TestedCFG) ;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        avoidPredInSeqFlag = 1;% avoid predict in-between the sequence of measurements
        
        
    end
    % update all sub-sequence at once
    xhat(:,i)=x_est;
    P(:,:,i) = P_est;
    %evaluate distance between prediction and actual for analysis of perfomance
    distCP2Bound(i) = sqrt((actualBoundary(1,i)-xhat(1,i))^2 + (actualBoundary(2,i)-xhat(2,i))^2);
    distCP2Obs(:,i) = sqrt((z(i,1,1:NumberOfObserver)-xhat(1,i)).^2 + (z(i,2,1:NumberOfObserver)-xhat(2,i)).^2);
    
    if NumberOfObserver>1 && distCP2Obs(2,i)<10
        warning('check scenario! (@time=%d)',i);
    end
end
% figure(1)
%
% legend(h,'Quantized Estimation (Covariances)','Sequence of Remote Measurements (Drone Trajectory)','Quantized Estimations (Expectations)');
% axis equal;
%plot for analysis
figure(2)
subplot(3,1,1);
hold on;
plot(1:N,xhat(1,1:N),'.-b');
plot(1:N,xhat(2,1:N),'--b');
plot(1:N,actualBoundary(1,1:N),'.-r');
plot(1:N,actualBoundary(2,1:N),'--r');

xlabel('time-step [s]','Interpreter','Latex');ylabel('pos [m]','Interpreter','Latex');
legend('x','y');
subplot(3,1,2);
hold on;
calcSpreadRate = sqrt(((xhat(1,2:N) - xhat(1,1:N-1))/delT).^2+((xhat(2,2:N) - xhat(2,1:N-1))/delT).^2);

plot(1:N,xhat(3,1:N),'.-r');
plot(1:N,xhat(4,1:N),'.-b');
plot(1:N,sqrt(xhat(3,1:N).^2+xhat(4,1:N).^2),'.-g');
plot(1:N-1,calcSpreadRate,'--k');

xlabel('time-step [s]','Interpreter','Latex');ylabel('v [m/s]','Interpreter','Latex');
legend('v_x','v_y','|v|','|v|_{calc}');

subplot(3,1,3);
hold on;
plot(2:N,distCP2Obs(1,2:N),'.-b');
plot(2:N,distCP2Bound(2:N),'--b');
plot(2:N,distObs2Bound(1,2:N),'.-r');
sigma_a = sqrt(Cov_a_b(1,1,2:N));
sigma_a = sigma_a(:,:);
plot(2:N,sigma_a,'.-g');

if QKF_flag
    %    plot(i,Ex_zA(1)-EzA,'.-r');
end
xlabel('time-step [s]','Interpreter','Latex');ylabel('$d [m]$','Interpreter','Latex');
legend('||{CP-DP}||','||{CP-Actual}||','||{DP-Actual}||','$\hat{\sigma}_{a}$');set(legend,'Interpreter','latex');

figure(3)
hold on;grid on;
imax=N;
% Define colors
cvalues=1:N;%(P(1,1,1:imax).*P(2,2,1:imax)).^4;
%// Sort and get their index to access the color array
[~,idx] = sort(cvalues);
%// Create colormap with the "jet" colors
colormap('jet');
ColorMap = jet(numel(cvalues));

% plot default 50% confidence interval
h(1) = error_ellipse( P([1 2],[1 2],1)*ellipse_factor,xhat([1 2],1),'style','--','Color',ColorMap(1,:));%PlotPdf(mu,var,'b');
h(2) = drawVehicle(z(1,1,1),z(1,2,1),heading,1,50);
if NumberOfObserver>1
    for j=2:NumberOfObserver
        drawVehicle(z(1,1,j),z(1,2,j),-heading,1,50);
    end
end
h(3) = plot(xhat(1,1),xhat(2,1),'sb');
len = 500;
brdLineX = actualBoundary(1,1)*ones(len,1);
brdLineY = -len/2:len/2-1;
%hline1 = plot(brdLineX,brdLineY,'r');
ax = gca;
h(4) = line(brdLineX+.06,brdLineY,...
    'LineWidth',4,...
    'Color',[.8 .0 .0],...
    'Parent',ax);
h(5) =   plot(actualBoundary(1,1),actualBoundary(2,1),'sr');
%hleg1 = legend(h,'Uncertianty ($\Sigma_{x,y}$)','Obseraver','Estimation ($E_{x,y}$)','Initial Boundary','Actual Boundary');
hleg1 = legend(h([1,3,4,5]),'Uncertianty ($\Sigma_{x,y}$)','Estimation ($E_{x,y}$)','Initial Boundary','Actual Boundary');

set(hleg1,'Interpreter','latex');

tit1 = title([TestedMainTitle, ' - Setup']);
% automatic
set(hleg1,'Location','best')
eval(['print -depsc  '  ['.\figures\',TestedCFG,num2str(testCase),'_',num2str(fig_ct)]]);


for j=2:20:N
    if ~all(eig(P([1 2],[1 2],j))>0)
        warning('The covariance matrix (P) must be positive definite (it has non-positive eigenvalues) (@time=%d)',j);
    else
        error_ellipse( P([1 2],[1 2],j)*ellipse_factor,xhat([1 2],j),'style','--','Color',ColorMap(j,:));
    end
end
error_ellipse( P([1 2],[1 2],end)*ellipse_factor,xhat([1 2],end),'style','--','Color',ColorMap(end,:));


plot(xhat(1,:),xhat(2,:),'sb');
plot(actualBoundary(1,2:20:end),actualBoundary(2,2:20:end),'sr');
for j=1:NumberOfObserver
    ind = z(:,3,j)==1;
    plot(z(ind,1,j),z(ind,2,j),'k>')
    plot(z(~ind,1,j),z(~ind,2,j),'g>')
    
end
%Plot complete actual periphery progress in time
for i = 1:N-1:N %start and end
    plot( [ GPxhat(end,1,i) ; GPxhat(:,1,i)],[ GPxhat(end,2,i); GPxhat(:,2,i)],'r');
    line([xhat(1,i) z(i,1,1)],[xhat(2,i) z(i,2,1)]);
end
colorbar;
caxis([min(cvalues) max(cvalues)]);
axis equal;
figure(4);
hold on;grid on;
plot(1:N,distCP2Obs(1,1:N),'-.b','LineWidth',2,'MarkerSize',12);%Dash-dotted line
plot(1:N,distCP2Bound(1:N),'-b','LineWidth',2,'MarkerSize',12);%Solid line
plot(1:N,distObs2Bound(1,1:N),'-r','LineWidth',2,'MarkerSize',12);%Solid line
plot(1:N,squeeze(sqrt(Cov_a_b(1,1,1:N))),':g','LineWidth',2,'MarkerSize',12);%Dotted line
plot(1:N,sqrt(sigma2_z(1:N)),'--k','LineWidth',2,'MarkerSize',12);%Dashed line

xlabel('time-step [s]','Interpreter','Latex');ylabel('$$d [m]$$','Interpreter','Latex');
legend('$$||{CP-DP}||$$','$$||{CP-Actual}||$$','$||{DP-Actual}||$','$$\hat{\sigma}_{a}$$','$$\sigma_{a}^{-}$$');  set(legend,'Interpreter','latex');
fig_ct = fig_ct+1;

tit1 = title([TestedMainTitle,TestedMainSubTitle]);
eval(['print -depsc  '  ['.\figures\',TestedCFG,num2str(testCase),'_',num2str(fig_ct)]]);

figure(5);
hold on;grid on;
plot(1:N,squeeze(sqrt(Cov_a_b(1,1,1:N))),':g','LineWidth',2,'MarkerSize',12);%Dotted line
plot(1:N,sqrt(sigma2_z(1:N)),'--k','LineWidth',2,'MarkerSize',12);%Dashed line
plot(1:N,squeeze(sqrt(P(1,1,1:N))),':r','LineWidth',2,'MarkerSize',12);%Dotted line
plot(1:N,squeeze(sqrt(P(2,2,1:N))),':b','LineWidth',2,'MarkerSize',12);%Dotted line

xlabel('time-step [s]','Interpreter','Latex');ylabel('$$\sigma [m]$$','Interpreter','Latex');
legend('$$\hat{\sigma}_{a}$$','$$\sigma_{a}^{-}$$','$$\sigma_{x}$$','$$\sigma_{y}$$');  set(legend,'Interpreter','latex');
fig_ct = fig_ct+1;

tit1 = title([TestedMainTitle,TestedMainSubTitle]);
eval(['print -depsc  '  ['.\figures\',TestedCFG,num2str(testCase),'_',num2str(fig_ct)]]);

figure(6);
hold on;grid on;
plot(1:N,xhat(1,1:N),':r','LineWidth',2,'MarkerSize',12);%Dotted line
%plot(2:N,xhat(1,2:N)-xhat(1,1:N-1),':r','LineWidth',2);%Dotted line
plot(1:N,xhat(2,1:N),'--b','LineWidth',2,'MarkerSize',12);%Dashed line
xlabel('time-step [s]','Interpreter','Latex');ylabel('CP Position [m]','Interpreter','Latex');
legend('Pos_x','Pos_y');
fig_ct = fig_ct+1;

tit1 = title([TestedMainTitle,TestedMainSubTitle]);
eval(['print -depsc  '  ['.\figures\',TestedCFG,num2str(testCase),'_',num2str(fig_ct)]]);

figure(7);
hold on;grid on;
plot(1:N,xhat(3,1:N),':r','LineWidth',3,'MarkerSize',12);%Dotted line
plot(1:N,xhat(4,1:N),'--b','LineWidth',3,'MarkerSize',12);%Dashed line
%plot(1:N-1,calcSpreadRate,'--k','LineWidth',3,'MarkerSize',12);

xlabel('time-step [s]','Interpreter','Latex');ylabel('CP Spread Rate [m/s]','Interpreter','Latex');
legend('V_x','V_y');%,'|v|_{calc}');
fig_ct = fig_ct+1;

tit1 = title([TestedMainTitle,TestedMainSubTitle]);
eval(['print -depsc  '  ['.\figures\',TestedCFG,num2str(testCase),'_',num2str(fig_ct)]]);

if testCase==8
    save mle.mat P xhat
elseif testCase==9
    save map.mat P xhat
elseif testCase==10
    save qkl.mat P xhat
end

if testCase==11
    save mlePred.mat P xhat
elseif testCase==12
    save mapPred.mat P xhat
elseif testCase==13
    save qklPred.mat P xhat
end

if testCase==15
    save mleFull.mat P xhat
elseif testCase==16
    save mapFull.mat P xhat
elseif testCase==17
    save qklFull.mat P xhat
end
% sumQ = sum(q) % determine Sum q which is Chiˆ2 on N d.o.f.
% r = xcorr(nu); % get autocorrealtion of innovation
if recFlag
    figure(8)
    hold on;grid on;
    xlim([-2030 2030]);
    ylim([-2030 2030])
    facecolors_red = [...
        [1, 0, 0];...    % left wing
        [1, 0, 0];...    % right wing
        [1, 0, 0];...     % tail
        ];
    facecolors_green = [...
        [0, 1, 0];...    % left wing
        [0, 1, 0];...    % right wing
        [0, 1, 0];...     % tail
        ];
    % Define colors
    cvalues=1:N;%(P(1,1,1:imax).*P(2,2,1:imax)).^4;
    %// Sort and get their index to access the color array
    [~,idx] = sort(cvalues);
    %// Create colormap with the "jet" colors
    colormap('jet');
    ColorMap = jet(numel(cvalues));
    %     colorbar;
    %     caxis([min(cvalues) max(cvalues)]);
    %     axis equal;
    
    % plot default 50% confidence interval
    h(1) = error_ellipse( P([1 2],[1 2],1)*ellipse_factor,xhat([1 2],1),'style','--','Color',ColorMap(1,:));%PlotPdf(mu,var,'b');
    
    h(3) = plot(xhat(1,1),xhat(2,1),'sb');
    h(4) = plot( [ GPxhat(end,1,1) ; GPxhat(:,1,1)],[ GPxhat(end,2,1); GPxhat(:,2,1)],'r');
    
    q = quiver(GPxhat(1:10:end,1,1),GPxhat(1:10:end,2,1),GPxhat(1:10:end,1,2),GPxhat(1:10:end,2,2),0.1);
    c = q.Color;
    q.Color = 'red';
    hleg1 = legend(h([1,3,4]),'Uncertianty ($\Sigma_{x,y}$)','Estimation ($E_{x,y}$)','Actual Boundary');
    set(hleg1,'Interpreter','latex');
    tit1 = title([TestedMainTitle, ' - Setup']);
    % automatic
    set(hleg1,'Location','best')
    %eval(['print -depsc  '  ['.\figures\',TestedCFG,num2str(testCase),'_',num2str(fig_ct)]]);
    for i=1:N
        hel=error_ellipse( P([1 2],[1 2],i)*ellipse_factor,xhat([1 2],i),'style','--','Color',ColorMap(i,:));
        plot(xhat(1,i),xhat(2,i),'sb');
        
        %Plot complete actual periphery progress in time
        hp=plot( [ GPxhat(end,1,i) ; GPxhat(:,1,i)],[ GPxhat(end,2,i); GPxhat(:,2,i)],'r');
        hl=line([xhat(1,i) z(i,1,1)],[xhat(2,i) z(i,2,1)]);
        for j=1:NumberOfObserver
            if i>1
                heading = atan2d(z(i,2,j)-z(i-1,2,j),z(i,1,j)-z(i-1,1,j)) + 90;
            else
                heading = 90;
            end
            if z(i,3,j)==1
                hz=plot(z(i,1,j),z(i,2,j),'k>');
                h(2) = drawVehicle(z(i,1,j),z(i,2,j),heading,1,50);
                set(h(2),'FaceVertexCData',facecolors_red);               
            else
                hz=plot(z(i,1,j),z(i,2,j),'g>');
                h(2) = drawVehicle(z(i,1,j),z(i,2,j),heading,1,50);
                set(h(2),'FaceVertexCData',facecolors_green);                
            end
        end
        drawnow;
        %pause (0.1);%Normal slow
        frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
        writeVideo(writerObj, frame);
        delete(hz);
        delete(hl);
        delete(hp);
        delete(hel);
        delete(h(2));
        
    end
    close(writerObj); % Saves the movie.
end



