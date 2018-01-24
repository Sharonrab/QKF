function UT_UpdateFitGaussEstimator(varargin)
if nargin==0
clear all;
close all;
testCase = 1;%1 to 5 MLE, 6,7,8 are MAP+MLE, 9 for QKF
elseif nargin==1    
testCase = varargin{1}; 
else
    error('Incorrect number of input arguments.')
end

recFlag =0;
DEBUG_FLAG =0;
REC_CNT =0;%the record counter is used to initiating the time step to record the resulted calculation
plottools('on') ;plotbrowser('off'); set(0,'DefaultFigureWindowStyle','docked')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Setup and initialization %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
QKF_flag = 0;%set the configuration to eighther quantization KF or else

N = 25;
delT = 1;
MaxAffectedDist = 2000;
MinAffectedDist = 0;
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

figure(1)
hold on
set(gca,'YDir','normal');
if recFlag
    %% Set up the movie.
    writerObj = VideoWriter(['.\figures\',TestedCFG,'Simulation',num2str(fig_ct),'.avi']); % Name it.
    writerObj.FrameRate = 2; % How many frames per second.
    open(writerObj);
end
set(gca,'XDir','normal');




% Q = sigma2Q * [ delT^3/3 delT^2/2
% delT^2/2 delT ];

F = [   1 0 delT 0
    0 1 0 delT
    0 0 1 0
    0 0 0 1 ];
if QKF_flag
    H = [1 0 0 0]; % measurement matrix
else
    H = [1 0 0 0;0 1 0 0;0 0 0 0;0 0 0 0]; % measurement matrix
end

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
%start with high uncertainty of the prior
%1000 meter of standard deviation and 2 meter per second in spread rate
P(:,:,1) = [1000^2 0 0 0; 0 1000^2 0 0; 0 0 MaxSpreadRate^2 0; 0 0 0 MaxSpreadRate^2];

if QKF_flag
    R = sigma2R * [ 1 ];
else
    R = sigma2R * [ 1 0;0 1];
end
x = zeros(4,N);% [ x y progress rate]
w = zeros(4,N);
Cov_a_b =[];
% z = zeros(2,N);% [ x y status] the status is quantization of the measurments -1 for inside the convex hull, +1 for outside of the convex hull
distCP2Obs = [];
distCP2Bound=[];
distObs2Bound= [];
v = zeros(2,N);
%number of cycles to avoid estimation after crossing actual boundary
avoidEstAftCrossDelay = 0;
%flag to enable limiter on calculated spread rate
limitCalcSpreadRateFlag = 0;
% propotional gain on the innovation of expected value
decayFactor = 0.1;% the size of innovation we allowed to the estimator given the observation
MAP_flag=0;
MLE_flag=0;
corr_corr_flag = 0;%enable correlation correction based on calculated ratio
% flag to indicate CP updated due to drone crossing actual boundary
updateOnCrossingFlag=0;
% flag to enable CP update feature
enableUpdateOnCrossingFlag=0;
%set Fire origin
origin = [0 ;0];
actualBoundary = [360;0];%after 30 minutes of fire start with nominal spread rate of 0.2 m/sec (0.2*30*60)
CP_Behind = actualBoundary(1)-180;% if the uncertianty of the spread rate is going down to 0.1 (0.1*30*60)
CP_Infront = actualBoundary(1)+180;% if the uncertianty of the spread rate is going up to 0.3 (0.3*30*60)
%common observation locations for first set of unit tests
drone_initial_location =[-1000 0];
distObs2Bound(1) = sqrt((drone_initial_location(1)-actualBoundary(1,1))^2 + (drone_initial_location(2)-actualBoundary(2,1))^2);

time_to_go = distObs2Bound(1)/(droneNominalSpeed - NominalSpreadRate);
z(1,:)= [drone_initial_location:droneNominalSpeed*delT:time_to_go*droneNominalSpeed+500];%approaching from the left
N= length(z);
z(2,:) = zeros(1,N);
%z1(3,:) = zeros(1,N);

switch testCase
    case 1
        %CASE I: drone on BURNED area approaching MOVING control point from the left
        MLE_flag = 1;      
        MAP_flag = 0;
        QKF_flag = 0;  
        TestedCFG = 'FGE_cdf_test';%fitted Gaussian Estimator
        TestedMainTitle = 'Fitted Gaussian Estimator (MLE)';
        TestedMainSubTitle = ' - Fixed Observer ( Predict-OFF, Estimation-ON, CP-Ahead) ';
        set(0,'defaultTextInterpreter','latex'); %trying to set the default
        fig_ct = 1;
        %%%%Initialize observation configuration
        OnBurnedArea=1;%set detected status
        NumberOfObserver = 1; %set the number of observer for interlacing to one sequential measurement
        z(3,:) = ones(1,N)*OnBurnedArea;
        z = z;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Uncertainty with measurments
        xhat = zeros(4,N); % state estimate
        xhat(1,1) = CP_Infront;%initial x position
        xhat(2,1) = 0;%initial y position
        xhat(3,1) = 0;%initial vx velocity
        xhat(4,1) = 0;%initial vy velocity
        plotDir = 'r-->';
        heading = pi/2;
        actualBorderVel = [NominalSpreadRate;0];
        avoidPredFlag = 0;
        avoidEstFlag = 1;
        avoidEstAftCrossFlag = 1;
        enableUpdateOnCrossingFlag =1;

    case 2
        %CASE II: drone on BURNED area approaching MOVING control point from the left
        MLE_flag = 1;        
        MAP_flag = 0;
        QKF_flag = 0;  
        TestedCFG = 'FGE_cdf_test';%fitted Gaussian Estimator
        TestedMainTitle = 'Fitted Gaussian Estimator (MLE)';
        TestedMainSubTitle = '- Predict-ON, Estimation-OFF, CP-Ahead';
        set(0,'defaultTextInterpreter','latex'); %trying to set the default
        fig_ct = 1;
        OnBurnedArea=1;%set detected status
        
        % Uncertainty with measurments
        xhat = zeros(4,N); % state estimate
        xhat(1,1) = CP_Infront;%initial x position
        xhat(2,1) = 0;%initial y position
        xhat(3,1) = NominalSpreadRate;%initial vx velocity
        xhat(4,1) = 0;%initial vy velocity
        plotDir = 'r-->';
        heading = pi/2;
        
        actualBorderVel = [NominalSpreadRate;0];
        avoidPredFlag = 0;
        avoidEstFlag = 1;
        avoidEstAftCrossFlag =1;
        enableUpdateOnCrossingFlag =1;

        %CASE III: drone on UNBURNED area approaching FIXED control point from the left
    case 3
        MLE_flag = 1;        
        MAP_flag = 0;
        QKF_flag = 0;  
        TestedCFG = 'FGE_cdf_test';%fitted Gaussian Estimator
        TestedMainTitle = 'Fitted Gaussian Estimator (MLE)';
        TestedMainSubTitle = ' - Predict-ON, Estimation-OFF, CP-Ahead, CrossingUpdate-On';
        set(0,'defaultTextInterpreter','latex'); %trying to set the default
        fig_ct = 1;
        OnBurnedArea=1;%set detected status
        % Uncertainty with measurments
        xhat = zeros(4,N); % state estimate
        xhat(1,1) = CP_Infront;%initial x position
        xhat(2,1) = 0;%initial y position
        xhat(3,1) = NominalSpreadRate;%initial vx velocity
        xhat(4,1) = 0;%initial vy velocity
        plotDir = 'r-->';
        heading = pi/2;
        
        actualBorderVel = [NominalSpreadRate;0];
        avoidPredFlag = 0;
        avoidEstFlag = 0;
        avoidEstAftCrossFlag = 0;
        limitCalcSpreadRateFlag =0;
        enableUpdateOnCrossingFlag =1;
        %CASE IV: drone on UNBURNED area approaching FIXED control point from the
        %right
    case 4
        MLE_flag = 1;  
        MAP_flag = 0;
        QKF_flag = 0;          
        TestedCFG = 'FGE_cdf_test';%fitted Gaussian Estimator
        TestedMainTitle = 'Fitted Gaussian Estimator (MLE)';
        TestedMainSubTitle = ' - Predict-ON, Estimation-ON, CP-Back, CrossingUpdate-On';
        set(0,'defaultTextInterpreter','latex'); %trying to set the default
        fig_ct = 1;
        OnBurnedArea=1;%set detected status
        % Uncertainty with measurments
        xhat = zeros(4,N); % state estimate
        xhat(1,1) = CP_Behind;%initial x position
        xhat(2,1) = 0;%initial y position
        xhat(3,1) = NominalSpreadRate;%initial vx velocity
        xhat(4,1) = 0;%initial vy velocity
        plotDir = 'r-->';
        heading = pi/2;
        
        actualBorderVel = [NominalSpreadRate;0];
        avoidPredFlag = 0;
        avoidEstFlag = 0;
        avoidEstAftCrossFlag = 0;
        enableUpdateOnCrossingFlag =1;

        %CASE V: drone on BURNED area approaching MOVING control point (Vx>0) from the
        %left
    case 5
        MLE_flag = 1;        
        MAP_flag = 0;
        QKF_flag = 0;  
        TestedCFG = 'FGE_cdf_test';%fitted Gaussian Estimator
        TestedMainTitle = 'Fitted Gaussian Estimator (MLE)';
        TestedMainSubTitle = ' - Predict-ON, Estimation-ON, CP-Back, CrossingUpdate-Delayed';
        set(0,'defaultTextInterpreter','latex'); %trying to set the default
        fig_ct = 1;
        OnBurnedArea=1;%set detected status
        % Uncertainty with measurments
        xhat = zeros(4,N); % state estimate
        xhat(1,1) = CP_Behind;%initial x position
        xhat(2,1) = 0;%initial y position
        xhat(3,1) = NominalSpreadRate;%initial vx velocity
        xhat(4,1) = 0;%initial vy velocity
        plotDir = 'r-->';
        heading = pi/2;
        
        actualBorderVel = [NominalSpreadRate;0];
        avoidPredFlag = 0;
        avoidEstFlag = 0;
        avoidEstAftCrossFlag = 0;
        avoidEstAftCrossDelay = 2;
        enableUpdateOnCrossingFlag =1;

        %CASE VI: drone on BURNED area approaching MOVING control point (Vx<0) from the
        %left
    case 6
        MLE_flag = 1;        
        MAP_flag = 1;  
        QKF_flag = 0;  
        TestedCFG = 'FGE_cdf_test';%fitted Gaussian Estimator
        TestedMainTitle = 'Fitted Gaussian Estimator (MAP)';
        TestedMainSubTitle = ' - Predict-ON, Estimation-ON, CP-Ahead, CrossingUpdate-On';
        set(0,'defaultTextInterpreter','latex'); %trying to set the default
        fig_ct = 1;
        OnBurnedArea=1;%set detected status
        % Uncertainty with measurments
        xhat = zeros(4,N); % state estimate
        xhat(1,1) = CP_Behind;%initial x position
        xhat(2,1) = 0;%initial y position
        xhat(3,1) = NominalSpreadRate;%initial vx velocity
        xhat(4,1) = 0;%initial vy velocity
        plotDir = 'r-->';
        heading = pi/2;
        
        actualBorderVel = [NominalSpreadRate;0];
        avoidPredFlag = 0;
        avoidEstFlag = 0;
        avoidEstAftCrossFlag = 0;
        avoidEstAftCrossDelay = 0;  
        enableUpdateOnCrossingFlag =1;
        
        decayFactor = 1;

        %CASE VII: drone on UNBURNED area approaching MOVING control point (Vx>0) from the
        %left
    case 7
        MLE_flag = 1;        
        MAP_flag = 1;  
        QKF_flag = 0;  
        TestedCFG = 'FGE_cdf_test';%fitted Gaussian Estimator
        TestedMainTitle = 'Fitted Gaussian Estimator (MAP)';
        TestedMainSubTitle = ' - Predict-ON, Estimation-ON, CP-Ahead, CrossingUpdate-On';
        set(0,'defaultTextInterpreter','latex'); %trying to set the default
        fig_ct = 1;
        OnBurnedArea=1;%set detected status
        % Uncertainty with measurments
        xhat = zeros(4,N); % state estimate
        xhat(1,1) = CP_Behind;%initial x position
        xhat(2,1) = 0;%initial y position
        xhat(3,1) = NominalSpreadRate;%initial vx velocity
        xhat(4,1) = 0;%initial vy velocity
        plotDir = 'r-->';
        heading = pi/2;
        
        actualBorderVel = [NominalSpreadRate;0];
        avoidPredFlag = 0;
        avoidEstFlag = 0;
        avoidEstAftCrossFlag = 0;
        avoidEstAftCrossDelay = 0;

        decayFactor = 0.7;
        enableUpdateOnCrossingFlag =1;
        
        %CASE IIX: drone on UNBURNED area approaching MOVING control point (Vx<0) from the
        %left
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    %%%% START test 3 solutions (MLE,MAP,QKF) with NO prediction step
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 8
        %Generate variance 
        MLE_flag = 1;        
        MAP_flag = 0;
        QKF_flag = 0; 
        decayFactor = 1;
        TestedCFG = 'MLE_Comp_test';%fitted Gaussian Estimator
        TestedMainTitle = 'MLE';
        TestedMainSubTitle = ' - Fixed Observer ( Predict-OFF, Estimation-ON, CP-Back) ';
        set(0,'defaultTextInterpreter','latex'); %trying to set the default
        fig_ct = 1;
        OnBurnedArea=1;%set detected status
        z(1,:)= -1000*ones(1,N);%fixed observer
        % Uncertainty with measurments
        xhat = zeros(4,N); % state estimate
        xhat(1,1) = CP_Behind;%initial x position
        xhat(2,1) = 0;%initial y position
        xhat(3,1) = 0;%initial vx velocity
        xhat(4,1) = 0;%initial vy velocity
        plotDir = 'r-->';
        heading = pi/2;
        
        actualBorderVel = [NominalSpreadRate;0];
        avoidPredFlag = 1;
        avoidEstFlag = 0;
        avoidEstAftCrossFlag = 0;
        %CASE IX: drone on UNBURNED area approaching MOVING control point (Vy<0) from the
        %left
    case 9
        %Generate variance for MAP comparision
        MLE_flag = 1;        
        MAP_flag = 1;
        QKF_flag = 0;  
        decayFactor = 1;
        TestedCFG = 'MAP_Comp_test';%fitted Gaussian Estimator
        TestedMainTitle = 'MAP';
        TestedMainSubTitle = ' - Fixed Observer ( - Predict-OFF, Estimation-ON, CP-Back) ';
        set(0,'defaultTextInterpreter','latex'); %trying to set the default
        fig_ct = 1;
        OnBurnedArea=1;%set detected status
        z(1,:)= -1000*ones(1,N);%fixed observer
        % Uncertainty with measurments
        xhat = zeros(4,N); % state estimate
        xhat(1,1) = CP_Behind;%initial x position
        xhat(2,1) = 0;%initial y position
        xhat(3,1) = 0;%initial vx velocity
        xhat(4,1) = 0;%initial vy velocity
        plotDir = 'r-->';
        heading = pi/2;
        
        actualBorderVel = [NominalSpreadRate;0];
        avoidPredFlag = 1;
        avoidEstFlag = 0;
        avoidEstAftCrossFlag = 0;
        %CASE X: drone on UNBURNED area approaching MOVING control point
        %(Vx>0) from the left
    case 10
        %Generate variance for QKF comparision
        MLE_flag = 0;        
        MAP_flag = 0;
        QKF_flag = 1;        
        decayFactor = 1;
        TestedCFG = 'QKF_Comp_test';%fitted Gaussian Estimator
        TestedMainTitle = 'QKF';
        TestedMainSubTitle = ' - Fixed Observer ( Predict-OFF, Estimation-ON, CP-Back) ';
        set(0,'defaultTextInterpreter','latex'); %trying to set the default
        fig_ct = 1;
        OnBurnedArea=1;%set detected status
        z(1,:)= -1000*ones(1,N);%fixed observer
        % Uncertainty with measurments
        xhat = zeros(4,N); % state estimate
        xhat(1,1) = CP_Behind;%initial x position
        xhat(2,1) = 0;%initial y position
        xhat(3,1) = 0;%initial vx velocity
        xhat(4,1) = 0;%initial vy velocity
        plotDir = 'r-->';
        heading = pi/2;
        
        actualBorderVel = [NominalSpreadRate;0];
        avoidPredFlag = 1;
        avoidEstFlag = 0;
        avoidEstAftCrossFlag = 0;
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    %%%% END test 3 solutions (MLE,MAP,QKF) with NO prediction step
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    %%%% START test 3 solutions (MLE,MAP,QKF) with prediction step
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 11
        %Generate variance 
        MLE_flag = 1;        
        MAP_flag = 0;
        QKF_flag = 0; 
        decayFactor = 1;
        TestedCFG = 'MLE_Comp_test';%fitted Gaussian Estimator
        TestedMainTitle = 'MLE';
        TestedMainSubTitle = ' - Fixed Observer ( Predict-ON, Estimation-ON, CP-Back, CrossingUpdate-On) ';
        set(0,'defaultTextInterpreter','latex'); %trying to set the default
        fig_ct = 1;
        OnBurnedArea=1;%set detected status
        z(1,:)= -1000*ones(1,N);%fixed observer
        % Uncertainty with measurments
        xhat = zeros(4,N); % state estimate
        xhat(1,1) = CP_Behind;%initial x position
        xhat(2,1) = 0;%initial y position
        xhat(3,1) = 0;%initial vx velocity
        xhat(4,1) = 0;%initial vy velocity
        plotDir = 'r-->';
        heading = pi/2;
        
        actualBorderVel = [NominalSpreadRate;0];
        avoidPredFlag = 0;
        avoidEstFlag = 0;
        avoidEstAftCrossFlag = 0;
        enableUpdateOnCrossingFlag =1;
        
        %CASE IX: drone on UNBURNED area approaching MOVING control point (Vy<0) from the
        %left
    case 12
        %Generate variance for MAP comparision
        MLE_flag = 1;        
        MAP_flag = 1;
        QKF_flag = 0;  
        decayFactor = 1;
        TestedCFG = 'MAP_Comp_test';%fitted Gaussian Estimator
        TestedMainTitle = 'MAP';
        TestedMainSubTitle = ' - Fixed Observer ( Predict-ON, Estimation-ON, CP-Back, CrossingUpdate-On) ';
        set(0,'defaultTextInterpreter','latex'); %trying to set the default
        fig_ct = 1;
        OnBurnedArea=1;%set detected status
        z(1,:)= -1000*ones(1,N);%fixed observer
        % Uncertainty with measurments
        xhat = zeros(4,N); % state estimate
        xhat(1,1) = CP_Behind;%initial x position
        xhat(2,1) = 0;%initial y position
        xhat(3,1) = 0;%initial vx velocity
        xhat(4,1) = 0;%initial vy velocity
        plotDir = 'r-->';
        heading = pi/2;
        
        actualBorderVel = [NominalSpreadRate;0];
        avoidPredFlag = 0;
        avoidEstFlag = 0;
        avoidEstAftCrossFlag = 0;
        enableUpdateOnCrossingFlag =1;
        
        %CASE X: drone on UNBURNED area approaching MOVING control point
        %(Vx>0) from the left
    case 13
        %Generate variance for QKF comparision
        MLE_flag = 0;        
        MAP_flag = 0;
        QKF_flag = 1;        
        decayFactor = 1;
        TestedCFG = 'QKF_Comp_test';%fitted Gaussian Estimator
        TestedMainTitle = 'QKF';
        TestedMainSubTitle = ' - Fixed Observer ( Predict-ON, Estimation-ON, CP-Back, CrossingUpdate-On) ';
        set(0,'defaultTextInterpreter','latex'); %trying to set the default
        fig_ct = 1;
        OnBurnedArea=1;%set detected status
        z(1,:)= -1000*ones(1,N);%fixed observer
        % Uncertainty with measurments
        xhat = zeros(4,N); % state estimate
        xhat(1,1) = CP_Behind;%initial x position
        xhat(2,1) = 0;%initial y position
        xhat(3,1) = 0;%initial vx velocity
        xhat(4,1) = 0;%initial vy velocity
        plotDir = 'r-->';
        heading = pi/2;
        
        actualBorderVel = [NominalSpreadRate;0];
        avoidPredFlag = 0;
        avoidEstFlag = 0;
        avoidEstAftCrossFlag = 0;
        enableUpdateOnCrossingFlag =1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    %%%% END test 3 solutions (MLE,MAP,QKF) with prediction step
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 case 14
        MLE_flag = 0;        
        MAP_flag = 0;  
        QKF_flag = 1;  
        TestedCFG = 'QKF_test';%fitted Gaussian Estimator
        TestedMainTitle = 'Fitted Gaussian Estimator (QKF)';
        TestedMainSubTitle = ' - Predict-ON, Estimation-ON, CP-Ahead, CrossingUpdate-On';
        set(0,'defaultTextInterpreter','latex'); %trying to set the default
        fig_ct = 1;
        OnBurnedArea=1;%set detected status
        % Uncertainty with measurments
        xhat = zeros(4,N); % state estimate
        xhat(1,1) = CP_Behind;%initial x position
        xhat(2,1) = 0;%initial y position
        xhat(3,1) = NominalSpreadRate;%initial vx velocity
        xhat(4,1) = 0;%initial vy velocity
        plotDir = 'r-->';
        heading = pi/2;
        
        actualBorderVel = [NominalSpreadRate;0];
        avoidPredFlag = 0;
        avoidEstFlag = 0;
        avoidEstAftCrossFlag = 0;
        avoidEstAftCrossDelay = 0;

        decayFactor = 0.7;
        enableUpdateOnCrossingFlag =1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    %%%% START test 3 solutions (MLE,MAP,QKF) with prediction step
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 15
   % corr_corr_flag=1;
        %Generate variance 
        MLE_flag = 1;        
        MAP_flag = 0;
        QKF_flag = 0; 
        decayFactor = 1;
        TestedCFG = 'MLE_Comp_test';%fitted Gaussian Estimator
        TestedMainTitle = 'MLE';
        TestedMainSubTitle = ' - Moving Observer ( Predict-ON, Estimation-ON, CP-Back, CrossingUpdate-On) ';
        set(0,'defaultTextInterpreter','latex'); %trying to set the default
        fig_ct = 1;
        OnBurnedArea=1;%set detected status
        % Uncertainty with measurments
        xhat = zeros(4,N); % state estimate
        xhat(1,1) = CP_Behind;%initial x position
        xhat(2,1) = 0;%initial y position
        xhat(3,1) = 0;%initial vx velocity
        xhat(4,1) = 0;%initial vy velocity
        plotDir = 'r-->';
        heading = pi/2;
        
        actualBorderVel = [NominalSpreadRate;0];
        avoidPredFlag = 0;
        avoidEstFlag = 0;
        avoidEstAftCrossFlag = 0;
        enableUpdateOnCrossingFlag =1;
        
        %CASE IX: drone on UNBURNED area approaching MOVING control point (Vy<0) from the
        %left
    case 16
        %Generate variance for MAP comparision
        MLE_flag = 1;        
        MAP_flag = 1;
        QKF_flag = 0;  
        decayFactor = 1;
        TestedCFG = 'MAP_Comp_test';%fitted Gaussian Estimator
        TestedMainTitle = 'MAP';
        TestedMainSubTitle = ' - Moving Observer ( Predict-ON, Estimation-ON, CP-Back, CrossingUpdate-On) ';
        set(0,'defaultTextInterpreter','latex'); %trying to set the default
        fig_ct = 1;
        OnBurnedArea=1;%set detected status
        % Uncertainty with measurments
        xhat = zeros(4,N); % state estimate
        xhat(1,1) = CP_Behind;%initial x position
        xhat(2,1) = 0;%initial y position
        xhat(3,1) = 0;%initial vx velocity
        xhat(4,1) = 0;%initial vy velocity
        plotDir = 'r-->';
        heading = pi/2;
        
        actualBorderVel = [NominalSpreadRate;0];
        avoidPredFlag = 0;
        avoidEstFlag = 0;
        avoidEstAftCrossFlag = 0;
        enableUpdateOnCrossingFlag =1;
        
        %CASE X: drone on UNBURNED area approaching MOVING control point
        %(Vx>0) from the left
    case 17
        %Generate variance for QKF comparision
        MLE_flag = 0;        
        MAP_flag = 0;
        QKF_flag = 1;        
        decayFactor = 1;
        TestedCFG = 'QKF_Comp_test';%fitted Gaussian Estimator
        TestedMainTitle = 'QKF';
        TestedMainSubTitle = ' - Moving Observer ( Predict-ON, Estimation-ON, CP-Back, CrossingUpdate-On) ';
        set(0,'defaultTextInterpreter','latex'); %trying to set the default
        fig_ct = 1;
        OnBurnedArea=1;%set detected status
        % Uncertainty with measurments
        xhat = zeros(4,N); % state estimate
        xhat(1,1) = CP_Behind;%initial x position
        xhat(2,1) = 0;%initial y position
        xhat(3,1) = 0;%initial vx velocity
        xhat(4,1) = 0;%initial vy velocity
        plotDir = 'r-->';
        heading = pi/2;
        
        actualBorderVel = [NominalSpreadRate;0];
        avoidPredFlag = 0;
        avoidEstFlag = 0;
        avoidEstAftCrossFlag = 0;
        enableUpdateOnCrossingFlag =1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    %%%% END test 3 solutions (MLE,MAP,QKF) with prediction step
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
     
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    %%%% START test solutions optimization for MAP
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    case 18
        %Generate variance for QKF comparision
        %R =R *1000;
        MLE_flag = 0;        
        MAP_flag = 0;
        QKF_flag = 1;        
        decayFactor = 1;
        TestedCFG = 'Optimization_test';%fitted Gaussian Estimator
        %TestedMainTitle = 'MAP';
        TestedMainTitle = 'QKF';
        z(1,:)= -1000*ones(1,N);%fixed observer
        z(2,:)= 0*ones(1,N);%fixed observer
        TestedMainSubTitle = ' - Stationary Observer ( Predict-OFF, Estimation-ON, CP-Back, CrossingUpdate-On) ';
        set(0,'defaultTextInterpreter','latex'); %trying to set the default
        fig_ct = 1;
        OnBurnedArea=1;%set detected status
        % Uncertainty with measurments
        xhat = zeros(4,N); % state estimate
        xhat(1,1) = CP_Behind;%initial x position
        xhat(2,1) = 0;%initial y position
        xhat(3,1) = 0;%initial vx velocity
        xhat(4,1) = 0;%initial vy velocity
        plotDir = 'r-->';
        heading = pi/2;
        
        actualBorderVel = [NominalSpreadRate;0];
        avoidPredFlag = 1;
        avoidEstFlag = 0;
        avoidEstAftCrossFlag = 0;
        enableUpdateOnCrossingFlag =1;
end

%z(1,:)= [-10:1:N-11];

%z(1,:)= [0:1:N-1];


h = [];
%initial plot parameters
distCP2Bound(1) = sqrt((actualBoundary(1,1)-xhat(1,1))^2 + (actualBoundary(2,1)-xhat(2,1))^2);
distCP2Obs(1) = sqrt((z(1,1)-xhat(1,1))^2 + (z(2,1)-xhat(2,1))^2);
distObs2Bound(1) = sqrt((z(1,1)-actualBoundary(1,1))^2 + (z(2,1)-actualBoundary(2,1))^2);
Cov_a_b = P(:,:,1);
sigZ(1) = P(1,1,1);
%running configuration
once=1;
passedBoundaryFlag = 0;
passedBoundaryTimetag = inf;
for i=2:N
    curMeas = z(1:2,i);%sequence of measurements from both side of the grid point
    %Update boundary
    actualBoundary(:,i) = actualBoundary(:,i-1) + actualBorderVel*delT;
    %Update tracking status
    distObs2Bound(i) = sqrt((curMeas(1)-actualBoundary(1,i))^2 + (curMeas(2)-actualBoundary(2,i))^2);
    boundTol = 0.1;
    if distObs2Bound(i) < boundTol && (~passedBoundaryFlag) % passing by the  boundary from the left
        OnBurnedArea= ~OnBurnedArea;%set detected status
        passedBoundaryFlag = 1;
        passedBoundaryTimetag = i;
    end
    %initiate timestep to record the resulted computation
    %     if i==2
    %         REC_CNT=2;
    %     else
    %         REC_CNT=0;
    %     end
    [xhat(:,i),P(:,:,i),Cov_a_b(:,:,i),sigZ(i),distCP2Obs(i)] = UpdateFGEstimator(i,delT,origin,xhat(:,i-1),curMeas,F,P(:,:,i-1),Q,H,R, OnBurnedArea,MLE_flag,MAP_flag,QKF_flag,avoidPredFlag,avoidEstFlag,limitCalcSpreadRateFlag,decayFactor,corr_corr_flag,MaxSpreadRate,MinSpreadRate,MaxAffectedDist,MinAffectedDist,MaxCP_Heading,MinCP_Heading,DEBUG_FLAG,REC_CNT,TestedCFG) ;
    %evaluate distance between prediction and actual for analysis of perfomance
    distCP2Bound(i) = sqrt((actualBoundary(1,i)-xhat(1,i))^2 + (actualBoundary(2,i)-xhat(2,i))^2);
    
    if enableUpdateOnCrossingFlag && passedBoundaryFlag &&~updateOnCrossingFlag
        xhat([1 2],i) = curMeas;
        updateOnCrossingFlag =1;%indicate that CP location has been updated
    end
    %Avoid estimation after passing by boundary for several cycles TBD
    if passedBoundaryFlag == 1 && i < (passedBoundaryTimetag + avoidEstAftCrossDelay)
        avoidEstFlag = 1;
    end
    % the following condition is to prevent estimation for couple cycles
    % after drone croses the boundary and only in case we want to support
    % that feature in the current unit test
    if (i > (passedBoundaryTimetag + avoidEstAftCrossDelay) ) && ~avoidEstAftCrossFlag
        avoidEstFlag = 0;
    end
    
    
    
    %Analyzed results
    %q(:,i) = nu(:,i)'*inv(S(:,i))*nu(:,i);
    % xhat(1,i) = xpred(1,1);
    % P(1,1) =Ppred(1,1);
    
    %Plot for Debug
    if DEBUG_FLAG
        
        figure(1)
        mu_x = xhat(1,i);
        var_x = P(1,1);
        mu_y = xhat(2,i);
        var_y = P(2,2);
        if once
            % plot default 50% confidence interval
            h(1) = error_ellipse( P([1 2],[1 2]),xhat([1 2],i),'style','--');%PlotPdf(mu,var,'b');
            h(2) = plot(curMeas(1),curMeas(2),plotDir);
            fadeInAndOut(h(2),i/N,1);
            
            h(2) = drawVehicle(curMeas(1),curMeas(2),heading,i,0.2)
            h(2) = plot(curMeas(1),curMeas(2),plotDir);
            
            h(3) = plot(xhat(1,i),xhat(2,i),'sb');
            len = 20;
            brdLineX = actualBoundary(1,1)*ones(len,1);
            brdLineY = -10:len/2-1;
            %hline1 = plot(brdLineX,brdLineY,'r');
            ax = gca;
            h(4) = line(brdLineX+.06,brdLineY,...
                'LineWidth',4,...
                'Color',[.8 .0 .0],...
                'Parent',ax);
            h(5) =   plot(actualBoundary(1,1),actualBoundary(2,1),'sr');
            
            once=0;
            leg1 = legend(h,'Estimations ($\Sigma_{x,y}$)','Sequence of Obseravtions',...
                'Estimations ($E_{x,y}$)','Initial Boundary','Actual Boundary');
            set(leg1,'Interpreter','latex');
            
            
        else
            %if ~mod(i,5)
            h(1) = error_ellipse( P([1 2],[1 2]),xhat([1 2],i),'style','--');%PlotPdf(mu,var,'b');
            %   fadeInAndOut(h(1),i/N,1);
            %end
            h(2) = plot(curMeas(1),curMeas(2),plotDir);
            
            fadeInAndOut(h(2),i/N,1);
            
            drawVehicle(curMeas(1),curMeas(2),heading,i,0.4);
            
            h(3) = plot(xhat(1,i),xhat(2,i),'sb');
            fadeInAndOut(h(3),i/N,1);
            plot(actualBoundary(1,i),actualBoundary(2,i),'sr');
            
            drawnow;
            
            if recFlag
                pause (0.1);%Normal slow
                frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
                writeVideo(writerObj, frame);
            end
        end
        %title(sprintf('p(x|z) mean (x,y): %.2f, %.2f,var (xx,yy): %.2f,  %.2f',mu_x,mu_y,var_x, var_y));
        xlabel('x');
        ylabel('y');
        
        figure(99)
        subplot(3,1,1);
        hold on;
        plot(i,xhat(1,i),'.-r');
        plot(i,xhat(2,i),'.-b');
        xlabel('time-step');ylabel('pos [m]','Interpreter','Latex');
        legend('x','y');
        subplot(3,1,2);
        hold on;
        plot(i,xhat(3,i),'.-r');
        plot(i,xhat(4,i),'.-b');
        plot(i,norm(xhat([3 4],i)),'.-g');
        xlabel('time-step');ylabel('v [m/s]','Interpreter','Latex');
        legend('v_x','v_y','|v|');
        
        subplot(3,1,3);
        hold on;
        plot(i,distCP2Obs(i),'.-b');
        
        plot(i,distCP2Bound(i),'.-r');
        plot(i,sqrt(Cov_a_b(1,1,i)),'.-g');
        if QKF_flag
            plot(i,Ex_zA(1)-EzA,'.-r');
        end
        xlabel('time-step');ylabel('$d [m] $','Interpreter','Latex');
        legend('||{CP-DP}||','||{CP-Actual}||','\sigma_{a}^{-}');
        
        
        %plot(cnt,,'.-b');
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
plot(2:N,distCP2Obs(2:N),'.-b');
plot(2:N,distCP2Bound(2:N),'--b');
plot(2:N,distObs2Bound(2:N),'.-r');
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
% plot default 50% confidence interval
h(1) = error_ellipse( P([1 2],[1 2],1)*4,xhat([1 2],1),'style','--');%PlotPdf(mu,var,'b');
h(2) = drawVehicle(z(1,1),z(2,1),heading,1,50)
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
hleg1 = legend(h,'Uncertianty ($\Sigma_{x,y}$)','Obseravtion','Estimation ($E_{x,y}$)','Initial Boundary','Actual Boundary');
set(hleg1,'Interpreter','latex');
tit1 = title([TestedMainTitle, ' - Setup']);
% automatic
set(hleg1,'Location','best')
eval(['print -depsc  '  ['.\figures\',TestedCFG,num2str(testCase),'_',num2str(fig_ct)]]);
for j=2:N
    error_ellipse( P([1 2],[1 2],j)*4,xhat([1 2],j),'style','--');
end
figure(4);
hold on;grid on;
plot(1:N,distCP2Obs(1:N),'-.b','LineWidth',2,'MarkerSize',12);%Dash-dotted line
plot(1:N,distCP2Bound(1:N),'-b','LineWidth',2,'MarkerSize',12);%Solid line
plot(1:N,distObs2Bound(1:N),'-r','LineWidth',2,'MarkerSize',12);%Solid line
plot(1:N,squeeze(sqrt(Cov_a_b(1,1,1:N))),':g','LineWidth',2,'MarkerSize',12);%Dotted line
plot(1:N,sqrt(sigZ(1:N)),'--k','LineWidth',2,'MarkerSize',12);%Dashed line

xlabel('time-step [s]','Interpreter','Latex');ylabel('$$d [m]$$','Interpreter','Latex');
legend('$$||{CP-DP}||$$','$$||{CP-Actual}||$$','$||{DP-Actual}||$','$$\hat{\sigma}_{a}$$','$$\sigma_{a}^{-}$$');  set(legend,'Interpreter','latex');
fig_ct = fig_ct+1;

tit1 = title([TestedMainTitle,TestedMainSubTitle]);
eval(['print -depsc  '  ['.\figures\',TestedCFG,num2str(testCase),'_',num2str(fig_ct)]]);

figure(5);
hold on;grid on; 
plot(1:N,squeeze(sqrt(Cov_a_b(1,1,1:N))),':g','LineWidth',2,'MarkerSize',12);%Dotted line
plot(1:N,sqrt(sigZ(1:N)),'--k','LineWidth',2,'MarkerSize',12);%Dashed line
plot(1:N,squeeze(sqrt(P(1,1,1:N))),':r','LineWidth',2,'MarkerSize',12);%Dotted line

xlabel('time-step [s]','Interpreter','Latex');ylabel('$$\sigma [m]$$','Interpreter','Latex');
legend('$$\hat{\sigma}_{a}$$','$$\sigma_{a}^{-}$$','$$\sigma_{x}$$');  set(legend,'Interpreter','latex');
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
plot(1:N-1,calcSpreadRate,'--k','LineWidth',3,'MarkerSize',12);

xlabel('time-step [s]','Interpreter','Latex');ylabel('CP Spread Rate [m/s]','Interpreter','Latex');
legend('V_x','V_y','|v|_{calc}');
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
    close(writerObj); % Saves the movie.
end


