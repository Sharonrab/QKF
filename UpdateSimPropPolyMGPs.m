function [GPxhat] = UpdateSimPropPolyMGPs(varargin)
if nargin==0
    clear all;
    close all;
    % T - simulation time duration
    T = 35;
    t=1;
elseif nargin==4
    %current update time-step
    t = varargin{1};
    % T - simulation time duration
    T = varargin{2};
    %GPxhat - array of grid points with full time sequence
    GPxhat = varargin{3};   
    % NominalSpreadRate - the  Nominal Spread Rate of the periphery
    NominalSpreadRate = varargin{4};     
elseif nargin==5
    %current update time-step
    t = varargin{1};
    % T - simulation time duration
    T = varargin{2};
    %GPxhat - array of grid points with full time sequence
    GPxhat = varargin{3};
    %GPxhat - array of grid points with full time sequence
    init_actualBoundary = varargin{4};  
    % NominalSpreadRate - the  Nominal Spread Rate of the periphery
    NominalSpreadRate = varargin{5}; 
else
    error('Incorrect number of input arguments.')
end
% debug mode to plot internal parmas
debugFlag=0;
recFlag=0;
% delT - time step
delT = 1;

% F - the system matrix
% F - the system matrix
F = [   1 0     delT    0
    0 1     0       delT
    0 0     1       0
    0 0     0       1 ];
% B - Control matrix
B = [0       0
    0       0
    1       0
    0       1];
%F = [1 0 1 0; 0 1 0 1; 0 0 1 0; 0 0 0 1];
% u - wind velocity
u = 0*[0.1 0.1]';
max_speed = NominalSpreadRate;
%Q - Model Uncertainty
Qfactor = 1; % process noise mult factor
sigmaQ = Qfactor*sqrt(0.1);
Q = sigmaQ^2 * [ 1 0 0 0% process noise covariance matrix
    0 1 0 0
    0 0 1 0
    0 0 0 1];
%configure scenario physical parameters
%center location of the phenomenon
origin = [0 0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Number of grid point
N_GP = 60;
N_state = 4;
% Pgp - Covariance per grid point
% Pgp= zeros(N_GP,N_state,N_state);
% P - Covariance
P = zeros(N_state,N_state);
% for i = 1: N_GP
%     Pgp(i,:,:) = 3*Q;
% end
% X ~ [ x y progress rate]
if t==1
    GPxhat = zeros(N_GP,N_state,T); % state estimate
    % Initialize all grid point around the origin
    for i = 1 : N_GP
        initialHeadOnAngle = 360/N_GP*(i-1)*pi/180;
        GPxhat(i,1,1) = init_actualBoundary;%initial x position
        GPxhat(i,2,1) = 0;%initial y position
        GPxhat(i,[1 2],1) = SphericalRotation2x2(GPxhat(i,[1 2],1)', initialHeadOnAngle);
        GPxhat(i,3,1) = NominalSpreadRate;%initial Vx velocity
        GPxhat(i,[3 4],1) = SphericalRotation2x2(GPxhat(i,[3 4],1)', initialHeadOnAngle);
    end
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    GPxhat(:,:,t) = SimPropPolyMGPs(N_GP,GPxhat(:,:,t-1),F,B,u,P,Q,origin,max_speed,debugFlag);
    if debugFlag
        figure(1)
        %         xlim([-2 35]);
        %         ylim([-15 35]);
        plot(GPxhat(:,1,t), GPxhat(:,2,t), 's');
        plot( [ GPxhat(N_GP,1,t) ; GPxhat(:,1,t)],[ GPxhat(N_GP,2,t); GPxhat(:,2,t)]);
        
        hold on;
    end
end
if debugFlag
    figure(1)
    plot( [ GPxhat(N_GP,1,t) ; GPxhat(:,1,t)],[ GPxhat(N_GP,2,t); GPxhat(:,2,t)]);
    
end
if recFlag
    pause (0.1);%Normal slow
    frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
    writeVideo(writerObj, frame);
end
if debugFlag
    
    figure(1)
    for i = 1:T
        plot( [ GPxhat(N_GP,1,i) ; GPxhat(:,1,i)],[ GPxhat(N_GP,2,i); GPxhat(:,2,i)]);
        
    end
    figure(2)
    hold on;
    for j = 1:N_GP
        speed_mag = sqrt(GPxhat(j,3,:).^2+GPxhat(j,4,:).^2);
        speed_mag = speed_mag(:);
        plot([1:T],speed_mag );
    end
    
end
if recFlag
    close(writerObj); % Saves the movie.
end

