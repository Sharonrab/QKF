function fadeInAndOut(h,fadeInt,fadeFlag)
% clear all
% close all
% % Plot something and make is white
% h = plot(sin(0:0.1:pi));
% set(h,'color',[1 1 1]);

%% Create a fade matrix
originalRGB = get(h,'color');%[0 0 1];
%Create some other matrices to the fade
% fadeInt = 50;
% if fadeFlag == 0 % fade in
% rgbFade = [linspace(1,originalRGB(1),fadeInt)',linspace(1,originalRGB(2),fadeInt)',linspace(1,originalRGB(3),fadeInt)'];
% else %fade out
% rgbFade = [linspace(0,originalRGB(1),fadeInt)',linspace(0,originalRGB(2),fadeInt)',linspace(0,originalRGB(3),fadeInt)'];
% end
rgbFade = [1 1 1] - ([1 1 1] - originalRGB) * fadeInt;
%% loop through

% for k = 1:fadeInt
   set(h,'color',rgbFade);%rgbFade(k,:));
%    pause(0.1);
% end