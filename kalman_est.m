% ECEN 946 Homework 4 Start
% author: Subharthi Banerjee
% starter code written by Dr Andrew Harms
% Create the output of an autoregressive (AR) model
%
% y(k) = - a_1 * y(k-1) - a_2 * y(k-2) + v(k)
% 
% where v(k) is a white Gaussian process.

% Simulate for 500 time steps
clear
clc
close all
iter = 500;

% Use these parameters for the AR model
a_1 = 0.7;
a_2 = 0.12;
R_k = 0.05;
phi_k = eye(2, 2); 
% State-space model
x   = [a_1; a_2];

% Driving Gaussian white noise sequence
v   = sqrt(R_k).*randn(iter,1);

% Simulate measurement sequence and measurement model
y   = zeros(1,iter);
H   = zeros(iter,2);
H(1,:) = [0   , 0]; y(1) =               v(1);
H(2,:) = [y(1), 0]; y(2) = -x(1).*y(1) + v(2);
for idk = 3:iter
  y(1,idk) = -x(1).*y(idk-1) - x(2).*y(idk-2) + v(idk);
  H(idk,:) = [y(idk-1), y(idk-2)];
end



%% HW4.14.  estimation starts from here ---
% according to the kalman filter we need to 
% calculate different parameters iteratively
% time starts from 1
disp("processing HW414 ..................");%....
x_hat = zeros(2, iter); % depends on a_1 and a_2
x_h = [0, 0]';
P = eye(2); % initialization
% start 
for idk = 1 : iter
    HH = H(idk, :);
    K = P * HH' ;%* inv(HH * P * HH' + R_k);
    R_k_eta = HH * P * HH'+ R_k;
    
    x_h = phi_k*x_h + phi_k * K * inv(R_k_eta) * (y(idk) - HH * x_h);
    x_hat(:, idk) = x_h;
    P = phi_k * P * phi_k' - phi_k * K * inv(R_k_eta) * K' * phi_k';
end
 figure
 % remember the coefficients will be negative to the estimate
h1 = plot(1:iter, -x_hat(1, :)); % true value of y is without noise?
hold on
h2 = plot(1:iter, -x_hat(2, :));
 set(h1                            , ...
            'LineWidth'       , 1.5           , ...
            'LineStyle'       , '--'            , ...
            'Marker'          , 'none'         , ...
            'MarkerSize'      , 1          , ...
            'MarkerEdgeColor' , [.4 .3 .1]  , ...
            'MarkerFaceColor' , 'none'  );
         set(h2                            , ...
            'LineWidth'       , 1.5           , ...
            'LineStyle'       , ':'            , ...
            'Marker'          , 'none'         , ...
            'MarkerSize'      , 1           , ...
            'MarkerEdgeColor' , [.7 .8 .8]  , ...
            'MarkerFaceColor' , 'none'  );
 xtt = get(gca, 'XTick');
 set(gca, 'FontSize', 12);
 ytt = get(gca, 'YTick');
 set(gca, 'FontSize', 12);

grid on
legend({'$\hat{a_1}$', '$\hat{a_2}$'}, 'Interpreter'...
    , 'latex', 'Location', 'Northeast');
xlabel('True values', 'Interpreter', 'latex')
ylabel('Time/Iterations', 'Interpreter', 'latex')
title('Estimated value as $\hat{x}$ and True value $y$',...
    'Interpreter', 'latex');
set(gcf, 'PaperPositionMode', 'auto');
    
 print('-depsc2', 'hw414.eps')
 
 


%% HW415 now the coefficients change with time

% $a_1(k) = 0.4 + 0.5\cos(3\pik/200)$
% $a_2(k) = 0.5 + 0.3\sin(2\pik/200)$

clc
close all
iter = 500;

% Use these parameters for the AR model

R_k = 0.05;
%phi_k = eye(2, 2); 
sig_k = eye(2);
% State-space model
%x   = [a_1; a_2];
% the value of PHI -- is it observability problem? I guess I needed 
% to work on finding the value of PHI
% other than that the simulation worked okay
phi_k = eye(2); % for testing that is what I found efficient
Q_k = 0.001; % for testing --- increasing this value increases error
% Driving Gaussian white noise sequence
v   = sqrt(R_k).*randn(iter,1);

% Simulate measurement sequence and measurement model
y   = zeros(1,iter);
H   = zeros(iter,2);
H(1,:) = [0   , 0]; y(1) =               v(1);
H(2,:) = [y(1), 0]; 
% for k = 2
a_1 = 0.4 + 0.5*cos(3 * pi *2/200);

y(2) = -a_1.*y(1) + v(2);
for idk = 3:iter
  a_1 = 0.4 + 0.5*cos(3 * pi *idk/200);
  a_2 = 0.5 + 0.3*sin(2 * pi *idk/200);
  x = [a_1; a_2];
  y(1,idk) = -x(1).*y(idk-1) - x(2).*y(idk-2) + v(idk);
  H(idk,:) = [y(idk-1), y(idk-2)];
end


%% HW4.15.  estimation starts from here ---
% according to the kalman filter we need to 
% calculate different parameters iteratively
% time starts from 1
disp("processing HW414 ..................");%....
x_hat = zeros(2, iter); % depends on a_1 and a_2
x_h = [0, 0]';
P = eye(2); % initialization
% start 
A1 = 0.4 + 0.5*cos(3 * pi *(1:iter)/200);
A2 = 0.5 + 0.3*sin(2 * pi *(1:iter)/200);
for idk = 1 : iter
    HH = H(idk, :);
    K = P * HH' ;%* inv(HH * P * HH' + R_k);
    R_k_eta = HH * P * HH'+ R_k;
    
    x_h = phi_k*x_h + phi_k * K * inv(R_k_eta) * (y(idk) - HH * x_h);
    x_hat(:, idk) = x_h;
    P = phi_k * P * phi_k' + sig_k * Q_k * sig_k'...
        - phi_k * K * inv(R_k_eta) * K' * phi_k';
end
 figure
 % remember the coefficients will be negative to the estimate
h1 = plot(1:iter, -x_hat(1, :)); % true value of y is without noise?
hold on
h2 = plot(1:iter, A1);
 set(h1                            , ...
            'LineWidth'       , 1.5           , ...
            'LineStyle'       , '--'            , ...
            'Marker'          , 'none'         , ...
            'MarkerSize'      , 1          , ...
            'MarkerEdgeColor' , [.4 .3 .1]  , ...
            'MarkerFaceColor' , 'none'  );
         set(h2                            , ...
            'LineWidth'       , 1.5           , ...
            'LineStyle'       , ':'            , ...
            'Marker'          , 'none'         , ...
            'MarkerSize'      , 1           , ...
            'MarkerEdgeColor' , [.7 .8 .8]  , ...
            'MarkerFaceColor' , 'none'  );
 xtt = get(gca, 'XTick');
 set(gca, 'FontSize', 12);
 ytt = get(gca, 'YTick');
 set(gca, 'FontSize', 12);

grid on
legend({'$\hat{a_1}$', '$a_1$'}, 'Interpreter'...
    , 'latex', 'Location', 'Northeast');
ylabel('Estimated values', 'Interpreter', 'latex')
xlabel('Time/Iterations', 'Interpreter', 'latex')
title('Estimated value as $\hat{x}$ and True value $y$',...
    'Interpreter', 'latex');
set(gcf, 'PaperPositionMode', 'auto');
    
 print('-depsc2', 'hw415a1.eps')
 
 
 figure
 % remember the coefficients will be negative to the estimate
h1 = plot(1:iter, -x_hat(2, :)); % true value of y is without noise?
hold on
h2 = plot(1:iter, A2);
 set(h1                            , ...
            'LineWidth'       , 1.5           , ...
            'LineStyle'       , '--'            , ...
            'Marker'          , 'none'         , ...
            'MarkerSize'      , 1          , ...
            'MarkerEdgeColor' , [.4 .3 .1]  , ...
            'MarkerFaceColor' , 'none'  );
         set(h2                            , ...
            'LineWidth'       , 1.5           , ...
            'LineStyle'       , ':'            , ...
            'Marker'          , 'none'         , ...
            'MarkerSize'      , 1           , ...
            'MarkerEdgeColor' , [.7 .8 .8]  , ...
            'MarkerFaceColor' , 'none'  );
 xtt = get(gca, 'XTick');
 set(gca, 'FontSize', 12);
 ytt = get(gca, 'YTick');
 set(gca, 'FontSize', 12);

grid on
legend({'$\hat{a_2}$', '$a_2$'}, 'Interpreter'...
    , 'latex', 'Location', 'Northeast');
xlabel('Iterations', 'Interpreter', 'latex')
ylabel('Estimated Values', 'Interpreter', 'latex')
title('Estimated value as $\hat{x}$ and True value $y$',...
    'Interpreter', 'latex');
set(gcf, 'PaperPositionMode', 'auto');
    
 print('-depsc2', 'hw415a2.eps')
 