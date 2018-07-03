
%% Initialization
clear ; close all; clc

%% ======================= Plotting =======================
fprintf('Plotting Data ...\n')
data1 = load('sim5.dat'); data2 = load('sim7.dat');
X1 = data1(:, 5); X2 = data2(:, 5);
y1_2 = data1(:, 2); y2_2 = data2(:, 2);
y1_3 = data1(:, 3); y2_3 = data2(:, 3);
y1_4 = data1(:, 4); y2_4 = data2(:, 4);


% Plot Data
figure; % open a new figure window
plot(X1, y1_2, 'rx', 'MarkerSize', 5)     
hold on 
plot(X2, y2_2, 'gx', 'MarkerSize', 5)  
grid on;
ylabel('R VALUE');             % Set the y−axis label
xlabel('Time elapsed...');     % Set the x−axis label

figure; % open a new figure window
plot(X1, y1_4, 'rx', 'MarkerSize', 5)     
hold on 
plot(X2, y2_4, 'gx', 'MarkerSize', 5)  
grid on;
ylabel('E VALUE');             % Set the y−axis label
xlabel('Time elapsed...');     % Set the x−axis label


fprintf('Program paused. Press enter to continue.\n');
pause;
