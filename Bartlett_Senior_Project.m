%%%%%%%%%% COVID-19 Senior Research Project %%%%%%%%%%%
% Simulating Earth's Field NMR Experiment 
% Author: Carson Bartlett 
% Date: May 7th, 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SINGLE ISOLATED PROTON IN A CONSTANT EXTERNAL MAGNETIC FIELD %%
%% Initializing Variables 
clear; clc; close all;
Be = 0.5E-4; %T
B_ext = [0, 0, Be]; % Earth's magnetic field in z direction
mu_p = 1; %(J/T)
gamma = 2.675E8/1E3; % (rad/ms*T)

w = gamma*norm(B_ext); % angular frequency (rad/ms)
freq = w/(2*pi); %precession frequency of proton about B field (1/ms)
T = 1/freq; %period of the proton's precession (ms)
fprintf('The Larmor Frequency is: %1$.4fms\n', T);

% This corresponds to its initial conditions
mu0 = (mu_p*[sin(pi/3)*cos(pi/2); sin(pi/3)*sin(pi/2); cos(pi/3)]);
%{
- By the definition of the cross-product:
  [a x b] = [(a2*b3 - a3*b2), (a3*b1 - a1*b3), (a1*b2 - a2*b1)]
  Use this to define state functions for solving for magnetic moments

- The equation of motion for a precessing magnetic moment in an external 
  magnetic field: du/dt = cross(gamma*u_vector, Bext_vector);
%}

%% Numerically Solving for Magnetic Moment Components
f = @(t,u) [gamma*(u(2)*B_ext(3) - u(3)*B_ext(2)); ...
            gamma*(u(3)*B_ext(1) - u(1)*B_ext(3)); ... 
            gamma*(u(1)*B_ext(2) - u(2)*B_ext(1))];
              
tspan = linspace(0, T, 1E2);

[t, u] = ode45(f, tspan, mu0);

%% Plotting Magnetic Moment in Space
figure
hold on;
curve1 = animatedline('LineWidth', 2);
set(gca, 'XLim', [-1 1], 'YLim', [-1 1], 'ZLim', [-0.5 1.5]);
view(43,24); % Sets the angle of the view to see in 3D
title('Single, Isolated Proton in a Constant Magnetic Field')

zlabel('\fontsize{14}z')
ylabel('\fontsize{14}y')
xlabel('\fontsize{14}x')
mu0_vector = plot3([0 mu0(1)], [0 mu0(2)], [0 mu0(3)]);
for i = 1:length(u)
    addpoints(curve1, u(i,1), u(i,2), u(i,3));
    head = scatter3(u(i,1), u(i,2), u(i,3),'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0 .75 .75]);
    % Magnetic Moment Vector as a Function of Time
    mu_vector = plot3([0 u(i,1)], [0 u(i,2)], [0 u(i,3)], 'Color', 'red');
    drawnow
    delete(mu_vector);
    delete(head);
end


%% BULK MAGNETIZATION IN A CONSTANT EXTERNAL MAGNETIC FIELD %%
%% Initializing Variables
M0 = 1;
T1 = 2500; %ms
T2 = 2500; %ms
tp = 3*T1; % Five T1 Time Constants 
%% Numerically Solving for Magnetization Components

% State Equations for Magnetization
fM = @(t, M) ...
    gamma*[((M(2)*B_ext(3)) - (M(3)*B_ext(2))); ...
           ((M(3)*B_ext(1)) - (M(1)*B_ext(3))); ...
           ((M(1)*B_ext(2)) - (M(2)*B_ext(1)))] ...
           - ([M(1), M(2), 0]/T2)' + ([0, 0, (M0 - M(3))]/T1)';
            
tt = linspace(0, tp, 1E3);
Mint = [0, M0, 0]; % Initial conditions for magnetization 
[t, v] = ode45(fM, tt, Mint);

%% Plotting Magnetization
figure
hold on;
curve2 = animatedline('LineWidth', 2);
set(gca, 'XLim', [-1 1], 'YLim', [-1 1], 'ZLim', [-1.5 1.5]);
view(43,24); %Sets the angle of the view to see in 3D

title('Bulk Magnetization in a Constant Magnetic Field')
zlabel('\fontsize{14}z')
ylabel('\fontsize{14}y')
xlabel('\fontsize{14}x')

curve2 = animatedline('LineWidth', 2);
for i = 1:length(v)
    addpoints(curve2, v(i,1), v(i,2), v(i,3));
    head2 = scatter3(v(i,1), v(i,2), v(i,3),'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0 .75 .75]);
    M_vector = plot3([0 v(i,1)], [0 v(i,2)], [0 v(i,3)], 'Color', 'red', 'LineWidth', 3);
    
    % Drawing Magnetization Component Vectors 
    M_vector_x = plot3([v(i,1) 2*v(i,1)], [v(i,2) v(i,2)], [v(i,3) v(i,3)], 'Color', 'blue', 'LineWidth', 2);
    M_vector_y = plot3([v(i,1) v(i,1)], [v(i,2) 2*v(i,2)], [v(i,3) v(i,3)], 'Color', 'green', 'LineWidth', 2);
    M_vector_z = plot3([v(i,1) v(i,1)], [v(i,2) v(i,2)], [v(i,3) 2*v(i,3)], 'Color', 'magenta', 'LineWidth', 2);
    
    drawnow
    delete(M_vector)
    delete(M_vector_x)
    delete(M_vector_y)
    delete(M_vector_z)
    delete(head2);
end  
hold off;


%% Plotting Transverse Magnetization 
figure;
plot(t, v(:,1));
hold on;

title('Transverse Magnetization')
M_x = zeros(1, length(v(:, 1)));
for i = 1:length(M_x)
    M_x(1, i) = sqrt(v(i,1)^2 + v(i,2)^2);
end
xlim([0 10])
ylabel('\fontsize{14}y')
xlabel('\fontsize{14}x')


% curve3 = animatedline('LineWidth', 2);
% set(gca, 'XLim', [0 20], 'YLim', [-1 1]);
% mag = @(x, y) sqrt(x^2 + y^2);
% M_perp = zeros(1, length(v));
% for i = 1:length(v)
%     M_perp(i) = mag(v(i,1), v(i,2));
%     addpoints(curve3, t(i), M_perp(i));
%     head3 = scatter(t(i), M_perp(i),'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0 .75 .75]);
%     drawnow
%     delete(head3);
% end

hold off


%% APP BUILDING WORKFLOW:
%{
    Add startup function:
    ->app.UIFigure in Component Browser 
    ->Callbacks
    ->Add StartupFcn Callback
    
    Sharing App:
    - Package App button under Designer Tab
    - Complete will create an installer file that can be shared

    Helper Function (Tutorial 12):
    - Functions in Editor Tab -> Private or Public if it needs to be defined
      for several figures 
    Ex. function plotgraph(app)
            x = 0:0.01: (360/180)*pi;
            y = sin(x)
            plot(app.UIAxes, y);

    - Use this in Startup Fnc or Button Callbacks

   Same Callback for Multiple Components (Tutorial 13): 
       x = 0:1:100
       m = app.mEditField.Value;
       y = m*x;
       plot(app.UIAxes)
       
       Add Existing Callback to Edit Field to add plot functionality to
       input field by pressing 'Enter' 
       - Go to component browser and right-click on edit field 
         => Callbacks
         => Select existing callback

  Multi-Variable Functions (Tutorial 13): 
      function plotgraph(app, m, c)
          x = 0:0.01: (360/180)*pi;
          y = sin(m*x + c)
          plot(app.UIAxes, y);

      - Add numeric fields 
      - In the callback for the plot button add two variables to the
        plotgraph() function we called in that callback
        
        m = app.mEditField.Value
        c = app.cEditField.Value
%}

