%% MATLAB code for gait analysis with laser range finder

% Code for automatic detection of gait phases achieved via the analysis
% of measurements taken with a laser range finder positioned in front of
% the person walking.
% Phases (classes):
%   1 - double support left forward
%   2 - right swing phase
%   3 - double support right forward
%   4 - left swing phase
%   5 - standing
% After manual labelling of the acquired data paired with the extracted
% features, we trained some models with the Matlab Toolbox Classification
% Learner. The best performing model was a KNN.
% We used part of the labelled dataset for testing while the majority of
% the data was used for training with k-fold cross validation. In
% particular, we used the first 100 examples for testing and the others for
% training.
% For the realtime implementation it is possible to export the C/C++ code
% of the model.

% The following code simulates the realtime data acquisition (read data round
% by round), features extraction and automatic classification/detection of gait
% phases. This file and the gait_analysis.m one, can be merged together in
% the C/C++ realtime application.

clc;
clearvars;

% filenames
% data_lab_adj/test_5_forward.txt   r_max = 148
% data_lab_adj/test_6_forward.txt   r_max = 108
% data_lab_adj/test_7_turn.txt      r_max = 172
% data_lab_adj/test_8_zigzag.txt    r_max = 170
% data_lab_adj/test_9_auto.txt      r_max = 148
% data_lab_adj/test_10_auto.txt     r_max = 142

%% Load data

r_max = 170;

fid = fopen('data_lab_adj/test_8_zigzag.txt', 'r');
format = '%d theta:%f Dist:%f Q:%d';
data = textscan(fid, format);
fclose(fid);

varNames = {'round', 'theta', 'distance', 'quality'};
d_test = table(data{1}, data{2}, data{3}, data{4}, 'VariableNames', varNames);

ppr = height(d_test)/r_max;  % average number of samples per turn, used to 
                             % compute speed (250*ppr = turn duration in
                             % microseconds)

% for supervised classification purpose
sup_samples = readtable('feat_labelled/feat_8.txt');
target = sup_samples.target;

% load trained model
load('Fine_KNN_model.mat')

% read images of gait phases to show in the plot
phase_img = {imread('phase_img/double_support_left.png'), ... 
             imread('phase_img/right_swing.png'), ...
             imread('phase_img/double_support_right.png'), ...
             imread('phase_img/left_swing.png'), ...
             imread('phase_img/standing.png')
            };
    
%% Features extraction, gait phases detection and plot

% The features extraction part works identically as in gait_analysis.m
% After features extraction there is the classification.
% This code structure is to simulate realtime behavior

% Open a figure
figure();

%Initialise plot
ax = axes;
ax.XLim = [-800 0];
ax.YLim = [-350 500];
hold on;
daspect([1 1 1])

% -- Init parameters ------------------------------------------------------

% area of interest = samples of point cloud to keep
THETA_MIN = 160;    % min angle of area of interest 
THETA_MAX = 200;    % max angle of area of interest
DIST_MAX = 900;    % max dist of area of interes

% parameters of search region - area in which we expect each single leg to be
THETA_INT = 8;  % Interval for leg search region
DIST_INT = 150;  % Interval for leg search region

% parameters for circumference fitting
R = 55;  % radius of fitting circumference
N_GRIDX = 10;   % number of x coordinates for grid search
N_GRIDY = 12;   % number of y coordinates for grid search
GRID_SPACE = 5;    % spacing between coordinates in grid
CIRC_PARAM = [R,N_GRIDX,N_GRIDY,GRID_SPACE];

% init features vectors - features derived from circumference fitting
% result
c_x = zeros(r_max,2);           % centers of legs - x coordinates
c_y = zeros(r_max,2);           % centers of legs - y coordinates
rel_pos_x = zeros(r_max,1);     % relative position of legs - x coordinate
rel_pos_y = zeros(r_max,1);     % relative position of legs - y coordinate
leg_dist = zeros(r_max,1);      % legs distance
leg_angle = zeros(r_max,1);     % angle between legs centers (the angle is the one
                                %  in the x-y plane where leg_angle = 0 when
                                %  the legs are aligned in the direction perpendicular
                                %  to the motion. e.g. people standing with paired
                                %  feet gives a null angle)
vel = zeros(r_max,2);           % velocity of the 2 legs
rel_vel = zeros(r_max,1);       % relative velocity of the legs

% class column (gait phase) - collected from labelled data for validation and test
target = [target;zeros(r_max-length(target),1)]; 

% features table intialization
sz = [r_max 13];
featNames = {'round', 'left_x', 'left_y', 'right_x', 'right_y', 'rel_x', 'rel_y', 'dist', 'angle','left_velocity','right_velocity','rel_velocity','target'};
% note: keep 'target' as last column
vt(1:length(featNames)) = "double";
feat = table('Size', sz, 'VariableNames', featNames, 'VariableTypes', vt);
feat.round(:) = linspace(1,r_max,r_max);

% description of gait phases to show in the plot
phase_desc = {'double support left forward', 'right swing phase', ...
              'double support right forward', 'left swing phase', ...
              'standing'};


% -- Main loop, work as would do in realtime ------------------------------

for r = 1:r_max
    
    % acquisition of 1 turn (grabScanData)
    theta = d_test.theta(d_test.round == r);        % angles of samples in point cloud
    dist = d_test.distance(d_test.round == r);      % distance of samples in point cloud
    qual = d_test.quality(d_test.round == r);       % quality of the measurements
    
    % mirroring correction
    theta = -1*theta+360;
    
    % remove samples out of range of interest
    in_range = (dist <= DIST_MAX & theta >= THETA_MIN & theta <= THETA_MAX);
    theta = theta(in_range);
    dist = dist(in_range);
    qual = qual(in_range);
    
    % clean data - remove samples with zero quality
    dist = dist(qual ~= 0);
    theta = theta(qual ~= 0);
    
    
    % ######### legs extraction #########
    
    % first round: legs supposed still and frontal wrt lidar, we use a
    %              threshold on theta to split the legs
    
    % other rounds: for each leg we extract the points of point cloud residing 
    %               in a search reagion identified by the previous center of mass
    %               +/- some intervals in distance and angle
    
    if r == 1
        % threshold - halfway between max and min angle
        th = (theta(1) + theta(end)) / 2;
        % left and right legs - polar coordinates
        ll = [theta(theta < th), dist(theta < th)];
        rl = [theta(theta > th), dist(theta > th)];
        % left and right legs - cartesian coordinates
        [ll_x,ll_y] = pol2cart(deg2rad(theta(theta < th)), dist(theta < th));
        [rl_x,rl_y] = pol2cart(deg2rad(theta(theta > th)), dist(theta > th));
    else
        % search regions
        interval_l = (theta < (cm_l(1) + THETA_INT) & theta > (cm_l(1) - THETA_INT) & dist < (cm_l(2) + DIST_INT) & dist > (cm_l(2) - DIST_INT));
        interval_r = (theta < (cm_r(1) + THETA_INT) & theta > (cm_r(1) - THETA_INT) & dist < (cm_r(2) + DIST_INT) & dist > (cm_r(2) - DIST_INT));
        % left and right legs - polar coordinates
        ll = [theta(interval_l) , dist(interval_l)];
        rl = [theta(interval_r) , dist(interval_r)];
        % left and right legs - cartesian coordinates
        [ll_x,ll_y] = pol2cart(deg2rad(theta(interval_l)) , dist(interval_l));
        [rl_x,rl_y] = pol2cart(deg2rad(theta(interval_r)) , dist(interval_r));
        % all leg points for complete point cloud plot (not needed in realtime scenario)
        dist = dist(interval_l | interval_r);
        theta = theta(interval_l | interval_r);
    end

    % ######### features extraction #########
    
    % centres of mass for the 2 legs (uniform weight) - polar
    cm_l = mean(ll,1);
    cm_r = mean(rl,1);
    
    % centres of mass for the 2 legs (uniform weight) - cartesian
    cm_x = [mean(ll_x,1); mean(rl_x,1)];
    cm_y = [mean(ll_y,1); mean(rl_y,1)];
    
    % circle fit
    [c_x(r,:),c_y(r,:),best_x,best_y] = circ_fit(ll_x,ll_y,rl_x,rl_y,cm_x,cm_y,CIRC_PARAM);
    
    % relative position - we assumed positive the position of forward left leg
    rel_pos_x(r) = c_x(r,1)-c_x(r,2);
    rel_pos_y(r) = c_y(r,1)-c_y(r,2);
    
    % distance between legs' centers
    leg_dist(r) = sqrt(rel_pos_x(r)^2+rel_pos_y(r)^2);
    
    % angle of line linking legs' centers
    leg_angle(r) = rad2deg(acot(rel_pos_y(r)/rel_pos_x(r)));

    % velocity of the single leg and relative velocity between legs - only x coordinate (advancement)   
    if r >= 2
        vel(r,:) = ((c_x(r,:) - c_x(r-1,:))*10^-3) / ((250*ppr)*10^-6);  
        rel_vel(r) = ((rel_pos_x(r) - rel_pos_x(r-1))*10^-3) / ((250*ppr)*10^-6);  
    end
    
    % feature vector
    feat(r,:) = {r, c_x(r,1),c_y(r,1),c_x(r,2),c_y(r,2),rel_pos_x(r),rel_pos_y(r),leg_dist(r),leg_angle(r),vel(r,1),vel(r,2),rel_vel(r),target(r)};
    
    
    % ######### CLASSIFICATION #########
    
    % Use the traine KNN model to predict the gait phase from the feature
    % vector of the current acquisition round
    gait_phase = Fine_KNN_model.predictFcn(feat(r,2:end-1));
    % Set options to be shown in plot in case of correct or wrong
    % prediction
    if gait_phase == feat{r,end}
        text_opt = 'green';
    else
        text_opt = 'red';
    end
    
    % ######### update and plot #########
    
    % plotting is just for visualization of points describing the legs,
    % legs circles and centers of mass. It allows to see the motion of the
    % legs. Not necessary in real time scenario.
    
    % complete legs point cloud
    [x,y] = pol2cart(deg2rad(theta),dist);
    
    % plot
    if r == 1
        p_xy = plot(ax,x,y,'.');
        cm_xy = plot(ax,cm_x,cm_y,'ro');
        best_xy = plot(ax,best_x,best_y,'r.');
        c_xy = plot(ax,c_x(r,:),c_y(r,:),'r--s');
    else    % plot update
        set(p_xy, 'YData', y, 'Xdata', x);
        set(cm_xy, 'YData', cm_y, 'Xdata', cm_x);
        set(best_xy, 'YData', best_y, 'Xdata', best_x);
        set(c_xy, 'YData', c_y(r,:), 'Xdata', c_x(r,:));
    end
    
    % Display round
    h = text(-780,480, ['Round: ', num2str(r)]);
    
    % Display description of the predicted gait phase
    phase = text(-780,450,phase_desc{gait_phase},'Color',text_opt);
    
    % Display image of the predicted gait phase
    imagesc('XData',[-125 -20],'YData',[480 270],'CData',phase_img{gait_phase});
    
    pause(.15);   % simulation step
    
    delete(h)
    delete(phase)
    
end
