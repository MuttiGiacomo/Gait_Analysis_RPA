function [c_x,c_y,best_x,best_y] = circ_fit(ll_x,ll_y,rl_x,rl_y,cm_x,cm_y,CIRC_PARAM)
    
    % -- DESCRIPTION ------------------------------------------------------

    % This function takes as input some information on the left and right
    % leg points and some parameters used to define a grid for the research
    % of the best fitting circumference for the legs' points.
    % The output is the center of the best fitting circumference and the
    % points of such circle. The points of the best circle are used only
    % for plot purposes, so that output can be suppressed in a realtime
    % scenario.
    
    % The best circumference is determined analyzing the error committed by
    % approximating the leg points with a circle of known radius centered
    % in a grid point.
    % The center of the best fitting circumference is found as the center
    % of mass of all the grid points where the mass is given by the inverse
    % of the error.

    % INPUTS:
    % ll_x: left leg points - x coordinates
    % ll_y: left leg points - y coordinates
    % rl_x: right leg points - x coordinates
    % rl_y: right leg points - y coordinates
    % cm_x: left leg center of mass - x coordinate
    % cm_y: left leg center of mass - y coordinate
    % CIRC_PARAM: parameters to define the grid and radius for best circumference fitting
    
    % OUTPUTS:
    % c_x: center of best fitting circumference - x coordinate
    % c_y: center of best fitting circumference - y coordinate
    % best_x: center of best fitting circumference - x coordinate
    % best_y: center of best fitting circumference - x coordinate
    
    % -- INITIALIZATION ---------------------------------------------------
    
    % CIRC_PARAM remapping
    [R,N_GRIDX,N_GRIDY,GRID_SPACE] = deal(CIRC_PARAM(1),CIRC_PARAM(2),CIRC_PARAM(3),CIRC_PARAM(4));
    
    % vector of angles for circle representation (only needed if interested in plotting)
    sim_theta = linspace(0,2*pi,360*2);
    
    % grid creation: set of (x,y) points used as centers of test circumferences
    % a grid for each leg based on the position of the center of mass
    l_gridx = ones(1,N_GRIDX)*cm_x(1) - (GRID_SPACE:GRID_SPACE:GRID_SPACE*(N_GRIDX));
    l_gridy = ones(1,N_GRIDY)*cm_y(1) + (-2*GRID_SPACE:GRID_SPACE:GRID_SPACE*(N_GRIDY-3));
    r_gridx = ones(1,N_GRIDX)*cm_x(2) - (GRID_SPACE:GRID_SPACE:GRID_SPACE*(N_GRIDX));
    r_gridy = ones(1,N_GRIDY)*cm_y(2) + (-GRID_SPACE*(N_GRIDY-3):GRID_SPACE:2*GRID_SPACE);
    
    % error vectors initialization
    l_error_ij = zeros(N_GRIDX,N_GRIDY);
    r_error_ij = zeros(N_GRIDX,N_GRIDY);
    
    % initialization of variables for computing the center
    l_x_weighted_pos_ij = zeros(N_GRIDX,N_GRIDY);
    l_y_weighted_pos_ij = zeros(N_GRIDX,N_GRIDY);
    r_x_weighted_pos_ij = zeros(N_GRIDX,N_GRIDY);
    r_y_weighted_pos_ij = zeros(N_GRIDX,N_GRIDY);
    l_x_weighted_sum = 0;
    l_y_weighted_sum = 0;
    l_inv_err_sum = 0;
    r_x_weighted_sum = 0;
    r_y_weighted_sum = 0;
    r_inv_err_sum = 0;
    
    % -- CORE -------------------------------------------------------------
    
    % for each grid point, compute the squared error obtained by
    % approximating the leg with a circumference of radius R centered in
    % the grid point
    for i=1:N_GRIDX
        for j=1:N_GRIDY
            for k = 1 : height(ll_x)     % possible to vectorize computation
                l_error_ij(i,j) = l_error_ij(i,j) + (R - sqrt((l_gridx(i)-ll_x(k))^2 + (l_gridy(j)-ll_y(k))^2))^2 ;
            end
            for k = 1 : height(rl_x)     % possible to vectorize computation
                r_error_ij(i,j) = r_error_ij(i,j) + (R - sqrt((r_gridx(i)-rl_x(k))^2 + (r_gridy(j)-rl_y(k))^2))^2 ;
            end
        end
    end
    
    % weightted sum 
    for i = 1:N_GRIDX
        for j = 1:N_GRIDY
            l_y_weighted_pos_ij(i,j) = (l_error_ij(i,j)^-1) * l_gridy(j);
            l_y_weighted_sum = l_y_weighted_sum + l_y_weighted_pos_ij(i,j);
            
            l_x_weighted_pos_ij(i,j) = (l_error_ij(i,j)^-1) * l_gridx(i);
            l_x_weighted_sum = l_x_weighted_sum + l_x_weighted_pos_ij(i,j);

            l_inv_err_sum = l_inv_err_sum + l_error_ij(i,j)^-1;


            r_y_weighted_pos_ij(i,j) = r_error_ij(i,j)^-1 * r_gridy(j);
            r_y_weighted_sum = r_y_weighted_sum + r_y_weighted_pos_ij(i,j);
            
            r_x_weighted_pos_ij(i,j) = r_error_ij(i,j)^-1 * r_gridx(i);
            r_x_weighted_sum = r_x_weighted_sum + r_x_weighted_pos_ij(i,j);

            r_inv_err_sum = r_inv_err_sum + r_error_ij(i,j)^-1;
        end
    end
    
    % center computation
    c_x_l = l_x_weighted_sum / l_inv_err_sum;
    c_y_l = l_y_weighted_sum / l_inv_err_sum;
    c_x_r = r_x_weighted_sum / r_inv_err_sum;
    c_y_r = r_y_weighted_sum / r_inv_err_sum;

    % find points of best fitting circles (for plot purposes only)
    l_best_x = R*cos(sim_theta) + c_x_l;
    l_best_y = R*sin(sim_theta) + c_y_l;
    r_best_x = R*cos(sim_theta) + c_x_r;
    r_best_y = R*sin(sim_theta) + c_y_r;
    
    % output
    c_x = [c_x_l,c_x_r];
    c_y = [c_y_l,c_y_r];
    best_x = [l_best_x,r_best_x];
    best_y = [l_best_y,r_best_y];
    
end

