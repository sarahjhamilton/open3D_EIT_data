% =========================================================================
% FILENAME: ACT5_box_eidors_sim.m
%
% USE: To solve the forward EIT problem to simulate voltage data for the
%           3D "box" setup used in the manuscript:
%
%       "Fast absolute 3d CGO-based electrical impedance tomography on
%          experimental tank data", Physiological Measurement 43 (2022),
%          no. 12, 124001 written by Sarah J Hamilton, P A Muller,
%          D Isaacson, V Kolehmainen, J Newell, O Rajabi Shishvan,
%          G Saulnier, and J Toivanen.
%
%
%      It is currently set to have lengths in cm and
%      conductivities in S/cm. There is a unit conversion variable to
%      adjust the conductivity values.
%
%
% Written by: Peter Muller, Villanova University May 28, 2021
%               and Sarah Hamilton, Marquette University
% Last Updated: PM, June 8, 2021 (Checks if EIDORS has been started and starts it if not)
% =========================================================================
clear;
%% ==================== User input ========================================
% Run EIDORS if it has not been started yet
if exist('eidors','dir')==0
    %     run 'YOUR PATH TO EIDORS HERE'    
end


ld_mesh             = 1; % If you want to just load the mesh in rather than recompute it
mesh_name           = 'ACT5_box_mesh';

% The following toggle is just for debugging
volt_solve_bg1      = 1; % 1 for running forward solver for the homogeneous 1S/m problem, 0 for no forward solve
volt_solve_hom      = 1; % 1 for running forward solver for the homogeneous problem, 0 for no forward solve
volt_solve_inhom    = 1; % 1 for running forward solver for inhomogeneous data

sv_hom_mesh         = 0; % Toggle 1 to save homogeneous mesh, 0 for no saving
sv_bg1_data         = 1; % Toggle 1 to save bg 1S/m voltages, 0 for no saving
sv_hom_data         = 1; % Toggle 1 to save inhomg voltages, 0 for no saving
sv_inhom_data       = 1; % Toggle 1 to save inhomg voltages, 0 for no saving

% Set electrode dimensions
elec_width          = 8; % width in cm, assumes the electrodes are square, but
%                   a different height could be accounted for with adjustments
% Show electrodes in show_fem_enhanced
show_electrodes     = false; % true shows electrodes as green, false doesn't show them
% Show mesh edges in show_fem_enhanced
show_mesh_edge      = false; % true shows mesh edges, false doesn't show them

% Set maximum edge length for the mesh for the different pieces
MaxEdge_dom         = 1; % max edge length for domain
MaxEdge_elecs       = 0.2; % max edge length for electrodes,
MaxEdge_incl        = .5; % max edge length for inclusions
% dom=1, elec=.2, incl=.5 yields ~180k nodes, ~864k elems, ~1.9k nodes/elec

cond_sc             = 10^-5; % Change accordingly. 10^-5 converts mS/m to S/cm
% Units will be converted to S/cm later using cond_sc
contact_imp         = .01; % Set the contact impedance. default in EIDORS is .01

% ------------------------------------------------------------------------
% Phantom/Example to simulate:
% ------------------------------------------------------------------------
% Example 1: Two targets example
sph_inc_loc         = [-2.25,-6.5,4;
    2.25,6.5,0];
sph_inc_rad         = [3;3];
sph_inc_cond        = [300;60];
bg_conduct          = 100;
ex_name             = 'Ex01'; % name of your example

% Convert conductivities to S/cm
bg_conduct          = bg_conduct*cond_sc; % background conductivity
inc_cond            = sph_inc_cond*cond_sc; % inclusion conductivities



% % Ellipsoid inclusion conductivities, set as [] for no ellipsoids
% ell_inc_cond        = []; % mS/m
%
% % Cylinder inclusion conductivities, set as [] for no cylinders
% cyl_inc_cond        = [];
%
% % Specify ellipsoidal center location. Each row should be the x,y,z (in cm)
% % of a single inclusion (put in same order as ellipsoid conductivities)
% ell_inc_loc         = [];
%
% % Specify ellipsodial inclusion semi-axes directions
% % x1 direction
% ell_inc_ax_a        = [];% Each row should be a vector defining the direction of the x1 semi-axis
%
% % x2 direction
% ell_inc_ax_b        = [];% Each row should be a vector defining the direction of the x2 semi-axis
%
% % x3 direction
% ell_inc_ax_c        = [];% Each row should be a vector defining the direction of the x3 semi-axis
%
% % Specify cylinder inclusion top/bottom center locations (put in same order as cylinder conductivities)
% cyl_inc_loc_top     = [];% Each row should be the x,y,z (in cm)
%
% cyl_inc_loc_bot     = []; % Each row should be the x,y,z (in cm)
%
% % Specify cylinder inclusion radii (in cm) in a vector (put in same
% % order as cylinder conductivities)
% cyl_inc_rad         = [];
%
% Specify current patterns to use for simulation
curr_type           = 'opt';
% CP Choices: 'haar','trig','opt',

% Specify where the electrode location information is in an Excel file
filename            = 'ElectrodeLayout.xlsx'; % Excel file's name with electrode locations
sheet               = 'Sheet1'; % Sheet name in Excel file
xLrange             = 'B3:D34';% Top left corner:Bottom right corner of data cells with electrode locations

% ===================== End User Input ====================================
%%
% Extract electrode centers from Excel file
elecs               = readmatrix(filename,'Sheet',sheet,'Range',xLrange);

%%%%% IMPORTANT NOTE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This mesh construction will place coordinates in [x1,x2,x3] order
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ld_mesh
    load([mesh_name '.mat'],'fmdl','body_geometry','electrode_geometry','mesh_info');
else
    % Identify the coordinate of each face of the box
    left_face       = min(elecs(:,2));
    bot_face        = min(elecs(:,3));
    back_face       = min(elecs(:,1));
    right_face      = max(elecs(:,2));
    top_face        = max(elecs(:,3));
    front_face      = max(elecs(:,1));
    back_left_bot   = [back_face,left_face,bot_face];
    front_right_top = [front_face,right_face,top_face];
    
    % Set up domain geometry (same for all inclusion types)
    body_geometry.ortho_brick.opposite_corner_a = back_left_bot;
    body_geometry.ortho_brick.opposite_corner_b = front_right_top;
    body_geometry.max_edge_length               = MaxEdge_dom;
    
    % Get electrode indices for each face
    top_elecs       = find(elecs(:,3)==top_face);%[5 6 13 14 21 22];
    bot_elecs       = find(elecs(:,3)==bot_face);%[9 10 17 18 25 26];
    left_elecs      = find(elecs(:,2)==left_face);%1:4;
    right_elecs     = find(elecs(:,2)==right_face);%29:32;
    front_elecs     = find(elecs(:,1)==front_face);%[11 12 19 20 27 28];
    back_elecs      = find(elecs(:,1)==back_face);%[7 8 15 16 23 24];
    
    % Define electrode corners (diametrically opposed) for each face
    % Electrodes need thickness, too (NOT SURE IF THIS IS CORRECTLY DONE)
    el_half_width   = elec_width/2;
    
    electrode_geometry      = cell(size(elecs,1),1);
    for num = 1:length(top_elecs)
        top_el_num          = top_elecs(num);
        elec_loc_a          = [elecs(top_el_num,1)-el_half_width,elecs(top_el_num,2)-el_half_width,elecs(top_el_num,3)];
        elec_loc_b          = [elecs(top_el_num,1)+el_half_width,elecs(top_el_num,2)+el_half_width,elecs(top_el_num,3)+.1];
        
        electrode_geometry{top_el_num}.ortho_brick.opposite_corner_a    = elec_loc_a;
        electrode_geometry{top_el_num}.ortho_brick.opposite_corner_b    = elec_loc_b;
        electrode_geometry{top_el_num}.max_edge_length                  = MaxEdge_elecs;
    end
    for num = 1:length(bot_elecs)
        bot_el_num          = bot_elecs(num);
        elec_loc_a          = [elecs(bot_el_num,1)-el_half_width,elecs(bot_el_num,2)-el_half_width,elecs(bot_el_num,3)-.1];
        elec_loc_b          = [elecs(bot_el_num,1)+el_half_width,elecs(bot_el_num,2)+el_half_width,elecs(bot_el_num,3)];
        
        electrode_geometry{bot_el_num}.ortho_brick.opposite_corner_a    = elec_loc_a;
        electrode_geometry{bot_el_num}.ortho_brick.opposite_corner_b    = elec_loc_b;
        electrode_geometry{bot_el_num}.max_edge_length                  = MaxEdge_elecs;
    end
    for num = 1:length(left_elecs)
        left_el_num         = left_elecs(num);
        elec_loc_a          = [elecs(left_el_num,1)-el_half_width,elecs(left_el_num,2)-.1,elecs(left_el_num,3)-el_half_width];
        elec_loc_b          = [elecs(left_el_num,1)+el_half_width,elecs(left_el_num,2),elecs(left_el_num,3)+el_half_width];
        
        electrode_geometry{left_el_num}.ortho_brick.opposite_corner_a   = elec_loc_a;
        electrode_geometry{left_el_num}.ortho_brick.opposite_corner_b   = elec_loc_b;
        electrode_geometry{left_el_num}.max_edge_length                 = MaxEdge_elecs;
    end
    for num = 1:length(right_elecs)
        right_el_num        = right_elecs(num);
        elec_loc_a          = [elecs(right_el_num,1)-el_half_width,elecs(right_el_num,2),elecs(right_el_num,3)-el_half_width];
        elec_loc_b          = [elecs(right_el_num,1)+el_half_width,elecs(right_el_num,2)+.1,elecs(right_el_num,3)+el_half_width];
        
        electrode_geometry{right_el_num}.ortho_brick.opposite_corner_a  = elec_loc_a;
        electrode_geometry{right_el_num}.ortho_brick.opposite_corner_b  = elec_loc_b;
        electrode_geometry{right_el_num}.max_edge_length                = MaxEdge_elecs;
    end
    for num = 1:length(back_elecs)
        back_el_num         = back_elecs(num);
        elec_loc_a          = [elecs(back_el_num,1)-.1,elecs(back_el_num,2)-el_half_width,elecs(back_el_num,3)-el_half_width];
        elec_loc_b          = [elecs(back_el_num,1),elecs(back_el_num,2)+el_half_width,elecs(back_el_num,3)+el_half_width];
        
        electrode_geometry{back_el_num}.ortho_brick.opposite_corner_a   = elec_loc_a;
        electrode_geometry{back_el_num}.ortho_brick.opposite_corner_b   = elec_loc_b;
        electrode_geometry{back_el_num}.max_edge_length                 = MaxEdge_elecs;
    end
    for num = 1:length(front_elecs)
        front_el_num        = front_elecs(num);
        elec_loc_a          = [elecs(front_el_num,1),elecs(front_el_num,2)-el_half_width,elecs(front_el_num,3)-el_half_width];
        elec_loc_b          = [elecs(front_el_num,1)+.1,elecs(front_el_num,2)+el_half_width,elecs(front_el_num,3)+el_half_width];
        
        electrode_geometry{front_el_num}.ortho_brick.opposite_corner_a  = elec_loc_a;
        electrode_geometry{front_el_num}.ortho_brick.opposite_corner_b  = elec_loc_b;
        electrode_geometry{front_el_num}.max_edge_length                = MaxEdge_elecs;
    end
    
    % For visualization only, create a mesh without electrodes, since they
    % obscure the interior
    % fmdl_no_elecs = ng_mk_geometric_models(body_geometry);
    % Create mesh with electrodes and refinements for inclusions
    fmdl = ng_mk_geometric_models(body_geometry,electrode_geometry);
    
    % Account for contact impedance
    for el_num=1:size(elecs,1)
        fmdl.electrode(el_num).z_contact=contact_imp;
    end
    mesh_info.electrodeCenters=elecs;
    mesh_info.electrodeWidth=[num2str(elec_width),'cm'];
    
    save('mesh_new.mat','fmdl','body_geometry','electrode_geometry','mesh_info');
end

lightGrey = 0.9*[1 1 1]; % It looks better if the lines are lighter
opts.edge.color=lightGrey;
opts.electrode.display.show=show_electrodes;
opts.edge.interior.show=show_mesh_edge;
%--------------------------------------------------------------------------
% Identify element center coordinates for placing inclusions.
%--------------------------------------------------------------------------
m_elem      = fmdl.elems;
m_nodes     = fmdl.nodes;

numNodes = size(fmdl.nodes,1);
numElems = size(fmdl.elems,1);

m_x         = zeros(numElems,1);
m_y         = m_x;
m_z         = m_x;

for q = 1:numElems
    m_x(q)  = mean(m_nodes(m_elem(q,:),1));
    m_y(q)  = mean(m_nodes(m_elem(q,:),2));
    m_z(q)  = mean(m_nodes(m_elem(q,:),3));
end
% -------------------------------------------------------------------------
% String for loading current patterns
curr_pat_mat    = ['curr_',curr_type,'.mat'];

% Extract the number of electrodes and number of current patterns
L           = size(elecs,1); % Number of electrodes
num_pats    = L-1; % Number of current patterns to be applied

% Load current patterns used by ACT5
load(curr_pat_mat); % From ACT5 should be saved as current_pattern
cur_pattern = current_patterns; clear current_patterns

% Set stimulation patterns
for pat = num_pats:-1:1 % Run backwards so size is determined on first step
    stim(pat).stimulation       = 'Amp'; % Quantity stimulated into electrodes
    stim(pat).meas_pattern      = eye(L); % Matrix representing the measurement patterns for stimulation i. Each column of this matrix represents the amplification of the signal at each electrode for a single measurement pattern. From: http://tudr.thapar.edu:8080/jspui/bitstream/10266/2720/4/2720.pdf
    stim(pat).stim_pattern      = cur_pattern(:,pat);% Vector of the stimulation quantity applied to each electrode during each pattern. From: http://tudr.thapar.edu:8080/jspui/bitstream/10266/2720/4/2720.pdf
end
%%
save_name       = [mesh_name,'_',ex_name];
% -------------------------------------------------------------------------
% Simulate homogeneous data:
% -------------------------------------------------------------------------
if volt_solve_bg1
    % Assign background conductivity of 1 S/m to homogeneous mesh
    img_1                       = mk_image(fmdl, 1000*cond_sc); % background of 1S/m
    cond_bg1                    = 100*cond_sc;
    % Add stimulation patterns to the forward model
    img_1.fwd_model.stimulation = stim;
    % Solve forward model (system_mat_fields.m line 58 says this uses CEM.)
    save_bg1_name               = ['volts_' save_name '_hom1Sm_w_CP',curr_type];
    fprintf('Starting BG 1 forward solver. \n')
    homog1_start                = tic;
    V_1                         = fwd_solve(img_1);
    homog1_solve_time           = toc(homog1_start);
    homog1_solve_time_disp      = [num2str(homog1_solve_time),' seconds'];
    if homog1_solve_time>60
        homog1_solve_time_disp  = [num2str(homog1_solve_time/60),' minutes'];
    end
    disp(['BG 1 forward solver done in ',homog1_solve_time_disp,'!'])
    % Reshape the voltages
    voltages                    = reshape(V_1.meas,L,num_pats);
    if sv_bg1_data
        save(save_bg1_name,'voltages','cur_pattern','cond_bg1');
    end
end

% -------------------------------------------------------------------------
% Setup to simulate inhomogeneous data
% -------------------------------------------------------------------------
condBg          = bg_conduct;
sigmaTrue       = condBg*ones(numElems,1);
numInclusions   = length(sph_inc_rad);

% Get ready to solve the homogenous background problem:
if volt_solve_hom
    img_homog       = mk_image(fmdl,sigmaTrue);
    
    % Add current patterns
    img_homog.fwd_model.stimulation=stim;
    
    fprintf('Starting homogeneous forward solver. \n')
    homog_start             = tic;
    V_hom                   = fwd_solve(img_homog);
    homog_solve_time        = toc(homog_start);
    homog_solve_time_disp   = [num2str(homog_solve_time),' seconds'];
    if homog_solve_time>60
        homog_solve_time_disp   = [num2str(homog_solve_time/60),' minutes'];
    end
    disp(['Homogeneous forward solver done in ',homog_solve_time_disp,'!'])
    % Reshape the voltages
    Vh=reshape(V_hom.meas,L,num_pats);
    
    if sv_hom_data        
        save(['volts_' save_name '_hom_w_CP' curr_type '.mat'],'Vh','cur_pattern','condBg');
    end
end

% Add spherical inclusions:
for ii = 1:numInclusions
    incInds = ((m_x - sph_inc_loc(ii,1)).^2 + ...
        (m_y - sph_inc_loc(ii,2)).^2 + ...
        (m_z - sph_inc_loc(ii,3)).^2) <  (sph_inc_rad(ii)^2);
    
    sigmaTrue(incInds) = inc_cond(ii);
end
% -------------------------------------------------------------------------
% Plot the inhomogeneous example:
% -------------------------------------------------------------------------
figure
img_inhomog     = mk_image(fmdl,sigmaTrue);
show_fem_enhanced(img_inhomog,opts)
view(50,30)
% -------------------------------------------------------------------------
% Solve for the voltages for the inhomogeneous example:
% -------------------------------------------------------------------------
if volt_solve_inhom
    % Add current patterns
    img_inhomog.fwd_model.stimulation=stim;
    
    % Begin forward solve
    fprintf('Starting inhomogeneous forward solver')
    inhomog_start           = tic;
    V_inhom                 = fwd_solve(img_inhomog);
    inhomog_solve_time      = toc(inhomog_start);
    inhomog_solve_time_disp = [num2str(inhomog_solve_time),' seconds'];
    if inhomog_solve_time>60
        inhomog_solve_time_disp=[num2str(inhomog_solve_time/60),' minutes'];
    end
    disp(['Inhomogeneous forward solver done in ',inhomog_solve_time_disp,' !'])
    Vi=reshape(V_inhom.meas,L,num_pats);
    
    if sv_inhom_data
        save(['volts_' save_name '_w_CP' curr_type '.mat'],'Vi','cur_pattern','condBg');
        % Sarah
        sigma_details = struct;
        sigma_details.centers = sph_inc_loc;
        sigma_details.radii = sph_inc_rad;
        sigma_details.cond = sph_inc_cond;
        sigma_details.condbg = condBg;
        sigma_details.truth = sigmaTrue;
        sigma_details.scFactor = cond_sc;
        save(['sigma_' save_name '.mat'],'sigma_details');
    end
end
%--------------------------------------------------------------------------
