%--------------------------------------------------------------------------
%                  FEM 1D Spring Element Practice Script
%--------------------------------------------------------------------------
% PURPOSE:
%   To practice the 1D spring finite element.
%
% DESCRIPTION:
%   This script reads FEM input data from a specified .txt file. It then 
%   performs finite element analysis based on the 1D spring element 
%   formulation. Finally, it generates several plots to visualize the 
%   results.
%
% VARIABLES:
%   inputname : Name of the .txt file containing the mesh data
%   uscale    : Displacement scaling parameter for post-processing
%   ndpnd     : Number of degrees of freedom (DOFs) per node
%   nnpe      : Number of nodes per element
%   ndim      : Number of spatial dimensions (e.g., 1 for 1D)
%   coord     : Nodal coordinates (nnd x 1) [m]
%   nnd       : Number of nodes
%   nel       : Number of elements
%   f         : Global force vector (nnd x 1) [N]
%   u         : Global displacement vector (nnd x 1) [m]
%   K         : Global stiffness matrix (nnd x nnd) [N/m]
%   Ke        : Elemental stiffness matrix (nnpe x nnpe) [N/m]
%   fr        : Reduced nodal force vector (nfreedof x 1) [N]
%   ur        : Reduced displacement vector (nfreedof x 1) [m]
%   Kr        : Reduced stiffness matrix (nfreedof x nfreedof) [N/m]
%   spcoef    : Spring coefficient [N/m]
%   coord_def : Nodal coordinates of the deformed configuration 
%               (nnd x 1) [m]
%
% AUTHOR:
%   Seung-Hoon Kang
%
%--------------------------------------------------------------------------

clear all; close all; clc

inputname = 'spring1d_input.txt'; % Input file name
uscale = 1.0; % Displacement scaling parameter (user-defined)
ndpnd  = 1; % Number of DOFs per node
nnpe   = 2; % Number of nodes per element
ndim   = 1; % Number of spatial dimensions

%% 1. Scan Input
fprintf('Scanning the input file ...\n')

% 1-1. Initialize variables
coord        = [];                              % Nodal coordinates [m]
ELEMENT      = struct('con', {}, 'PROPid', {}); % Element data
PROP         = struct('spcoef', {});            % Property data
FORCE_NODE   = struct('NODEid', {}, 'val', {}); % Nodal Force data
BC_PRECRIBED = struct('NODEid', {}, 'val', {}); % Prescribed boundary 
                                                % condition data
current_section = ''; % To track which part of the file we are in

% 1-2. Open the file
fid = fopen(inputname, 'r');

while 1
    % Read one line from the file and remove leading/trailing whitespace
    line = strtrim(fgetl(fid));

    % Skip empty lines and comment lines
    if isempty(line) || startsWith(line, '%')
        continue;
    end

    % 1-3. Check for section keywords
    if startsWith(line, '*')
        % Switch to the new section based on the keyword
        if strcmpi(line, '*NODE')
            current_section = 'NODE';
        elseif strcmpi(line, '*ELEMENT')
            current_section = 'ELEMENT';
        elseif strcmpi(line, '*PROP')
            current_section = 'PROP';
        elseif strcmpi(line, '*FORCE_NODE')
            current_section = 'FORCE_NODE';     
        elseif strcmpi(line, '*BC_PRESCRIBED')
            current_section = 'BC_PRESCRIBED';
        elseif strcmpi(line, '*EOF')
            % Stop reading at the End-Of-File marker
            break;
        end
        % Skip to the next line after processing a keyword
        continue;
    end

    % 1-4. Process data based on the current section
    switch current_section
        case 'NODE' % Scan nodal coordinates
            node_scan = sscanf(line, '%f')';
            coord = [coord; node_scan];
        case 'ELEMENT' % Scan element information
            elem_scan = sscanf(line, '%d')';
            ELEMENT(end+1).con = elem_scan(1:nnpe); % Connectivity
            ELEMENT(end).PROPid = elem_scan(nnpe+1); % PROP ID
        case 'PROP' % Scan material and physical properties
            prop_scan = sscanf(line, '%f')';
            PROP(end+1).spcoef = prop_scan; % Spring coefficient [N/m]
        case 'FORCE_NODE' % Scan nodal force data
            force_scan = sscanf(line, '%f')';
            FORCE_NODE(end+1).NODEid = round(force_scan(1)); % NODE ID
            FORCE_NODE(end).val = force_scan(2); % Value
        case 'BC_PRESCRIBED' % Scan prescribed boundary condition data
            bc_scan = sscanf(line, '%f')';
            BC_PRECRIBED(end+1).NODEid = round(bc_scan(1)); % NODE ID
            BC_PRECRIBED(end).val = bc_scan(2); % Value    
    end
end

nnd = size(coord,1); % Number of nodes
nel = length(ELEMENT); % Number of elements

% 1-5. Display the input data (finite element mesh)
fprintf('Done input scan.\n');
fprintf('- Number of nodes: %d\n', nnd)
fprintf('- Number of Elements: %d\n', nel)

%% 2. Numerical Solver
fprintf('Performing finite element analysis ...\n')
% 2-1. Initialize matrices and vectors
K = sparse(nnd,nnd); % Global stiffness matrix (in sparse storage format)

f = zeros(nnd,1); % Global force vector
for i = 1:length(FORCE_NODE)
    NODEid = FORCE_NODE(i).NODEid;
    val = FORCE_NODE(i).val;
    f(NODEid) = f(NODEid) + val; % Add nodal force
end

u = zeros(nnd,1); % Global displacement vector
bcdof = []; % ID set of constrained DOFs
for i = 1:length(BC_PRECRIBED)
    NODEid = BC_PRECRIBED(i).NODEid;
    val = BC_PRECRIBED(i).val;
    u(NODEid) = val; % Assign value of prescribed boundary condition 
    bcdof(end+1) = NODEid; % Append DOF ID.
end

% 2-2. Assemble stiffness matrix
for iel = 1:nel % Loop for each element
    cone = ELEMENT(iel).con; % Connectivity of current element
    PROPid = ELEMENT(iel).PROPid; % PROP ID of current element
    spcoef = PROP(PROPid).spcoef; % Spring coefficient of current element;
    Ke = spcoef*[1 -1;-1 1]; % Elemental stiffness matrix
    K(cone,cone) = K(cone,cone) + Ke; % Add elemental stiffness matrix
end

% 2-3. Obtain reduced matrices and vector (boundary condition)
freedof = setdiff(1:nnd, bcdof); % ID set of free DOFs
Kr = K(freedof,freedof); % Reduced stiffness matrix
fr = f(freedof) - K(freedof,bcdof)*u(bcdof); % Reduced force vector

% 2-4. Compute reduced displacement vector
ur = Kr\fr;

% 2-5. Obtain global displacement and force vectors
u(freedof) = ur;
f(bcdof) = K(bcdof,:)*u;

fprintf('Done finite element analysis.\n')

%% 3. Display the numerical results
fprintf('Displaying the numerical results ...\n')

fprintf('- Force vector\n')
fprintf('%14.3e\n',f)
fprintf('- Displacement vector\n')
fprintf('%14.3e\n',u)

coordy = zeros(nnd,1); % Y-coordinates for visualization
con = zeros(nel,nnpe); % Elemental connectivity matrix for visualization
for iel = 1:nel
    con(iel,:) = ELEMENT(iel).con; % Assign elemental connectivity
end

% 3-1. Display the input data (finite element mesh)
figure(1)
patch('Faces', con, 'Vertices', [coord, coordy], 'LineWidth', 3, ...
    'Marker', 'o', 'MarkerSize', 10);
set(gca,'DataAspectRatio',[1 1 1])
title('FE mesh')
set(gca, 'Fontsize', 16)
xlabel('X [m]')
ylabel('Y [m]')

% 3-2. Display the X-displacement on the deformed configuration
%   coord_def: Nodal coordinates of the deformed configuration [m]
coord_def = zeros(nnd, ndim);
coord_def(:) = coord + uscale*u; % x <- x0+u

figure(2)
patch('Faces', con, 'Vertices', [coord_def, coordy], 'LineWidth', 3, ...
    'FaceVertexCData', u, 'Edgecolor', 'interp', 'Marker', 'o', ...
    'MarkerSize', 10);
set(gca,'DataAspectRatio',[1 1 1])
title('X-displacement [m]')
set(gca, 'Fontsize', 16)
xlabel('X [m]')
ylabel('Y [m]')
colorbar


fprintf('Done displaying.\n')

