%--------------------------------------------------------------------------
%                FEM 2D 2-Node Bar Element Practice Script
%--------------------------------------------------------------------------
% PURPOSE:
%   To practice the 2D 2-node bar finite element.
%
% DESCRIPTION:
%   This script reads FEM input data from a specified .txt file. It then 
%   performs finite element analysis based on the 2D 2-node bar element 
%   formulation. Finally, it generates several plots to visualize the 
%   results.
%
% VARIABLES:
%   inputname : Name of the .txt file containing the mesh data
%   uscale    : Displacement scaling parameter for post-processing
%   ndpnd     : Number of degrees of freedom (DOFs) per node
%   nnpe      : Number of nodes per element
%   ndpe      : Number of DOFs per element
%   nnpl      : Number of nodes per line
%   ndim      : Number of spatial dimensions (e.g., 2 for 2D)
%   coord     : Nodal coordinates (nnd x 2) [m]
%   nnd       : Number of nodes
%   nel       : Number of elements
%   ndof      : Number of DOFs
%   f         : Global force vector (nnd x 1) [N]
%   u         : Global displacement vector (ndof x 1) [m]
%   K         : Global stiffness matrix (ndof x ndof) [N/m]
%   Ke        : Elemental stiffness matrix (ndpe x ndpe) [N/m]
%   fr        : Reduced nodal force vector (nfreedof x 1) [N]
%   ur        : Reduced displacement vector (nfreedof x 1) [m]
%   Kr        : Reduced stiffness matrix (nfreedof x nfreedof) [N/m]
%   E         : Elastic modulus [Pa]
%   A         : Cross-sectional area [m^2]
%   siga      : Element-wise axial stress
%   coord_def : Nodal coordinates of the deformed configuration 
%               (nnd x 2) [m]
%
% AUTHOR:
%   Seung-Hoon Kang
%
%--------------------------------------------------------------------------

clear all; close all; clc

inputname = 'CSAL05_bar2d2_input.txt'; % Input file name
uscale = 50.0; % Displacement scaling parameter (user-defined)
ndpnd  = 2; % Number of DOFs per node
nnpe   = 2; % Number of nodes per element
nnpl   = 2; % Number of nodes per line
ndim   = 2; % Number of spatial dimensions

%% 1. Scan Input
fprintf('Scanning the input file ...\n')

% 1-1. Initialize variables
coord         = [];                              % Nodal coordinates [m]
ELEMENT       = struct('con', {}, 'PROPid', {}); % Element data
PROP          = struct('E', {}, 'A', {});        % Property data
FORCE_NODE    = struct('NODEid', {}, 'val', {}); % Nodal force data
FORCE_LINE    = struct('con', {}, 'val', {});    % Line force data
% Prescribed boundary condition data
BC_PRESCRIBED = struct('NODEid', {}, 'val', {}, 'type', {}); 
                                                 
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
        elseif strcmpi(line, '*FORCE_LINE')
            current_section = 'FORCE_LINE';      
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
            PROP(end+1).E = prop_scan(1); % Elastic modulus [N/m]
            PROP(end).A = prop_scan(2); % Cross-sectional area [m^2]
        case 'FORCE_NODE' % Scan nodal force data
            force_scan = sscanf(line, '%f')';
            FORCE_NODE(end+1).NODEid = round(force_scan(1)); % NODE ID
            FORCE_NODE(end).val = force_scan(2:3); % Value [N]
        case 'FORCE_LINE' % Scan line force data
            force_scan = sscanf(line, '%f')';
            % Connectivity
            FORCE_LINE(end+1).con = round(force_scan(1:nnpl));
            FORCE_LINE(end).val = force_scan(nnpl+1); % Value [N/m]
        case 'BC_PRESCRIBED' % Scan prescribed boundary condition data
            bc_scan = sscanf(line, '%f')';
            BC_PRESCRIBED(end+1).NODEid = round(bc_scan(1)); % NODE ID
            BC_PRESCRIBED(end).val  = bc_scan(2); % Value [m]
            % Type (1:ux, 2:uy)
            BC_PRESCRIBED(end).type = round(bc_scan(3)); 
            if ne(BC_PRESCRIBED(end).type, 1) && ...
                    ne(BC_PRESCRIBED(end).type, 2)
                error("Invalid type for BC_PRESCRIBED")
            end
    end
end

nnd  = size(coord,1); % Number of nodes
nel  = length(ELEMENT); % Number of elements
ndof = ndpnd*nnd; % Number of total DOFs

% 1-5. Display the input data (finite element mesh)
fprintf('Done input scan.\n');
fprintf('- Number of nodes: %d\n', nnd)
fprintf('- Number of elements: %d\n', nel)
fprintf('- Number of total DOFs: %d\n', ndof)

%% 2. Numerical Solver
fprintf('Performing finite element analysis ...\n')
% 2-1. Initialize matrices and vectors
K = sparse(ndof,ndof); % Global stiffness matrix (in sparse storage format)

f = zeros(ndof,1); % Global force vector
for i = 1:length(FORCE_NODE)
    NODEid = FORCE_NODE(i).NODEid;
    val = FORCE_NODE(i).val;
    f(2*NODEid-1) = f(2*NODEid-1) + val(1); % Add nodal force (fx)
    f(2*NODEid)   = f(2*NODEid)   + val(2); % Add nodal force (fy)
end

for i = 1:length(FORCE_LINE)
    cone = FORCE_LINE(i).con;
    val  = FORCE_LINE(i).val;
    lvec = coord(cone(2),:)-coord(cone(1),:); % Line vector
    f(2*cone-1) = f(2*cone-1) + 0.5*val*lvec(1); % Add line force (fx)
    f(2*cone)   = f(2*cone)   + 0.5*val*lvec(2); % Add line force (fy)
end

u = zeros(ndof,1); % Global displacement vector
bcdof = []; % ID set of constrained DOFs
for i = 1:length(BC_PRESCRIBED)
    NODEid = BC_PRESCRIBED(i).NODEid;
    val = BC_PRESCRIBED(i).val;
    type = BC_PRESCRIBED(i).type;
    if eq(type,1)     % BC type: ux
        u(2*NODEid-1) = val; % Assign value of prescribed BC 
        bcdof(end+1)  = 2*NODEid-1; % Append DOF ID.
    elseif eq(type,2) % BC type: uy
        u(2*NODEid)   = val; % Assign value of prescribed BC 
        bcdof(end+1)  = 2*NODEid;   % Append DOF ID.
    end
end

% 2-2. Assemble stiffness matrix
for iel = 1:nel % Loop for each element
    cone = ELEMENT(iel).con; % Connectivity of current element
    coorde = coord(cone,:); % Coordinate of current element
    PROPid = ELEMENT(iel).PROPid; % PROP ID of current element

    E = PROP(PROPid).E; % Elastic modulus of current element;
    A = PROP(PROPid).A; % Cross-sectional area of current element;

    lvec = coorde(2,:)-coorde(1,:); % Line vector
    L = norm(lvec); % Length of current element;
    c = lvec(1)/L;  % Cosine value
    s = lvec(2)/L;  % Sine value

    % Transformation matrix
    T = [ c  s  0  0
         -s  c  0  0
          0  0  c  s
          0  0 -s  c];

    % Elemental stiffness matrix (local coordinates)
    Keloc = (E*A/L)*[ 1  0 -1  0
                      0  0  0  0
                     -1  0  1  0
                      0  0  0  0];

    % Elemental stiffness matrix (global coordinates)
    Ke = T'*Keloc*T; 

    % Elemental DOF connectivity
    conde = [2*cone(1)-1 2*cone(1) 2*cone(2)-1 2*cone(2)];
    % Add elemental stiffness matrix
    K(conde,conde) = K(conde,conde) + Ke;
end

% 2-3. Obtain reduced matrices and vector (boundary condition)
freedof = setdiff(1:ndof, bcdof); % ID set of free DOFs
Kr = K(freedof,freedof); % Reduced stiffness matrix
fr = f(freedof) - K(freedof,bcdof)*u(bcdof); % Reduced force vector

% 2-4. Compute reduced displacement vector
ur = Kr\fr;

% 2-5. Obtain global displacement and force vectors
u(freedof) = ur;
f(bcdof) = K(bcdof,:)*u;

fx = f(1:2:ndof-1); fy = f(2:2:ndof); % Directional force
ux = u(1:2:ndof-1); uy = u(2:2:ndof); % Directional displacement

% 2-6. Obtain axial stress
siga = zeros(nel,1); % Element-wise axial stress
for iel = 1:nel % Loop for each element
    cone = ELEMENT(iel).con; % Connectivity of current element
    coorde = coord(cone,:); % Coordinate of current element
    PROPid = ELEMENT(iel).PROPid; % PROP ID of current element

    E = PROP(PROPid).E; % Elastic modulus of current element;

    lvec = coorde(2,:)-coorde(1,:); % Line vector
    L = norm(lvec); % Length of current element;
    c = lvec(1)/L;  % Cosine value
    s = lvec(2)/L;  % Sine value

    % Transformation matrix
    T = [c s 0 0
         0 0 c s];

    % Elemental DOF connectivity
    conde = [2*cone(1)-1 2*cone(1) 2*cone(2)-1 2*cone(2)];
    % Elemental displacement (global coordinates)
    ue = u(conde);

    % Elemental axial displacement (local coordinates)
    uxeloc = T*ue;
    epsxeloc = (uxeloc(2)-uxeloc(1))/L; % Elemental axial strain
    sigxeloc = E*epsxeloc; % Elemental axial stress

    siga(iel) = sigxeloc;
end

fprintf('Done finite element analysis.\n')

%% 3. Display the numerical results
fprintf('Displaying the numerical results ...\n')

fprintf('- Force\n')
fprintf('   Node ID            Fx            Fy\n')
for ind = 1:nnd
    fprintf('%10i%14.3e%14.3e\n',ind,fx(ind),fy(ind))
end
fprintf('- Displacement\n')
fprintf('   Node ID            Ux            Uy\n')
for ind = 1:nnd
    fprintf('%10i%14.3e%14.3e\n',ind,ux(ind),uy(ind))
end
fprintf('- Axial stress\n')
fprintf('  Elem. ID            Sa\n')
for iel = 1:nel
    fprintf('%10i%14.3e\n',iel,siga(iel))
end

con = zeros(nel,nnpe); % Elemental connectivity matrix for visualization
for iel = 1:nel
    con(iel,:) = ELEMENT(iel).con; % Assign elemental connectivity
end

% 3-1. Display the input data (finite element mesh)
figure(1)
patch('Faces', con, 'Vertices', coord, 'LineWidth', 3, ...
    'Marker', 'o', 'MarkerSize', 10);
set(gca,'DataAspectRatio',[1 1 1])
title('FE mesh')
set(gca, 'Fontsize', 16)
xlabel('X [m]')
ylabel('Y [m]')

% 3-2. Display the X-displacement on the deformed configuration
%   coord_def: Nodal coordinates of the deformed configuration [m]
coord_def = zeros(nnd, ndim);
coord_def(:,1) = coord(:,1) + uscale*ux; % x <- x0+ux
coord_def(:,2) = coord(:,2) + uscale*uy; % y <- x0+uy

figure(2)
patch('Faces', con, 'Vertices', coord_def, 'LineWidth', 3, ...
    'FaceVertexCData', ux, 'Edgecolor', 'interp', 'Marker', 'o', ...
    'MarkerSize', 10);
set(gca,'DataAspectRatio',[1 1 1])
title('X-displacement [m]')
set(gca, 'Fontsize', 16)
xlabel('X [m]')
ylabel('Y [m]')
colorbar
colormap jet

% 3-3. Display the Y-displacement on the deformed configuration
figure(3)
patch('Faces', con, 'Vertices', coord_def, 'LineWidth', 3, ...
    'FaceVertexCData', uy, 'Edgecolor', 'interp', 'Marker', 'o', ...
    'MarkerSize', 10);
set(gca,'DataAspectRatio',[1 1 1])
title('Y-displacement [m]')
set(gca, 'Fontsize', 16)
xlabel('X [m]')
ylabel('Y [m]')
colorbar
colormap jet

% 3-4. Display the total displacement on the deformed configuration
coord_def_xsym = coord_def;
coord_def_xsym(:,1) = -coord_def(:,1);

figure(4)
patch('Faces', con, 'Vertices', coord_def, 'LineWidth', 3, ...
    'FaceVertexCData', sqrt(ux.^2+uy.^2), 'Edgecolor', 'interp', ...
    'Marker', 'o', 'MarkerSize', 10);
hold on
patch('Faces', con, 'Vertices', coord_def_xsym, 'LineWidth', 3, ...
    'FaceVertexCData', sqrt(ux.^2+uy.^2), 'Edgecolor', 'interp', ...
    'Marker', 'o', 'MarkerSize', 10);
set(gca,'DataAspectRatio',[1 1 1])
title('Total displacement [m]')
set(gca, 'Fontsize', 16)
xlabel('X [m]')
ylabel('Y [m]')
colorbar
colormap jet


fprintf('Done displaying.\n')
