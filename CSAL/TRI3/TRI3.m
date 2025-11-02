%--------------------------------------------------------------------------
%            FEM 2D 3-Node Triangular (TRI3) Element Practice Script
%--------------------------------------------------------------------------
% PURPOSE:
%   To practice the 2D 3-node triangular (TRI3) finite element.
%
% DESCRIPTION:
%   This script reads FEM input data from a specified .txt file. It then 
%   performs finite element analysis based on the 2D TRI3 element 
%   formulation. Finally, it generates several plots to visualize the 
%   results.
%
% VARIABLES:
%   inputname : Name of the .txt file containing the mesh data
%   uscale    : Displacement scaling parameter for post-processing
%   ndpnd     : Number of degrees of freedom (DOFs) per node (2)
%   nnpe      : Number of nodes per element (3 for TRI3)
%   ndpe      : Number of DOFs per element (6 for TRI3)
%   nnpl      : Number of nodes per line (for traction)
%   ndim      : Number of spatial dimensions (2 for TRI3)
%   coord     : Nodal coordinates (nnd x 2) [m]
%   nnd       : Number of nodes
%   nel       : Number of elements
%   ndof      : Number of DOFs
%   f         : Global force vector (ndof x 1) [N]
%   u         : Global displacement vector (ndof x 1) [m]
%   K         : Global stiffness matrix (ndof x ndof) [N/m]
%   Ke        : Elemental stiffness matrix (ndpe x ndpe) [N/m]
%   fr        : Reduced nodal force vector (nfreedof x 1) [N]
%   ur        : Reduced displacement vector (nfreedof x 1) [m]
%   Kr        : Reduced stiffness matrix (nfreedof x nfreedof) [N/m]
%   E         : Elastic modulus [Pa]
%   nu        : Poisson's ratio
%   sigx      : Node-wise x-directional normal stress [Pa]
%   sigy      : Node-wise y-directional normal stress [Pa]
%   sigz      : Node-wise z-directional normal stress (plain strain) [Pa]
%   tauxy     : Node-wise xy-directional shear stress [Pa]
%   sigv      : Node-wise von-Mises stress [Pa]
%   coord_def : Nodal coordinates of the deformed configuration 
%               (nnd x 2) [m]
%
% AUTHOR:
%   Seung-Hoon Kang
%
%--------------------------------------------------------------------------

clear all; close all; clc

inputname = 'CSAL08_TRI3_input.txt'; % Input file name
uscale = 1.0; % Displacement scaling parameter (user-defined)
ndpnd  = 2; % Number of DOFs per node
nnpe   = 3; % Number of nodes per element
nnpl   = 2; % Number of nodes per line
ndim   = 2; % Number of spatial dimensions

%% 1. Scan Input
fprintf('Scanning the input file ...\n')

% 1-1. Initialize variables
coord          = [];                              % Nodal coordinates [m]
ELEMENT        = struct('con', {}, 'PROPid', {}); % Element data
% Property data
PROP           = struct('E', {}, 'nu', {}, 't', {}, 'type', {});
FORCE_NODE     = struct('NODEid', {}, 'val', {}); % Nodal force data
FORCE_TRACTION = struct('con', {}, 'val', {});    % Traction force data
FORCE_BODY     = struct('con', {}, 'val', {});    % Body force data
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
        elseif strcmpi(line, '*FORCE_TRACTION')
            current_section = 'FORCE_TRACTION';
        elseif strcmpi(line, '*FORCE_BODY')
            current_section = 'FORCE_BODY';    
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
            PROP(end+1).E  = prop_scan(1); % Elastic modulus [Pa]
            PROP(end).nu   = prop_scan(2); % Poisson's ratio
            PROP(end).t    = prop_scan(3); % Thickness [m]
            % Type (1:plane stress, 2:plane strain)
            PROP(end).type = round(prop_scan(4));
            if ne(PROP(end).type, 1) && ne(PROP(end).type, 2)
                error("Invalid type for PROP")
            end
        case 'FORCE_NODE' % Scan nodal force data
            force_scan = sscanf(line, '%f')';
            FORCE_NODE(end+1).NODEid = round(force_scan(1)); % NODE ID
            FORCE_NODE(end).val = force_scan(2:3); % Value [N]
        case 'FORCE_TRACTION' % Scan traction force data
            force_scan = sscanf(line, '%f')';
            % Connectivity
            FORCE_TRACTION(end+1).con = round(force_scan(1:nnpl));
            FORCE_TRACTION(end).val = force_scan(nnpl+1:nnpl+2); % Value [Pa]
        case 'FORCE_BODY' % Scan body force data
            force_scan = sscanf(line, '%f')';
            % Connectivity
            FORCE_BODY(end+1).con = round(force_scan(1:nnpe));
            FORCE_BODY(end).val = force_scan(nnpe+1:nnpe+2); % Value [N/m^3]    
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

for i = 1:length(FORCE_TRACTION)
    cone = FORCE_TRACTION(i).con; % Connectivity of current line segment
    found = false; % Flag to check if matching element is found
    for iel = 1:nel
        if (all(ismember(cone, ELEMENT(iel).con)))
            PROPid = ELEMENT(iel).PROPid; % PROP ID of current element
            t = PROP(PROPid).t; % Thickness of current element
            FORCE_TRACTION(i).t = t; % Thickness of current line segment
            found = true;
            break;
        end
    end
    if ~found
        error('No matching element found for FORCE_TRACTION on nodes %d-%d', ...
            cone(1), cone(2));
    end
end

for i = 1:length(FORCE_BODY)
    cone = FORCE_BODY(i).con; % Connectivity of current body element
    found = false; % Flag to check if matching element is found
    for iel = 1:nel
        if (all(ismember(cone, ELEMENT(iel).con)))
            PROPid = ELEMENT(iel).PROPid; % PROP ID of current element
            t = PROP(PROPid).t; % Thickness of current element
            FORCE_BODY(i).t = t; % Thickness of current body element
            found = true;
            break;
        end
    end
    if ~found
         error('No matching element found for FORCE_BODY on nodes %d-%d-%d', ...
             cone(1), cone(2), cone(3));
    end
end

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

for i = 1:length(FORCE_TRACTION)
    cone = FORCE_TRACTION(i).con; % Connectivity of current line segment
    val  = FORCE_TRACTION(i).val; % val(1): Tangential traction force (qt)
                                  % val(2): Normal traction force (qn)
    t    = FORCE_TRACTION(i).t;   % Thickness of current line segment

    lvec = coord(cone(2),:)-coord(cone(1),:); % Line vector
    L = norm(lvec); % Length of current line segment;

    c = lvec(1)/L;  % Cosine value
    s = lvec(2)/L;  % Sine value
    Tx = c*val(1) - s*val(2); % Tx = c*qt - s*qn
    Ty = s*val(1) + c*val(2); % Ty = s*qt + c*qn

    f(2*cone-1) = f(2*cone-1) + L*t*Tx/2.0; % Add traction force (fx)
    f(2*cone)   = f(2*cone)   + L*t*Ty/2.0; % Add traction force (fy)
end

for i = 1:length(FORCE_BODY)
    cone = FORCE_BODY(i).con; % Connectivity of current body element
    val  = FORCE_BODY(i).val; % val(1): x-directional body force (bx)
                              % val(2): y-directional body force (by)
    t    = FORCE_BODY(i).t;   % Thickness of current body element   

    coorde = coord(cone,:);  % Coordinate of current body element
    [A] = Area_TRI3(coorde); % Area of current body element;

    f(2*cone-1) = f(2*cone-1) + A*t*val(1)/3.0; % Add body force (fx)
    f(2*cone)   = f(2*cone)   + A*t*val(2)/3.0; % Add body force (fy)
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

    E    = PROP(PROPid).E;    % Elastic modulus of current element
    nu   = PROP(PROPid).nu;   % Poisson's ratio of current element
    t    = PROP(PROPid).t;    % Thickness of current element
    type = PROP(PROPid).type; % 2D stress state type of current element

    [A] = Area_TRI3(coorde); % Area of current element;
    [B] = Bmat_TRI3(coorde); % Strain-displacement matrix 
    [C] = Cmat_2D(E, nu, type); % Constitutive matrix
    
    % Elemental stiffness matrix
    Ke = B'*C*B*A*t;

    % Elemental DOF connectivity
    conde = [2*cone(1)-1 2*cone(1) 2*cone(2)-1 2*cone(2) 2*cone(3)-1 2*cone(3)];
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

% 2-6. Obtain stress
sigx  = zeros(nnd,1); % x-directional normal stress
sigy  = zeros(nnd,1); % y-directional normal stress
sigz  = zeros(nnd,1); % z-directional normal stress (plain strain)
tauxy = zeros(nnd,1); % xy-directional shear stress
sigv  = zeros(nnd,1); % von-Mises stress

for iel = 1:nel % Loop for each element
    cone = ELEMENT(iel).con; % Connectivity of current element
    coorde = coord(cone,:); % Coordinate of current element
    PROPid = ELEMENT(iel).PROPid; % PROP ID of current element

    E    = PROP(PROPid).E;    % Elastic modulus of current element
    nu   = PROP(PROPid).nu;   % Poisson's ratio of current element
    t    = PROP(PROPid).t;    % Thickness of current element
    type = PROP(PROPid).type; % 2D stress state type of current element

    [A] = Area_TRI3(coorde); % Area of current element;
    [B] = Bmat_TRI3(coorde); % Strain-displacement matrix 
    [C] = Cmat_2D(E, nu, type); % Constitutive matrix

    % Elemental DOF connectivity
    conde = [2*cone(1)-1 2*cone(1) 2*cone(2)-1 2*cone(2) 2*cone(3)-1 2*cone(3)];
    
    ue = u(conde); % Elemental displacement
    epse = B*ue;   % Elemental strain (epsx, epsy, gammaxy)
    sige = C*epse; % Elemental stress (sigx, sigy, tauxy)

    % Nodal stress averaging
    sigx(cone)  = sigx(cone)  + sige(1);
    sigy(cone)  = sigy(cone)  + sige(2);
    tauxy(cone) = tauxy(cone) + sige(3);

    % Elemental z-directional normal strain for plane strain state
    if eq(type,2)
        sigze = E*nu/(1+nu)/(1-2*nu)*(epse(1)+epse(2));
        % Nodal stress averaging
        sigz(cone) = sigz(cone) + sigze;
    end
end

% Nodal stress averaging
nepn = zeros(nnd,1); % Number of elements per node
for iel = 1:nel
    cone = ELEMENT(iel).con; % Connectivity of current element
    nepn(cone) = nepn(cone)+1;
end
sigx  = sigx./nepn;
sigy  = sigy./nepn;
sigz  = sigz./nepn;
tauxy = tauxy./nepn;

% Compute von-Mises stress
for ind = 1:nnd
    sigv(ind) = sqrt(0.5*((sigx(ind)-sigy(ind))^2 ...
        + (sigy(ind)-sigz(ind))^2 + (sigz(ind)-sigx(ind))^2) ...
        + 3*tauxy(ind)^2);
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
fprintf('- von-Mises stress\n')
fprintf('   Node ID            Sv\n')
for ind = 1:nnd
    fprintf('%10i%14.3e\n',iel,sigv(ind))
end

con = zeros(nel,nnpe); % Elemental connectivity matrix for visualization
for iel = 1:nel
    con(iel,:) = ELEMENT(iel).con; % Assign elemental connectivity
end

% 3-1. Display the input data (finite element mesh)
figure(1)
patch('Faces', con, 'Vertices', coord, 'LineWidth', 0.5, ...
    'Facecolor', 'c');
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
patch('Faces', con, 'Vertices', coord_def, 'FaceVertexCData', ux, ...
    'Facecolor', 'interp');
set(gca,'DataAspectRatio',[1 1 1])
title('X-displacement [m]')
set(gca, 'Fontsize', 16)
xlabel('X [m]')
ylabel('Y [m]')
colorbar
colormap jet

% 3-3. Display the Y-displacement on the deformed configuration
figure(3)
patch('Faces', con, 'Vertices', coord_def, 'FaceVertexCData', uy, ...
    'Facecolor', 'interp');
set(gca,'DataAspectRatio',[1 1 1])
title('Y-displacement [m]')
set(gca, 'Fontsize', 16)
xlabel('X [m]')
ylabel('Y [m]')
colorbar
colormap jet

% 3-4. Display the total displacement on the deformed configuration
figure(4)
patch('Faces', con, 'Vertices', coord_def, 'FaceVertexCData', ...
    sqrt(ux.^2+uy.^2), 'Facecolor', 'interp');
set(gca,'DataAspectRatio',[1 1 1])
title('Total displacement [m]')
set(gca, 'Fontsize', 16)
xlabel('X [m]')
ylabel('Y [m]')
colorbar
colormap jet

% 3-5. Display the von-Mises stress on the deformed configuration
figure(5)
patch('Faces', con, 'Vertices', coord_def, 'FaceVertexCData', sigv, ...
    'Facecolor', 'interp');
set(gca,'DataAspectRatio',[1 1 1])
title('von-Mises stress [Pa]')
set(gca, 'Fontsize', 16)
xlabel('X [m]')
ylabel('Y [m]')
colorbar
colormap jet

fprintf('Done displaying.\n')


function [A] = Area_TRI3(coorde)
% coorde : Nodal coordinate data of 3-node triangle (3 x 2)
% A      : Area of 3-node triangle

x = coorde(:,1);
y = coorde(:,2);

x21  = x(2)-x(1);
x13  = x(1)-x(3);
y12  = y(1)-y(2);
y31  = y(3)-y(1);

A = 0.5*abs(x21*y31-x13*y12);
end


function [B] = Bmat_TRI3(coorde)
% coorde : Nodal coordinate data of 3-node triangle (3 x 2)
% A      : Area of 3-node triangle
% B      : Strain-displacement matrix of 3-node triangle

x = coorde(:,1);
y = coorde(:,2);

x21  = x(2) - x(1);
x13  = x(1) - x(3);
x32  = x(3) - x(2);
y12  = y(1) - y(2);
y31  = y(3) - y(1);
y23  = y(2) - y(3);

A = 0.5*abs(x21*y31-x13*y12);
B = 1/(2*A)*[y23   0 y31   0 y12   0
               0 x32   0 x13   0 x21
             x32 y23 x13 y31 x21 y12];
end


function [C] = Cmat_2D(E, nu, type)
% E    : Elastic modulus
% nu   : Poisson's ratio
% type : 2D stress state type
% C    : 2D constitutive matrix

if eq(type,1)      % Plane stress
    C = E/(1-nu^2)*[ 1 nu 0
                    nu  1 0
                     0  0 (1-nu)/2];
elseif eq(type,2)  % Plane strain
    C = E/(1+nu)/(1-2*nu)*[1-nu   nu 0
                             nu 1-nu 0
                              0    0 (1-2*nu)/2];
end

end