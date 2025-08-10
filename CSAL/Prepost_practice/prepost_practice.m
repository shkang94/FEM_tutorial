%--------------------------------------------------------------------------
%            FEM Pre- and Post-Processing Practice Script
%--------------------------------------------------------------------------
% PURPOSE:
%   To practice the pre-processing (reading input) and post-processing
%   (visualizing results) stages of a Finite Element Method (FEM) workflow.
%
% DESCRIPTION:
%   This script reads FEM mesh data (nodal coordinates and element 
%   connectivity) from a specified .txt file. It then loads pre-computed 
%   results from a specified .mat file. Finally, it generates several plots 
%   to visualize the results.
%
%   Note: This script does not contain the numerical solver itself.
%
% VARIABLES:
%   inputname : Name of the .txt file containing the mesh data
%   uscale    : Displacement scaling parameter for post-processing
%   ndpnd     : Number of degrees of freedom (DOFs) per node
%   nnpe      : Number of nodes per element
%   ndim      : Number of spatial dimensions (e.g., 2 for 2D)
%   coord     : Nodal coordinates (nnd x ndim) [m]
%   con       : Elemental connectivity (nel x nnpel)
%   nnd       : Number of nodes
%   nel       : Number of elements
%   ux        : Nodal x-displacement (nnd x 1) [m]
%   uy        : Nodal y-displacement (nnd x 1) [m]
%   sigv      : Nodal von-Mises stress (nnd x 1) [Pa]
%   coord_def : Nodal coordinates of the deformed configuration 
%               (nnd x ndim) [m]
%
% AUTHOR:
%   Seung-Hoon Kang (seunghoon.kang94@gmail.com)
%
% LAST MODIFIED:
%   2025-08-10
%
%--------------------------------------------------------------------------

clear all; close all; clc

inputname = 'prepost_practice_input.txt'; % Input file name
uscale = 1.0; % Displacement scaling parameter (user-defined)
ndpnd  = 2; % Number of DOFs per node
nnpe   = 4; % Number of nodes per element (Q4)
ndim   = 2; % Number of spatial dimensions

%% 1. Scan Input
fprintf('Scanning the input file ...\n')
% 1-1. Initialize variables
coord = []; % Nodal coordinates [m]
con   = []; % Elemental connectivity
current_section = ''; % To track which part of the file we are in

% 1-2. Open the file
fid = fopen(inputname, 'r');

while 1
    % Read one line from the file and remove leading/trailing whitespace
    line = strtrim(fgetl(fid));

    % Skip empty lines
    if isempty(line)
        continue;
    end

    % 1-3. Check for section keywords
    if startsWith(line, '*')
        % Switch to the new section based on the keyword
        if strcmpi(line, '*NODE')
            current_section = 'NODE';
        elseif strcmpi(line, '*ELEMENT')
            current_section = 'ELEMENT';
        elseif strcmpi(line, '*EOF')
            % Stop reading at the End-Of-File marker
            break;
        end
        % Skip to the next line after processing a keyword
        continue;
    end

    % 1-4. Process data based on the current section
    switch current_section
        case 'NODE'
            % Convert the line of text into a row of numbers (coordinates)
            node_scan = sscanf(line, '%f')';
            coord = [coord; node_scan]; % Append to coord
        case 'ELEMENT'
            % Convert the line of text into a row of numbers (connectivity)
            elem_scan = sscanf(line, '%d')';
            con = [con; elem_scan]; % Append to con
    end
end

nnd = size(coord,1); % Number of nodes
nel = size(con,1);   % Number of elements

% 1-5. Display the input data (finite element mesh)
fprintf('Done input scan.\n');
fprintf('- Number of nodes: %d\n', nnd)
fprintf('- Number of Elements: %d\n', nel)

%% 2. Numerical Solver
% Load pre-computed solver results.
%   ux: Nodal x-displacement [m]
%   uy: Nodal y-displacement [m]
%   sigv: Nodal von-Mises stress [Pa]
load('prepost_practice_results.mat')

%% 3. Display the numerical results
fprintf('Displaying the numerical results ...\n')

% 3-1. Display the input data (finite element mesh)
figure(1)
patch('Faces',con,'Vertices',coord,'Facecolor','c','LineWidth',1);
set(gca,'DataAspectRatio',[1 1 1])
title('FE mesh')
set(gca, 'Fontsize', 16)
xlabel('X [m]')
ylabel('Y [m]')

% 3-2. Display the deformed configuration
%   coord_def: Nodal coordinates of the deformed configuration [m]
coord_def = zeros(nnd, ndim);
coord_def(:,1) = coord(:,1) + uscale*ux; % x <- x0+ux
coord_def(:,2) = coord(:,2) + uscale*uy; % y <- y0+uy

figure(2)
patch('Faces',con,'Vertices',coord_def,'Facecolor','c', 'LineWidth',1);
set(gca,'DataAspectRatio',[1 1 1])
title('Deformed configuration')
set(gca, 'Fontsize', 16)
xlabel('X [m]')
ylabel('Y [m]')

% 3-3. Display the y-displacement on the deformed configuration
figure(3)
patch('Faces',con,'Vertices',coord_def,'FaceVertexCData',uy,...
    'Facecolor','interp','LineWidth',1);
set(gca,'DataAspectRatio',[1 1 1])
title('Y-displacement [m]')
set(gca, 'Fontsize', 16)
xlabel('X [m]')
ylabel('Y [m]')
colorbar

% 3-4. Display the von-Mises stress on the deformed configuration
figure(4)
patch('Faces',con,'Vertices',coord_def,'FaceVertexCData',sigv,...
    'Facecolor','interp','LineWidth',1);
set(gca,'DataAspectRatio',[1 1 1])
title('von-Mises stress [Pa]')
set(gca, 'Fontsize', 16)
xlabel('X [m]')
ylabel('Y [m]')
colorbar

fprintf('Done displaying.\n')