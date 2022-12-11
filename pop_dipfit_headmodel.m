% pop_dipfit_headmodel() - generate headmodel from MRI
%
% Usage: 
%  >> OUTEEG = pop_dipfit_headmodel( INEEG ); 
%
% Inputs:
%   INEEG     - input dataset
%
% Outputs:
%   OUTEEG      output dataset
%
% Authors: Arnaud Delorme, SCCN, 2022

% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [OUTEEG, com] = pop_dipfit_headmodel( EEG, mriFile, varargin )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<1
  help pop_dipfit_headmodel;
  return
end

if nargin < 2
    str = strvcat('Fisrt, select the subject''s anatomical T1 MRI.', ...
                  'This tool only work if fiducials are defined is', ...
                  'an associated JSON file. See the DIPFIT tutorial', ...
                  'https://eeglab.org/tutorials/09_source/Custom_head_model.html');

    questdlg2(str, 'Subject''s MRI', 'Continue', 'Continue');
    [filename, filepath] = uigetfile('*', 'Select the MRI file');
    mriFile = fullfile(filepath, filename);

%     commandload = [ '[filename, filepath] = uigetfile(''*'', ''Select a text file'');' ...
%                     'if filename ~=0,' ...
%                     '   set(findobj(''parent'', gcbf, ''tag'', ''mrifile''), ''string'', [ filepath filename ]);' ...
%                     'end' ...
%                     'clear filename filepath tagtest;' ];
%     
%     geometry = { [2 1] [2 1] [0.8 0.3 1.5] [2.05 0.26 .75] [2.05 0.26 .75] [2.05 0.26 .75] ...
%                  [2.05 0.26 .75] [2.05 0.26 .75] [2.05 0.26 .75] [2.05 0.26 .75] [2 1] };
%     if isstruct(EEG.dipfit.mrifile)
%         mristring = 'set in DIPFIT settings';
%         mrienable = 'off';
%     else
%         mristring = EEG.dipfit.mrifile;
%         mrienable = 'on';
%     end
%     uilist = { { 'style' 'text' 'string' 'MRI file' } ...
%                { 'style' 'edit' 'string' mriFile } ...
%                { 'style' 'text' 'string' 'Select EEG (3 surfaces) or MEG (1 surface)' } ...
%                { 'style' 'edit' 'string' '' } ...
%                { 'style' 'text' 'string' 'Background image' } ...
%                { 'style' 'pushbutton' 'string' '...' 'enable' mrienable 'callback' commandload } ...
%                { 'style' 'edit' 'string' mristring 'enable' mrienable 'tag' 'mrifile' } ...
%                { 'style' 'text' 'string' 'Plot summary mode' } ...
%                { 'style' 'checkbox' 'string' '' } {} ...
%                { 'style' 'text' 'string' 'Plot edges' } ...
%                { 'style' 'checkbox' 'string' '' } {} ...
%                { 'style' 'text' 'string' 'Plot closest MRI slide' } ...
%                { 'style' 'checkbox' 'string' '' } {} ...
%                { 'style' 'text' 'string' 'Plot dipole''s 2-D projections' } ...
%                { 'style' 'checkbox' 'string' '' } {} ...
%                { 'style' 'text' 'string' 'Plot projection lines' } ...
%                { 'style' 'checkbox' 'string' '' 'value' 0 } {} ...
%                { 'style' 'text' 'string' 'Make all dipoles point out' } ...
%                { 'style' 'checkbox' 'string' '' } {} ...
%                { 'style' 'text' 'string' 'Normalized dipole length' } ...
%                { 'style' 'checkbox' 'string' '' 'value' 1 } {} ...
%                { 'style' 'text' 'string' 'Additionnal dipplot() options' } ...
%                { 'style' 'edit' 'string' '' } };
%      
% 	result = inputgui( geometry, uilist, 'pophelp(''pop_dipplot'')', 'Plot dipoles - pop_dipplot');
% 	if length(result) == 0 return; end

end

g.plotfiducial = [];
g.plotmesh = [];

mri = ft_read_mri(mriFile);
[mriPath,mriFile] = fileparts(mriFile);
[~,mriFile] = fileparts(mriFile);
fileSideCar = fullfile(mriPath, [ mriFile '.json' ]);

% read JSON file
if exist(fileSideCar, 'file')
    disp('JSON sidecar file found');
    fid = fopen(fileSideCar, 'r');
    if fid == -1
        error('Cannot open file %s', fileSideCar)
    end
    raw = fread(fid,inf);
    fclose(fid);
    coordinates = jsondecode(char(raw'));
    
    if isfield(coordinates, 'AnatomicalLandmarkCoordinates')
        if ~isfield(coordinates.AnatomicalLandmarkCoordinates, 'Nasion') || ...
                ~isfield(coordinates.AnatomicalLandmarkCoordinates, 'LPA') || ...
                ~isfield(coordinates.AnatomicalLandmarkCoordinates, 'RPA')
            error('Some anatomical coordinate missing')
        end

        cfg              = [];
        cfg.method       = 'fiducial';
        % this information has been obtained from the .json associated with the anatomical image
        cfg.fiducial.nas = coordinates.AnatomicalLandmarkCoordinates.Nasion';
        cfg.fiducial.lpa = coordinates.AnatomicalLandmarkCoordinates.LPA';
        cfg.fiducial.rpa = coordinates.AnatomicalLandmarkCoordinates.RPA';
        fiducials = [...
            cfg.fiducial.nas;
            cfg.fiducial.lpa;
            cfg.fiducial.rpa];

        % plot fiducial
        if ~isempty(g.plotfiducial)
            ind = strmatch(lower(g.plotfiducial), {'nasion' 'lpa' 'rpa'} );
            if isempty(ind)
                error('Fiducial not found')
            end
            cfg = [];
            cfg.locationcoordinates = 'voxel'; % treat the location as voxel coordinates
            cfg.location = fiducials(2,:);
            mri2 = mri;
            mri2.transform = eye(4);
            mri2.transform(:,4) = 1;
            ft_sourceplot(cfg, mri2);
        end

        mri = ft_volumerealign(cfg, mri);        

    else
        disp('JSON sidecar file found, but it does not contain the AnatomicalLandmarkCoordinates field');
    end
else
    error('Fiducials must be provided')
end

%mri  = ft_convert_coordsys(mri, 'acpc');
% extract volume
meshes        = {'brain','skull','scalp'};
cfg           = [];
cfg.output    = meshes;
segmentedmri  = ft_volumesegment(cfg, mri);

% convert fiducials
fiducials = segmentedmri.transform * [ fiducials ones(size(fiducials,1),1)]';
fiducials = fiducials(1:3,:)';
chanlocs = [];
chanlocs.labels = 'Nasion';
chanlocs.X      = fiducials(1,1);
chanlocs.Y      = fiducials(1,2);
chanlocs.Z      = fiducials(1,3);
chanlocs(end+1).labels = 'LPA';
chanlocs(end).X = fiducials(2,1);
chanlocs(end).Y = fiducials(2,2);
chanlocs(end).Z = fiducials(2,3);
chanlocs(end+1).labels = 'RPA';
chanlocs(end).X = fiducials(3,1);
chanlocs(end).Y = fiducials(3,2);
chanlocs(end).Z = fiducials(3,3);
chanlocs = convertlocs(chanlocs, 'cart2all');

if 0 % single shell MEG
    segmentedmri2 = segmentedmri;
    segmentedmri2 = rmfield(segmentedmri2, 'skull');
    segmentedmri2 = rmfield(segmentedmri2, 'scalp');
    cfg        = [];
    cfg.method = 'singleshell';
    headmodel  = ft_prepare_headmodel(cfg, segmentedmri2);
end

% plot mesh on MRI
if ~isempty(g.plotmesh)
    ind = strmatch(lower(g.plotmesh), meshes, 'exact' );
    if isempty(ind)
        error('Mesh not found')
    end
    cfg = [];
    cfg.intersectmesh = headmodel.bnd(ind);
    ft_sourceplot(cfg, mri);
end

if 0 % optional
    % prepare mesh
    cfg=[];
    cfg.tissue={'brain','skull','scalp'};
    cfg.numvertices = [3000 2000 1000];
    bnd=ft_prepare_mesh(cfg,segmentedmri);

    % create head model
    cfg        = [];
    cfg.method ='bemcp'; % You can also specify 'openmeeg', 'bemcp', or another method.
    headmodel  = ft_prepare_headmodel(cfg, bnd);
else
    % create head model
    cfg        = [];
    cfg.method ='bemcp'; % You can also specify 'openmeeg', 'bemcp', or another method.
    headmodel  = ft_prepare_headmodel(cfg, segmentedmri);
end

% save data in DIPFIT structure
EEG.dipfit.mrifile  = mri;
EEG.dipfit.hdmfile  = headmodel;
EEG.dipfit.chanfile = chanlocs;
EEG.dipfit.coordformat = 'MNI';


