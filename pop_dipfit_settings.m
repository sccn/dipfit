% pop_dipfit_settings() - select global settings for dipole fitting through a pop up window
%
% Usage:
%   >> OUTEEG = pop_dipfit_settings ( INEEG ); % pop up window
%   >> OUTEEG = pop_dipfit_settings ( INEEG, 'key1', 'val1', 'key2', 'val2' ... )
%
% Inputs:
%   INEEG	input dataset
%
% Optional inputs:
%   'hdmfile'  - [string] file containing a head model compatible with
%                the Fieldtrip dipolefitting() function ("vol" entry)
%   'mrifile'  - [string] file containing an anatomical MR head image. 
%                The MRI must be normalized to the MNI brain. See the .mat 
%                files used by the sphere and boundary element models
%                (For instance, select the sphere model and study 'EEG.dipfit'). 
%                If SPM2 software is installed, dipfit will be able to read 
%                most MRI file formats for plotting purposes (.mnc files, etc...). 
%                To plot dipoles in a subject MRI, first normalize the MRI 
%                to the MNI brain using SPM2.
%   'coordformat' - ['MNI'|'Spherical'] Coordinates returned by the selected
%                head model. May be MNI coordinates or spherical coordinates.
%                (For spherical coordinates, the head radius is assumed to be 85 mm.
%   'chanfile' - [string] template channel locations file. (This function will
%                check whether your channel locations file is compatible with 
%                your selected head model).
%   'chansel'  - [integer vector] indices of channels to use for dipole fitting. 
%                {default: all}
%   'electrodes' - [integer array] indices of channels to include
%                  in the dipole model. {default: all}
%   'model'      - ['standardBEM'|'standardBESA'] use one of the standard
%                  models. {default: all}
%   'coord_transform' - [float array] Talairach transformation matrix for
%                       aligning the dataset channel locations to the selected 
%                       head model. You may also use the string
%                       'warpfiducials' to 
%
% Outputs:
%   OUTEEG	output dataset
%
% Author: Arnaud Delorme, SCCN, La Jolla 2003-
%         Robert Oostenveld, SMI/FCDC, Nijmegen 2003

% MEG flag:
%   'gradfile' - [string] file containing gradiometer locations
%                ("gradfile" parameter in Fieldtrip dipolefitting() function)

% SMI, University Aalborg, Denmark http://www.smi.auc.dk/
% FC Donders Centre, University Nijmegen, the Netherlands http://www.fcdonders.kun.nl

% Copyright (C) 2003 arno@salk.edu, Arnaud Delorme, SCCN, La Jolla 2003-2005
%
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
function [OUTEEG, com] = pop_dipfit_settings ( EEG, varargin )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1
   help pop_dipfit_settings;
   return;
end

if ~plugin_askinstall('Fieldtrip-lite', 'ft_dipolefitting'), return; end

OUTEEG = EEG;
com = '';

% get the default values and filenames
dipfitdefs;

if nargin < 2

    if ischar(EEG)
        userdat = get(gcf, 'userdata');
        current   = userdat.current_model;
        templates = userdat.template_models;
        chanfile  = userdat.chanfile;
        valmodel  = get(findobj(gcf, 'tag', 'listmodels'), 'value');
        
        if isequal(EEG, 'setmodel') % redraw
            [hdmEnable,  hdmString ] = checkStructure(current.hdmfile);
            [mriEnable,  mriString ] = checkStructure(current.mrifile);
            [chanEnable, chanString] = checkStructure(current.chanfile);
            coordformat = current.coordformat;
            set(findobj(gcf, 'tag', 'model'  ), 'string', hdmString, 'enable', hdmEnable);
            set(findobj(gcf, 'tag', 'mri'    ), 'string', mriString, 'enable', mriEnable);
            set(findobj(gcf, 'tag', 'chans'  ), 'string', chanString, 'enable', chanEnable);
            set(findobj(gcf, 'tag', 'coord'), 'value' , fastif(strcmpi(current.coordformat,'MNI'),2, ...
                fastif(strcmpi(current.coordformat,'Spherical'),1,3)));
            set(findobj(gcf, 'tag', 'modelBrowse'  ), 'enable', hdmEnable);
            set(findobj(gcf, 'tag', 'mriBrowse'    ), 'enable', mriEnable);
            set(findobj(gcf, 'tag', 'chansBrowse'  ), 'enable', chanEnable);

            [allkeywordstrue, transform] = lookupchantemplate(chanfile, templates(valmodel).coord_transform);

            if allkeywordstrue
                set(findobj(gcf, 'tag', 'coregtext'), 'string', char(vararg2str({ transform })));
                if isempty(transform)
                    set(findobj(gcf, 'tag', 'coregcheckbox'), 'value', 1);
                else set(findobj(gcf, 'tag', 'coregcheckbox'), 'value', 0);
                end
            else
                set(findobj(gcf, 'tag', 'coregtext'), 'string', '');
            end

        elseif isequal(EEG, 'changemodel') % redraw
            userdat.current_model = userdat.template_models(valmodel);
            set(gcf, 'userdata', userdat);
            pop_dipfit_settings('setmodel');
            
        elseif isequal(EEG, 'coregister') % redraw
            coregObj = findobj(gcbf, 'tag', 'coregtext');
            tmptransf = get( coregObj, 'string');
            [~, tmptransf] = coregister(userdat.chanlocs, current.chanfile, 'mesh', current.hdmfile, ...
                    'transform', str2num(tmptransf), 'chaninfo1', userdat.chaninfo, 'helpmsg', 'on');
            if ~isempty(tmptransf)
                set( coregObj, 'string', num2str(tmptransf)); 
            end
            
        end
            return;
    end
    
    % detect DIPFIT1.0x structure
    % ---------------------------
    if isfield(EEG(1).dipfit, 'vol')
        str = [ 'Dipole information structure from DIPFIT v1.02 detected.' ...
                'Keep or erase the old dipole information including dipole locations? ' ...
                'In either case, a new dipole model can be constructed.' ];
        
        tmpButtonName=questdlg2( strmultiline(str, 60), 'Old DIPFIT structure', 'Keep', 'Erase', 'Keep');
        if strcmpi(tmpButtonName, 'Keep'), return; end    

    end  
    
    % determine default model
    % -----------------------
    dipfitdefs; % contains template_model
    valmodel    = 2; % default model now MNI
    userdata    = [];
    if isfield(EEG(1).dipfit, 'hdmfile') && ~isempty(EEG(1).dipfit.hdmfile)
        OUTEEG.dipfit = EEG(1).dipfit;
        if isequal(EEG(1).dipfit.hdmfile, template_models(1).hdmfile)
            valmodel = 1;
        elseif isequal(EEG(1).dipfit.hdmfile, template_models(2).hdmfile)
            valmodel = 2;
        else
            valmodel = 3;
        end
        nocoreg = false;
    else
        if isfield(EEG(1).chaninfo, 'filename')
            if ~isempty(findstr(lower(EEG(1).chaninfo.filename), 'standard-10-5-cap385')), valmodel = 1; end
            if ~isempty(findstr(lower(EEG(1).chaninfo.filename), 'standard_1005')),        valmodel = 2; end
        end
        OUTEEG.dipfit.hdmfile     = template_models(valmodel).hdmfile;
        OUTEEG.dipfit.mrifile     = template_models(valmodel).mrifile;
        OUTEEG.dipfit.chanfile    = template_models(valmodel).chanfile;
        OUTEEG.dipfit.coordformat = template_models(valmodel).coordformat;
        nocoreg = ~template_models(valmodel).coregval;
    end
   
    % coordinate format
    if isequal(OUTEEG.dipfit.coordformat, 'spherical')
        valcoord = 1;
    elseif isequal(OUTEEG.dipfit.coordformat, 'MNI')
        valcoord = 2;
    else
        valcoord = 3;
    end
    geomvert = [1 1 1 1 1 1 1 1];
    geomhorz = {
        [1 2] 
        [0.2 1 1.3 0.5 0.5 ]
        [0.2 1 1.3 0.9 0.1 ]
        [0.2 1 1.3 0.5 0.5 ]
        [0.2 1 1.3 0.5 0.5 ]
        [1]
        [1 1 0.5 0.5 ]
        [1 1 0.5 0.5 ]
        };
    
    % define each individual graphical user element
    comhelp1 = [ 'warndlg2(strvcat(''The two default head models are in standard_BEM and standard_BESA'',' ...
                 ''' sub-folders in the DIPFIT2 plugin folder, and may be modified there.''), ''Model type'');' ];
    comhelp3 = [ 'warndlg2(strvcat(''Any MR image normalized to the MNI brain model may be used for plotting'',' ...
                 '''(see the DIPFIT 2.0 tutorial for more information)''), ''Model type'');' ];
    comhelp2 = [ 'warndlg2(strvcat(''The template location file associated with the head model'',' ...
                 '''you are using must be entered (see tutorial).''), ''Template location file'');' ];
    commandload1 = [ '[filename, filepath] = uigetfile(''*'', ''Select a text file'');' ...
                    'if filename ~=0,' ...
                    '   set(findobj(''parent'', gcbf, ''tag'', ''model''), ''string'', [ filepath filename ]);' ...
                    'end;' ...
                    'clear filename filepath tagtest;' ];    
    commandload2 = [ '[filename, filepath] = uigetfile(''*'', ''Select a text file'');' ...
                    'if filename ~=0,' ...
                    '   set(findobj(''parent'', gcbf, ''tag'', ''meg''), ''string'', [ filepath filename ]);' ...
                    'end;' ...
                    'clear filename filepath tagtest;' ];    
    commandload3 = [ '[filename, filepath] = uigetfile(''*'', ''Select a text file'');' ...
                    'if filename ~=0,' ...
                    '   set(findobj(''parent'', gcbf, ''tag'', ''mri''), ''string'', [ filepath filename ]);' ...
                    'end;' ...
                    'clear filename filepath tagtest;' ];    
    cb_selectcoreg = 'pop_dipfit_settings(''coregister'');';
    cb_selectelectrodes = 'tmplocs = EEG(1).chanlocs; set(findobj(gcbf, ''tag'', ''elec''), ''string'', int2str(pop_chansel({tmplocs.labels}))); clear tmplocs;';
    changemodel = 'pop_dipfit_settings(''changemodel'');';
    
    templatenames = { template_models.name };
    elements  = { ...
        { 'style' 'text'        'string'  'Select a head model' 'fontweight' 'bold' } ...
        { 'style' 'popupmenu'     'string'  strvcat(templatenames{:}) ... 
                                'callback' changemodel 'value' valmodel 'tag' 'listmodels' } ...
        { } { 'style' 'text'        'string' 'Output coordinates' } ...
        { 'style' 'popupmenu'   'string' 'spherical (head radius 85 mm)|MNI|Custom' 'tag' 'coord' ...
          'value' valcoord  'userdata' 'editable' 'enable' 'off'} { } { } ...
        { } { 'style' 'text'        'string' 'Head model file' 'tag' 'headstr' } ...
        { 'style' 'edit'        'string' '_______' 'tag' 'model' 'horizontalalignment' 'right' 'userdata' 'editable' 'enable' 'on'} ...
        { 'style' 'pushbutton'  'string' 'Browse'  'tag' 'modelBrowse'   'callback' commandload1       'userdata' 'editable' 'enable' 'off' } ...
        { 'style' 'pushbutton'  'string' 'Help'      'callback' comhelp1 } ...
        { } { 'style' 'text'        'string' 'Associated MRI file for plotting' 'tag' 'mristr' } ...
        { 'style' 'edit'        'string' '_______' 'tag' 'mri' 'horizontalalignment' 'right' 'enable' 'on' } ...
        { 'style' 'pushbutton'  'string' 'Browse'  'tag' 'mriBrowse'     'callback' commandload3 } ...
        { 'style' 'pushbutton'  'string' 'Help'         'callback' comhelp3 } ...
        { } { 'style' 'text'        'string' 'Associated coord. landmarks if any', 'tag', 'chanstr' } ...
        { 'style' 'edit'        'string' '_______' 'tag' 'chans'  'horizontalalignment' 'right' 'userdata' 'editable' 'enable' 'on'} ...
        { 'style' 'pushbutton'  'string' 'Browse'  'tag' 'chansBrowse' 'callback' commandload2  'userdata' 'editable' 'enable' 'on'} ...
        { 'style' 'pushbutton'  'string' 'Help'         'callback' comhelp2 } ...
        { } ...
        { 'style' 'text'        'string' 'Matrix to align chan. locs. with head model' 'userdata' 'coreg' } ...
        { 'style' 'edit'        'string' '' 'tag' 'coregtext' 'userdata' 'coreg' } ...
        { 'style' 'pushbutton'  'string' 'Co-register' 'fontweight' 'bold' 'tag' 'manualcoreg' 'callback' cb_selectcoreg } ... 
        { 'style' 'checkbox'    'string' 'No Co-Reg.'    'tag' 'coregcheckbox' 'value' nocoreg  'userdata' 'coreg' } ... 
        { 'style' 'text'        'string' 'Channels to omit from dipole fitting' } ...
        { 'style' 'edit'        'string' ''             'tag' 'elec' } ...
        { 'style' 'pushbutton'  'string' '...' 'callback' cb_selectelectrodes } { } ...
                };
    
    % plot GUI and protect parameters
    % -------------------------------
    userdata.template_models  = template_models;
    userdata.current_model    = OUTEEG.dipfit;
    userdata.chanlocs         = EEG(1).chanlocs;
    userdata.chaninfo         = EEG(1).chaninfo;
    if isfield(EEG(1).chaninfo, 'filename')
         userdata.chanfile         = lower(EEG(1).chaninfo.filename);
    else userdata.chanfile         = '';
    end
    optiongui = { 'geometry', geomhorz, 'uilist', elements, 'helpcom', 'pophelp(''pop_dipfit_settings'')', ...
                  'title', 'Dipole fit settings - pop_dipfit_settings()', ...
                  'userdata', userdata, 'geomvert', geomvert 'eval' 'pop_dipfit_settings(''setmodel'');' };
	[result, ~, ~, outstruct] = inputgui( optiongui{:});
    if isempty(result), return; end
    
    if test_wrong_parameters(outstruct)
    	return;
    end

    % decode GUI inputs
    % -----------------
    options = {};
    if ~strcmpi(result{3}, 'Matlab structure')
        options = { options{:} 'hdmfile' result{3} };
    end
    if ~strcmpi(result{4}, 'Matlab structure')
        options = { options{:} 'mrifile' result{4} };
    end
    if ~strcmpi(result{5}, 'Matlab structure')
        options = { options{:} 'chanfile' result{5} };
    end
    if valcoord == 1
        options = { options{:} 'coordformat' 'Spherical' };
    elseif valcoord == 2
        options = { options{:} 'coordformat' 'MNI' };
    else
        if isfield(EEG.dipfit, 'coordformat') && ~isempty(EEG.dipfit.coordformat)
            options = { options{:} 'coordformat' EEG.dipfit.coordformat };
        else
            options = { options{:} 'coordformat' 'other' };
        end
    end
    if ~result{7}, options = { options{:} 'coord_transform' str2num(result{6}) }; end
    options = { options{:} 'chansel' setdiff(1:EEG(1).nbchan, str2num(result{8})) };

else
    options = varargin;
end

% process multiple datasets
% -------------------------
if length(EEG) > 1
    if nargin < 2
        [ OUTEEG, com ] = eeg_eval( 'pop_dipfit_settings', EEG, 'warning', 'on', 'params', options );
    else
        [ OUTEEG, com ] = eeg_eval( 'pop_dipfit_settings', EEG, 'params', options );
    end
    return;
end

g = finputcheck(options, { 'hdmfile'  'string'    []         '';
                                 'mrifile'  ''    []         '';
                                 'chanfile' ''    []         '';
                                 'chanfile' ''    []         '';
                                 'chansel'  'integer'   []         1:EEG.nbchan;
                                 'electrodes' 'integer'   []         [];
                                 'coord_transform' ''     []         [];
                                 'plotalignment' 'string'     {'on' 'off'}         'off';
                                 'model'      'string'    {'standardBEM' 'standardBESA' ''}         '';
                                 'coordformat' 'string'    { 'MNI','spherical','CTF' } 'MNI' });
if isstr(g), error(g); end

dipfitdefs;
OUTEEG.dipfit = EEG.dipfit;
if ~isempty(g.hdmfile)
    OUTEEG.dipfit.hdmfile = g.hdmfile;
end
if ~isempty(g.mrifile)
    OUTEEG.dipfit.mrifile = g.mrifile;
end
if ~isempty(g.chanfile)
    OUTEEG.dipfit.chanfile = g.chanfile;
end
if ~isempty(g.coordformat)
    OUTEEG.dipfit.coordformat = g.coordformat;
end
if strcmpi(g.model, 'standardBESA')
    OUTEEG.dipfit.hdmfile     = template_models(1).hdmfile;
    OUTEEG.dipfit.mrifile     = template_models(1).mrifile;
    OUTEEG.dipfit.chanfile    = template_models(1).chanfile;
    OUTEEG.dipfit.coordformat = template_models(1).coordformat;
elseif strcmpi(g.model, 'standardBEM')
    OUTEEG.dipfit.hdmfile     = template_models(2).hdmfile;
    OUTEEG.dipfit.mrifile     = template_models(2).mrifile;
    OUTEEG.dipfit.chanfile    = template_models(2).chanfile;
    OUTEEG.dipfit.coordformat = template_models(2).coordformat;
end
OUTEEG.dipfit.chansel     = g.chansel;
if ischar(g.coord_transform)
    if isequal(g.coord_transform, 'warpfiducials')
        [~,coord_transform] = coregister(EEG.chaninfo.nodatchans, OUTEEG.dipfit.chanfile, 'warp', 'auto', 'manual', 'off');
        OUTEEG.dipfit.coord_transform = coord_transform;
    end
else
    OUTEEG.dipfit.coord_transform = g.coord_transform;
end
if ~isempty(g.electrodes), OUTEEG.dipfit.chansel = g.electrodes; end

% removing channels with no coordinates
% -------------------------------------
[~, ~, ~, ~, indices] = readlocs(EEG.chanlocs);
if length(indices) < length(EEG.chanlocs)
    disp('Warning: Channels without coordinates removed from dipole fitting!');
    OUTEEG.dipfit.chansel = intersect( OUTEEG.dipfit.chansel, indices);
end

if isfield(OUTEEG.chanlocs, 'type') && ~isempty(strfind(OUTEEG.chanlocs(1).type, 'meg'))

    % Standard head model
    comp = eeglab2fieldtrip(OUTEEG, 'componentanalysis', 'dipfit');
    if ischar(OUTEEG.dipfit.hdmfile) && ~isempty(strfind(OUTEEG.dipfit.hdmfile, 'dipfit')) && ~isempty(strfind(OUTEEG.dipfit.hdmfile, 'standard_vol'))
        disp('Coordinate transform for MEG detected - no MRI provided, always try to use the subject''s MRI')

        if ~isempty(OUTEEG.dipfit.coord_transform)
            % Rational: we change the MRI transformation matrix to match
            % the final sensor space based on the provided transformation
            % matrix (usually computed using fiducials)

            transform_mat = OUTEEG.dipfit.coord_transform;
            if isfield(OUTEEG.chaninfo, 'originalnosedir')
                if strcmpi(OUTEEG.chaninfo.originalnosedir, '+Y')
                    transform_mat(6) = transform_mat(6)+pi/2;
                end
            end
            tra = traditionaldipfit( transform_mat );
            tra = pinv(tra);

            % change MRI coordinate system
            mri = ft_read_mri(OUTEEG.dipfit.mrifile);
            mri.transform = tra * mri.transform;
            mri.coordsys  = comp.grad.coordsys; % target coordinate system
            mri.unit      = comp.grad.unit;
            OUTEEG.dipfit.mrifile = mri;

            % change head model coordinate system
            hdm = ft_read_headmodel(OUTEEG.dipfit.hdmfile);
            hdm.bnd(3).pos = tra * [ hdm.bnd(3).pos ones(length(hdm.bnd(3).pos),1)]';
            hdm.bnd(3).pos = hdm.bnd(3).pos(1:3,:)';
            hdm.unit = comp.grad.unit;
            OUTEEG.dipfit.hdmfile = hdm;
        else
            % Rational: we perform automated coordinate transformation
            % based on the MRI and MEG coordinate space coordinates

            mri = ft_read_mri(OUTEEG.dipfit.mrifile);
            mri.coordsys = 'acpc';
            OUTEEG.dipfit.mrifile  = ft_convert_coordsys(mri, comp.grad.coordsys); % better to use the fiducials above

            hdm = ft_read_headmodel(OUTEEG.dipfit.hdmfile);
            hdm.coordsys = 'acpc';
            OUTEEG.dipfit.hdmfile  = ft_convert_coordsys(hdm, comp.grad.coordsys); % better to use the fiducials above

        end
        cfg = [];
        cfg.method='singleshell';
        OUTEEG.dipfit.hdmfile = ft_prepare_headmodel(cfg, OUTEEG.dipfit.hdmfile.bnd(3)); % use brain surface only
        OUTEEG.dipfit.coord_transform = []; % all contained in the MRI transformation and model coordinate

    else
        if ~isempty(OUTEEG.dipfit.coord_transform)
            error([ 'Coordinate transform cannot be used with custom head model.' 10 'The model is assumed to be aligned with the sensors.' ])
        end
        OUTEEG.dipfit.mrifile = ft_read_mri(OUTEEG.dipfit.mrifile);
        OUTEEG.dipfit.hdmfile = ft_read_headmodel(OUTEEG.dipfit.hdmfile);
    end

    % ft_convert_coordsys assumes the same units
    if isequal(comp.grad.unit, 'cm') && isequal(OUTEEG.dipfit.mrifile.unit, 'mm')
        OUTEEG.dipfit.mrifile.transform = traditionaldipfit([0 0 0 0 0 0 1/10 1/10 1/10])*OUTEEG.dipfit.mrifile.transform;
        OUTEEG.dipfit.mrifile.unit      = 'cm';
        OUTEEG.dipfit.hdmfile.bnd.pos   = OUTEEG.dipfit.hdmfile.bnd.pos/10;
        OUTEEG.dipfit.hdmfile.unit      = 'cm';
    end

    % check alignment
    if strcmpi(g.plotalignment, 'on')
        if ischar(OUTEEG.dipfit.hdmfile)
            hdm = ft_read_headmodel(OUTEEG.dipfit.hdmfile);
        else
            hdm = OUTEEG.dipfit.hdmfile;
        end
        if isequal(hdm.unit, 'mm') && isequal(comp.grad.unit, 'cm')
            hdm.bnd.pos = hdm.bnd.pos/10;
        end
        figure; ft_plot_mesh(comp.grad.chanpos(1:2:end,:), 'vertexindex', true);
        ft_plot_mesh(hdm.bnd(1), 'facealpha', 0.5)
    end
else
    % check alignment
    if strcmpi(g.plotalignment, 'on')
        if ischar(OUTEEG.dipfit.hdmfile)
            hdm = ft_read_headmodel(OUTEEG.dipfit.hdmfile);
        else
            hdm = OUTEEG.dipfit.hdmfile;
        end
        figure; 
        ft_plot_mesh(comp.elec.chanpos(1:2:end,:), 'vertexindex', true);
        ft_plot_mesh(hdm.bnd(1), 'facealpha', 0.5)
    end
end

% checking electrode configuration
% --------------------------------
if 0
    disp('Checking the electrode configuration');
    tmpchan          = readlocs(OUTEEG.dipfit.chanfile);
    [tmp1 ind1 ind2] = intersect( lower({ tmpchan.labels }), lower({ OUTEEG.chanlocs.labels }));
    if isempty(tmp1)
        disp('No channel labels in common found between template and dataset channels');
        if ~isempty(findstr(OUTEEG.dipfit.hdmfile, 'BESA'))
            disp('Use the channel editor to fit a head sphere to your channel locations.');
            disp('Check for inconsistency in dipole info.');
        else
            disp('Results using standard BEM model are INACCURATE when the chan locations are not on the head surface!');
        end
    else % common channels: performing best transformation
        TMP = OUTEEG;
        elec1 = eeglab2fieldtrip(TMP, 'elec');
        elec1 = elec1.elec;
        TMP.chanlocs = tmpchan;
        elec2 = eeglab2fieldtrip(TMP, 'elec');
        elec2 = elec2.elec;
        cfg.elec     = elec1;
        cfg.template = elec2;
        cfg.method   = 'warp';
        elec3 = electrodenormalize(cfg);
    
        % convert back to EEGLAB format
        OUTEEG.chanlocs = struct( 'labels', elec3.label, ...
                              'X'     , mat2cell(elec3.elecpos(:,1)'), ...
                              'Y'     , mat2cell(elec3.elecpos(:,2)'), ...
                              'Z'     , mat2cell(elec3.elecpos(:,3)') );
        OUTEEG.chanlocs = convertlocs(OUTEEG.chanlocs, 'cart2all');
    end
end

com = sprintf('EEG = pop_dipfit_settings( EEG, %s);', vararg2str(options));

% test for structure
% ------------------
function [enableFlag, str] = checkStructure(file)
    if isstruct(file)
        str = 'Matlab structure';
        enableFlag = 'off';
    else
        str = file;
        enableFlag = 'on';
    end

% test for wrong parameters
% -------------------------
function bool = test_wrong_parameters(result)

    coreg1 = result.coregtext;
    coreg2 = result.coregcheckbox; 
    meg    = result.coord; 
    
    bool = 0;
    if meg == 3, return; end
    if coreg2 == 0 && isempty(coreg1)
         bool = 1; warndlg2(strvcat('INCORRECT SETTINGS - SELECT THIS MENU AGAIN', ...
                                    'You must co-register your channel locations', ...
                                    'with the head model (Press button, "Co-register".', ...
                                    'and follow instructions); To bypass co-registration,', ...
                                    'check the checkbox " No Co-Reg".'), 'Error');
    end
