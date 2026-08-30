%% plot_bdm_simulation
% Reproduce matlab/images/bdm_simulation.png from the coupled BDM output.
%
% The script plots radial water activity for selected aerosol bins and
% labels the figure according to the diffusion treatment in namelist.in.
% It plots bins
%   7, 13, 20, 27, 33, 40, 47, 53, 60.
%
% The script may be launched from any working directory.  By default it
% first runs main.exe using namelist.in, then reads the
% resulting NetCDF file and overwrites matlab/images/bdm_simulation.png.

clearvars;
close all;
clc;

%% User controls
run_model = false;             % false -> only replot an existing NetCDF file
bins_to_plot = [7 13 20; ...
                27 33 40; ...
                47 53 60];

%% Resolve paths from this script, not from MATLAB's current directory
script_path = mfilename('fullpath');
script_dir  = fileparts(script_path);
tmp_dir  = fileparts('/tmp/');
root_dir    = fileparts(script_dir);

main_namelist    = fullfile(root_dir,'namelist.in');
model_exe        = fullfile(root_dir,'main.exe');
ncfile           = fullfile(tmp_dir,'output3.nc');
image_dir        = fullfile(script_dir,'images');
pngfile          = fullfile(image_dir,'bdm_simulation.png');

if ~exist(image_dir,'dir')
    mkdir(image_dir);
end

%% Run the model using the current main namelist
if run_model
    if ~isfile(model_exe)
        error(['Cannot find compiled model:\n  %s\n\n' ...
               'Compile the model first, then rerun this script.'],model_exe);
    end
    if ~isfile(main_namelist)
        error('Cannot find main namelist: %s',main_namelist);
    end

    old_dir = pwd;
    restore_dir = onCleanup(@() cd(old_dir));
    cd(root_dir);

    cmd = sprintf('"%s" "%s"',model_exe,main_namelist);
    fprintf('Running: %s\n',cmd);
    [status,command_output] = system(cmd);
    fprintf('%s',command_output);
    if status ~= 0
        error('BDM run failed with exit status %d.',status);
    end
end

if ~isfile(ncfile)
    error(['Cannot find BDM NetCDF output:\n  %s\n' ...
           'Set run_model=true or point ncfile at an existing output file.'],ncfile);
end

%% Read output and put dimensions into a known MATLAB order
% c: kp x ncomp x nbins x times
% r: kp x nbins x times
c = read_netcdf_in_order(ncfile,'c',{'kp','ncomp','nbins','times'});
r = read_netcdf_in_order(ncfile,'r',{'kp','nbins','times'});
time = double(ncread(ncfile,'time'));
time = time(:);

nt = min([numel(time), size(c,4), size(r,3)]);
time = time(1:nt);
c = c(:,:,:,1:nt);
r = r(:,:,1:nt);

if size(c,2) < 2
    error('Variable c must contain at least water and solute components.');
end
if max(bins_to_plot,[],'all') > size(c,3)
    error('Requested bin %d, but output contains only %d bins.', ...
          max(bins_to_plot,[],'all'),size(c,3));
end

%% Water activity
% For the supplied BMM example there is one solute with van't Hoff factor
% nu_core1 = 3.  This is the same expression used by icenucleation_diff:
%
%     a_w = c_water / (c_water + nu * c_solute)
%
% Read nu from bmm/namelist.diffusion so the plotting example remains tied
% to the model setup rather than silently hard-coding it here.
bmm_namelist = fullfile(root_dir,'bmm','namelist.diffusion');
nu = read_namelist_number(bmm_namelist,'nu_core1',3.0);

cw = c(:,1,:,:);
cs = c(:,2,:,:);
denominator = cw + nu.*cs;
activity = cw ./ denominator;
activity(denominator <= 0) = NaN;

%% Diffusion treatment used by this simulation
% diffusion_type = 0: use the constant coefficient nmd%d_coeff from MBD.
% diffusion_type = 1: use the diffusion-coefficient parameterisation file.
diffusion_type = round(read_namelist_number(main_namelist,'diffusion_type',NaN));

if ~isfinite(diffusion_type)
    error('Could not read diffusion_type from %s.',main_namelist);
end

fprintf('diffusion_type = %d (from %s)\n',diffusion_type,main_namelist);

switch diffusion_type
    case 0
        diff_namelist = fullfile(root_dir,'mbd','namelist.in');
        D = read_namelist_number(diff_namelist,'nmd%d_coeff',NaN);
        if ~isfinite(D)
            error('Could not read nmd%%d_coeff from %s.',diff_namelist);
        end
        simulation_title = sprintf('simulation for constant D=%g m^2 s^{-1}',D);
        fprintf('Diffusion treatment: constant D = %g m^2 s^-1\n',D);

    case 1
        dc_filename = read_namelist_string(main_namelist,'diffusion_coeff_file', ...
                                           'namelist.diff_coeffs');
        if is_absolute_path(dc_filename)
            dc_namelist = dc_filename;
        else
            dc_namelist = fullfile(root_dir,dc_filename);
        end

        param = round(read_namelist_number(dc_namelist,'param',NaN));
        compound = round(read_namelist_number(dc_namelist,'compound',NaN));
        [param_name,compound_name] = diffusion_setup_name(param,compound);

        if isempty(compound_name)
            simulation_title = sprintf('simulation using %s',param_name);
        else
            simulation_title = sprintf('simulation using %s (%s)', ...
                                       param_name,compound_name);
        end
        fprintf('Diffusion treatment: %s; file = %s\n', ...
                simulation_title(18:end),dc_namelist);

    otherwise
        simulation_title = sprintf('simulation with diffusion type %d',diffusion_type);
        warning('Unknown diffusion_type=%d.',diffusion_type);
end

%% Plot
fig = figure('Color','w','Position',[50 50 1680 945]);
tlo = tiledlayout(fig,3,3,'TileSpacing','compact','Padding','compact');
colormap(fig,parula(256));

for ip = 1:numel(bins_to_plot)
    ibin = bins_to_plot(ip);
    ax = nexttile(tlo,ip);

    rr = squeeze(r(:,ibin,:));
    aa = squeeze(activity(:,:,ibin,:));
    aa = reshape(aa,size(rr));

    % Do not colour cells outside the particle.  The diffusion solver stores
    % zero concentrations there, so their activity is undefined.
    aa(~isfinite(aa) | rr <= 0) = NaN;

    tt = repmat(time.',size(rr,1),1);

    surf(ax,tt,rr,zeros(size(rr)),aa, ...
        'EdgeColor','none','FaceColor','flat');
    view(ax,2);
    set(ax,'YScale','log','Layer','top');
    box(ax,'on');

    xlabel(ax,'time (s)');
    ylabel(ax,'radius (m)');
    xlim(ax,[time(1) time(end)]);

    % Match the radial range of the moving-boundary grid used by the example.
    positive_r = rr(isfinite(rr) & rr > 0);
    if ~isempty(positive_r)
        ylim(ax,[max(1e-9,min(positive_r)) 5e-2]);
    else
        ylim(ax,[1e-9 5e-2]);
    end

    text(ax,0.10,0.80,sprintf('bin number %d',ibin), ...
        'Units','normalized','HorizontalAlignment','left', ...
        'VerticalAlignment','middle');

    cb = colorbar(ax);
    cb.Label.String = 'activity of water';

    % Let each panel use its own data range, as in the reference figure.
    aval = aa(isfinite(aa));
    if ~isempty(aval)
        amin = min(aval);
        amax = max(aval);
        if amax > amin
            clim(ax,[amin amax]);
        end
    end
end

sgtitle(tlo,simulation_title,'FontWeight','bold');

%% Save in the same repository location and with the same filename
try
    exportgraphics(fig,pngfile,'Resolution',150);
catch
    % Fallback for older MATLAB releases without exportgraphics.
    print(fig,pngfile,'-dpng','-r150');
end

fprintf('Wrote %s\n',pngfile);


%% Local helper functions
function A = read_netcdf_in_order(ncfile,varname,target_dims)
%READ_NETCDF_IN_ORDER Read a variable and permute it by NetCDF dimension name.
% This avoids depending on the Fortran/C/MATLAB dimension-order convention.

    info = ncinfo(ncfile,varname);
    source_dims = {info.Dimensions.Name};
    A = double(ncread(ncfile,varname));

    % Retain all dimensions even if MATLAB has collapsed a singleton one.
    source_size = info.Size;
    if numel(source_size) < numel(source_dims)
        source_size(end+1:numel(source_dims)) = 1;
    end
    if numel(A) == prod(source_size)
        A = reshape(A,source_size);
    end

    perm = zeros(1,numel(target_dims));
    for k = 1:numel(target_dims)
        idx = find(strcmp(source_dims,target_dims{k}),1);
        if isempty(idx)
            error(['Variable "%s" in %s does not have dimension "%s".\n' ...
                   'Found dimensions: %s'], ...
                   varname,ncfile,target_dims{k},strjoin(source_dims,', '));
        end
        perm(k) = idx;
    end

    if numel(source_dims) ~= numel(target_dims)
        error('Unexpected number of dimensions for variable "%s".',varname);
    end

    A = permute(A,perm);
end


function value = read_namelist_number(filename,key,default_value)
%READ_NAMELIST_NUMBER Extract an active numeric namelist assignment.
% Handles Fortran forms such as 1.e-17, 1.0e-17, .5e-3 and 1.d-17.
% Text after ! is stripped before matching, so commented examples cannot
% override the live namelist setting.

    value = default_value;
    if ~isfile(filename)
        warning('Could not find %s; using %g for %s.',filename,default_value,key);
        return;
    end

    txt = strip_fortran_comments(fileread(filename));
    escaped_key = regexptranslate('escape',key);
    number_pattern = '[-+]?(?:(?:\d+(?:\.\d*)?)|(?:\.\d+))(?:[eEdD][-+]?\d+)?';
    expression = ['(?m)^\s*' escaped_key ...
                  '\s*(?:\([^\)]*\))?\s*=\s*(' number_pattern ')'];
    tok = regexp(txt,expression,'tokens','once');
    if isempty(tok)
        warning('Could not read %s from %s; using %g.',key,filename,default_value);
        return;
    end

    number_string = regexprep(tok{1},'[dD]','e');
    parsed = str2double(number_string);
    if isfinite(parsed)
        value = parsed;
    else
        warning('Could not parse %s from %s; using %g.',key,filename,default_value);
    end
end


function value = read_namelist_string(filename,key,default_value)
%READ_NAMELIST_STRING Extract an active quoted string namelist assignment.

    value = default_value;
    if ~isfile(filename)
        warning('Could not find %s; using "%s" for %s.', ...
                filename,default_value,key);
        return;
    end

    txt = strip_fortran_comments(fileread(filename));
    escaped_key = regexptranslate('escape',key);
    expression = ['(?m)^\s*' escaped_key '\s*=\s*[''\"]([^''\"]+)[''\"]'];
    tok = regexp(txt,expression,'tokens','once');
    if ~isempty(tok)
        value = strtrim(tok{1});
    else
        warning('Could not read %s from %s; using "%s".', ...
                key,filename,default_value);
    end
end


function txt = strip_fortran_comments(txt)
%STRIP_FORTRAN_COMMENTS Remove ! comments line-by-line.
% This is sufficient for these namelists, where ! is not used inside paths.

    lines = regexp(txt,'\r\n|\n|\r','split');
    for i = 1:numel(lines)
        bang = strfind(lines{i},'!');
        if ~isempty(bang)
            lines{i} = lines{i}(1:bang(1)-1);
        end
    end
    txt = strjoin(lines,newline);
end


function tf = is_absolute_path(filename)
%IS_ABSOLUTE_PATH True for Unix or Windows absolute paths.
    tf = startsWith(filename,filesep) || ...
         ~isempty(regexp(filename,'^[A-Za-z]:[\\/]','once'));
end


function [param_name,compound_name] = diffusion_setup_name(param,compound)
%DIFFUSION_SETUP_NAME Human-readable names matching namelist.diff_coeffs.

    names = {'constant','Darken','Vignes','Lienhard2014','Lienhard2015', ...
             'Price2014','Price2015','Price2016','Zobrist2011','Shiraiwa2013'};

    if isfinite(param) && param >= 1 && param <= numel(names)
        param_name = names{param};
    else
        param_name = sprintf('parameterisation %g',param);
    end

    compound_name = '';
    if ~isfinite(compound)
        return;
    end

    switch param
        case 4
            compounds = {'citric acid'};
        case 5
            compounds = {'levoglucosan','levoglucosan/NH4HSO4','raffinose', ...
                         '3-MBTCA','alpha-pinene','sucrose','citric acid', ...
                         'shikimic acid'};
        case 6
            compounds = {'sucrose','levoglucosan','MgSO4','raffinose'};
        case 7
            compounds = {'alpha-pinene'};
        case 8
            compounds = {'sucrose','water'};
        case 9
            compounds = {'sucrose'};
        case 10
            compounds = {'alpha-pinene'};
        otherwise
            compounds = {};
    end

    if compound >= 1 && compound <= numel(compounds)
        compound_name = compounds{compound};
    elseif ~isempty(compounds)
        compound_name = sprintf('compound %d',compound);
    end
end
