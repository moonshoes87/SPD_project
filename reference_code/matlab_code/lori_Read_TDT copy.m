clear vars *;
tic


% Method codes
% 0 - Undefined
% 1 - Coe
% 2 - IZZI
% 3 - Aitken
% 4 - Thellier


DataDIR='./Data/';

%folders=dir(strcat(DataDIR, 'Etna*'));
%files1 = dir(strcat(DataDIR, folders(1).name, '/*tdt'));

files=dir(strcat(DataDIR, '*.tdt'));
n_files=length(files);


% Load the AARM tensor data for applying the anistropy correction
Tensors=tdfread(strcat(DataDIR, 'AARM_Tensors.dat') );

% Load the NLT factor data for applying the non-linear TRM correction
NLT=tdfread(strcat(DataDIR, 'NLT_factors.dat') );

Specimens = struct();

%% The main loop

for f=1:n_files
    
    %% Read the TDT files
    flags=struct;
    Data=[];
    
    % Open and read the file
    FID=fopen(strcat(DataDIR, files(f).name),'r');
    files(f).name
    H2=textscan(FID, '%f', 'HeaderLines', 1); % Read the second headerline, skipping the first header
    Flab=H2{1,1}(1); % Lab field
    Az=H2{1,1}(2); % Core azimuth
    Pl=-H2{1,1}(3); % Core plunge
    
    [input]=textscan(FID, '%s %f %f %f %f'); % Read the data into an array
    
    temperature=num2str(input{1,2}); % second column is temperature
    Int=input{1,3}; % Intensity
    Dec=input{1,4}; % Declination
    Inc=input{1,5}; % Inclination
    
    fclose(FID); % close the input file
    
    % Organize measurement types and temperature steps
    % Note - The files i work with use the 150.2 convention, not the 150.12 style that is sometimes used
    Temps=floor(input{1,2}); % get the temperatures without the treatments
    Treatment=NaN(length(Int), 1);
    
    for i=1:length(Int)
        splt=regexp(temperature(i,:), '\.', 'split'); %split the temperature at the decimal point, i.e., 150.2 becomes 150 and 2
        if length(splt) < 2
            % Catch cases where no 0.0 is added, i.e., 150 instead of 150.0
            Treatment(i)=0;
        else
            Treatment(i)=str2double(splt{1,2});
        end
    end
    
    % Get the number of steps
    UT=unique(Temps); % Unique temperatures
    points=length(UT); % No of unique heating steps
    
    
    % Convert to x, y, z and put into Mvec
    Mvec=NaN(length(Int), 3); %Create an empty matrix
    [Mvec(:,1), Mvec(:,2), Mvec(:,3)]=dir2cart(Dec, Inc, Int);
    
    % a is the same as spec.pars['zijdblock'][0]
    a=[Temps(1), Dec(1), Inc(1), Int(1)];
    b=[Temps(2), Dec(2), Inc(2), Int(2)];
    
   
   
    
    
    %% Setup the flags for the expected fields, anisotropy and NLT, methods, etc
    % This sets up the various flags and meta data needed. This data is
    % needed for SPD, but is not included in the TDT files
    % Intensity method
    % Blab orientation
    % Expected B
    % Rotaion information - may not be entirely correct.
    
    % Experimental flag
    Exp_Flag=1; % default is the Coe/Aitken/IZZI protocols
    Meth=0; % Catch any unassigned samples
    
    % Lab Field orientation - Default assumes negative z-axis
    F_orient=[0,0,-1];
    % lj for debugging
    F_orient = [0,0,1];
    %lj
    
    
    flags.Anis=0;
    s_tensor=[];
    
    flags.NLT=0;
    NLT_hat=[];
    
    % NRM rotation flags, default is no rotation with NRM along negative z
    flags.Rot=0;
    %     Az=[];
    %     Pl=[];
    N_orient=[];
    
    
    if (strcmp(files(f).name(1:2), 'LV')) % Paterson Lascar
        % -Z-axis
        F_exp=24.0;
        Meth=1;
        flags.Rot=1;
        [N_orient(1), N_orient(2), N_orient(3)]=dir2cart(385.5, -18.7, 1);
    elseif (strcmp(files(f).name(1:2), 'vm')) % Muxworthy Vesuvius
        % Z-axis
        F_orient=[0,0,1];
        F_exp=44;
        Meth=1;
    elseif (strcmp(files(f).name(1:1), 'p')) % Muxworthy
        % Z-axis
        F_exp=45;
        F_orient=[0,0,1];
        Meth=1;
    elseif (strcmp(files(f).name(1:4), 'Mux_') )
        % Z-axis
        F_exp=100;
        F_orient=[0,0,1];
        flags.Rot=1;
        N_orient=F_orient;
    elseif (strcmp(files(f).name(1:4), 'Kras') )
        flags.Rot=1;
        Meth=1;
        if (strcmp(files(f).name(1:7), 'Krasa_M') )
            % -Z-axis
            F_exp=25;
            F_orient=[0,0,1];
            N_orient=[1/3, -2/3, -2/3];
        elseif (strcmp(files(f).name(1:7), 'Krasa_W') )
            % -Z-axis
            F_exp=60;
            N_orient=F_orient;
        else
            disp(files(f).name(1:end))
            error('IntAnalysis:F_exp', 'Unrecogonised sample name from Krasa et al., 2003.');
        end
        
    elseif (strcmp(files(f).name(1:1), 'M')) % Paterson MSH
        % -Z-axis
        F_exp=55.6;
        flags.Rot=1;
        [N_orient(1), N_orient(2), N_orient(3)]=dir2cart(17.1, 67.5, 1);
        Meth=0; % This is a mix of Coe and IZZI
    elseif (strcmp(files(f).name(1:2), 'TS'))
        % X-axis
        F_exp=45.7;
        F_orient=[1,0,0];
        flags.Rot=1;
        [N_orient(1), N_orient(2), N_orient(3)]=dir2cart(45.0, -4.0, 1);
        Meth=1;
    elseif (strcmp(files(f).name(1:2), 'SW'))
        % X-axis
        F_exp=46.0;
        F_orient=[1,0,0];
        flags.Rot=1;
        [N_orient(1), N_orient(2), N_orient(3)]=dir2cart(44.9, -4.8, 1);
        Meth=1;
    elseif (strcmp(files(f).name(1:3), 'ET1'))
        % -Z-axis
        F_exp=43.3;
        Meth=1;
    elseif (strcmp(files(f).name(1:3), 'ET2'))
        % -Z-axis
        F_exp=44.1;
        Meth=1;
    elseif (strcmp(files(f).name(1:3), 'ET3'))
        % -Z-axis
        F_exp=44.2;
        Meth=1;
    elseif (strcmp(files(f).name(1:2), 'A-') || strcmp(files(f).name(1:2), 'B-') || strcmp(files(f).name(1:2), 'C-'))
        % -Z-axis
        F_exp=36.2;
        Az=Az-90;
        Meth=1;
    elseif (strcmp(files(f).name(1:3), 'HEL') || strcmp(files(f).name(1:3), 'BR0'))
        % +X-axis
        F_exp=49.6;
        F_orient=[1,0,0];
        Meth=1;
    elseif (strcmp(files(f).name(1:3), 'SBG'))
        % +Z-axis
        F_exp=37.0;
        F_orient=[0,0,1];
        Meth=1;
    elseif (strcmp(files(f).name(1:2), 'rs'))
        % -Z-axis
        F_orient=[0,0,-1];
        flags.Anis=1;
        Meth=2;
        if (strcmp(files(f).name(1:4), 'rs25'))
            F_exp=30.0;
            flags.NLT=1;
        elseif(strcmp(files(f).name(1:4), 'rs26'))
            F_exp=60.0;
            flags.NLT=1;
        elseif(strcmp(files(f).name(1:4), 'rs27'))
            F_exp=90.0;
        else
            disp(files(f).name(1:end))
            error('IntAnalysis:F_exp', 'Unrecogonised sample name: Shaar et al., 2010.');
        end
    elseif (strcmp(files(f).name(1:2), 'AL'))
        % +Z-axis
        F_orient=[0,0,1];
        F_exp=35.8;
        Meth=1;
    elseif (strcmp(files(f).name(1:4), 'm428'))
        % +z-axis
        F_orient=[0, 0, 1];
        flags.Anis=1;
        flags.Rot=1;
        N_orient=F_orient;
        F_exp=25;
        Meth=1;
    elseif (strcmp(files(f).name(1:2), 'kf'))
        % Tanaka
        F_orient=[0,0,1];
        F_exp=52.1;
        Meth=1;
        [N_orient(1), N_orient(2), N_orient(3)]=dir2cart(339.9, 76.3, 1); % From IGRF
    else
        disp(files(f).name(1:end))
        error('IntAnalysis:F_exp', 'Unrecogonised sample name.');
    end
    
    
    % Get the ansiotropy tensor and the NLT data
    if flags.Anis==1
        
        match=zeros(length(Tensors.ID),1);
        for i=1:length(Tensors.ID)
            match(i)=strcmp(deblank(Tensors.ID(i,1:end)), deblank(files(f).name(1:end-4)));
        end
        
        if sum(match)==0
            error('PInt:Anis', 'Missing anisotropy tensors for %s', files(f).name(1:end-4));
        end
        
        A_ind=find(match==1);
        s_tensor=[Tensors.x1(A_ind), Tensors.x2(A_ind), Tensors.x3(A_ind), Tensors.x4(A_ind), Tensors.x5(A_ind), Tensors.x6(A_ind)];
    end
    
    
    if flags.NLT==1
        
        match=zeros(length(NLT.ID),1);
        for i=1:length(NLT.ID)
            match(i)=strcmp(deblank(NLT.ID(i,1:end)), deblank(files(f).name(1:end-4)));
        end
        
        if sum(match)==0
            error('PInt:NLT', 'Missing non-linear TRM data for %s', files(f).name(1:end-4));
        end
        
        NLT_ind=find(match==1);
        NLT_hat=[NLT.alpha(NLT_ind), NLT.beta(NLT_ind) ];
    end
    
    
    % Get the method for undefined samples
    if Meth==0
        
        
        % Assign Axis to temp variable and remove first NRM and checks
        Axis_T=Treatment(2:end);
        Axis_T(Axis_T==2 | Axis_T==3 | Axis_T==4)=[];
        
        if sum(Axis_T==5)>0
            error('PInt_Params:Meth_Check', 'Thellier-Thellier method detected!');
        end
        
        DA=diff(Axis_T); % Take the approx. deriv.
        
        if sum(DA(1:2:end)) == length(Axis_T)/2
            % Coe protocol
            Meth=1;
        elseif sum(DA(1:2:end)) == -length(Axis_T)/2
            % Aitken protocol
            Meth=3;
        elseif sum(DA(1:2:end))==-1 || sum(DA(1:2:end))==0 || sum(DA(1:2:end))==1
            % IZZI protocol
            Meth=2;
        else
            error('PInt_Params:Meth_Check', 'Unrecognized method detected!');
        end
    end
    
    
    %% Call GetPintParams
    % Get the params from point 1 to 7
    
    max_tmax = length(Temps(Treatment==0));
    
    
    %for i=1:max_tmax:
    
    
    Params=GetPintParams_v6(Mvec, Temps, Treatment, 1, 7, F_orient, Flab, flags.Rot, Az, Pl, N_orient, flags.Anis, s_tensor, flags.NLT, NLT_hat, .1);
    %Params=GetPintParams_v5c(Mvec, Temps, Treatment, 1, 7, F_orient, Flab, flags.Rot, Az, Pl, N_orient, flags.Anis, s_tensor, flags.NLT, NLT_hat);
           %GetPintParams_v6(Mvec, Temps, Treatment, start_pt, end_pt, Blab_orient, Blab, NRM_rot_flag, Az, Pl, ChRM, A_corr, s_tensor, NLT_corr, NLT_hat, beta_T)

    
    % norming by NRM
    NRM = Params.Ypts(1);
    Params.NRM = NRM;
    Params.Xpts = Params.Xpts / NRM;
    Params.Ypts = Params.Ypts / NRM;
    Params.Y_int = Params.Y_int / NRM;
    Params.X_int = Params.X_int / NRM;
    Params.VDS = Params.VDS / NRM;
    Params.x_prime = Params.x_prime / NRM;
    Params.y_prime = Params.y_prime / NRM;
    %extra
    Params.x_err = Params.x_err / NRM;
    Params.y_err = Params.y_err / NRM;
    %Params.x_tag = Params.x_tag / NRM;
    %Params.y_tag = Params.y_tag / NRM;
    
    % end extra
    Params.delta_x_prime = Params.delta_x_prime / NRM;
    Params.delta_y_prime = Params.delta_y_prime / NRM;
    %Params.lj_delta_x_prime = abs(Params.x_prime(end) - Params.x_prime(1));
    %Params.lj_delta_y_prime = abs(Params.y_prime(end) - Params.y_prime(1));
    
    
    
    % norming ends


    
    
    %
    split = strsplit(files(f).name, '.');
    spec_name = split{1}; % removes .tdt from spec_name
    Params.s = spec_name;
    spec_name = strrep(spec_name, '-', '_'); % - in name will cause matlab error
    spec_name = strcat(spec_name, '_', '1', '_', '7');
    Specimens.(spec_name) = Params;

    %end
    
    
    
    
    
end

clear FID;
fclose all;

toc

% lori

%lj how to make FRAC the same as in my code
% basically, NRMvec is normed and one addition vector diff is included
% compared to Greig's code
%NRMint = Int(Treatment==0);
%NRMvec = Params.NRMvec(:,2:end);
%zdata = NRMvec(:,:) / NRM;
%seg_min = 1;
%seg_max = 7;
%FRAC=sum(sqrt(sum((diff(zdata(seg_min:seg_max,1:end)).^2),2)))/ Params.VDS;
%FRAC;
%vector_diffs = sqrt(sum((diff(zdata(seg_min:seg_max,1:end)).^2),2));
%lj  that's it    


%for i=1:length(NRMint)
%    NRMint(i)
%    z = NRMvec(i,:) / NRMint(i)
%    zdata = [zdata; z]
%end


% Greig's names for everything
basic_reqd = {'s', 'n', 'missing start', 'missing end', 'Tmin', 'Tmax'};
arai_reqd = {'x_err', 'y_err', 'b', 'sigma_b', 'Banc', 'sigma_B', 'Y_int', 'X_int', 'VDS', 'delta_x_prime', 'delta_y_prime', 'f', 'f_vds', 'frac', 'beta', 'missing g', 'GAP-MAX', 'qual', 'w', 'k', 'SSE', 'SCAT', 'R_corr', 'R_det', 'Z', 'Z_star', 'IZZI_MD'};
directional_reqd = {'Dec_F', 'Dec_A', 'Inc_F', 'Inc_A', 'MAD_free', 'MAD_anc', 'alpha', 'missing theta', 'DANG', 'NRM_dev', 'gamma' };
ptrm_reqd = {'n_pTRM', 'check', 'dCK', 'DRAT', 'maxDev', 'CDRAT', 'CDRAT_prime', 'DRATS', 'DRATS_prime', 'mean_DRAT', 'mean_DRAT_prime', 'mean_DEV', 'mean_DEV_prime', 'dpal'};
tail_reqd = {'n_tail', 'DRATtail', 'dTR', 'MDvds'};
additivity_reqd = {'n_add', 'dAC'};
all_reqd = {'s', 'n', 'Tmin', 'Tmax', 'b', 'sigma_b', 'Blab', 'Banc', 'sigma_B', 'Y_int', 'X_int', 'VDS', 'delta_x_prime', 'delta_y_prime', 'f', 'f_vds', 'FRAC', 'beta', 'gap', 'GAP_MAX', 'qual', 'w', 'k', 'SSE', 'SCAT', 'R_corr', 'R_det', 'Z', 'Z_star', 'IZZI_MD', 'Dec_F', 'Dec_A', 'Inc_F', 'Inc_A', 'MAD_free', 'MAD_anc', 'alpha', 'Theta', 'DANG', 'NRM_dev', 'gamma', 'n_pTRM', 'check', 'dCK', 'DRAT', 'maxDEV', 'CDRAT', 'CDRAT_prime', 'DRATS', 'DRATS_prime', 'mean_DRAT', 'mean_DRAT_prime', 'mean_DEV', 'mean_DEV_prime', 'dpal', 'n_tail', 'DRATtail', 'dTR', 'MDvds', 'n_add', 'dAC'};

% my names for everything
basic_stats = {'s', 'specimen_n', 'start', 'end', 'tmin', 'tmax'};
arai_plot_stats = {'specimen_b', 'specimen_b_sigma', 'B_anc', 'B_anc_sigma', 'specimen_YT', 'specimen_XT', 'specimen_vds', 'x_prime', 'y_prime', 'delta_x_prime', 'delta_y_prime', 'specimen_f', 'specimen_fvds', 'FRAC', 'specimen_b_beta', 'specimen_g', 'GAP-MAX', 'specimen_q', 'specimen_w', 'specimen_k', 'SSE', 'SCAT', 'R_corr2', 'R_det2', 'Z', 'Zstar', 'IZZI_MD'};
directional_stats = {'Dec_Free', 'Dec_Anc', 'Inc_Free', 'Inc_Anc', 'MAD_Free', 'MAD_Anc', 'alpha', 'theta', 'DANG', 'NRM_dev', 'gamma'};
ptrm_stats = {'n_ptrm', 'max_ptrm_check_percent', 'delta_CK', 'DRAT', 'max_DEV', 'CDRAT', 'CDRAT_prime', 'DRATS', 'DRATS_prime', 'mean_DRAT', 'mean_DRAT_prime', 'mean_DEV', 'mean_DEV_prime', 'delta_pal'};
tail_stats = {'n_tail', 'DRAT_tail', 'delta_TR', 'MD_VDS'};
additivity_stats = {'n_add', 'delta_AC'};
all_stats = {'s', 'specimen_n', 'tmin', 'tmax', 'specimen_b', 'specimen_b_sigma', 'B_lab', 'B_anc', 'B_anc_sigma', 'specimen_YT', 'specimen_XT', 'specimen_vds',  'delta_x_prime', 'delta_y_prime', 'specimen_f', 'specimen_fvds', 'FRAC', 'specimen_b_beta', 'specimen_g', 'GAP-MAX', 'specimen_q', 'specimen_w', 'specimen_k', 'SSE', 'SCAT', 'R_corr2', 'R_det2', 'Z', 'Zstar', 'IZZI_MD', 'Dec_Free', 'Dec_Anc', 'Inc_Free', 'Inc_Anc', 'MAD_Free', 'MAD_Anc', 'alpha', 'theta', 'DANG', 'NRM_dev', 'gamma', 'n_ptrm', 'max_ptrm_check_percent', 'delta_CK', 'DRAT', 'max_DEV', 'CDRAT', 'CDRAT_prime', 'DRATS', 'DRATS_prime', 'mean_DRAT', 'mean_DRAT_prime', 'mean_DEV', 'mean_DEV_prime', 'delta_pal', 'n_tail', 'DRAT_tail', 'delta_TR', 'MD_VDS', 'n_add', 'delta_AC'}; 


outfile = fopen('Bowles_etal_2006_matlab.out', 'w');


% prints column names at top of outfile
for i=1:length(all_stats);
    fprintf(outfile, all_stats{i});
    fprintf(outfile, '\t');
    %disp(all_stats{i})
end
fprintf(outfile, '\n')


o = struct();

spec_names = fieldnames(Specimens);

for i=1:length(spec_names);

    for j=1:length(all_reqd);
        %disp(all_reqd(j))
        %class(all_reqd{j})
        %disp(all_stats(j))
        field_name = all_reqd{j};
        if field_name == 's';
            spec_names{i};
            field_name;
            o.(field_name) = spec_names{i};
        elseif length(field_name) > 6 & field_name(1:7) == 'missing';
            %disp('MISSING');
            %disp(field_name);
            o.(field_name) = 'missing';
        else
            %disp('valid field name');
            %disp(field_name);
            %disp(Params.(field_name));
            %o.(field_name) = Params.(field_name); % fixing here
            o.(field_name) = Specimens.(spec_names{i}).(field_name);
        
        end
        %disp(Params.(all_reqd{j}))
    
    end

    length(fieldnames(o));
    length(all_reqd);
    for k=1:length(all_reqd);
        fieldname = all_reqd{k};
        fprintf(outfile, [fieldname ': ']);
        fprintf(outfile, num2str(o.(fieldname)));
        fprintf(outfile, '\t');
    end
    fprintf(outfile, '\n');

end


fclose(outfile);