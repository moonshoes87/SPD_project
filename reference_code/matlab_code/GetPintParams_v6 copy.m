function [Params] = GetPintParams_v6(Mvec, Temps, Treatment, start_pt, end_pt, Blab_orient, Blab, NRM_rot_flag, Az, Pl, ChRM, A_corr, s_tensor, NLT_corr, NLT_hat, beta_T)
% function to return all paleointensity parameters for a sample
% v5 combines all output into a structure for easily data handling. It also
% now includes all custom subfunctions to increase portability.
%
%
% LAST UPDATE 09-11 Feb 2014
%
% 09-11 Feb 2014 - GAP
% 1) Updated version to v6
% 2) Corrected R_det to use the segment and not all of the points
% 3) Corrected the denominator for GAP_MAX
% 4) Changed theta to use the free-floating PCA fit and not the anchored
% 5) Corrected calculation for VDS, which affects FRAC and GAP_MAX
% 7) Added Params.Plot_orth Params.Plot_line to output the NRM vector and best-fit Arai line for plotting
% 8) Added Params.x_prime, Params.y_prime, Params.Delta_x_prime, and Params.Delta_y_prime, to replace Px, Py, TRM_len, and dyt, respectively.
% 9) Added beta_T for the SCAT box as an input variable with a default value of 0.1
% 10) Changed T_orient to Blab_orient
% 11) Added PearsonCorr2 function to get the square of the Pearson Correlation (R_corr). This removes a dependency on the inbuilt MATLAB function
% 12) Brought more variable names inline with SPD names
% 13) Updated the IZZI_MD calculation to follow the pseudo-code outline in SPD (logic rearrangement, no calculation change)
% 14) Removed Params.alpha(Params.alpha>90)=Params.alpha-90 for alpha, alpha_prime, DANG, and alpha_TRM. Legacy from early version of the models 
%     when direction of PCA was not constrained in the PmagPCA routine
% 15) Changed Ints(j) to Params.Ypts(j) in the calculation of CRM(j) for CRM_R. This SPD notation.
% 16) Made Params.Mean_DRAT calculation explicit (removed the use of CDRAT).
% 17) Added Params.dpal_ratio (= log(corrected_slope/measured_slope)) and Params.dpal_signed (signed dpal) to test later - I think dpal is biased
% 18) Change the anisotropy correction to use the free-floating PCA fit instead of the vector mean
% 19) Simplified the calculation of tstar(j) (part of dt*) in the last IF segment - it included terms that cancelled
% 20) Corrected parameters that use a normalizer that can be negative. Now they use the absolute value of the normalizer (f, Z, Zstar, NRM_dev, dCK, dt_star, dAC).
% 21) Added additional comment lines that are added for the SPD.m version of this code
% 22) Corrected several spelling mistakes in the comments and updated the descriptions for the output statistics
% 23) v6 of the code will be the basis for the version to made publicly available (as SPD.m). The public version will remove extra functions and
%     statistics that are not in SPD (e.g., tests for common slopes/elevations, commented code, etc.). The plotting outputs will also be removed.
%
% 23/24 Nov. 2013 - GAP
% 1) Added Params.pTRM_Lines to output the lines for plotting pTRM checks
% on Arai plots
% 2) Clarified orientation convention in subfunction dirot (rotates from core to geographic coords)
%
% 25 Sept. 2013 - GAP
% 1) Added Params.pTRM_sign, the sign of the maximum pTRM check stats
% 2) Added Params.CpTRM_sign, the sign of the cumulative pTRM check stats
% 3) Added Params.tail_sign, the sign of the maximum pTRM tail check stats
%
% 24 Sept. 2013 - GAP
% 1)Corrected an error in the handling of Thellier data, due to additivity checks not being added to the 
% data handling routine.The error caused a crash and did not influence any stats.
% 2) Rearranged the pTRM check data handling for Thellier data to make it more efficient.
% 3) Added a missing line in Thellier data handling for the calculation of check(%)
%
% 28-31 Aug. 2013 - GAP
% 1) Added rounding of the stats. Placed in a single code block at the end to allow easy adjustments
% 3) Updated terminology to fit with SPD. 
% 3) Corrected variable name "corr_Params.slope" to "corr_slope" lines 963/4
% 4) Change Params.MDpoints to output the measure NRM for the tail checks
% 5) Added Params.ADpoints to output additivity data and overhauled/corrected/updated the additivity check calculations
% 6) Added Params.PLOT, an array containing the data required to create an Arai plot
% 7) Updated pTRM checks such that they are only included if both the check temp and the peak experiment temp are <= Tmax
% 8) Updated IZZI_MD to correct the calculation of the the ZI curve length
%

%% Input

% Mvec          - a (n x 3) matrix of magnetization values when n is the total number of steps, which includes NRM, TRM, and all check measurements (i.e., the raw data measurements)
% Temps         - a (n x 1) vector of the temperature of each treatment
% Treatment     - a (n x 1) vector of integer values describing the treatment type at each step, this follows the convention of the ThellierTool
%                 (0=NRM demag, 1=TRM remag, 2=pTRM check, 3=pTRM tail check, 4=additivity check, 5=inverse TRM step). If any Treatment is set to '5' all data are treated as a Thellier experiment
% start_pt      - the temperature of the start point for the Arai plot best-fit line - ONLY THE INDEX IS CODED, NOT THE TEMPERATURE. THIS IS BETTER FOR TESTING
% end_pt        - the temperature of the end point for the Arai plot best-fit line - ONLY THE INDEX IS CODED, NOT THE TEMPERATURE. THIS IS BETTER FOR TESTING
% Blab_orient      - a (1 x 3) unit vector containing the x, y, z values for a know Blab orientation.
% Blab          - the strength of the laboratory field in muT
% NRM_rot_flag  - flag for directional rotation from core coords to stratigraphic coords (0 or []=no rotation, 1=apply rotation)
% Az            - the sample azimuth in degrees
% Pl            - the sample plunge in degrees - angle from horizontal, positive down
% ChRM          - a (1 x 2) or a (1 x 3) vector describing an independent measure of the direction of Banc (can be a known direction).
%                 If ChRM is a (1 x 2) vector it is assumed to be a Dec/Inc direction and converted to a Cartesian vector.
%                 If rotation is applied ChRM is assumed to already be in the final coordinate system
% A_corr        - flag for anisotropy correction (0=no correction, 1=apply correction)
% s_tensor      - a (1 x 6) vector with the six unique elements that describe an anisotropy tensor
% NLT_corr      - flag for non-linear TRM correction (0=no correction, 1=apply correction)
% NLT_hat       - a (1 x 2) vector containing the coefficients of the best-fit hyperbolic tangent function
% beta_T        - a (1 x 1) scalar for the beta threshold for defining the SCAT box, if missing the default value is beta_T=0.1

%% Output
% Output is a MATLAB structure with fields listed below. Use Params.Xpts to access Xpts, for example.
% See the Standard Paleointensity Definitions for full details of the statistics.
%
% Xpts         -   TRM points on the Arai plot
% Ypts         -   NRM points on the Arai plot
% N            -   number of points for the best-fit line on the Arai plot
% nmax         - the total number of NRM-TRM points on the Arai plot
% Seg_Ends     -   the indices for the start and end of the best-fit line [2x1]
% b            -   the slope of the best-fit line
% sigma_b      -   the standard error of the slope
% beta         -   sigma_b/b
% X_int        -   the intercept of the best-fit line of the TRM axis
% Y_int        -   the intercept of the best-fit line of the NRM axis
% x_prime      -   the Arai plot TRM points projected onto the best-fit line
% y_prime      -   the Arai plot NRM points projected onto the best-fit line
% Delta_x_prime -   the TRM length of the best-fit line
% Delta_y_prime -   the NRM length of the best-fit line
% VDS          -   the vector difference sum of the NRM vector
% f            -   the fraction of NRM used for the best-fit line
% f_vds        -   the NRM fraction normalized by the VDS
% FRAC         -   the NRM fraction with full vector calculation (Shaar & Tauxe, 2013; G-Cubed)
% gap          -   the gap factor
% GAP_MAX      -   the maximum gap (Shaar & Tauxe, 2013; G-Cubed)
% qual         -   the quality factor
% w            -   the weighting factor of Prevot et al. (1985; JGR)
% R_corr      -   Linear correlation coefficient (Pearson correlation)
% R_det        -   Coefficient of determination of the SMA linear model fit
% Line_Len     -   the length of the best-fit line
% f_vds        -   the fraction of vector difference sum NRM used for the best-fit line
% k            -   the curvature of the Arai plot following Paterson (2012; JGR)
% SSE          -   the fit of the circle used to determine curvature (Paterson, 2012; JGR)

% MAD_anc      -   the maximum angular deviation of the anchored PCA directional fit
% MAD_free     -   the maximum angular deviation of the free-floating PCA directional fit
% alpha        -   the angle between the anchored and free-floating PCA directional fits
% alpha_prime  -   the angle between the anchored PCA directional fit and the true NRM direction (assumed to be well known)
% alpha_TRM    -   the angle between the applied field and the acquire TRM direction (determined as an anchored PCA fit to the TRM vector)
% DANG         -   the deviation angle (Tauxe & Staudigel, 2004; G-Cubed)
% Theta        -   the angle between the applied field and the NRM direction (determined as a free-floating PCA fit to the TRM vector)
% a95          -   the alpha95 of the Fisher mean of the NRM direction of the best-fit segment
% NRM_dev      -   the intensity deviation of the free-floating principal component from the origin, normalized by Y_int (Tanaka & Kobayashi, 2003; EPS)
% CRM_R        -   the potential CRM% as defined by Coe et al. (1984; JGR)

% BY_Z       -   the zigzag parameter of Ben-Yosef et al. (2008; JGR)
% Z          -   the zigzag parameter of Yu & Tauxe (2005; G-Cubed)
% Z_star     -   the zigzag parameter of Yu (2012; JGR)
% IZZI_MD    -   the zigzag parameter of Shaar et al. (2011, EPSL)
% MD_area    -   the unsigned Arai plot area normalized by the length of the best-fit line (following the triangle method of Shaar et al. (2011, EPSL))

% N_alt           -   the number of pTRM checks used to the maximum temperature of the best-fit segment
% PCpoints        -   (N_alt x 3) array. First column is temperature the check is to (Ti), second is the temperature the check is from (Tj), third is total pTRM gained at each pTRM check (= pTRM_check_i,j in SPD)
% check           -   the maximum pTRM difference when normalized by pTRM acquired at each check step
% dCK             -   pTRM check normalized by the total TRM (i.e., X_int) (Leonhardt et al., 2004; G-Cubed)
% DRAT            -   pTRM check normalized by the length of the best-fit line (Selkin & Tauxe, 2000; Phil Trans R Soc London)
% maxDEV          -  maximum pTRM difference normalized by the length of the TRM segment used for the best-fit slope (Blanco et al., 2012; PEPI)
% CDRAT           -   cumulative DRAT (Kissel & Laj, 2004; PEPI)
% CDRAT_prime     -   cumulative absolute DRAT
% DRATS           -   cumulative pTRM check normalized by the maximum pTRM of the best-fit line segment (REF - Tauxe...)
% DRATS_prime     -   cumulative absolute pTRM check normalized by the maximum pTRM of the best-fit line segment
% mean_DRAT       -   CDRAT divided by number of checks
% mean_DRAT_prime -   the average DRAT (Herro-Bervera & Valet, 2009; EPSL)
% mean_DEV        -   average pTRM difference normalized by the length of the TRM segment used for the best-fit slope (Blanco et al., 2012; PEPI)
% mean_DEV_prime  -   average absolute pTRM difference normalized by the length of the TRM segment used for the best-fit slope
% dpal            -   cumulative check of Leonhardt et al. (2004; G-Cubed)

% N_MD       -   the number of tail checks used to the maximum temperature of the best-fit segment
% MDpoints   -   (N_MD x 2) array. First column is temperature, second is the remanence remaining after the demagnetization step for a pTRM tail check (= tail_check_i in SPD)
% dTR        -   tail check normalized by the total NRM (i.e., Y_int) (Leonhardt et al., 2004; G-Cubed)
% DRATtail   -   tail check normalized by the length of the best-fit line (Biggin et al., 2007; EPSL)
% MDvds      -   tail check normalized by the vector difference sum corrected NRM (REF - Tauxe...)
% dt_star    -   pTRM tail after correction for angular dependence (Leonhardt et al., 2004; G-Cubed)

% N_AC       - the number of additivity checks used to the maximum temperature of the best-fit segment
% ADpoints   - (N_AC x 3) array. First 2 columns are lower and upper temperatures for remaining pTRM, third is the remanence remaining after a repeat demagnetization additivity check (= Mrem in SPD)
% d_AC       - maximum additivity check normalized by the total TRM (i.e., X_int) (Leonhardt et al., 2004; G-Cubed)

% SCAT            - SCAT parameters of Shaar & Tauxe (2013). N.B. the beta threshold is hard-wired to 0.1
% com_slope_pval  - the probability of a common slope when comparing the NRM-TRM slope, with that defined by pTRM and pTRM tail checks (Warton et al., 2006; Biol. Rev.)
% com_elev_pval   - the probability of a common elevation (intercept) when comparing the NRM-TRM slope, with that defined by pTRM and pTRM tail checks (Warton et al., 2006; Biol. Rev.)

% Hanc         -   the unit vector in the direction of Banc as calculated by Selkin et al. (2000; EPSL) for the correction of anisotropy
% anis_scale   -   factor used to scale the TRM vectors to correct for the effects of anisotropy
% IMG_flag     -   flag to determine if non-linear TRM correction returns a complex number (1 if true, 0 if false)


%% Initial input variable checks

if size(Mvec,2)~=3
    error('GetPintParams:Input', 'Input magnetization data must be a [n x 3] magnetization vector')
end

% Mvec, Temps, and Treatment should all have the same length
if ~isequal( size(Mvec,1), length(Temps), length(Treatment) )
    error('GetPintParams:Input', 'Input magnetization matrix, temperature vector, and Treatment vector must be the same length')
end

%   start undefined  OR end undefined OR start @ neg temp OR end @ neg temp OR end comes before start
if isempty(start_pt) || isempty(end_pt) || start_pt < 0 || end_pt < 0 || end_pt <= start_pt
    error('GetPintParams:Input', 'Start and end points for best-fit analysis must be properly defined')
end


if isempty(Az) || isempty(Pl)
    NRM_rot_flag=0;
    if     NRM_rot_flag==1
        warning('GetPintParams:Input', 'No sample azimuth or plunge defined: no rotation applied');
    end
end


if ~isempty(ChRM)
    if length(ChRM)==2 % Assume to be Dec, Inc
        tmp_var=ChRM./norm(ChRM);
        [ChRM(1), ChRM(2), ChRM(3)]=dir2cart(tmp_var(1), tmp_var(2), 1);
    elseif length(ChRM)==3
        % a 3D unit vector, but normalize just in case it is not a unit vector
        ChRM=ChRM./norm(ChRM);
    else
        warning('GetPintParams:Input', 'Unrecognised ChRM format. Statistics that require ChRM will be ignored')
        ChRM=[];
    end
end


if isempty(A_corr) % No anisotropy correction
    A_corr=0;
    s_tensor=[];
end

if A_corr==1 && isempty(s_tensor)
    A_corr=0;
    warning('GetPintParams:Input', 'Anisotropy correction is desired, but no tensor is provided. No correction will be applied')
end


if isempty(NLT_corr) % No non-linear TRM correction
    NLT_corr=0;
    NLT_hat=[];
end

if NLT_corr==1 && isempty(NLT_hat)
    NLT_corr=0;
    warning('GetPintParams:Input', 'Non-linear TRM correction is desired, but no linearity estimate is provided. No correction will be applied')
end

if isempty(beta_T)
    beta_T=0.1;
end

% Set the experimental flag for separating the data
% Exp_Flag=0 for Coe, Aitken, and IZZI
% Exp_Flag>0 for Thellier
Exp_Flag=sum(Treatment==5);

%% Input data handling
% Separate the data into NRM, TRM, pTRM and tail measurements etc
% For the TRM, NRM, and tail check  matrices, the first column is the temperature
% For pTRM and additivity checks, the first column is the temperature the check is TO, the second is the temperature the check is FROM

if Exp_Flag==0
    % Coe/Aitken/IZZI Experiment
    % Find the NRM steps
    NRMvec=[];
    
    % solves differences between model and real data - for model data the
    % initial NRM has Axis=1 (generates the initial NRM), but for real data Axis=0
    if Treatment(1) ~=0
        NRMvec(1,:)=[0, Mvec(1,:)];
    end
    NRMvec=[NRMvec; Temps(Treatment==0), Mvec(Treatment==0,:)]; %#ok<*AGROW>
    tail_vec=[Temps(Treatment==3), Mvec(Treatment==3, :)];
    
    
    % Calculate the tail differences
    %     MD_vec=tail_vec;
    %     MD_vec(:,2:end)=NaN;
    MD_scalar=NaN(length(tail_vec),2);
    for n=1:size(tail_vec,1)
        %         MD_vec(n, 2:end)=tail_vec(n, 2:end) - NRMvec( NRMvec(:,1)==tail_vec(n,1), 2:end );
        MD_scalar(n,1)=tail_vec(n,1);
        MD_scalar(n,2)=sqrt(sum( tail_vec(n, 2:end).^2, 2)) - sqrt(sum( NRMvec( NRMvec(:,1)==tail_vec(n,1), 2:end ).^2, 2) );
    end
    
    
    TRMvec=[0,0,0,0]; % Set the first step to zeros
    pCheck=[]; % the full vector check
    SCAT_check=[];% a matrix for the SCAT parameter
    ALT_vec=[]; % the pTRM check difference
    ALT_scalar=[]; % the pTRM check difference as a scalar
    check_pct=[]; % the pTRM difference as a percentage of the pTRM at that step
    ADD_vec=[]; % the matrix for additivity check vectors
    ADD_scalar=[];
    for n=2:length(Treatment)
        
        if Treatment(n)==1 %TRM step
            TRMvec=[TRMvec; Temps(n), Mvec(n,:)-NRMvec(NRMvec(:,1)==Temps(n),2:end)];
        end
        
        if Treatment(n)==2 % pTRM check
            % In this experiment they are always performed after a
            % demagnetization experiment Axis==0 || ==3
            if Treatment(n-1) ~=0 && Treatment(n-1)~=3
                disp(['Temperature step ', num2str(Temps(n))])
                error('Unsupported experiment')
            end
            % Since it is a repeat TRM we can search TRMvec for the correct
            % temperature
            pCheck_vec=Mvec(n,:)-Mvec(n-1,:);
            pCheck=[pCheck; pCheck_vec];
            SCAT_check=[SCAT_check;Temps(n), Temps(n-1), sqrt(sum(pCheck_vec.^2, 2))]; 
            ALT_vec=[ALT_vec; Temps(n), Temps(n-1), pCheck_vec-TRMvec(TRMvec(:,1)==Temps(n),2:end)];
            ALT_scalar=[ALT_scalar; Temps(n), Temps(n-1), sqrt(sum(pCheck_vec.^2, 2)) - sqrt(sum( TRMvec(TRMvec(:,1)==Temps(n),2:end).^2, 2) )];
            check_pct =[check_pct; Temps(n), Temps(n-1), (sqrt(sum(pCheck_vec.^2, 2)) - sqrt(sum( TRMvec(TRMvec(:,1)==Temps(n),2:end).^2, 2) )) / sqrt(sum( TRMvec(TRMvec(:,1)==Temps(n),2:end).^2, 2))];
        end
        
        
        if Treatment(n)==4 % additivity check
            % Additivity check is a repeat demag and must always follow an
            % in-field step
            
            switch Treatment(n-1)
                
                case 1 % Previous step was a TRM step
                    %                                      total vec  -     previous NRM
                    ADD_vec=[ADD_vec; Temps(n), Temps(n-1), Mvec(n,:) - Mvec(Temps==Temps(n-1) & Treatment==0, :)]; % This observed Mrem
                    ADD_scalar=[ADD_scalar; Temps(n), Temps(n-1), sqrt(sum(Mvec(n,:).^2,2)) - sqrt(sum( Mvec(Temps==Temps(n-1) & Treatment==0, :).^2,2 ) ) ];
                    
                case 2 % Previous step was a pTRM check step
                    
                    if Treatment(n-2)==0
                        ADD_vec=[ADD_vec; Temps(n), Temps(n-2), Mvec(n,:) - Mvec(n-2, :)]; % This observed Mrem
                        ADD_scalar=[ADD_scalar; Temps(n), Temps(n-2), sqrt(sum(Mvec(n,:).^2,2)) - sqrt(sum( Mvec(n-2, :).^2,2 ) ) ];
                    elseif Treatment(n-2)==1 || Treatment(n-2)==3
                        %                                      total vec  -     previous NRM
                        ADD_vec=[ADD_vec; Temps(n), Temps(n-2), Mvec(n,:) - Mvec(Temps==Temps(n-2) & Treatment==0, :)]; % This observed Mrem
                        ADD_scalar=[ADD_scalar; Temps(n), Temps(n-2), sqrt(sum(Mvec(n,:).^2,2)) - sqrt(sum( Mvec(Temps==Temps(n-2) & Treatment==0, :).^2,2 ) ) ];
                    else
                        warning('GetPintParams:Additivity', 'Ignoring additivity check: Unsupported additivity check sequence')
                    end
                    
                    
                case 4
                    
                    if Treatment(n-2)==1
                        ADD_vec=[ADD_vec; Temps(n), Temps(n-2),  Mvec(n,:) - Mvec(Temps==Temps(n-2) & Treatment==0, :)]; % This observed Mrem;
                        ADD_scalar=[ADD_scalar; Temps(n), Temps(n-2), sqrt(sum(Mvec(n,:).^2,2)) - sqrt(sum( Mvec(Temps==Temps(n-2) & Treatment==0, :).^2,2 ) ) ];

                    else
                        warning('GetPintParams:Additivity', 'Ignoring additivity check following an additivity check: this sequence is not yet supported')
                    end
                    
                otherwise
                    disp(['Temperature step ', num2str(Temps(n))])
                    error('GetPintParams:Additivity', 'Unsupported additivity check sequence')
                    
            end
            
        end
        
    end
    
else
    % Thellier-Thellier Experiment
    
    Dirvec=[Temps(Treatment==1), Mvec(Treatment==1, :)];
    Dirvec(1,:)=[]; % Remove the first row, which is the initial NRM
    
    Invvec=[Temps(Treatment==5), Mvec(Treatment==5, :)];
    
    NRMvec=[];
    NRMvec(1,:)=[0, Mvec(1,:)];
    NRMvec=[NRMvec; Dirvec(:,1), (Dirvec(:,2:end) + Invvec(:,2:end))./2];
    
    TRMvec=[0,0,0,0]; %Set the first step to zeros
    TRMvec=[TRMvec; Dirvec(:,1), (Dirvec(:,2:end) - Invvec(:,2:end))./2];
    
    
    % Find the pTRM tail checks
    tail_vec=[Temps(Treatment==3), Mvec(Treatment==3, :)];
    % Calculate the tail differences
    %     MD_vec=tail_vec;
    %     MD_vec(:,2:end)=NaN;
    MD_scalar=NaN(length(tail_vec),2);
    ADD_vec=[]; % the matrix for additivity check vectors
    ADD_scalar=[];

    for n=1:size(tail_vec, 1)
        %         MD_vec(n, 2:end)=tail_vec(n, 2:end) - NRMvec( NRMvec(:,1)==tail_vec(n,1), 2:end );
        MD_scalar(n,1)=tail_vec(n,1);
        MD_scalar(n,2)=sqrt(sum( tail_vec(n, 2:end).^2, 2)) - sqrt(sum( NRMvec( NRMvec(:,1)==tail_vec(n,1), 2:end ).^2, 2) );
    end
    
    pCheck=[];% NaN(nAlt,3); % the full vector check
    SCAT_check=[];% a matrix for the SCAT parameter
    ALT_vec=[]; % the pTRM check difference
    ALT_scalar=[]; % the pTRM check difference as a scalar
    check_pct=[]; % the pTRM difference as a percentage of the pTRM at that step
    for n=2:length(Treatment)
        
        if Treatment(n)==2 % pTRM check
            % In this experiment they are always performed after an inverse
            % TRM acquisition Axis==5 or a pTRM tail check Axis==3
            % Since it is a repeat TRM we can search TRMvec for the correct
            % temp
            if Treatment(n-1)==5
                pCheck_vec=(Mvec(n,:)-Mvec(n-1,:))./2;
            elseif Treatment(n-1)==3
                pCheck_vec=(Mvec(n,:)-Mvec(n-1,:));
            else
                error('GetPintParams:Thellier_pTRM', 'Should not be here!\nTreatments: %f and %f\nTemperature: %f and %f', Treatment(n), Treatment(n-1), Temps(n), Temps(n-1));
            end
            
            ALT_vec=[ALT_vec; Temps(n), Temps(n-1), pCheck_vec-TRMvec(TRMvec(:,1)==Temps(n),2:end)];
            ALT_scalar=[ALT_scalar; Temps(n), Temps(n-1), sqrt(sum(pCheck_vec.^2, 2)) - sqrt(sum( TRMvec(TRMvec(:,1)==Temps(n),2:end).^2, 2) )];
            check_pct= [check_pct; Temps(n), Temps(n-1), (sqrt(sum(pCheck_vec.^2, 2)) - sqrt(sum( TRMvec(TRMvec(:,1)==Temps(n),2:end).^2, 2) )) / sqrt(sum( TRMvec(TRMvec(:,1)==Temps(n),2:end).^2, 2))];

            pCheck=[pCheck; pCheck_vec];
            SCAT_check=[SCAT_check; Temps(n), Temps(n-1), sqrt(sum(pCheck_vec.^2, 2))]; % a matrix for the SCAT parameter
        elseif Treatment(n)==4 % additivity check
            warning('GetPintParams:Additivity', 'Additivity checks not yet supported for Thellier-Thellier routine  - if required contact GAP')
        else
            % Skip over the direct and inverse steps, which are already taken care of
        end
            
        
    end
    
    
end


Params.Xpts=sqrt(sum(TRMvec(:,2:end).^2,2));
Params.Ypts=sqrt(sum(NRMvec(:,2:end).^2,2));
Params.nmax=length(Params.Xpts);
Params.Temp_steps=TRMvec(:,1);

Params.Blab=Blab;

% Output the various check points
if isempty(pCheck)
    Params.PCpoints=[];
else
    Params.PCpoints=[ALT_scalar(:,1), ALT_scalar(:,2), sqrt(sum(pCheck.^2,2))];
end

if isempty(tail_vec)
    Params.MDpoints=[];
else
    Params.MDpoints=[tail_vec(:,1), sqrt(sum(tail_vec(:,2:end).^2,2))];
end

if isempty(ADD_vec)
    Params.ADpoints=[];
else
    Params.ADpoints=[ADD_vec(:,1:2), sqrt( sum(ADD_vec(:,3:end).^2,2) )];
end


% Define the indices of the best-fit segment
seg_min=start_pt;
seg_max=end_pt;
seg=(seg_min:seg_max);

%% Anisotropy correction
Params.Anis_c=NaN;
Params.Hanc=NaN(1,3);
if A_corr==1
    % Follows the method of Veitch et al(1984; Arch. Sci., 37, 359-373) as recommended by 
    % Paterson (2013; Geophys. J. Int.; 193, 684-710, doi: 10.1093/gji/ggt033)   
    
    % Find the anisotropy corrected NRM direction
    % Get the NRM Decs and Incs then PCA the mean unit vector  (mX, mY, mZ)
    [Ds, Is, Int]=cart2dir(NRMvec(seg,2), NRMvec(seg,3), NRMvec(seg,4));
    [mD, mI]=PmagPCA(Ds, Is, Int, 'free');
    [mX, mY, mZ]=dir2cart(mD, mI, 1); % Mhat_ChRM
    
    A=Anis_mat(s_tensor); % The anisotropy tensor
    
    Params.Hanc=(A\[mX, mY, mZ]')';
    Params.Hanc=Params.Hanc./norm(Params.Hanc); % Unit vector in the direction of the ancient field
    
    Manc=(A*Params.Hanc')';
    Mlab=(A*Blab_orient')';
    
    Params.Anis_c=norm(Mlab)/norm(Manc);
end

%% Arai stats

Params.n=length(seg);
Params.Seg_Ends=[seg_min, seg_max];

% Max and min temperatures
Params.Tmin=TRMvec(seg_min,1);
Params.Tmax=TRMvec(seg_max,1);

X_seg=Params.Xpts(seg);
Y_seg=Params.Ypts(seg);
Params.xbar=mean(X_seg);
Params.ybar=mean(Y_seg);
U=detrend(X_seg,0); % (Xi-Xbar)
V=detrend(Y_seg,0); % (Yi-Ybar)

%lj 
Params.x_err = U;
Params.y_err = V;
Params.NRMvec = NRMvec;
Params.TRMvec = TRMvec;
%lj


% Get the paleointensity estimate
Params.b=sign(sum(U.*V))*std(Y_seg)/std(X_seg);
if A_corr==1 && NLT_corr==0
    Params.Banc = Blab * abs( Params.Anis_c*Params.b );
elseif A_corr==0 && NLT_corr==1
    Params.Banc = real( atanh( abs(Params.b) * tanh(NLT_hat(2)*Blab) ) / NLT_hat(2) );
elseif A_corr==1 && NLT_corr==1
    Params.Banc = real( atanh( abs(Params.Anis_c*Params.b) * tanh(NLT_hat(2)*Blab) ) / NLT_hat(2) );
else
    Params.Banc = abs(Params.b)*Blab;
end
    

% Params.IMG_flag=NaN; % flag to check imaginary fixes
% if ~isempty(NLT_hat)
%     if abs(Params.b)*tanh(NLT_hat(2)*Blab) >=1
%         Params.IMG_flag=1;
%     end
% end


Params.sigma_b=sqrt( (2*sum(V.^2)-2*(Params.b)*sum(U.*V)) / ( (Params.n-2)*sum(U.^2)) );
Params.beta=abs(Params.sigma_b/Params.b);
Params.sigma_B=Blab*Params.sigma_b;

Params.Y_int=mean(Y_seg)-Params.b*mean(X_seg);
Params.X_int=-Params.Y_int/Params.b;

% Project the data onto the best-fit line
Rev_x=(Params.Ypts-Params.Y_int)./Params.b; % The points reflected about the bes-fit line
Rev_y=Params.b.*Params.Xpts+Params.Y_int;
Params.x_prime=(Params.Xpts+Rev_x)./2; % Average the both sets to get the projected points
Params.y_prime=(Params.Ypts+Rev_y)./2;

% Get the TRM, NRM, and line lengths
Params.Delta_x_prime=abs( max(Params.x_prime(seg))-min(Params.x_prime(seg)) );
Params.Delta_y_prime=abs( max(Params.y_prime(seg))-min(Params.y_prime(seg)) );
%lj
Params.delta_y_prime = Params.Delta_y_prime;
Params.delta_x_prime = Params.Delta_x_prime;
%lj
Params.Line_Len=sqrt(Params.Delta_x_prime^2 + Params.Delta_y_prime^2);

% Get the VDS and fraction related stuff
Params.VDS=sum(sqrt(sum((diff(NRMvec(:,2:end)).^2),2)))+sqrt(sum(NRMvec(end,2:end).^2));
sumdy=sum( diff( Params.y_prime(seg) ).^2);

Params.f=abs(Params.Delta_y_prime/Params.Y_int);
Params.f_vds=abs(Params.Delta_y_prime/Params.VDS);
Params.FRAC=sum( sqrt( sum( (diff(NRMvec(seg_min:seg_max,2:end)).^2), 2 ) ) ) /Params.VDS;
Params.gap=1-(sumdy/(Params.Delta_y_prime^2));
Params.GAP_MAX=max( sqrt( sum( diff(NRMvec(seg_min:seg_max,2:end)).^2, 2 ) ) ) / sum(sqrt(sum((diff(NRMvec(seg_min:seg_max,2:end)).^2),2)));
Params.qual=Params.f*Params.gap/Params.beta;
Params.w=Params.qual/sqrt(Params.n-2);
Params.R_corr=PearsonCorr2(X_seg, Y_seg);
Params.R_det=1 - (sum((Y_seg-Params.y_prime(seg)).^2) / sum((Y_seg-Params.ybar).^2) );


% Curvature
kparams=AraiCurvature(Params.Xpts, Params.Ypts);
Params.k=kparams(1);
Params.SSE=kparams(4);

% curvature using only the best-fit segment
% kparams=AraiCurvature(Params.Xpts(seg), Params.Ypts(seg));
% Params.kprime=kparams(1);
% Params.SSEprime=kparams(4);

%% Directional stats
[Decs, Incs, Ints]=cart2dir(NRMvec(seg,2), NRMvec(seg,3), NRMvec(seg,4));
[TDecs, TIncs, TInts]=cart2dir(TRMvec(seg,2), TRMvec(seg,3), TRMvec(seg,4));

% Rotate the NRM directions - Only used for alpha_prime and CRM_R
if NRM_rot_flag==1
    [Rot_Decs, Rot_Incs]=dirot(Decs, Incs, Az, Pl);
else
    Rot_Decs=Decs;
    Rot_Incs=Incs;
end

if sum(isnan(Incs))+sum(isnan(Decs))>0 %PmagPCA can't handle NaN or inf, so skip - only affects models with no noise
    Params.MAD_anc=0;
    Params.MAD_free=0;
    Params.alpha=0;
    Params.alpha_prime=0;
    Params.DANG=0;
    Params.NRM_dev=0;
else
    [Params.Dec_A, Params.Inc_A, Params.MAD_anc]=PmagPCA(Decs, Incs, Ints, 'anc');
    [Params.Dec_F, Params.Inc_F, Params.MAD_free]=PmagPCA(Decs, Incs, Ints, 'free');
    
    Params.alpha=calc_angle([Params.Dec_A, Params.Inc_A], [Params.Dec_F, Params.Inc_F]);
%     Params.alpha(Params.alpha>90)=Params.alpha-90; % Take the smaller angle
    
    if isempty(ChRM)
        Params.alpha_prime=NaN;
%         Params.Rot_Dec_A, Params.Rot_Inc_A, Params.MAD_anc
    else
        % determine the PCA fit on the rotated directions
        [Rot_Dec_A, Rot_Inc_A, ~]=PmagPCA(Rot_Decs, Rot_Incs, Ints, 'anc');
        [NRM_dec, NRM_inc]=cart2dir(ChRM(1), ChRM(2), ChRM(3));
        Params.alpha_prime=calc_angle([Rot_Dec_A, Rot_Inc_A], [NRM_dec, NRM_inc]);
%         Params.alpha_prime(Params.alpha_prime>90)=Params.alpha_prime-90; % Take the smaller angle
    end
    
    
    % Calculate DANG
    dirfit=NaN(1,3);
    [dirfit(1), dirfit(2), dirfit(3)]=dir2cart(Params.Dec_F, Params.Inc_F, 1);
    Centre=mean(NRMvec(seg, 2:end)); %Centre of mass
    Params.DANG=rad2deg( atan2(norm(cross(dirfit, Centre)), dot(dirfit, Centre)) );
%     Params.DANG(Params.DANG>90)=Params.DANG-90; % Take the smaller angle
    
    Params.NRM_dev=(norm(Centre)*sin(deg2rad(Params.DANG))) / abs(Params.Y_int) * 100;
    
end
% Get the Fisher Mean and stats [Mdec, Minc, k, a95, R]
[~, ~, ~, Params.a95, ~]=FisherMeanDir(Decs, Incs);


% Do some directional stats on the TRM data

%PmagPCA can't handle NaN or inf, so remove - THIS IS FOR THE STOCHASTIC MODELS ONLY - only affects models with no noise
Bad_data=unique([find(isnan(TIncs)), find(isnan(TDecs)), find(isinf(TIncs)), find(isinf(TDecs)) ]);

TIncs(Bad_data)=[];
TDecs(Bad_data)=[];
TInts(Bad_data)=[];

if length(TIncs)< 3 %PmagPCA can't handle NaN or inf, so skip - should only affect models with no or very low noise
    Params.alpha_TRM=NaN;
else
    [Dec_TA, Inc_TA]=PmagPCA(flipud(TDecs), flipud(TIncs), flipud(TInts), 'anc'); % flip the data up/down (reverses order) so it behaves like demag data
    [TRM_dec, TRM_inc]=cart2dir(Blab_orient(1), Blab_orient(2), Blab_orient(3));
    Params.alpha_TRM=calc_angle([Dec_TA, Inc_TA], [TRM_dec, TRM_inc]);
%     Params.alpha_TRM(Params.alpha_TRM>90)=Params.alpha_TRM-90; % Take the smaller angle
end


% AnisotroParams.y_prime check
Params.gamma=rad2deg( atan2(norm(cross(TRMvec(seg_max,2:end), Blab_orient)), dot(TRMvec(seg_max,2:end), Blab_orient)) );


% Angle between measured NRM and Blab
[NRMhat(1), NRMhat(2), NRMhat(3)]=dir2cart(Params.Dec_F, Params.Inc_F);
Params.Theta=rad2deg( atan2(norm(cross(NRMhat, Blab_orient)), dot(NRMhat, Blab_orient)) );
%lj -- old version from v5c
Params.NRMhat = NRMhat;
Params.Blab_orient = Blab_orient;
%seems to be the same, but renamed T_orient
%Params.Theta=rad2deg( atan2(norm(cross(NRMhat, T_orient)), dot(NRMhat, T_orient)) );
%lj

% Coe et al. (1984) CRM parameter
if isempty(ChRM)
    % Need the definition of ChRM
    Params.CRM_R=NaN;
else
    
    if NRM_rot_flag==1 % We need to also rotate the Blab vector into geographic coords
        [tmp_D, tmp_I]=cart2dir(Blab_orient(1), Blab_orient(2), Blab_orient(3));
        [tmp_D, tmp_I]=dirot(tmp_D, tmp_I, Az, Pl);
        [tmp_O(1),tmp_O(2),tmp_O(3)]=dir2cart(tmp_D, tmp_I, 1);
        phi2=( atan2(norm(cross(ChRM, tmp_O)), dot(ChRM, tmp_O)) );
    else
        phi2=( atan2(norm(cross(ChRM, Blab_orient)), dot(ChRM, Blab_orient)) );
    end
    
    [fit_vec(:,1), fit_vec(:,2), fit_vec(:,3)]=dir2cart(Rot_Decs, Rot_Incs, 1); % Even if we don't rotate Rot_Decs/Incs contains the unrotated directions
    CRM=NaN(size(fit_vec,2),1); % create an empty vector
    
    for j=1:1:size(fit_vec,2)
        phi1=( atan2(norm(cross(fit_vec(j,:), ChRM)), dot(fit_vec(j,:), ChRM)) );
        CRM(j)=Params.Ypts(j)*sin(phi1)/sin((phi2));
    end
    
    Params.CRM_R=100*max(CRM)/(Params.Delta_x_prime);
    
end


%% Zig-Zag
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ben-Yosef (2008; JGR) method %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Based on PmagParams.y_prime implementation
%

% Preallocate variable
Frat=0;
Trat=0; %#ok<NASGU>

% Directions
% Get the Fisher stats of the alternating steps
[~, ~, ~, ~, R]=FisherMeanDir(Decs(1:end), Incs(1:end)); % For all the steps
[D1, I1, ~, ~, R1]=FisherMeanDir(Decs(1:2:end), Incs(1:2:end)); % For odd steps
[D2, I2, ~, ~, R2]=FisherMeanDir(Decs(2:2:end), Incs(2:2:end)); % For even steps

Params.BY_Z=NaN;
if length(Decs(1:2:end)) > 2 && length(Decs(2:2:end)) > 2
    Params.BY_Z=0;
    
    if calc_angle([D1, I1], [D2, I2])  > 3 % check that the angle between the mean directions is greater that 3 degrees
        F=(Params.n-2) * (R1 + R2 - R) / ( Params.n - R1 - R2 );
        Frat=F/finv(1-0.05, 2, 2*(Params.n-2));
        
        if Frat > 1
            Params.BY_Z=Frat;
        end
    end
    
end

% Slopes
dy=diff(Y_seg);
dx=diff(X_seg);
b_izzi=dy./dx;
b1=b_izzi(1:2:end);
b2=b_izzi(2:2:end);

% Suppress noise ala Tauxe
r1=sqrt( dy(1:2:end).^2 + dx(1:2:end).^2 );
r2=sqrt( dy(2:2:end).^2 + dx(2:2:end).^2 );

b1(r1<=0.1*Params.VDS)=[];
b2(r2<=0.1*Params.VDS)=[];

if length(b1) > 2 && length(b2) > 2
    
    if abs(atan(mean(b1))-atan(mean(b2))) > deg2rad(3) % If the angle between the mean Params.slopes is greater than 3 degrees
        [~, ~, ~, tb]=ttest2(b1, b2, 0.05);
        Trat=tb.tstat/tinv(1-0.05, Params.n-2);
        if Trat > 1 && Trat > Frat
            Params.BY_Z=Trat;
        end
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%
% Yu (2012; JGR) method %
%%%%%%%%%%%%%%%%%%%%%%%%%

bi=(Params.Y_int-Params.Ypts)./Params.Xpts;
bi(1)=0;
dummy_var=abs( (bi-abs(Params.b)).*Params.Xpts);

Params.Z=sum( dummy_var(seg_min:seg_max)./abs(Params.X_int)); % Yu & Tauxe (2005; G-Cubed)
Params.Z_star=sum(100.*sum( dummy_var(seg_min:seg_max))./abs(Params.Y_int))/(Params.n-1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shaar et al. (2011, EPSL) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tri_area=NaN(length(Params.Xpts)-3,1); % The triangle areas
Zsign=NaN(length(Params.Xpts)-3,1); % The sign of the areas
ZI_steps=NaN(length(Params.Xpts)-3,1); % The ZI step for the ZI line length

% create a vector of treatments and remove checks
Zig_Treatment=Treatment(2:end); % Ignore the first point as does Shaar
Zig_Treatment(Zig_Treatment==2)=[]; % remove all checks
Zig_Treatment(Zig_Treatment==3)=[];
Zig_Treatment(Zig_Treatment==4)=[];

% Zig_Treatment should have an even number of elements, i.e., it must have a
% demag and a remag
if mod(numel(Zig_Treatment), 2) ~=0
    error('GetPintParams:ZigZag', 'Zigzag measurement must be even');
end

Zig_Treatment=reshape(Zig_Treatment, 2, numel(Zig_Treatment)/2)';

IZZI_flag=0;
if sum(Zig_Treatment(:,1)==0) > 0  && sum(Zig_Treatment(:,2)==0) > 0
    % IZZI experiment
    IZZI_flag=1;
end

% Shaar uses the normalized Arai points and excludes the NRM point
Xn=Params.Xpts(2:end)./Params.Ypts(1);
Yn=Params.Ypts(2:end)./Params.Ypts(1);
% Xn=Params.Xpts(2:end);
% Yn=Params.Ypts(2:end);

for j=1:1:length(Xn)-2
    
    % Determine the best-fit line of the endpoints (of the 3 used for the triangle)and use the
    % intercept of this best-fit line (EP_fit(2) and the intercept of the line with the same slope, 
    % but through the middle point (Mid_int) to determine which points are higher
    EP_fit=polyfit(Xn([j,j+2]), Yn([j,j+2]), 1); % calculate the endpoint fit
    
    a1=EP_fit(2);
    a2=Yn(j+1)-EP_fit(1)*Xn(j+1); % Intercept when the line is passed through the midpoint
    
    if IZZI_flag==1
        
        % determine what step midpoint (j+1) is
        if Zig_Treatment(j+1,1)==1 && Zig_Treatment(j+1,2)==0
            % Mid point is IZ step
            
            % Assign surrounding points as ZI steps
            ZI_steps(j)=1;
            ZI_steps(j+2)=1;
            
            if a1 < a2 % Midpoint is above
                Zsign(j)=-1;
            else % midpoint is above
                Zsign(j)=1;
            end
            
        elseif Zig_Treatment(j+1,1)==0 && Zig_Treatment(j+1,2)==1
            % Midpoint is ZI step
            
            % Assign midpoint as a ZI step
            ZI_steps(j+1)=1;
            
            if a1 < a2 % Midpoint is above
                Zsign(j)=1;
            else % midpoint is above
                Zsign(j)=-1;
            end
            
        else
            error('GetPintParams:ZigZag', 'Unexpected IZZI sequence');
        end
        
    else
        ZI_steps(j)=1;
        
        Zsign(j)=1;
    end
    
    Tri_area(j)=polyarea(Xn(j:j+2), Yn(j:j+2));
end

ZI_inds=find(ZI_steps); % find the non-zero elements

ZI_len=sum( sqrt(diff(Xn(ZI_inds)).^2 + diff(Yn(ZI_inds)).^2) );

% ZI_ends=[Xn(ZI_inds([1,end])), Yn(ZI_inds([1,end]))];
% ZI_len=sqrt( diff(ZI_ends(:,1))^2 + diff(ZI_ends(:,2))^2 );


Params.MD_Area=sum(Tri_area)/Params.Line_Len;
Params.IZZI_MD=sum(Zsign.*Tri_area)/ZI_len;

%% SCAT

sigma_T=beta_T*abs(Params.b);
b1=Params.b - 2*sigma_T;
b2=Params.b + 2*sigma_T;

% determine the intercepts
a1=Params.ybar-b1*Params.xbar; % upper
a2=Params.ybar-b2*Params.xbar; % lower
a3=-a2/b2; % the upper 
a4=-a1/b1; % and lower x-axis intercepts

C1=[0, a2]; % lower left corner
C2=[0, a1]; % upper left corner
C3=[a3, 0]; % upper right corner
C4=[a4, 0]; % lower right corner
Params.SCAT_BOX=[C1; C2; C3; C4; C1];

% find the corners of the SCAT box - METHOD DESCRIBED IN THE PAPER, BUT WAS
% NOT USED IN THE CODING OF THELLIER_GUI
% C1=[X_seg(1), (-a2/a4)*X_seg(1)+a2 ]; % lower left corner
% C2=[X_seg(1), (-a1/a3)*X_seg(1)+a1 ]; % upper left corner
% C3=[X_seg(end), (-a1/a3)*X_seg(end)+a1]; % upper right corner
% C4=[X_seg(end), (-a2/a4)*X_seg(end)+a2 ]; % lower right corner
% if C4(2) < 0 % the SCAT box goes outside the region
%     % Set the y-coord of the 4th corner to zero
%     C4(2)=0;
%     
%     % Add an extra point
%     C5=[a4, 0];
%     Params.SCAT_BOX=[C1; C2; C3; C4; C5; C1];
%     
% else
%     Params.SCAT_BOX=[C1; C2; C3; C4; C1];
% end

Check_points=[]; % the x-, y-coords of all the checks within the SCAT range

% find the pTRM checks
if ~isempty(SCAT_check)
    tmp_x=SCAT_check(SCAT_check(:,1)<=Params.Tmax & SCAT_check(:,1)>=Params.Tmin & SCAT_check(:,2)<=Params.Tmax, 3); % Check performed within range, and not from above range
    tmp_T=SCAT_check(SCAT_check(:,1)<=Params.Tmax & SCAT_check(:,1)>=Params.Tmin & SCAT_check(:,2)<=Params.Tmax, 1);
    tmp_inds=[];
    for j=1:length(tmp_T)
        tmp_inds(j)=find(NRMvec(:,1)==tmp_T(j));
    end
else 
    tmp_x=[];
    tmp_inds=[];    
end

if ~isempty(tmp_x) && ~isempty(tmp_inds)
    tmp_y=sqrt(sum( NRMvec( tmp_inds, 2:end ).^2, 2));
    Check_points=[Check_points; [tmp_x, tmp_y]];
end
clear vars tmp_x tmp_y tmp_inds;

% find the tail checks
if ~isempty(tail_vec)
    tmp_y=sqrt(sum(tail_vec(tail_vec(:,1)<= Params.Tmax & tail_vec(:,1)>=Params.Tmin, 2:end).^2,2));
    tmp_T=tail_vec(tail_vec(:,1)<= Params.Tmax & tail_vec(:,1)>=Params.Tmin, 1);
    tmp_inds=[];
    for j=1:length(tmp_T)
        tmp_inds(j)=find(TRMvec(:,1)==tmp_T(j));
    end
else 
    tmp_y=[];
    tmp_inds=[];    
end
if ~isempty(tmp_y) && ~isempty(tmp_inds)
    tmp_x=sqrt(sum( TRMvec( tmp_inds, 2:end ).^2, 2));
    Check_points=[Check_points; [tmp_x, tmp_y]];
end
clear vars tmp_x tmp_y tmp_inds;

% Create an array with the points to test
Params.SCAT_points=[X_seg, Y_seg; Check_points]; % Add the TRM-NRM Arai plot points

IN=inpolygon(Params.SCAT_points(:,1), Params.SCAT_points(:,2), Params.SCAT_BOX(:,1), Params.SCAT_BOX(:,2));
Params.SCAT=floor(sum(IN)/length(IN)); % The ratio ranges from 0 to 1, the floor command rounds down to nearest integer (i.e., rounds to 0 or 1)

%% Common slope and common elevation tests
% Based on the algorithms of Warton et al. (2006; Biol. Rev.)
Params.Check_points=Check_points;
N_check_points=size(Check_points, 1);

Params.check_slope=NaN;
Params.check_int=NaN;
Params.com_slope_pval=NaN;
Params.com_elev_pval=NaN;

if ~isempty(Check_points) % do the test
    % descriptive parameters
    cU=detrend(Check_points(:,1),0);
    cV=detrend(Check_points(:,2),0);
    barX=[Params.xbar, mean(Check_points(:,1))];
    barY=[Params.ybar, mean(Check_points(:,2))];
    varX=[var(X_seg), var(Check_points(:,1))];
    varY=[var(Y_seg), var(Check_points(:,2))];
    varXY=[sum(U.*V)/(Params.n-1) , sum(cU.*cV)/(N_check_points-1)];
    Ns=[Params.n, N_check_points];
    
    Params.check_slope=sign(sum(cU.*cV))*std(Check_points(:,2))/std(Check_points(:,1));
    Params.check_int=barY(2)-Params.check_slope*barX(2);
    
    % make initial guess - the average of the two slopes
    bhat=mean([Params.b, Params.check_slope]);
    b_com=fminsearch(@(x) abs(common_slope(x, varX, varY, varXY, Ns)), bhat, optimset('TolX', 1e-10, 'TolFun', 1e-10, 'MaxIter', 1e3, 'Display', 'off'));
    
    % determine the correlation between the fitted and residual axes
    resid_axis=Y_seg-b_com.*X_seg;
    fit_axis=Y_seg+b_com.*X_seg;
    rrf(1)=corr(resid_axis, fit_axis);
    
    resid_axis=Check_points(:,2)-b_com.*Check_points(:,1);
    fit_axis=Check_points(:,2)+b_com.*Check_points(:,1);
    rrf(2)=corr(resid_axis, fit_axis);
    
    test_val=-sum((Ns-2.5).*log(1-rrf.^2));
    Params.com_slope_pval=1-chi2cdf(test_val, 1); % probability that the slopes are different
    
    % Common elevation
    % Pearson correlations of the data used to determine the slopes
    rxy(1)=corr(X_seg, Y_seg);
    rxy(2)=corr(Check_points(:,1), Check_points(:,2));
    
    % The variances on the individual slope estimates
    varB=(1./(Ns-2)) .* (varY./varX) .* (1-rxy.^2);
    
    % variance of the common slope (b_com)
    varB_com=1/(sum(1./varB));
    
    % variance of residuals
    varRes=(Ns-1)./(Ns-2) .* (varY - 2.*b_com.*varXY + (b_com.^2).*varX);
    
    tmp_X=[barX; barX];
    
    var_AS = diag(varRes./Ns) + (varB_com *  tmp_X.*tmp_X');
    
    % intercepts using the common slope
    Ahats=[barY(1)-b_com*barX(1), barY(2)-b_com*barX(2)];
    
    L=[ones(1,1), -eye(1)];
    % stat=(L*Ahats')' * inv(L*var_AS*L') * (L*Ahats');
    stat=((L*Ahats')' * (L*Ahats'))/(L*var_AS*L');
    
    Params.com_elev_pval=1-chi2cdf(stat, 1);
end

%% pTRM checks

Params.n_pTRM=0;
Params.check=NaN;
Params.dCK=NaN;
Params.DRAT=NaN;
Params.maxDEV=NaN;

Params.CDRAT=NaN;
Params.CDRAT_prime=NaN;
Params.DRATS=NaN;
Params.DRATS_prime=NaN;
Params.mean_DRAT=NaN;
Params.mean_DRAT_prime=NaN;
Params.mean_DEV=NaN;
Params.mean_DEV_prime=NaN;

Params.dpal=NaN;

Params.pTRM_sign=NaN;
Params.CpTRM_sign=NaN;


% The first two columns in check_pct, ALT_scalar, and ALT_vec are structured in the same way
% (:,1) - Temperature the check is to
% (:,2) - Temperature the check is from


if ~isempty(ALT_scalar)
    
    % Both temperatures must be less than the maximum
    pTRM_checks=ALT_scalar(ALT_scalar(:,1)<=Params.Tmax & ALT_scalar(:,2)<=Params.Tmax, 3);
    
    Params.n_pTRM=length(pTRM_checks);
    
    
    Params.check=100*max(abs(check_pct(check_pct(:,1)<=Params.Tmax & check_pct(:,2)<=Params.Tmax, 3)));
    Params.dCK=100*max(abs(pTRM_checks/Params.X_int));
    Params.DRAT=100*max(abs(pTRM_checks/Params.Line_Len));
    Params.maxDEV=100*max(abs(pTRM_checks/Params.Delta_x_prime));

    Params.CDRAT=abs( 100*sum(pTRM_checks)/Params.Line_Len );
    Params.CDRAT_prime=abs( 100*sum(abs(pTRM_checks))/Params.Line_Len );

    Params.DRATS=abs( 100*sum(pTRM_checks)/Params.Xpts(seg(end)) );
    Params.DRATS_prime=abs( 100*sum(abs(pTRM_checks))/Params.Xpts(seg(end)) );

    Params.mean_DRAT=100*abs(mean(pTRM_checks./Params.Line_Len));
    Params.mean_DRAT_prime=100*mean(abs(pTRM_checks./Params.Line_Len));
    
    Params.mean_DEV=100*abs(mean((pTRM_checks./Params.Delta_x_prime)));
    Params.mean_DEV_prime=100*mean(abs(pTRM_checks./Params.Delta_x_prime));    
    
    Params.pTRM_sign=sign(pTRM_checks(abs(pTRM_checks)==max(abs(pTRM_checks))));
    Params.CpTRM_sign=sign(sum(pTRM_checks));

    
    % dpal
    % the cumulative sum of the vector difference between the original TRM
    % and the repeat measurement
    
    % The matrix for the cumulative sum needs to be padded with zeros where
    % no pTRM checks were performed
    to_sum=zeros(length(TRMvec),3); %create empty matrix of zeros
    for j=1:size(ALT_vec,1)
        ind= TRMvec(:,1)==ALT_vec(j,1);
        to_sum(ind,:)=ALT_vec(j,3:end);
    end
    
    
    dPal_sum=cumsum(to_sum);
    corr_TRM=zeros(length(TRMvec),1);
    corr_TRM(1)=0;
    for j=2:1:size(TRMvec,1)
        corr_TRM(j)=sqrt(sum( (TRMvec(j,2:end)+dPal_sum(j-1,:)).^2, 2 ) );
    end
    
    Xcorr=corr_TRM(seg);
    Ucorr=detrend(Xcorr,0);
    corr_slope=sign(sum(Ucorr.*V))*std(Y_seg)/std(Xcorr);
    Params.dpal=abs(100*(Params.b-corr_slope)/Params.b);
    
    Params.dpal_signed=(100*(Params.b-corr_slope)/Params.b);
    Params.dpal_ratio=log(corr_slope/Params.b);
    
    
end

% Catch the cases with no checks
if Params.n_pTRM==0
    Params.check=NaN;
    Params.dCK=NaN;
    Params.DRAT=NaN;
    Params.maxDEV=NaN;
    
    Params.CDRAT=NaN;
    Params.CDRAT_prime=NaN;
    Params.DRATS=NaN;
    Params.DRATS_prime=NaN;
    Params.mean_DRAT=NaN;
    Params.mean_DRAT_prime=NaN;
    Params.mean_DEV=NaN;
    Params.mean_DEV_prime=NaN;
    
    Params.dpal=NaN;
    Params.dpal_signed=NaN;
    Params.dpal_ratio=NaN;
    
    Params.pTRM_sign=NaN;
    Params.CpTRM_sign=NaN;

end

%% pTRM tail checks

Params.n_tail=0;

if ~isempty(MD_scalar)
    % do the checks
    
    tail_checks=MD_scalar(MD_scalar(:,1)<=Params.Tmax, 2);
    
    Params.n_tail=length(tail_checks);
    
    Params.dTR=max(abs(100*tail_checks/Params.Y_int));
    Params.DRATtail=max(abs(100*tail_checks/Params.Line_Len));
    Params.MDvds=max(abs(100*tail_checks/Params.VDS));
    
    Params.tail_sign=sign(tail_checks(abs(tail_checks)==max(abs(tail_checks))));

    
    % dt*
    tstar=NaN(1,Params.nmax); % Assign the vector
    
    for j=1:1:Params.nmax
        
        if sum((NRMvec(j,1)==tail_vec(:,1))) <= 0
            % No check at this temperature step so skip
            
        elseif sum((NRMvec(j,1)==tail_vec(:,1))) > 1
            % should not be true - multiple tail checks to the same
            % temperature
            error('GetPintParams:tails', 'Too many tail checks?');
            
        else % We have a tail check
            tind=find(tail_vec(:,1)==NRMvec(j,1)); % the index of the tail check
            
            MDx=tail_vec(tind,2);
            MDy=tail_vec(tind,3);
            MDz=tail_vec(tind,4);
            
            % This is more accurate than dot product when theta is small
            theta_dt = atan2(norm(cross(NRMvec(j,2:end)./Params.Ypts(j),Blab_orient)),dot(NRMvec(j,2:end)./Params.Ypts(j), Blab_orient));
            
            % Define horizontal and vertical according to Blab_orient, such that Blab_orient is always "vertical"
            if abs(Blab_orient(1)/norm(Blab_orient))==1 % Blab is along x
                dH=sqrt(sum(NRMvec(j,2:3).^2))-sqrt(MDy^2+MDz^2);
                dZ=NRMvec(j,2)-MDx;
                [~, F_inc]=cart2dir(Blab_orient(2), Blab_orient(3), Blab_orient(1));
                [~, N_inc]=cart2dir(NRMvec(j,3), NRMvec(j,4), NRMvec(j,2));
                inc_diff=F_inc-N_inc;
            elseif abs(Blab_orient(2)/norm(Blab_orient))==1 % Blab is along y
                dH=sqrt(sum(NRMvec(j,2,4).^2))-sqrt(MDx^2+MDz^2);
                dZ=NRMvec(j,3)-MDy;
                [~, F_inc]=cart2dir(Blab_orient(2), Blab_orient(1), Blab_orient(2));
                [~, N_inc]=cart2dir(NRMvec(j,4), NRMvec(j,2), NRMvec(j,3));
                inc_diff=F_inc-N_inc;
            elseif abs(Blab_orient(3)/norm(Blab_orient))==1 % Blab is along z
                dH=sqrt(sum(NRMvec(j,2:3).^2,2))-sqrt(MDx^2+MDy^2);
                dZ=NRMvec(j,4)-MDz;
                [~, F_inc]=cart2dir(Blab_orient(1), Blab_orient(2), Blab_orient(3));
                [~, N_inc]=cart2dir(NRMvec(j,2), NRMvec(j,3), NRMvec(j,4));
                inc_diff=F_inc-N_inc;
            else
                error('GetPintParams:dt_star', 'Blab is expected to be along either x, y, or z. If you see this contact GAP.');
                
                % THIS IN NOT FULLY TESTED YET              
                % In this case Blab_orient is at an angle to x, y and z
                % to calculate dt* we transform coordinates such that Blab_orient
                % is along z, and we rotate the NRM and MD vector
%                 Rot_params=vrrotvec(Blab_orient./norm(Blab_orient), [0,0,1]); % The parameters needed to rotate the field vector to z
%                 
%                 % rotate the NRM
%                 if mod(Rot_params(4), pi/2)==0 % fields along x/y/z should be dealt with above
%                     % no need to rotate
%                     error('GetPintParams:dt_star', 'Unexpected field angle');
%                 else
%                     disp(['in the rot  ', num2str( rad2deg(Rot_params(4)) )])
%                     NRMrot=Rotate_Vec(NRMvec(j,2:4), Rot_params(4), Rot_params(1:3));
%                     MDrot=Rotate_Vec(tail_vec(tind,2:4), Rot_params(4), Rot_params(1:3));
%                 end
%                 
%                 MDx=MDrot(1);
%                 MDy=MDrot(2);
%                 MDz=MDrot(3);
%                 
%                 dH=sqrt(sum(NRMrot(1:2).^2,2))-sqrt(MDx^2+MDy^2);
%                 dZ=NRMrot(3)-MDz;
%                 [~, F_inc]=cart2dir(Blab_orient(1), Blab_orient(2), Blab_orient(3));
%                 [~, N_inc]=cart2dir(NRMrot(1), NRMrot(2), NRMrot(3));
%                 inc_diff=F_inc-N_inc;
                
            end
            
            
            %             dH=sqrt(sum(NRMvec(j,2:3).^2,2))-sqrt(MDx.^2+MDy.^2);
            %             dZ=NRMvec(j,4)-MDz;
            B=dH./(tan(theta_dt));
            
            %             [~, F_inc]=cart2dir(Blab_orient(1), Blab_orient(2), Blab_orient(3));
            %             [~, N_inc]=cart2dir(NRMvec(j,2), NRMvec(j,3), NRMvec(j,4));
            %             inc_diff=F_inc-N_inc; % Inclination difference in degrees
            
            % Roman Leonhardt's dt* implementation - as of v4.2
            if (floor(theta_dt*1000) < 2968 && floor(theta_dt*1000) > 175) %% TT Version 4.1
                if (inc_diff > 0)
                    tstar(j) = (-dZ + B) * abs(Params.b) * 100.0/abs(Params.Y_int); %% sign dependent  diff (new in Vers. 1.8)
                else
                    tstar(j) = (dZ - B) * abs(Params.b) * 100.0/abs(Params.Y_int); %% sign dependent  diff (new in Vers. 1.8)
                end
            else
                if (floor(theta_dt*1000) <= 175)
                    tstar(j) = 0;
                elseif (floor(theta_dt*1000) >= 2968) %% TT Version 4.1
                    tstar(j) =-dZ*100.0/(abs(Params.X_int)+abs(Params.Y_int));% -minmax.maxmagTH/(minmax.maxmagTH+minmax.maxmagPT)*dZ * 100.0/minmax.maxmagTH;
                end
            end
            
        end
    end
    
   
    Params.dt_star=max(tstar(NRMvec(:,1)<=Params.Tmax));
    
    if Params.dt_star<0
        Params.dt_star=0;
    end
    
end

% Catch the case where no pTRM tail checks are used in the segment or at all
if Params.n_tail==0
    Params.dTR=NaN;
    Params.DRATtail=NaN;
    Params.MDvds=NaN;
    Params.dt_star=NaN;
    Params.tail_sign=NaN;
end

%% Additivity checks

Params.n_add=NaN;
Params.dAC=NaN;

if ~isempty(ADD_vec)
    % The structure of the Array
    %  ADD_vec(:,1) - lower temperature of the remaining pTRM
    %  ADD_vec(:,2) - upper temperature of the remaining pTRM
    %  ADD_vec(:,3) - x component
    %  ADD_vec(:,4) - y component
    %  ADD_vec(:,5) - z component
    
    % We want to add the pTRM acquired at ADD_vec(:,1), i.e., TRMvec(TRMvec(:,1)==ADD_vec(:,1), :)
    
    tmp_N=size(ADD_vec, 1);
    AC=NaN(tmp_N, 1);
    for j=1:tmp_N
        
%         % The previously observed pTRM, i.e., pTRM(Tk, T0)
%         TRM1= TRMvec(Params.Temp_steps==ADD_vec(j,1), 2:end);
%         % pTRM(Ti, T0)
%         TRM2=TRMvec(Params.Temp_steps==ADD_vec(j,2), 2:end);
%         % Mrem - the observed pTRM remaining, i.e., pTRM(Ti,Tk)
%         Mrem=ADD_vec(j,3:end);
%         
%         % For SD grains Mrem = TRM2 - TRM1 (= 0)
%         % Additivity check is the difference between the LHS and RHS of the above
%         % Done by vector difference and then calculate the length
%         AC(j) = sqrt( sum( (Mrem - (TRM2 - TRM1)).^2, 2 ) );
        
        
        % Calculate it by scalar difference - this is consistent with pTRM and pTRM tail checks
        AC(j) = (ADD_scalar(j,3) - (Params.Xpts(Params.Temp_steps==ADD_vec(j,2)) - Params.Xpts(Params.Temp_steps==ADD_vec(j,1))) );
        
        
    end
    
    AC=AC(ADD_vec(:,1) <= Params.Tmax & ADD_vec(:,2) <= Params.Tmax, 1); % Select only those within selected temperature range
    
%         AC=AC(ADD_vec(:,1) < Params.Tmax, 1); % I believe that this is the selection used by the ThellierTool

    
    Params.dAC=100*max(abs(AC))/abs(Params.X_int);
    
    Params.n_add=length(AC);
    
end


%% Create an output for plotting

% Params.PLOT=NaN(Params.nmax, 5);
pTRM_plot=NaN(Params.nmax, 1);
tail_plot=NaN(Params.nmax, 1);
Line_pts=cell(1,4);

for i=1:Params.nmax
    if ~isempty(Params.PCpoints) && ~isempty(Params.PCpoints(Params.PCpoints(:,1)==Params.Temp_steps(i),3))
        pTRM_plot(i)=Params.PCpoints(Params.PCpoints(:,1)==Params.Temp_steps(i),3);
        
        % Get the points for the lines
        % Horizontal line
        From_temp=Params.PCpoints(Params.PCpoints(:,1)==Params.Temp_steps(i), 2); % The temp that the check was from
        
        Hy1=Params.Ypts(Params.Temp_steps==From_temp);
        Hy2=Hy1;
        
        Hx1=Params.Xpts(Params.Temp_steps==From_temp);
        Hx2=pTRM_plot(i);
        
        Vx1=Hx2;
        Vx2=Hx2;
        
        Vy1=Hy1;
        Vy2=Params.Ypts(Params.Temp_steps==Params.Temp_steps(i));
        
        % Vertical line
        
        Line_pts=[Line_pts;  num2cell( [Hx1, Hy1, Vx1, Vy1; Hx2, Hy2, Vx2, Vy2]); cell(1,4) ];
        
    end
    
    if ~isempty(Params.MDpoints) && ~isempty(Params.MDpoints(Params.MDpoints(:,1)==Params.Temp_steps(i),2))
        tail_plot(i)=Params.MDpoints(Params.MDpoints(:,1)==Params.Temp_steps(i),2);
    end
        
end
Params.Plot_Arai=[Params.Temp_steps, Params.Xpts, Params.Ypts, pTRM_plot, tail_plot];
Params.Plot_Line=[Params.x_prime([Params.Seg_Ends]), Params.y_prime([Params.Seg_Ends])];

Params.Plot_pTRM_Lines=Line_pts;

Params.Plot_orth=NRMvec;

% Get the Cartesian coords of the best fits
[x_a, y_a, z_a]=dir2cart(Params.Dec_A, Params.Inc_A, 1); % Get the unit vector for the anchored PCA fit
Params.Plot_PCA_anc=repmat([x_a, y_a, z_a], Params.n, 1).*repmat(Params.Ypts(seg), 1,3); % Replicate n times and scale by the selected segment
% Params.Plot_PCA_anc=[0, 0, 0; x_a, y_a, z_a].*max(Params.Ypts); % Replicate n times and scale by the selected segment

[x_f, y_f, z_f]=dir2cart(Params.Dec_F, Params.Inc_F, 1); % Get the unit vector for the free-floating PCA fit
Params.Plot_PCA_free=repmat([x_f, y_f, z_f], Params.n, 1).*repmat(Params.Ypts(seg), 1,3); % Replicate n times and scale by the selected segment
% Params.Plot_PCA_anc=[0, 0, 0; x_f, y_f, z_f].*max(Params.Ypts); % Replicate n times and scale by the selected segment



%% Round the stats to the recommended SPD precision

% Arai plot
Params.b=round(1000*Params.b)/1000;
Params.sigma_b=round(1000*Params.sigma_b)/1000;
Params.Blab=round(10*Params.Blab)/10;
Params.Banc=round(10*Params.Banc)/10;
Params.sigma_B=round(10*Params.sigma_B)/10;
Params.f=round(1000*Params.f)/1000;
Params.f_vds=round(1000*Params.f_vds)/1000;
Params.FRAC=round(1000*Params.FRAC)/1000;
Params.beta=round(1000*Params.beta)/1000;
Params.gap=round(1000*Params.gap)/1000;
Params.GAP_MAX=round(1000*Params.GAP_MAX)/1000;
Params.qual=round(10*Params.qual)/10;
Params.w=round(10*Params.w)/10;
Params.k=round(1000*Params.k)/1000;
Params.SSE=round(1000*Params.SSE)/1000;
Params.R_corr=round(1000*Params.R_corr)/1000;
Params.R_det=round(1000*Params.R_det)/1000;
Params.Z=round(10*Params.Z)/10;
Params.Z_star=round(10*Params.Z_star)/10;
Params.IZZI_MD=round(1000*Params.IZZI_MD)/1000;

% Directional
Params.Dec_A=round(10*Params.Dec_A)/10;
Params.Inc_A=round(10*Params.Inc_A)/10;
Params.Dec_F=round(10*Params.Dec_F)/10;
Params.Inc_F=round(10*Params.Inc_F)/10;

Params.MAD_anc=round(10*Params.MAD_anc)/10;
Params.MAD_free=round(10*Params.MAD_free)/10;
Params.alpha=round(10*Params.alpha)/10;
Params.alpha_prime=round(10*Params.alpha_prime)/10;

Params.Theta=round(10*Params.Theta)/10;
Params.DANG=round(10*Params.DANG)/10;
Params.NRM_dev=round(10*Params.NRM_dev)/10;
Params.gamma=round(10*Params.gamma)/10;
Params.CRM_R=round(10*Params.CRM_R)/10;


% pTRM checks
Params.check=round(10*Params.check)/10;
Params.dCK=round(10*Params.dCK)/10;
Params.DRAT=round(10*Params.DRAT)/10;
Params.maxDEV=round(10*Params.maxDEV)/10;

Params.CDRAT=round(10*Params.CDRAT)/10;
Params.CDRAT_prime=round(10*Params.CDRAT_prime)/10;
Params.DRATS=round(10*Params.DRATS)/10;
Params.DRATS_prime=round(10*Params.DRATS_prime)/10;
Params.mean_DRAT=round(10*Params.mean_DRAT)/10;
Params.mean_DRAT_prime=round(10*Params.mean_DRAT_prime)/10;
Params.mean_DEV=round(10*Params.mean_DEV)/10;
Params.mean_DEV_prime=round(10*Params.mean_DEV_prime)/10;
Params.dpal=round(10*Params.dpal)/10;


% tail checks
Params.DRATtail=round(10*Params.DRATtail)/10;
Params.dTR=round(10*Params.dTR)/10;
Params.MDvds=round(10*Params.MDvds)/10;
Params.dt_star=round(10*Params.dt_star)/10;


% Additivity checks
Params.dAC=round(10*Params.dAC)/10;


% Anis stats
Params.Anis_c=round(1000*Params.Anis_c)/1000;


% NLT stats


end


%% Required functions


function [R2]=PearsonCorr2(X, Y)
%
% function to determine the Pearson linear correlation between two input
% vector, X and Y

Xd=detrend(X, 0); % (x-xbar)
Yd=detrend(Y, 0); % (y-ybar)

R2 = sum((Xd.*Yd))^2 ./ ( sum(Xd.^2).*sum(Yd.^2) );

end

function [InRadians]=deg2rad(InDegrees)
%
% Convert an angle in degrees to radians

InRadians = (pi/180) .* InDegrees;

end

function [InDegrees]=rad2deg(InRadians)
%
% Convert an angle in radians to degree

InDegrees = (180/pi) .* InRadians;

end

function [A]=Anis_mat(s)
%
% Build the anisotroParams.y_prime tensor

A(1,1)=s(1);
A(2,2)=s(2);
A(3,3)=s(3);
A(2,1)=s(4);
A(3,2)=s(5);
A(3,1)=s(6);
A(1,2)=A(2,1);
A(1,3)=A(3,1);
A(2,3)=A(3,2);

end % AnisotroParams.y_prime matrix

function [theta]=calc_angle(Dir1, Dir2)
%
% Function to calculate the angle between two paleomagnetic directions
%

[x, y, z]=dir2cart(Dir1(1), Dir1(2));
Dir1_cart=[x, y, z];
Dir1_cart=Dir1_cart./sqrt(sum(Dir1_cart.^2));

[x, y, z]=dir2cart(Dir2(1), Dir2(2));
Dir2_cart=[x, y, z];
Dir2_cart=Dir2_cart./sqrt(sum(Dir2_cart.^2));

theta = atan2(norm(cross(Dir1_cart,Dir2_cart)),dot(Dir1_cart, Dir2_cart));
theta=rad2deg(theta);


end % Angle between two directions

function return_val=common_slope(bhat, varX, varY, varXY, Ns)
%
% Minimization Function to determine the probability of a common slope
% Derived from Warton et al. (2006), Bivariate line-fitting methods for
% allomettry, Biol. Rev., 81, 259-291, doi: 10.17/S1464793106007007
%

Sr=(Ns-1)./(Ns-2) .* (varY - 2.*bhat.*varXY + (bhat.^2).*varX);
Sf=(Ns-1)./(Ns-2) .* (varY + 2.*bhat.*varXY + (bhat.^2).*varX);
Srf=(Ns-1)./(Ns-2) .* (varY - (bhat.^2).*varX);

return_val=sum(Ns.*( 1./(Sr) + 1./(Sf) ).*Srf.^2);

end % common slope minimization function

function [parameters]=AraiCurvature(x,y)
%
% Function for calculating the radius of the best fit circle to a set of 
% x-y coordinates.
% Paterson, G. A., (2011), A simple test for the presence of multidomain
% behaviour during paleointensity experiments, J. Geophys. Res., doi: 10.1029/2011JB008369
%
% parameters(1) = k
% parameters(2) = a
% parameters(3) = b
% parameters(4) = SSE (goodness of fit)

% Reshape vectors for suitable input
x=reshape(x, length(x), 1);
y=reshape(y, length(y), 1);

% Normalize vectors
x=x./max(x);
y=y./max(y);

% Provide the initial estimate
E1=TaubinSVD([x,y]);

% Determine the iterative solution
E2=LMA([x,y], E1);

estimates=[E2(3), E2(1), E2(2)];

% Define the function to be minimized and calculate the SSE
func=@(v) sum((sqrt((x-v(2)).^2+(y-v(3)).^2)-v(1)).^2);
SSE=func(estimates);

if (E2(1)<=mean(x) && E2(2)<=mean(y))
    k=-1/E2(3);
else
    k=1/E2(3);
end

parameters=[k; E2(1); E2(2); SSE];

end % Arai plot curvature

function Par = TaubinSVD(XY)

%--------------------------------------------------------------------------
%  
%     Algebraic circle fit by Taubin
%      G. Taubin, "Estimation Of Planar Curves, Surfaces And Nonplanar
%                  Space Curves Defined By Implicit Equations, With 
%                  Applications To Edge And Range Image Segmentation",
%      IEEE Trans. PAMI, Vol. 13, pages 1115-1138, (1991)
%
%     Input:  XY(n,2) is the array of coordinates of n points x(i)=XY(i,1), y(i)=XY(i,2)
%
%     Output: Par = [a b R] is the fitting circle:
%                           center (a,b) and radius R
%
%     Note: this is a version optimized for stability, not for speed
%
%--------------------------------------------------------------------------

centroid = mean(XY);   % the centroid of the data set

X = XY(:,1) - centroid(1);  %  centering data
Y = XY(:,2) - centroid(2);  %  centering data
Z = X.*X + Y.*Y;
Zmean = mean(Z);
Z0 = (Z-Zmean)/(2*sqrt(Zmean));
ZXY = [Z0 X Y];
[U,S,V]=svd(ZXY,0); %#ok<ASGLU>
A = V(:,3);
A(1) = A(1)/(2*sqrt(Zmean));
A = [A ; -Zmean*A(1)];
Par = [-(A(2:3))'/A(1)/2+centroid , sqrt(A(2)*A(2)+A(3)*A(3)-4*A(1)*A(4))/abs(A(1))/2];

end   %  TaubinSVD

function Par = LMA(XY,ParIni)

%--------------------------------------------------------------------------
%  
%     Geometric circle fit (minimizing orthogonal distances)
%     based on the Levenberg-Marquardt scheme in the 
%     "algebraic parameters" A,B,C,D  with constraint B*B+C*C-4*A*D=1
%        N. Chernov and C. Lesort, "Least squares fitting of circles",
%        J. Math. Imag. Vision, Vol. 23, 239-251 (2005)
%
%     Input:  XY(n,2) is the array of coordinates of n points x(i)=XY(i,1), y(i)=XY(i,2)
%             ParIni = [a b R] is the initial guess (supplied by user)
%
%     Output: Par = [a b R] is the fitting circle:
%                           center (a,b) and radius R
%
%--------------------------------------------------------------------------

factorUp=10;  factorDown=0.04;
lambda0=0.01;  epsilon=0.000001;  
IterMAX = 50;  AdjustMax = 20;
Xshift=0;  Yshift=0;  dX=1;  dY=0;

n = size(XY,1);      % number of data points

%     starting with the given initial guess

anew = ParIni(1) + Xshift;
bnew = ParIni(2) + Yshift;

Anew = 1/(2*ParIni(3));
aabb = anew*anew + bnew*bnew;
Fnew = (aabb - ParIni(3)*ParIni(3))*Anew;
Tnew = acos(-anew/sqrt(aabb));
if (bnew > 0) 
    Tnew = 2*pi - Tnew;
end
%        if 1+4*Anew*Fnew < epsilon
%            fprintf(1,' +++ violation:  %f\n',1+4*Anew*Fnew);
%        end

VarNew = VarCircle(XY,ParIni);

%     initializing lambda and iter

lambda = lambda0;  finish = 0;

for iter=1:IterMAX

    Aold = Anew;  Fold = Fnew;  Told = Tnew;  VarOld = VarNew;

    H = sqrt(1+4*Aold*Fold);
    aold = -H*cos(Told)/(Aold+Aold) - Xshift;
    bold = -H*sin(Told)/(Aold+Aold) - Yshift;
    Rold = 1/abs(Aold+Aold);
%    fprintf(1,'%2d  (%f, %f)  %f  %.8f\n',iter,aold,bold,Rold,sqrt(VarOld));

%           computing moments

    DD = 1 + 4*Aold*Fold;
    D = sqrt(DD);
    CT = cos(Told);
    ST = sin(Told);

    H11=0; H12=0; H13=0; H22=0; H23=0; H33=0; F1=0; F2=0; F3=0;
    
    for i=1:n
        Xi = XY(i,1) + Xshift;
        Yi = XY(i,2) + Yshift;
        Zi = Xi*Xi + Yi*Yi;
        Ui = Xi*CT + Yi*ST;
        Vi =-Xi*ST + Yi*CT;

        ADF = Aold*Zi + D*Ui + Fold;
        SQ = sqrt(4*Aold*ADF + 1);
        DEN = SQ + 1;
        Gi = 2*ADF/DEN;
        FACT = 2/DEN*(1 - Aold*Gi/SQ);
        DGDAi = FACT*(Zi + 2*Fold*Ui/D) - Gi*Gi/SQ;
        DGDFi = FACT*(2*Aold*Ui/D + 1);
        DGDTi = FACT*D*Vi;

        H11 = H11 + DGDAi*DGDAi;
        H12 = H12 + DGDAi*DGDFi;
        H13 = H13 + DGDAi*DGDTi;
        H22 = H22 + DGDFi*DGDFi;
        H23 = H23 + DGDFi*DGDTi;
        H33 = H33 + DGDTi*DGDTi;

        F1 = F1 + Gi*DGDAi;
        F2 = F2 + Gi*DGDFi;
        F3 = F3 + Gi*DGDTi;
    end

    for adjust=1:AdjustMax

%             Cholesly decomposition

        G11 = sqrt(H11 + lambda);
        G12 = H12/G11;
        G13 = H13/G11;
        G22 = sqrt(H22 + lambda - G12*G12);
        G23 = (H23 - G12*G13)/G22;
        G33 = sqrt(H33 + lambda - G13*G13 - G23*G23);

        D1 = F1/G11;
        D2 = (F2 - G12*D1)/G22;
        D3 = (F3 - G13*D1 - G23*D2)/G33;

        dT = D3/G33;
        dF = (D2 - G23*dT)/G22;
        dA = (D1 - G12*dF - G13*dT)/G11;

%            updating the parameters

        Anew = Aold - dA;
        Fnew = Fold - dF;
        Tnew = Told - dT;
%        fprintf(1,'%2d   %.8f\n',iter,lambda);    

        if (1+4*Anew*Fnew < epsilon && lambda>1)
%            fprintf(1,'     violation:  %f\n',1+4*Anew*Fnew);
            Xshift = Xshift + dX;
            Yshift = Yshift + dY;

            H = sqrt(1+4*Aold*Fold);
            aTemp = -H*cos(Told)/(Aold+Aold) + dX;
            bTemp = -H*sin(Told)/(Aold+Aold) + dY;
            rTemp = 1/abs(Aold+Aold);

            Anew = 1/(rTemp + rTemp);
            aabb = aTemp*aTemp + bTemp*bTemp;
            Fnew = (aabb - rTemp*rTemp)*Anew;
            Tnew = acos(-aTemp/sqrt(aabb));
            if bTemp > 0 
               Tnew = 2*pi - Tnew;
            end
            VarNew = VarOld;
            break;
        end

        if 1+4*Anew*Fnew < epsilon
           lambda = lambda * factorUp;
           continue;
        end

        DD = 1 + 4*Anew*Fnew;
        D = sqrt(DD);
        CT = cos(Tnew);
        ST = sin(Tnew);

        GG = 0;

        for i=1:n
            Xi = XY(i,1) + Xshift;
            Yi = XY(i,2) + Yshift;
            Zi = Xi*Xi + Yi*Yi;
            Ui = Xi*CT + Yi*ST;

            ADF = Anew*Zi + D*Ui + Fnew;
            SQ = sqrt(4*Anew*ADF + 1);
            DEN = SQ + 1;
            Gi = 2*ADF/DEN;
            GG = GG + Gi*Gi;
        end
    
        VarNew = GG/(n-3);

        H = sqrt(1+4*Anew*Fnew);
        anew = -H*cos(Tnew)/(Anew+Anew) - Xshift;
        bnew = -H*sin(Tnew)/(Anew+Anew) - Yshift;
        Rnew = 1/abs(Anew+Anew);
       
%             checking if improvement is gained

        if VarNew <= VarOld      %   yes, improvement
           progress = (abs(anew-aold) + abs(bnew-bold) + abs(Rnew-Rold))/(Rnew+Rold);
           if progress < epsilon
              Aold = Anew;
              Fold = Fnew;
              Told = Tnew;
              VarOld = VarNew; %#ok<NASGU>
              finish = 1;
              break;
           end
           lambda = lambda * factorDown;
           break;
        else                     %   no improvement
           lambda = lambda * factorUp;
           continue;
        end
    end
    if finish == 1
       break;
    end
end

H = sqrt(1+4*Aold*Fold);
Par(1) = -H*cos(Told)/(Aold+Aold) - Xshift;
Par(2) = -H*sin(Told)/(Aold+Aold) - Yshift;
Par(3) = 1/abs(Aold+Aold);

end  % LMA

function Var = VarCircle(XY,Par)

%--------------------------------------------------------------------------
%  
%             computing the sample variance of distances 
%             from data points (XY) to the circle Par = [a b R]
%
%--------------------------------------------------------------------------
 

n = size(XY,1);      % number of data points
Dx = XY(:,1) - Par(1);  Dy = XY(:,2) - Par(2);
D = sqrt(Dx.*Dx + Dy.*Dy) - Par(3);

Var = D'*D/(n-3);

end  %  VarCircle

function [vector]=Rotate_Vec(vec, angle, rotaxis) %#ok<DEFNU>
%
% Function to rotate an input 3D vector (vec) by a desired angle (angle)
% around a rotation axis given by rotaxis
%

% Data checking
if nargin < 3
    error('Rotate_Vector:TooFewInputs', 'At least 3 inputs are required.');
end

s=size(rotaxis);
if s(1)~=length(angle)
    error('Rotate_Vector:VectorLength','Rotaxis and angle vectors must be the same length.');
end

% Convert rotation axis into a unit vector - in case it is not
vec_len=sqrt(sum(rotaxis.^2, 2));
for i=1:1:3
    rotaxis(:,i)=rotaxis(:,i)./vec_len;
end

vector(length(angle), 3)=zeros;
cAng=cos(angle);

% Expand vec to the size of rotaxis if vec is 1x3
s=size(vec);
if s(1)==1
    vec=vec(ones(length(angle),1), :);
end


for j=1:1:length(angle)
    if angle(j)==0
        vector(j,:)=vec(j,:);
    else
        rotaxis(j,:)=rotaxis(j,:)./norm(rotaxis(j,:));
        vector(j,:)=vec(j,:).*cAng(j) + dot(vec(j,:),rotaxis(j,:))*(1-cAng(j)).*rotaxis(j,:) + cross(rotaxis(j,:),vec(j,:)).*sin(angle(j));
    end
end

end % 3-D vector rotation

function [x, y, z]=dir2cart(Dec, Inc, Mag)
%
% Converts a paleomagnetic direction to Cartesian coordinates
%

if nargin <3
    Mag=1;
end

Dec=deg2rad(Dec);
Inc=deg2rad(Inc);

x=cos(Dec).*cos(Inc).*Mag;
y=sin(Dec).*cos(Inc).*Mag;
z=sin(Inc).*Mag;

end % convert a direction to a Cartesian vector

function [Dec, Inc, R]=cart2dir(x, y, z)
%
% Convert Cartesian coordinates to declination/inclination
% 

if nargin <3
    error('cart2dir:Error', 'Needs 3 input coordinates')
end


R=sqrt((x).^2+(y).^2+(z).^2);

if (R==0)
    warning('cart2dir:Error', 'R==0')

else
    
    
    Dec=rad2deg(atan2(y,x) );
    Inc=rad2deg(asin(z./R) );
    
    Dec(Dec<0)=Dec(Dec<0)+360;
    Dec(Dec>360)=Dec(Dec>360)-360;
end

Dec(R==0)=NaN;
Inc(R==0)=NaN;

end % convert a Cartesian vector to a direction

function [Dgeo, Igeo]=dirot(Dec, Inc, Az, Pl)

% Converts a direction to geographic coordinates using az,pl as azimuth and
% plunge (inclination) of Z direction
% Based on 2G software which uses the strike for the calculations 
% Strike = Az +90


[x,y,z]=dir2cart(Dec, Inc, 1);

str_rad=deg2rad(Az+90); % Here add 90 to get the strike
pl_rad=deg2rad(Pl);

xp=(x.*sin(pl_rad) + z.*cos(pl_rad)).*sin(str_rad) + y.*cos(str_rad);
yp=-(x.*sin(pl_rad) + z.*cos(pl_rad)).*cos(str_rad) + y.*sin(str_rad);
zp=-x.*cos(pl_rad) + z.*sin(pl_rad);


[Dgeo, Igeo]=cart2dir(xp, yp, zp);


end % rotate data from core to geographic coords

function [Mdec, Minc, k, a95, R]=FisherMeanDir(Dec, Inc, RFlag)


if nargin < 3
    RFlag=0;
end

if RFlag==1 % Flip Reverse directions
    Inc(Dec<270 & Dec >90)=-Inc(Dec<270 & Dec >90);
    Dec(Dec<270 & Dec >90)=Dec(Dec<270 & Dec >90)+180;
    Dec(Dec<0)=Dec(Dec<0)+360;
    Dec(Dec>360)=Dec(Dec>360)-360;
end
% Dec(Dec> 90 & Deac <270)=Dec(Dec> 90 & Deac <270)+180

[X, Y, Z]=dir2cart(Dec, Inc);
N=length(X);

% xbar=mean(X);
% ybar=mean(Y);
% zbar=mean(Z);


xsum=sum(X);
ysum=sum(Y);
zsum=sum(Z);

R2=xsum.^2+ysum.^2+zsum.^2;
R=sqrt(R2);

x=xsum./R;
y=ysum./R;
z=zsum./R;

Mdec=rad2deg(atan2(y,x));
Minc=rad2deg(asin(z));

Mdec(Mdec<0)=Mdec(Mdec<0)+360;
Mdec(Mdec>360)=Mdec(Mdec>360)-360;

k=(N-1)/(N-R);

pwr=1/(N-1);
bracket=((1/0.05).^pwr)-1;
OB=(N-R)/R;
Calpha=1-(OB*bracket);
alpha=acos(Calpha);
a95=rad2deg(alpha);

end % get the Fisher mean

function [Mdec, Minc, MAD]=PmagPCA(Dec, Inc, Int, type)
%
% Function to determine the best-fit line through demagnetization data
% Based on PCA methods of Kirschvink (1980)
% Note: This only fits lines, not planes
%

% Variable checking
if nargin < 4
    type='free'; % Default fitting method
end

% Setup input
[X, Y, Z]=dir2cart(Dec, Inc, Int);
input=[X,Y,Z];

if strcmpi(type, 'free') % Free fit
    bars=[mean(X), mean(Y), mean(Z)];
elseif strcmpi(type, 'anc') || strcmpi(type, 'anchored') % Anchored fit
    bars=[0, 0, 0];
elseif strcmpi(type, 'origin')  % Free fit with the origin
    input=[input; 0, 0, 0];
    bars=[mean(X), mean(Y), mean(Z)];
else
    error('PCADir:Fit_Type', 'Fit type unrecognised');
end

n=size(input,1);
x=zeros(n, 3);
for j=1:1:n
    x(j,:)=input(j,:)-bars;
end

% Perform PCA
% based on PmagParams.y_prime

T=Tmatrix(x);
[V,tau]=eig(T);
tau=diag(tau); % take the diagonal to get the three taus
max_tau_ind= tau==max(tau); % find the index of the maximum tau
v1=V(:, max_tau_ind)';

tau=tau./sum(tau);% rescale tau to sum-to-one
tau=sort(tau, 'descend'); % sort into descending order so that tau(1) is the max


% reference vector for defining the direction of principal component
P(1,:)=input(1,:);
P(2,:)=input(end,:);

reference=P(1,:)-P(2,:);
ref_dot=sum(v1.*reference);

ref_dot(ref_dot<-1)=-1;
ref_dot(ref_dot>1)=1;

if acos(ref_dot) > pi/2
    v1=-v1;
end


if tau(2)+tau(3) <=eps % MAD is too small to calculate - set to zero
    MAD=0;
else
    MAD=atand(  sqrt( tau(2)+tau(3) ) / sqrt(tau(1))  );
end

[Mdec, Minc]=cart2dir(v1(1), v1(2), v1(3));


end % PCA for a paleomag direction - lines only, not planes

function T=Tmatrix(X)
%
% function to create the orientation matrix for Pmag data
%

T(1,1)=sum(X(:,1).*X(:,1));
T(2,1)=sum(X(:,2).*X(:,1));
T(3,1)=sum(X(:,3).*X(:,1));

T(1,2)=sum(X(:,1).*X(:,2));
T(2,2)=sum(X(:,2).*X(:,2));
T(3,2)=sum(X(:,3).*X(:,2));

T(1,3)=sum(X(:,1).*X(:,3));
T(2,3)=sum(X(:,2).*X(:,3));
T(3,3)=sum(X(:,3).*X(:,3));

end % build the orientation matrix for PCA


