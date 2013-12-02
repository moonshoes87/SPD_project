function [Params] = GetPintParams_v5b(Mvec, Temps, Treatment, start_pt, end_pt, T_orient, Blab, NRM_rot_flag, Az, Pl, ChRM, A_corr, s_tensor, NLT_corr, NLT_hat)
%
%% Input
  % n = 16
  Mvec = [1, 2, 3; 1, 2, 3; 1, 2, 3; 1, 2, 3; 1, 2, 3; 1, 2, 3; 1, 2, 3; 1, 2, 3]
  temps = [10, 10, 20, 20, 30, 30, 40, 40, 50, 50, 60, 60, 70, 70, 80, 80]
  Treatment = [0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3] % possibly should exclude 5s
  start_pt = 20
  end_pt = 80
  T_orient = [x, y, z]  % unit vector for a known Blab orientation.  not sure what numbers are appropriate
  Blab = 0.00004  %scalar value.  strength in uT
  NRM_rot_flat = 0 % 0 for no rotation, 1 for apply rotation
  Az = 30 % sample azimuth in degrees (no idea)
  Pl = 50 % the sample plunge in degrees
  ChRM = [25, 87] %  - a (1 x 2) or a (1 x 3) vector describing an independent measure of the direction of Banc (can be a known direction). If ChRM is a (1 x 2) vector it is assumed to be a Dec/Inc direction and converted to a Cartesian vector.
  A_corr = 0 % - flag for anisotropy correction (0=no correction, 1=apply correction)
  s_tensor = [0, 0, 0, 0, 0, 0] %(no clue) a (1 x 6) vector with the six unique elements that describe an anisotropy tensor
  NLT_corr = 0 %      - flag for non-linear TRM correction (0=no correction, 1=apply correction)
  NLT_hat = [1, 2] % (no clue)       - a (1 x 2) vector containing the coefficients of the best-fit hyperbolic tanget funtion

  
  
  % Mvec          - a (n x 3) matrix of magnetization values when n is the total number of steps, which includes NRM, TRM, and all check measurements (i.e., the raw data measurements)
  % Temps         - a (n x 1) vector of the temperature of each treatment
  % Treatment     - a (n x 1) vector of integer values describing the treatment type at each step, this follows the convetion of the ThellierTool
  %                 (0=NRM demag, 1=TRM remag, 2=pTRM check, 3=pTRM tail check, 4=additivity check, 5=inverse TRM step). If any Treatment is set to '5' all data are treated as a Thellier experiment
  % start_pt      - the temperature of the start point for the Arai plot best-fit line - ONLY THE INDEX IS CODED, NOT THE TEMPERATURE. THIS IS BETTER FOR TESTING
  % end_pt        - the temperature of the end point for the Arai plot best-fit line - ONLY THE INDEX IS CODED, NOT THE TEMPERATURE. THIS IS BETTER FOR TESTING
  % T_orient      - a (1 x 3) unit vector containing the x, y, z values for a know Blab orientation.
% Blab          - the strength of the laboratry field in muT
  % NRM_rot_flag  - flag for directional rotation from core coords to stratigraphic coords (0 or []=no rotation, 1=apply rotation)
% Az            - the sample azimuth in degrees
% Pl            - the sample plunge in degrees
    % ChRM          - a (1 x 2) or a (1 x 3) vector describing an independent measure of the direction of Banc (can be a known direction).
    %                 If ChRM is a (1 x 2) vector it is assumed to be a Dec/Inc direction and converted to a Cartesian vector.
%                 If rotation is applied ChRM is assumed to already be in the final coordinate system
  % A_corr        - flag for anisotropy correction (0=no correction, 1=apply correction)
  % s_tensor      - a (1 x 6) vector with the six unique elements that describe an anisotropy tensor
  % NLT_corr      - flag for non-linear TRM correction (0=no correction, 1=apply correction)
  % NLT_hat       - a (1 x 2) vector containing the coefficients of the best-fit hyperbolic tanget funtion
