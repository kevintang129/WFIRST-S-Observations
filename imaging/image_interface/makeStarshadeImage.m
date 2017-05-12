function [ efDefectImg, lambdaIn, vecPetalArray ] = makeStarshadeImage( Nx, delta_lambda, r_planet, psi_planet, opt, pupil_file )
% makeStarshadeImage
% Eric Cady, 4/25/17: first complete version
% Sergi Hildebrandt, 4/29/17: modification of the interface, pupil handling, ...
%
% Include RGB colors
% A sample program which creates a locus of edge point from an apodization
% function and propagates the resulting field to a telescope focal plane.

% Size of the pupil data in pixels (square)
  if ~exist( 'Nx', 'var' )
  Nx = 128 ;
  end

% Step of wavelength to consider
  if ~exist( 'delta_lambda', 'var' )
  dlt_lmbd = 50 ; % nm
  else
  dlt_lmbd = delta_lambda ;
  end

% Separation of the planet from the center of the pointing
% This is temporary, meant for tests only
  if ~exist( 'r_planet' )
  r_plnt = 0 ; % mas
  else
  r_plnt = floor( r_planet ) ;
  end

% Angle of the planet with respect the horizontal axis
  if ~exist( 'psi_planet', 'var' )
  psi_plnt = 0 ; % degrees
  else
  psi_plnt = floor( psi_planet ) ;
  end

% Replace by a file in FITS format with an Nx x Nx array if you want a specific pupil  
  if ~exist( 'pupil_file', 'var' )
  ppl_fl = './in/pupil_D1Kpix_2048.fits' ;
  else
  ppl_fl = pupil_file ;
  end

  if strcmp( ppl_fl, '0' )
  ppl_fl = '' ;
  end

% Settings for saving fields
useSave = opt.save_all ; % 1 save, 0 don't
  if isfield( opt, 'save_path' )
  savePath = opt.save_path ;
  else
  savePath = './out/';
  end

  if ~isdir( savePath ), system( [ 'mkdir -p ' savePath ] ) ; end
saveFilename = sprintf( 'starshade_out_Nx_%i_pix_dl_%inm_dr_%i_mas_psi_%i_deg', Nx, dlt_lmbd, r_plnt, psi_plnt ) ;
  if numel( ppl_fl ) == 0
  saveFilename = sprintf( '%s_ideal', saveFilename ) ;
  end

% Directory where to store the figures
useSaveImg = 1; % 1 save, 0 don't
savePathImg = './fig/';
  if ~isdir( savePathImg ), system( [ 'mkdir -p ' savePathImg ] ) ; end


%---------------------------
% Step 1: Load up starshade
%---------------------------
  if isfield( opt, 'developer' )
  units_image
  else
  units
  end

% Load occulter definition (change this for your file structure)
occulterName = 'NW2';
folderpath = './in/' ; %'../Occulter design/Profiles/';
load([folderpath occulterName '.mat']); % Load up the comparison occulter

% Simulation parameters
lambdaIn = 300 * nm + dlt_lmbd * nm * ( 0 : 1 : 700 / dlt_lmbd ) ;
n_lmbd = numel( lambdaIn ) ;
disp( sprintf( 'Considering %i wavelengths', numel( lambdaIn ) ) ) 
imagePlaneDiamInMAS = 2000; % Diameter of image plane in milliarcseconds
Nxi = 400; % Number of pixels across image plane 
% Build/load
if isempty( ppl_fl )
    secondarySize = 0; % If you're fine with a generic circular pupil with a secondary (but no struts), set this to the fraction of the radius the secondard covers (e.g. 0.2)
    if isfield( opt, 'developer' )
    pupil = makePupil_image( Nx, Nx, 1, secondarySize, 0, 0 ) ;
    else
    pupil = makePupil(Nx, Nx, 1, secondarySize, 0, 0);
    end
else
    pupil = fitsread( ppl_fl );
    n_ppl = sqrt( numel( pupil ) ) ;
    % Check
    disp( sprintf( 'Reference pupil size is %ix%i pixels.', n_ppl, n_ppl ) ) 
    % Reducing the size of the pupil grid (fast<1s)
      if Nx ~= n_ppl
      % finding where the pupil starts
      q_ppl_1 = min( find( pupil ~= 0 ) ) ;
      q_ppl_2 = max( find( pupil ~= 0 ) ) ;
      % Column where the pupil data starts and ends
      clmn_1 = floor( q_ppl_1 / n_ppl ) ;
      clmn_2 = floor( q_ppl_2 / n_ppl ) ;
      % Actual resize (method bicubic, some low pass filter is applied first with 11 pix box over the 2048, fine)
      %% Asssuming symmertic pupil when cutting it out
      pupil = resizem( pupil( clmn_1 : clmn_2, clmn_1 : clmn_2 ), Nx / ( clmn_2 - clmn_1 + 1 ), 'bicubic' ) ; 
      % Sharping the edges
      pupil( find( pupil < .7 ) ) = 0 ;
      pupil( find( pupil ~= 0 ) ) = 1 ;
      disp( sprintf( 'New pupil size is reduced to %ix%i pixels.', Nx, Nx ) )
      end
% Exception for checks
%pupil = fitsread( 'pupil_D1Kpix_64.fits' );
end
dbstop if error
%make_an_error
% Lateral offsets in meters
deltaX = 0;
deltaY = 0;

% -------------------------------
% Step 2: Build edge
% -------------------------------

% Build edge; this function is overkill for an unaberrated edge but will do
% the job
tic
  if isfield( opt, 'developer' )
  vecPetalArray = createErrorProfileV1_image(r, Profile, occulterDiameter, petalLength, numPetals, {});
  else
  vecPetalArray = createErrorProfileV1(r, Profile, occulterDiameter, petalLength, numPetals, {});
  end
t = toc ;
  if isfield( opt, 'developer' )
  createErr_lbl = 'createErrorProfileV1_image' ;
  else
  createErr_lbl = 'createErrorProfileV1' ;
  end
disp( sprintf( '%s took %3.2f seconds', createErr_lbl, t ) )


tic
tmpxVals = [];
tmpyVals = [];
tmpzVals = [];%**
for j = 1:numPetals
    tmpxVals = [tmpxVals vecPetalArray{j}{1}(1, :)]; %#ok<AGROW>
    tmpyVals = [tmpyVals vecPetalArray{j}{1}(2, :)]; %#ok<AGROW>
    tmpzVals = [tmpzVals vecPetalArray{j}{1}(3, :)]; %#ok<AGROW> %**
end
xVals = [tmpxVals tmpxVals(1)];
yVals = [tmpyVals tmpyVals(1)];
zVals = [tmpzVals tmpzVals(1)];
% xVals, yVals, zVals give the 3D edge locus
% Under most circumstances zVals will be all zeros
t = toc ;
disp( sprintf( 'xyzVals took %3.2f seconds', t ) )

%--------------------------
% Step 3: Compute field at 
%--------------------------

% Propagate to telescope aperture plane with line integral
tic
  if isfield( opt, 'developer' )
  UTotL = bdwf_image(xVals, yVals, zVals, Z, lambdaIn, telescopeDiameter/Nx, Nx, r_plnt * mas, psi_plnt * degree, deltaX, deltaY);
  else
  UTotL = bdwf(xVals, yVals, zVals, Z, lambdaIn, telescopeDiameter/Nx, Nx, r_plnt * mas, psi_plnt * degree, deltaX, deltaY);
  end
t = toc;
  if isfield( opt, 'developer' )
  bdwf_lbl = 'bdwf_image' ;
  else
  bdwf_lbl = 'bdwf' ;
  end

disp( sprintf('%s took %3.2f seconds, %3.2f per wavelength bin', bdwf_lbl, t, t / n_lmbd ) )

tic
  if isfield( opt, 'developer' )
  [ ~, ~, impeak ] = mft_image( 1, Nxi, pupil ) ; % Use 1 L/D square
  else
  [ ~, ~, impeak ] = mft( 1, Nxi, pupil ) ; % Use 1 L/D square
  end
peak = max(abs(impeak(:)));

masPerPixel = imagePlaneDiamInMAS/Nxi;
efDefectImg = zeros(Nxi, Nxi, length(lambdaIn));
for ll = 1:length(lambdaIn)
    UTot = UTotL(:,:,ll);
    lambda = lambdaIn(ll);
    disp(['Wavelength: ' num2str(lambda*1e9) 'nm'])
    
    D = telescopeDiameter;
    imagePlaneDiameterInMAS = Nxi*masPerPixel;
    imagePlaneDiameterInLambdaOverD = imagePlaneDiameterInMAS*mas*D/lambda;
   
      if isfield( opt, 'developer' ) 
      [Xout, Yout, imagePlane] = mft_image(imagePlaneDiameterInLambdaOverD, Nxi, UTot.*pupil);
      else
      [Xout, Yout, imagePlane] = mft(imagePlaneDiameterInLambdaOverD, Nxi, UTot.*pupil);
      end
    imagePlane = imagePlane/peak; % Normalize to unocculted on-axis peak = 1
    
    efDefectImg(:,:,ll) = imagePlane;
end
t = toc ;
disp( sprintf( 'efDefectImg took %3.2f seconds', t ) )

if useSave == 1
    save([savePath saveFilename '.mat' ], 'efDefectImg', 'lambdaIn', 'vecPetalArray')
end

dbstop if error
%make_an_error

