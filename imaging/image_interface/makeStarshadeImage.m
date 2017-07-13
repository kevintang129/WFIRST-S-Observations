function [ efDefectImg, lambdaIn, vecPetalArray ] = makeStarshadeImage( opt_in )
% makeStarshadeImage
% A sample program which creates a locus of edge point from an apodization
% function and propagates the resulting field to a telescope focal plane.
% History:
% 4/25/17: first complete version, Eric Cady (JPL/Caltech)
% 4/29/17: modification of the interface, pupil handling, ..., Sergi Hildebrandt (JPL/Caltech)
%
% TBD:
% 1) Include RGB colors

% Check some options are being passed
  if ~exist( 'opt_in', 'var' )
  disp( '(makeStarshadeImage) Provide some options: opt.x=y ; makeStarshadeImage( opt ) ; returning ...' )
  return
  end

% Get default options (look inside the function for specific definitions)
opt = get_default_options( opt_in ) ;

% Developer version? Here is just this
  if ( opt.developer ), units_image ; else, units ; end

% Main parameters
Nx = opt.Nx_pupil_pix ;
dlt_lmbd = opt.delta_lambda_nm ;
r_src = opt.r_source_mas ;
psi_src = opt.psi_source_deg ;
ppl_fl = opt.pupil_file ;

% Settings for saving fields
useSave = opt.save ; % 1 save, 0 don't
savePath = opt.save_path ;
  if ~isdir( savePath ), system( [ 'mkdir -p ' savePath ] ) ; end
% Skipping the simulation if it is saved and does not need to be re-done
  if ( ~opt.redo ) && ( exist( [ savePath opt.save_filename '.mat' ] ) == 2 )
  disp( sprintf( '(makeStarshadeImage) Simulation %s exists. Skipping.', opt.save_filename ) )
  load( [ savePath '/' opt.save_filename '.mat' ] )
  return
  end
  if ( opt.redo ) && ( exist( [ savePath '/' opt.save_filename '.mat' ] ) == 2 )
  disp( sprintf( '(makeStarshadeImage) Simulation %s exists, but re-doing it.', opt.save_filename ) )
  end

%---------------------------
% Step 1: Load up starshade
%---------------------------
  if ( opt.developer )
  units_image
  else
  units
  end

% Load occulter definition (change this for your file structure)
folderpath = './in/' ; %'../Occulter design/Profiles/';
occulterName = [ folderpath opt.occulter_name '.mat' ] ;
load( occulterName ); % Load up the comparison occulter
disp( sprintf( '(makeStarshadeImage) Occulter name: %s', occulterName ) )

% Simulation parameters (lambdaRange comes from the occulter file)
  if ~exist( 'lambdaRange', 'var' )
  disp( sprintf( '(makeStarshadeImage) lambdaRange variable not present in the occulter file %s. Returning.', occulterName ) )
  dbstop if error
  make_an_error
  end
% More than enough
lambdaIn = lambdaRange( 1 ) + dlt_lmbd * nm * ( 0 : 1 : 100 / dlt_lmbd ) ;
% Those who are acceptable within the range of the occulter
q = find( lambdaIn <= lambdaRange( end ) ) ;
lambdaIn = lambdaIn( q ) ;
n_lmbd = numel( lambdaIn ) ;
disp( sprintf( 'Considering %i wavelengths', numel( lambdaIn ) ) ) 
% Build/load
if strcmp( ppl_fl, '0' )
    secondarySize = 0; % If you're fine with a generic circular pupil with a secondary (but no struts), set this to the fraction of the radius the secondard covers (e.g. 0.2)
    if ( opt.developer )
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
  if ( opt.developer )
  vecPetalArray = createErrorProfileV1_image(r, Profile, occulterDiameter, petalLength, numPetals, {});
  else
  vecPetalArray = createErrorProfileV1(r, Profile, occulterDiameter, petalLength, numPetals, {});
  end
t = toc ;
  if ( opt.developer )
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
  if ( opt.developer )
  UTotL = bdwf_image(xVals, yVals, zVals, Z, lambdaIn, telescopeDiameter/Nx, Nx, r_src * mas, psi_src * degree, deltaX, deltaY);
  else
  UTotL = bdwf(xVals, yVals, zVals, Z, lambdaIn, telescopeDiameter/Nx, Nx, r_src * mas, psi_src * degree, deltaX, deltaY);
  end
t = toc;
  if ( opt.developer )
  bdwf_lbl = 'bdwf_image' ;
  else
  bdwf_lbl = 'bdwf' ;
  end

disp( sprintf('%s took %3.2f seconds, %3.2f per wavelength bin', bdwf_lbl, t, t / n_lmbd ) )
  if useSave == 1
  % These fields do not belong to this computation and may vary later on
    if isfield( opt, 'Nx_image_pix' ), opt = rmfield( opt, 'Nx_image_pix' ) ; end
    if isfield( opt, 'Nx_img' ), opt = rmfield( opt, 'Nx_img' ) ; end
    if isfield( opt, 'diam_img_mas' ), opt = rmfield( opt, 'diam_img_mas' ) ; end
  save([savePath '/' opt.save_filename '.mat' ], 'lambdaIn', 'opt', 'pupil', 'telescopeDiameter', 'UTotL' )
  end

dbstop if error
%make_an_error

