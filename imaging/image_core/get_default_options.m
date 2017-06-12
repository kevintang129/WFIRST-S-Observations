function opt = get_default_options( opt )
% Function to define the default options for Starshade imaging
% History:
% 05/12/17: created, Sergi Hildebrandt (JPL/Caltech)


% Developer version or Eric Cady's original version
  if ~isfield( opt, 'developer' )
  opt.developer = 0 ;
  end

% Switch to control some work that needs to be re-done or may be skipped if it already exists
  if ~isfield( opt, 'redo' )
  opt.redo = 0 ;
  end

% Occulter Name
  if ~isfield( opt, 'occulter_name' )
  opt.occulter_name = 'NI2' ;
  end

% Size of the pupil data in pixels (square)
  if ~isfield( opt, 'Nx_pupil_pix' )
  opt.Nx_pupil_pix = 64 ;
  end

% Number of pixels across image plane
  if ~isfield( opt, 'Nx_image_pix' )
  opt.Nx_img = 400 ;
  else
  opt.Nx_img = opt.Nx_image_pix ;
  end

% Diameter of image plane in milliarcseconds
  if ~isfield( opt, 'diam_image_mas' )
  opt.diam_img_mas = 2000 ;
  else
  opt,diam_img_mas = opt.diam_image_mas ;
  end

% Step of wavelength to consider
  if ~isfield( opt, 'delta_lambda_nm' )
  opt.delta_lambda_nm = 50 ; % nm
  end

% X/Y positions on the focal plane (instead of r_source_mas and psi_source_deg)
% Matlab: column major. Therefore in an image, X/Y are to be understood as -y/x in usual Cartesain convention (e.g., X=1, Y=0 means (0,-1) in sual x/y axis)
% Convention ox X/Y with respect to r/Psi as in bdwf core function
% s1 = sin(psi1);
% c1 = cos(psi1);
% s2 = sin(psi2);
% c2 = cos(psi2);


% Checks of consistency
  if isfield( opt, 'r_source_mas' ) && ~isfield( opt, 'psi_source_deg' )
  disp( '(get_default_options) r_source_mas set, but psi_source_deg not. Inconsistent. Returning.' )
  return
  end
  if isfield( opt, 'psi_source_deg' ) && ~isfield( opt, 'r_source_mas' )
  disp( '(get_default_options) psi_source_deg set, but r_source_mas not. Inconsistent. Returning.' )
  return
  end

% Separation of the source from the center of the pointing
  if ~isfield( opt, 'r_source_mas' )
  opt.r_source_mas = 0 ; % mas
  end

% Angle of the source with respect the horizontal axis
  if ~isfield( opt, 'psi_source_deg' )
  opt.psi_source_deg = 0 ; % degrees
  end

% Consistency check
  if isfield( opt, 'y_source_mas' ) && ~isfield( opt, 'x_source_mas' )
  disp( '(get_default_options) If opt.y_source_mas is set, opt.x_source_mas must also be set. Returning.' )
  return
  end
  if isfield( opt, 'x_source_mas' ) && ~isfield( opt, 'y_source_mas' )
  disp( '(get_default_options) If opt.x_source_mas is set, opt.y_source_mas must also be set. Returning.' )
  return
  end

  if isfield( opt, 'x_source_mas' )
  % Consistency check
    if ~isfield( opt, 'y_source_mas' )
    disp( '(get_default_options) If opt.x_source_mas is set, opt.y_source_mas must also be set. Returning.' ) 
    return
    end
  % Transform into r and psi
  %% Check for consistency if r_source_mas and psi_source_deg exist:
    if isfield( opt, 'r_source_mas' ) && isfield( opt, 'psi_source_deg' ) && ( sqrt( opt.x_source_mas^2 + opt.y_source_mas^2 ) )
      if ( ( sqrt( opt.x_source_mas^2 + opt.y_source_mas^2 ) ~= opt.r_source_mas ) || ( atan2( opt.y_source_mas, opt.x_source_mas ) * 180 / pi ~= opt.psi_source_deg ) ) 
      r_source_mas_tmp = sqrt( opt.x_source_mas^2 + opt.y_source_mas^2 ) ;
      psi_source_deg_tmp = atan2( opt.y_source_mas, opt.x_source_mas ) * 180 / pi ; % deg
      disp( sprintf( '(get_default_options) Changing the values of r_source_mas and psi_source_deg from %3.3f, %3.3f to %3.3f, %3.3f', opt.r_source_mas, opt.psi_source_deg, r_source_mas_tmp, psi_source_deg_tmp ) )
      opt.r_source_mas = r_source_mas_tmp ;
      opt.psi_source_deg = psi_source_deg_tmp ;
      end
    end
  end
  
  %% If they are not fields, create them
    if ~isfield( opt, 'r_source_mas' ) 
    opt.r_source_mas = sqrt( opt.x_source_mas^2 + opt.y_source_mas^2 ) ;
    end
    %%% Recall matlab column major convention and atan2(y,x)
    if ~isfield( opt, 'psi_source_deg' )
      if opt.r_source_mas % if the radius is zero, the angle does not matter.
      opt.psi_source_deg = atan2( opt.y_source_mas, opt.x_source_mas ) * 180 / pi ; % deg
      else
      opt.psi_source_deg = 0 ;
      end
    end

% For ease identification of the pixel on the image
opt.x_source_mas = opt.r_source_mas * cos( opt.psi_source_deg * pi / 180 ) ;
opt.y_source_mas = opt.r_source_mas * sin( opt.psi_source_deg * pi / 180 ) ;

% Replace by a file in FITS format with an Nx x Nx array if you want a specific pupil
  if ~isfield( opt, 'pupil_file' )
  opt.pupil_file = './in/pupil_D1Kpix_2048.fits' ;
  end

% Saving all the output results and images (0=No, 1=Yes)
  if ~isfield( opt, 'save_all' )
  opt.save_all = 0 ;
  end

% Saving only output results
  if ~isfield( opt, 'save' )
  opt.save = 0 ;
  end
% Saving figures
  if ~isfield( opt, 'save_fig' )
  opt.save_fig = 0 ;
  end

% If all saved, then:
  if opt.save_all
  opt.save = 1 ;
  opt.save_fig = 1 ;
  end

% paths to save the results
  if ~isfield( opt, 'save_path' )
  opt.save_path = './out' ;
    if ( opt.developer )
    opt.save_path = './out_dev/' ;
    end
  end

% Saving the figures
  if ~isfield( opt, 'save_path_fig' )
  opt.save_path_fig = './fig' ;
    if ( opt.developer )
    opt.save_path_fig = './fig_dev/' ;
    end
  end

% For the interpolation analysis
  if ~isfield( opt, 'star' )
  opt.star = 0 ;
  end

  if ~isfield( opt, 'planet' )
  opt.planet = 0 ;
  end

  if ~isfield( opt, 'super_resolution' )
  opt.super_resolution.res = 1 ;
  opt.super_resolution.interp_method = 'linear' ;
  end

  if ~isfield( opt.super_resolution, 'res' )
  opt.super_resolution.res = 1 ;
  end

  if ~isfield( opt.super_resolution, 'interp_method' )
  opt.super_resolution.interp_method = 'linear' ;
  end

  if ~isfield( opt, 'low_resolution' )
  opt.low_resolution.res = 2 ;
  opt.low_resolution.interp_method = 'linear' ;
  end

  if ~isfield( opt.low_resolution, 'res' )
  opt.low_resolution.res = 2 ;
  end

  if ~isfield( opt.low_resolution, 'interp_method' )
  opt.low_resolution.interp_method = 'linear' ;
  end

  % Number of exact simulations to derive interpolation results: for instance, 2*N+1, -N, -N+1, ..., -1, 0, 1, ..., N-1, N
  if ~isfield( opt, 'n_basis_interpolation' )
  opt.n_basis_interpolation = 5 ;
  end

  % Step in mas for a series of simulations
  if ~isfield( opt, 'step_mas' )
  opt.step_mas = 5 ;
  end
 


