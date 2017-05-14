function opt = get_default_options( opt )
% Function to define the default options for Starshade imaging
% History:
% 05/12/17: created, Sergi Hildebrandt (JPL/Caltech)


% Developer version or Eric Cady's original version
  if ~isfield( opt, 'developer' )
  opt.developer = 0 ;
  end

% Size of the pupil data in pixels (square)
  if ~isfield( opt, 'Nx_pupil_pix' )
  opt.Nx_pupil_pix = 128 ;
  end

% Step of wavelength to consider
  if ~isfield( opt, 'delta_lambda_nm' )
  opt.delta_lambda_nm = 50 ; % nm
  end

% Separation of the source from the center of the pointing
  if ~isfield( opt, 'r_source_mas' )
  opt.r_source_mas = 0 ; % mas
  end

% Angle of the source with respect the horizontal axis
  if ~isfield( opt, 'psi_source_deg' )
  opt.psi_source_deg = 0 ; % degrees
  end

% Replace by a file in FITS format with an Nx x Nx array if you want a specific pupil
  if ~isfield( opt, 'pupil_file' )
  opt.pupil_file = './in/pupil_D1Kpix_2048.fits' ;
  end

% Saving the output
  if ~isfield( opt, 'save_path' )
  opt.save_path = './out' ;
  end

% Saving the figures
  if ~isfield( opt, 'save_fig' )
  opt.save_fig = './fig' ;
  end

