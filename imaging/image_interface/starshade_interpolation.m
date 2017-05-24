function starshade_interpolation( opt_in )
% Function that tests the results of imaging starshade exactly and with spatical/wavelength interpolation methods
% History
% 05/14/17: first version. Sergi Hildebrandt (JPL/Caltech)
% 05/22/17: added the analysis of the simulations to evaluate the goodness of interpolation

% Get default options
opt = get_default_options( opt_in ) ;

% 1) Create the reference simulations or study interpolation
  if ( opt.create_reference )
  opt.save = 1 ;
  % Star simulation only or with planet
    if ( opt.star )
    % Studying some relative positions
    %% Small displacements for the star (pointing errors) centered at the center
    opt.x_source_mas_array = [ -5 : 1 : 5 ] ; % mas
    opt.y_source_mas_array = [ -5 : 1 : 5 ] ; % mas
    end
    %% Medium displacements for the planet centered at some position
    if ( opt.planet )
    opt.x_source_mas_array = [ opt.x_planet_mas - 20 : 5 : opt.x_planet_mas + 20 ] ;
    opt.y_source_mas_array = [ opt.y_planet_mas - 20 : 5 : opt.y_planet_mas + 20 ] ;
    end

  n_sm_x = numel( opt.x_source_mas_array ) ;
  n_sm_y = numel( opt.y_source_mas_array ) ;
  disp( sprintf( '(starshade_interpolation) Considering %i reference simulations', n_sm_x * n_sm_y ) )
  opt_sm = opt ;
    for i_sm_x = 1 : n_sm_x   
    disp( sprintf( '(starshade_interpolation) Considering the simulation %i/%i', i_sm_x * n_sm_y, n_sm_x * n_sm_y ) )
    opt_sm.x_source_mas = opt.x_source_mas_array( i_sm_x ) ;
      for i_sm_y = 1 : n_sm_y
      opt_sm.y_source_mas = opt.y_source_mas_array( i_sm_y ) ;
      makeStarshadeImage( opt_sm ) ;
      end % i_sm_y
    end % i_sm_x
  end % Create the reference simulations

% 2) Analyze the results
%% 2.1) Exploration (non-optimal in terms of speed and memory, but most accurate)
  if ( 1 )
    % Collecting all the data
  tic
    for i_sm_x = 1 : n_sm_x 
    opt_sm = opt ;
    opt_sm.x_source_mas = opt.x_source_mas_array( i_sm_x ) ;
      for i_sm_y = 1 : n_sm_y
      opt_sm.y_source_mas = opt.y_source_mas_array( i_sm_y ) ;
      opt_sm = get_default_options( opt_sm ) ;
      % Same naming as in makeStarshadeImage
      %% Main parameters
      Nx = opt_sm.Nx_pupil_pix ;
      dlt_lmbd = opt_sm.delta_lambda_nm ;
      r_src = opt_sm.r_source_mas ;
      psi_src = opt_sm.psi_source_deg ;
      ppl_fl = opt_sm.pupil_file ;
      savePath = opt_sm.save_path ;
      saveFilename = sprintf( 'starshade_out_Nx_%i_pix_dl_%inm_dr_%3.1f_mas_psi_%3.1f_deg', Nx, dlt_lmbd, r_src, psi_src ) ;
        if strcmp( ppl_fl, '0' ) == 1, saveFilename = sprintf( '%s_ideal', saveFilename ) ; end
      load( [ savePath '/' saveFilename '.mat' ] )
      % Working with the intensity only
      tmp = abs( efDefectImg( :, :, : ) ).^2 ;
      %% Number of wavelengths considered
      n_lmbd = numel( lambdaIn ) ;
        % First the grid needs some increased resolution
        %% For each wavelength
        for i_lmbd = 1 : n_lmbd
        tmp_spr_rsl( :, : ) = interp2( squeeze( tmp( :, :, i_lmbd ) ), opt.super_resolution.res, opt.super_resolution.interp_method ) ;
        IntDefectImg( :, :, i_lmbd, i_sm_y, i_sm_x ) = tmp_spr_rsl( :, : ) ;
        end
      end % i_sm_y
    end % i_sm_x
 % Size of the super-resolution simulations
 n_2 = size( IntDefectImg, 1 ) ;
 n_1 = size( IntDefectImg, 2 ) ;
 % Mesh where final interpolation will be performed
 [ xq, yq ] = meshgrid( (1 : n_sm_x), (1 : n_sm_y) ) ;
 disp( sprintf( '(starshade_interpolation) Time spent preparing the interpolation data %2.1f seconds', toc ) )
   % Estimating the error performed doing an interpolation
   tic
   %% For each wavelength
     for i_lmbd = 1 : n_lmbd
     % Selecting a subset of simulations as the interpolation basis
     x_basis_arry = ( 1 : opt.low_resolution.res : n_sm_x ) ;
     n_x_basis = numel( x_basis_arry ) ;
     y_basis_arry = ( 1 : opt.low_resolution.res : n_sm_y ) ;
     n_y_basis = numel( y_basis_arry ) ;
     % Redefinition
     [ x_basis, y_basis ] = meshgrid( x_basis_arry, y_basis_arry ) ;
     % Simulation basis
     sim_basis = reshape( squeeze( IntDefectImg( :, :, i_lmbd, y_basis_arry, x_basis_arry ) ), n_2 * n_1, n_y_basis, n_x_basis ) ;
     % this is very slow
     intrp_sim = NaN( n_1 * n_2, n_sm_y, n_sm_x ) ;
     % Going case by case, to center the region of importance
       for i_x = 1 : n_x_basis
         for i_y = 1 : n_y_basis
         % Peak
         tmp = squeeze( IntDefectImg( :, :, i_lmbd, i_x, i_y ) ) ;
         q_pk = find( tmp == max( tmp( : ) ) ) ;
         x_pk = mod( q_pk, size( tmp_spr_rsl, 1 ) ) ;
         y_pk = floor( q_pk / size( tmp_spr_rsl, 1 ) ) ;
           for i_2 = x_pk - 20 : x_pk + 20 
             for i_1 = y_pk - 20 : y_pk + 20
             intrp_sim( n_1 * i_1 + i_2 - n_1, :, : ) = interp2( x_basis, y_basis, squeeze( sim_basis( n_1 * i_1 + i_2 - n_1, :, : ) ), xq, yq, opt.low_resolution.interp_method );
             end
           end
         end
       disp( sprintf( '(starshade_interpolation) Interpolating data: %2.1f/100 done', ( i_2 / n_2 * 100 ) ) )
       end
     disp( sprintf( '(starshade_interpolation) Time spent interpolating data %2.1f seconds', toc ) )
     % Re-shaping the result
     intrp_sim_2d( :, :, i_lmbd, 1 : n_sm_y, 1 : n_sm_x ) = reshape( intrp_sim, n_2, n_1, n_sm_y, n_sm_x ) ;
     end % i_lmbd
dbstop if error
make_an_error
  end % Exploration method

% Some figures
  if ( 0 )
  % Wavelength
  i_lmbd = 1 ;
  % Interpolation basis (odd double indices)
  i_1 = 1 ; i_2 = 1 ;
  % Predictive result (any other case)
  i_1 = 2 ; i_2 = 2 ;
 
  setwinsize(gcf,900,650)
  tmp_1 = squeeze( IntDefectImg( 428:452, 428:452, i_lmbd, i_1, i_2 ) ) ;
  tmp_2 = squeeze( intrp_sim_2d( 428:452, 428:452, i_lmbd, i_1, i_2 ) ) ;
  subplot( 2, 2, 1 )
  imagesc( tmp_1 ) ; title( 'Exact simulation', 'FontSize', 14 ) ; colorbar
  subplot( 2, 2, 2 )
  imagesc( tmp_2 ) ; title( 'Interpolated Simulation', 'FontSize', 14 ) ; colorbar
  subplot( 2, 2, 3 )
  imagesc( tmp_1 - tmp_2 ) ; title( 'Difference exact minus interpolated', 'FontSize', 14 ) ; colorbar ; 
  subplot( 2, 2, 4 )
  imagesc( tmp_1 / max( tmp_1( : ) ) - tmp_2 / max( tmp_2( : ) ) ) ; title( 'Difference exact minus interpolated with peak normalized', 'FontSize', 14 ) ; colorbar ;
  suptitle( sprintf( 'Occulter: %s, lambda=%3.0f nm', opt.occulter_name, lambdaIn( i_lmbd ) * 1e9 ) ) 


  end


  


