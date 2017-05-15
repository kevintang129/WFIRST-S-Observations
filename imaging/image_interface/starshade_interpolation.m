function starshade_interpolation( opt_in )
% Function that tests the results of imaging starshade exactly and with spatical/wavelength interpolation methods
% History
% 05/14/17: first version. Sergi Hildebrandt (JPL/Caltech)

% Get default options
opt = get_default_options( opt_in ) ;

% Create the reference simulations or study interpolation
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
    disp( sprintf( '(starshade_interpolation) Considering the simulation %i/%i', i_sm_x, n_sm_x * n_sm_y ) )
    opt_sm.x_source_mas = opt.x_source_mas_array( i_sm_x ) ;
      for i_sm_y = 1 : n_sm_y
      opt_sm.y_source_mas = opt.y_source_mas_array( i_sm_y ) ;
      makeStarshadeImage( opt_sm ) ;
      end % i_sm_y
    end % i_sm_x
  end % Create the reference simulations
  


