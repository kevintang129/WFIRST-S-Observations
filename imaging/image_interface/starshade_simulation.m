% NB: include imagescnan in the distribution

function starshade_simulation( opt )

% Simple interface to test some simulation setup
% opt has the ingredients for the simulation
% opt.r_star position of the star in mas, usually on-axis (=0)
% opt.psi_star angle position of the star with respect the on-axis in degrees
% opt.r_planet is an array of planet orbital radius in mas
% opt.phase_planet is an array with the orbital phases of the planets in degrees
% opt.contrast_planet is an array of contrast values for each planet
% opt.save sets whether the final simulation will be stored (1/0=yes/no)
% opt.save_all sets whether each planet simulation will be stored (1/0=yes/no)
% opt.save_fig sets whether to produce/store a figure of the simulation (1/0=yes/no) 
% opt.save_fig_all sets whether to produce/stroe a figure of each element of the simulation before coadding them (1/0=yes/no)
% History:
% Sergi R. Hildebrandt 05/01/17: first version

  if isfield( opt, 'developer' )
  units_image
  else
  units
  end

n_plnt = numel( opt.r_planet ) 

% Simulation
%% Star
[ efDefectImg_str lambdaIn  vecPetalArray ] = makeStarshadeImage( 64, 500, opt.r_star, opt.psi_star, opt ) ;
IntDefectImg_sim = abs( efDefectImg_str ).^2 ; % This is correct, it does compute the intensity for each wavelength
dbstop if error
%make_an_error
%% Accumulating the results for each planet
  for i_plnt = 1 : n_plnt
  clear opt_tmp 
  opt_tmp.r_planet = opt.r_planet( i_plnt ) ;
  [ efDefectImg_tmp ] = makeStarshadeImage( 64, 500, opt.r_planet( i_plnt ), opt.phase_planet( i_plnt ), opt ) ;
  IntDefectImg_sim = IntDefectImg_sim + opt.contrast_planet( i_plnt ) * abs( efDefectImg_tmp ).^2 ;
  end

% Saving the simulation
savePath = './out/';
% Temporary ... FIXME
saveFilename = 'starshade_star_planet_sim.mat' ; 
if opt.save == 1
    save([savePath saveFilename '.mat' ], 'IntDefectImg_sim', 'lambdaIn', 'vecPetalArray')
end


% A figure
% Non-linear transformation of the original image
abs_img = squeeze( IntDefectImg_sim( :, :, 2 ) ) ;
fct_scl = 0. ;
nw_img = log10( abs_img - fct_scl * min( abs_img( : ) ) ) ;
%nw_img = abs_img ;
set(0,'defaultlinelinewidth',1.0);
set(0,'DefaultAxesFontSize',14);
figure( 1 )
clf; setwinsize(gcf,600,500)
imagescnan( nw_img, [ log10( min( opt.contrast_planet ) ) - 3, log10( max( opt.contrast_planet ) ) + 1 ] ) ;
%imagescnan( nw_img, [ ( min( opt.contrast_planet ) ) / 1000, 10*( max( opt.contrast_planet ) ) ] ) ;
colorbar
suptitle( sprintf( 'Starshade simulation: Intensity at %3.0f nm', lambdaIn( 2 ) / nm ) )
title( sprintf( 'Contrasts: %1.1e, %1.1e, %1.1e, %1.1e', opt.contrast_planet( 1 ), opt.contrast_planet( 2 ), opt.contrast_planet( 3 ), opt.contrast_planet( 4 ) ), 'FontSize', 16 )  


dbstop if error
make_an_error

