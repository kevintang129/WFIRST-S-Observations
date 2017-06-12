function starshade_translation( opt )


% For plotting (temporary)
  if ~exist( 'i_src_plt', 'var' )
  i_src_plt = 2 ;
  end


% Some sources
r_src = [ 400, 350, 300.5, 300, 200.5, 200.1, 200 ] ; % mas
psi_src = [ 0, 0, 0, 0, 0, 0, 0 ] ; % deg
n_src=  numel( r_src ) ;
  if numel( psi_src ) ~= n_src
  disp( '(starshade_translation) Number of elements for radius and angle is different. Stopped.' )
  return
  end

opt.developer = 1 ;
% Only one wavelength
opt.delta_lambda_nm = 200 ;
opt.redo = 0 ;
opt.save = 1 ;
% Only one wavelength
opt.delta_lambda_nm = 200 ;
opt = get_default_options( opt ) ;
  for i_src = 1 : n_src
  opt.r_source_mas = r_src( i_src ) ;
  opt.psi_source_deg = psi_src( i_src ) ;
  makeStarshadeImage( opt ) ;
  end

% Loading each simulation (notation from makeStarshadeImage)
dr_out = 'out_dev/' ;
% FOr files created before the new parameter was introduced
Nx_img_tmp = opt.Nx_img ;
diam_img_mas_tmp = opt.diam_img_mas ;
  for i_src = 1 : n_src
  fl_out = sprintf( 'starshade_out_Nx_%i_pix_dl_%inm_dr_%3.1f_mas_psi_%3.1f_deg', opt.Nx_pupil_pix, opt.delta_lambda_nm, r_src( i_src ), psi_src( i_src ) ) ;
    if opt.Nx_img ~= 400, fl_out = sprintf( '%s_Nx_img_%04i', fl_out, opt.Nx_img ) ; end
    if opt.diam_img_mas ~= 2000, fl_out = sprintf( '%s_diam_%04i', fl_out, opt.diam_img_mas ) ; end
  load( [ fl_out '.mat' ] ) 
  opt.Nx_img = Nx_img_tmp ;
  opt.diam_img_mas = diam_img_mas_tmp ;
  int_src( i_src, :, : ) = abs( efDefectImg ).^2 ;
  end

% Testing the translation

% Simple circular shift (assuming the source positions correspond to integer values of pixel numbers)
%% Reference source is the first one
fct_img_1 = opt.Nx_img / 400 ;
mas_pr_px = opt.diam_img_mas / opt.Nx_img ;

  for i_src = 2 : n_src
  idx_circ = ceil( ( r_src( 1 ) - r_src( i_src ) ) / mas_pr_px ) ;
  src_dff( i_src - 1, :, : ) = circshift( squeeze( int_src( i_src, :, : ) ), idx_circ ) - squeeze( int_src( 1, :, : ) ) ; 
  end

% Test
n_tmp = 30 * fct_img_1 ;
a1 = opt.Nx_img / 2 + r_src( 1 ) * sin( 180 / pi * psi_src( 1 ) ) / mas_pr_px - 30 * fct_img_1 ;
a2 = opt.Nx_img / 2 + r_src( 1 ) * sin( 180 / pi * psi_src( 1 ) ) / mas_pr_px + 30 * fct_img_1 ;
b1 = opt.Nx_img / 2 + r_src( 1 ) * cos( 180 / pi * psi_src( 1 ) ) / mas_pr_px - 30 * fct_img_1 ;
b2 = opt.Nx_img / 2 + r_src( 1 ) * cos( 180 / pi * psi_src( 1 ) ) / mas_pr_px + 30 * fct_img_1 ;
  for i_src = 2 : n_src
    for i_tmp = 1 : n_tmp
    idx_circ = ceil( ( r_src( 1 ) - r_src( i_src ) + ( i_tmp - n_tmp / 2 ) ) / mas_pr_px ) ;
    tmp = abs( circshift( squeeze( int_src( i_src, :, : ) ), idx_circ ) - squeeze( int_src( 1, :, : ) ) ) ;
    tmp = tmp( b1 : b2, a1 : a2 ) ;
    src_dff_circ( i_src - 1, i_tmp ) = sum( tmp( : ) ) ;
    end
  end

close all
hold all
figure( 1 )
  for i_src = 1 : n_src - 1
  plot( squeeze( src_dff_circ( i_src, : ) ), '+-' )
  xlim( [ n_tmp / 2 - round( mas_pr_px ), n_tmp / 2  ] )
  lgnd{ i_src } = sprintf( 'r_{src}=%3.2f mas', r_src( i_src + 1 ) ) ;
  end
legend( lgnd ) 
hold off


% Some plotting
figure( 2 )
opt.log10 = 0 ;
  if ~isfield( opt, 'log10' )
  opt.log10 = 0 ;
  end

  if ( opt.log10 )
  imagesc( log10( abs( squeeze( src_dff( i_src_plt - 1, :, : ) ) ) ) ) ; h = colorbar ; ylabel( h, 'Log10(Intensity)', 'FontSize', 16 ) ;
  else
  imagesc( ( abs( squeeze( src_dff( i_src_plt - 1, :, : ) ) ) ) ) ; h = colorbar ; ylabel( h, 'Intensity (Linear)', 'FontSize', 16 ) ;
  end
xlim( [ opt.Nx_img / 2 + r_src( 1 ) * sin( 180 / pi * psi_src( 1 ) ) / mas_pr_px - 50 * fct_img_1, opt.Nx_img / 2 + r_src( 1 ) * sin( 180 / pi * psi_src( 1 ) ) / mas_pr_px + 50 * fct_img_1 ] ) ;
ylim( [ opt.Nx_img / 2 + r_src( 1 ) * cos( 180 / pi * psi_src( 1 ) ) / mas_pr_px - 50 * fct_img_1, opt.Nx_img / 2 + r_src( 1 ) * cos( 180 / pi * psi_src( 1 ) ) / mas_pr_px + 50 * fct_img_1 ] ) ;
title( sprintf( 'r_{src}=%3.2f mas', r_src( i_src_plt ) ) )

dbstop if error
make_an_error

