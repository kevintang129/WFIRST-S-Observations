function starshade_translation( opt )


% For plotting (temporary)
  if ~exist( 'i_src_plt', 'var' )
  i_src_plt = 3 ;
  end


% Some sources
r_src = [ 400, 305, 302.5, 302, 300, 202.5, 202, 200, 152.5, 150.5, 150, 102.5, 101, 100 ] ; % 350, 304, 303, 301, 300.5, 200.5, 200.1, 200 ] ; % mas
psi_src = [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ] ; % deg
%r_src = [ 400, 302.5, 300, 200, 102.5, 101, 100 ] ;
%psi_src = [ 90, 90, 90, 90, 90, 90, 90 ] ;
%r_src = [ 400, sqrt( 200^2 + 200^2 ), 200, 200, sqrt( 151^2 + 17^2 ) ] ;
%psi_src = [ 0, 45, 0, 90, 180 / pi * atan2( 17, 151 ) ] ;
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
opt = get_default_options( opt ) ;
  for i_src = 1 : n_src
  opt.r_source_mas = r_src( i_src ) ;
  opt.psi_source_deg = psi_src( i_src ) ;
  makeStarshadeImage( opt ) ;
  end

% Loading each simulation (notation from makeStarshadeImage)
dr_out = 'out_dev/' ;
% For files created before the new parameter was introduced
Nx_img_tmp = opt.Nx_img ;
diam_img_mas_tmp = opt.diam_img_mas ;
  for i_src = 1 : n_src
  fl_out = sprintf( 'starshade_out_Nx_%i_pix_dl_%inm_dr_%3.1f_mas_psi_%3.1f_deg', opt.Nx_pupil_pix, opt.delta_lambda_nm, r_src( i_src ), psi_src( i_src ) ) ;
    if opt.Nx_img ~= 400, fl_out = sprintf( '%s_Nx_img_%04i', fl_out, opt.Nx_img ) ; end
    if opt.diam_img_mas ~= 2000, fl_out = sprintf( '%s_diam_%04i', fl_out, opt.diam_img_mas ) ; end
  load( [ fl_out '.mat' ] ) 
  opt.Nx_img = Nx_img_tmp ;
  opt.diam_img_mas = diam_img_mas_tmp ;
  n_lmbd = size( efDefectImg, 3 ) ;
    for i_lmbd = 1 : n_lmbd
    int_src( i_src, i_lmbd, :, : ) = abs( squeeze( efDefectImg( :, :, i_lmbd ) ) ).^2 ;
    end
  end

% Testing the translation

% Simple circular shift (assuming the source positions correspond to integer values of pixel numbers)
%% Reference source is the first one
fct_img_1 = opt.Nx_img / 400 ;
mas_pr_px = opt.diam_img_mas / opt.Nx_img ;
clear src_dff
  for i_src_rf = 1 : n_src ;
    for i_src = 1 : n_src
    idx_circ_x = ceil( ( r_src( i_src_rf ) * cos( pi / 180 * psi_src( i_src_rf ) ) - r_src( i_src ) * cos( pi / 180 * psi_src( i_src ) ) )  / mas_pr_px ) ;
    idx_circ_y = ceil( ( r_src( i_src_rf ) * sin( pi / 180 * psi_src( i_src_rf ) ) - r_src( i_src ) * sin( pi / 180 * psi_src( i_src ) ) ) / mas_pr_px ) ;
      for i_lmbd =1 : n_lmbd
      a1 = opt.Nx_img / 2 + r_src( i_src_rf ) * sin( pi / 180 * psi_src( i_src_rf ) ) / mas_pr_px - 50 * fct_img_1 ;
      a2 = opt.Nx_img / 2 + r_src( i_src_rf ) * sin( pi / 180 * psi_src( i_src_rf ) ) / mas_pr_px + 50 * fct_img_1 ; 
      b1 = opt.Nx_img / 2 + r_src( i_src_rf ) * cos( pi / 180 * psi_src( i_src_rf ) ) / mas_pr_px - 50 * fct_img_1 ;
      b2 = opt.Nx_img / 2 + r_src( i_src_rf ) * cos( pi / 180 * psi_src( i_src_rf ) ) / mas_pr_px + 50 * fct_img_1 ;
      tmp = circshift( squeeze( int_src( i_src, i_lmbd, :, : ) ), [ idx_circ_x idx_circ_y ] ) - squeeze( int_src( i_src_rf, i_lmbd, :, : ) ) ;
      src_dff( i_src_rf, i_src, i_lmbd, :, : ) = tmp( round( b1 ) : round( b2 ),  round( a1 ) : round( a2 ) ) ;
      end
    end
  end

% Reduced to the standard size
%% Two reductions before & after
clear src_dff_0
  for i_src_rf = 1 : n_src ;
    for i_src = 1 : n_src
    idx_circ_x = ceil( ( r_src( i_src_rf ) * cos( pi / 180 * psi_src( i_src_rf ) ) - r_src( i_src ) * cos( pi / 180 * psi_src( i_src ) ) ) / mas_pr_px / fct_img_1 ) ;
    idx_circ_y = ceil( ( r_src( i_src_rf ) * sin( pi / 180 * psi_src( i_src_rf ) ) - r_src( i_src ) * sin( pi / 180 * psi_src( i_src ) ) ) / mas_pr_px / fct_img_1 ) ;
      for i_lmbd =1 : n_lmbd
      a1 = round( opt.Nx_img / 2 / fct_img_1 + r_src( i_src_rf ) * sin( pi / 180 * psi_src( i_src_rf ) ) / mas_pr_px / fct_img_1 - 50 ) ;
      a2 = round( opt.Nx_img / 2 / fct_img_1 + r_src( i_src_rf ) * sin( pi / 180 * psi_src( i_src_rf ) ) / mas_pr_px / fct_img_1 + 50 ) ;
      b1 = round( opt.Nx_img / 2 / fct_img_1 + r_src( i_src_rf ) * cos( pi / 180 * psi_src( i_src_rf ) ) / mas_pr_px / fct_img_1 - 50 ) ;
      b2 = round( opt.Nx_img / 2 / fct_img_1 + r_src( i_src_rf ) * cos( pi / 180 * psi_src( i_src_rf ) ) / mas_pr_px / fct_img_1 + 50 ) ;
      tmp = circshift( imresize( squeeze( int_src( i_src, i_lmbd, :, : ) ), 1 / fct_img_1 ), [ idx_circ_x idx_circ_y ] ) - imresize( squeeze( int_src( i_src_rf, i_lmbd, :, : ) ), 1 / fct_img_1 ) ;
      % Before
      src_dff_0( 1, i_src_rf, i_src, i_lmbd, 1 : b2 - b1 + 1, 1 : a2 - a1 + 1 ) = tmp( b1 : b2, a1 : a2 ) ;
      % after
      tmp = imresize( squeeze( src_dff( i_src_rf, i_src, i_lmbd, :, : ) ), 1 / fct_img_1 ) ;
      src_dff_0( 2, i_src_rf, i_src, i_lmbd, 1 : size( tmp, 1 ), 1 : size( tmp, 2 ) ) = tmp;
      end
    end
  end

% Some plotting
n_src  = numel( r_src ) ;
opt.log10 = 0 ;
%% Reduced resolution
lbl_rd = { 'before shifting', 'after shifting' } ;
if ( 0 )
  for i_rd = 1 : 2
    for i_src_rf = 1 : n_src
    figure( ( i_rd - 1 ) * n_src + i_src_rf )
    setwinsize( gcf, 1300, 750 )
    clf
    i_lmbd = 1 ;
      for i_src_plt = 1 : n_src
      subplot( 4, 4, i_src_plt )
        if ( opt.log10 )
        imagesc( log10( abs( squeeze( src_dff_0( i_rd, i_src_rf, i_src_plt, i_lmbd, :, : ) ) ) ) ) ; h = colorbar ; ylabel( h, 'Log10(Intensity)', 'FontSize', 16 ) ;
        else
        imagesc( ( abs( squeeze( src_dff_0( i_rd, i_src_rf, i_src_plt, i_lmbd, :, : ) ) ) ) ) ; h = colorbar ; ylabel( h, 'Intensity (Linear)', 'FontSize', 16 ) ;
        end
      set( gca, 'xtick', [], 'ytick', [] )
      title( sprintf( 'r_{src}=%3.2f mas, r_{ref}=%3.2f', r_src( i_src_plt ), r_src( i_src_rf ) ) )
      end
    suptitle( sprintf( 'Reduction %s', lbl_rd{ i_rd } ) ) ;
    end
  end % i_rd
end

%% Input resolution
  for i_src_rf = 1 : n_src
  figure( 2 * n_src + i_src_rf )
  setwinsize( gcf, 1300, 750 )
  clf
  i_lmbd = 1 ;
    for i_src_plt = 1 : numel( r_src )
    subplot( 4, 4, i_src_plt )
      if ( opt.log10 )
      imagesc( log10( abs( squeeze( src_dff( i_src_rf, i_src_plt, i_lmbd, :, : ) ) ) ) ) ; h = colorbar ; ylabel( h, 'Log10(Intensity)', 'FontSize', 16 ) ;
      else
      imagesc( ( abs( squeeze( src_dff( i_src_rf, i_src_plt, i_lmbd, :, : ) ) ) ) ) ; h = colorbar ; ylabel( h, 'Intensity (Linear)', 'FontSize', 16 ) ;
      end
    set( gca, 'xtick', [], 'ytick', [] ) 
    title( sprintf( 'r_{src}=%3.2f mas, r_{ref}=%3.2f', r_src( i_src_plt ), r_src( i_src_rf ) ) )
    end
suptitle( sprintf( 'Initial resolution %ix%i times the input one', fct_img_1, fct_img_1 ) )
  end

dbstop if error
make_an_error

