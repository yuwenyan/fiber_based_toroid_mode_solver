
function L_det=calc_L_det(lambda_guess,fiber_geom,fiber_modes,fem,m_guess,scale_fac)
% disp('debug calc_T_det');
    L_mat=calc_L_mat(fiber_geom,fiber_modes,fem,lambda_guess,m_guess,scale_fac);
    L_det=real(det(L_mat));
    
%     [s_mat,v_mat,d_mat]=svd(L_mat);
%     L_det=real(det(v_mat));
%     disp('out of svd')
    
%     figure(1);clf; mesh(real(v_mat));
%  max(max(v_mat))
%  min(min(v_mat))
%  diag(v_mat)
%  disp('press any key to continue...');
%  pause