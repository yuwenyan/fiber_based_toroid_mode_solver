% toroid mode solver based on normal fiber modes
% BPM direct expansion
% 2013.06

clear all;
close all;

tic;

% constant
c_m_per_sec=299792458;
u0=4*pi*1e-7;
e0=8.8542e-12;

% fiber parameter
fiber_geom = toroid_structure;
val_m=fiber_geom.M_guess;
num_toroid_mode=6;
scale_fac=1;

toroid_comp_str={'Erad','Eaxi','Eazi','Hrad','Haxi','Hazi'};
fiber_comp_str={'Ex','Ey','Ez','Hx','Hy','Hz'};

m_result=val_m; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mode_calculation=num_toroid_mode;
res_lambda_record=zeros(1,mode_calculation);
res_Q_record=zeros(1,mode_calculation);
fiber_res=fiber_geom;
Resonance=cell(1,mode_calculation);
for mode_order_count=1:mode_calculation
    
    mode_No=mode_order_count;
    
    max_iteration=30;
    count_iteration=0;
    
    lambda_diff=1;
    m_diff=1;
    limit_lambda_diff=0.02e-9;
    limit_m_diff=0.00001;
    
    res_lambda_vec= zeros(1,max_iteration);
    res_lambda_error= zeros(1,max_iteration);
    res_m_error= zeros(1,max_iteration);
    
    
    while (count_iteration<max_iteration)&&(lambda_diff>limit_lambda_diff)
        
        disp(['Toroid mode order: ' num2str(mode_order_count) '   Scan round: ' num2str(count_iteration+1)]);
        disp(['Start computing fiber modes near ' num2str(round(fiber_res.lambda*1e9)) ' nm']);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [fiber_modes,fem_fib]=fiber_fem_solver(fiber_res);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fiber_res.mode_count=length(fiber_modes.Neff);
        disp(['Total ' int2str(fiber_res.mode_count) ' modes found']);
        grid_pts=length(fem_fib.mesh.p(1,:));
        disp(['Total FEM mesh points: ' int2str(grid_pts) ]);
        
        M_mat=calc_M_matrix(fiber_res,fiber_modes,fem_fib);
        res_k0=(2*pi/fiber_res.lambda);
        res_M_mat=res_k0*fiber_res.toroid_R*M_mat;
        % [V_mat,D_mat]=eigs(res_T_mat,1,'LR');
        % [V_mat,D_mat]=eig(res_T_mat);
        % D_vec=diag(D_mat);
        
        [V_matrix,D_mat]=eig(res_M_mat);
        D_array_real=real(diag(D_mat));
        D_array_imag=imag(diag(D_mat));
        [sol,sol_ind]=sort(D_array_real,'descend');
        
        D_indx=sol_ind(mode_No);
        D_val_real=D_array_real(D_indx);
        D_val_imag=D_array_imag(D_indx);
        V_mat=V_matrix(:,D_indx);
        
        
        res_m=m_result; % m
        res_lambda_temp=fiber_res.lambda.*D_val_real./res_m;
        
        lambda_diff=abs(res_lambda_temp-fiber_res.lambda);
        m_diff=abs(D_val_real-res_m);
        res_lambda_vec(count_iteration+1)= res_lambda_temp;
        res_lambda_error(count_iteration+1)=lambda_diff;
        
        fiber_res.lambda=res_lambda_temp;
        count_iteration=count_iteration+1;
        
        disp(['Result: lambda = ' num2str(res_lambda_temp*1e9) ' nm.']);
        
    end % while
    
    time_temp=toc;
    disp(['Elapsed time is ' num2str(round(time_temp/60)) ' mins.']);
    if (count_iteration==max_iteration);
        disp('Maxium iterations!');
    end
    
    %     x_axis=(1:count_iteration);
    %
    %     figure;
    %     [ax,h1,h2]=plotyy(x_axis,res_lambda_vec,x_axis,res_lambda_error);
    %     set(ax(2),'YScale','log');
    %
    %     ylabel(ax(2),'lambda error');
    %     figure;
    %     [m_ax,m_h1,m_h2]=plotyy(x_axis,res_lambda_vec,x_axis,res_m_error);
    %     set(m_ax(2),'YScale','log');
    %     ylabel(m_ax(2),'m error');
    
    res_lambda_record(mode_order_count)=res_lambda_temp;
    res_Q_record(mode_order_count)=D_val_real/D_val_imag;
    
    num_of_eigens=length(D_val_real);
    Field_temp=zeros(6,size(fem_fib.mesh.p,2));
    
    for count_eig=1:num_of_eigens
        Resonance{mode_order_count}.Field=zeros(size(Field_temp));
        Resonance{mode_order_count}.lambda_resonance=res_lambda_temp;
        Resonance{mode_order_count}.m=res_m;
        Resonance{mode_order_count}.V=V_mat;
        V=Resonance{mode_order_count}.V;
        for count=1:size(V,1);
            for count_field=1:length(fiber_comp_str)
                Field_temp(count_field,:)=postinterp(fem_fib,fiber_comp_str{count_field},fem_fib.mesh.p,'solnum',fiber_modes.Neff_indx(count));
            end %for count_field
            Resonance{mode_order_count}.Field=Resonance{mode_order_count}.Field+V(count)*Field_temp;
        end %for count
        
        fiber_fields_bend=Resonance{mode_order_count}.Field;
        
        file_label='Toroid_';
        fig_id=100+mode_order_count;
        figid=plotfields_mode(fem_fib,fiber_fields_bend,fig_id,file_label);
        
    end % count_eig
    
end % mode_order_count

time=toc;
