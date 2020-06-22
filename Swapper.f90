!****************************************************************************
!
!   Module: Swapper
!
!   PURPOSE: Swap parameters and output to file
!
!****************************************************************************
Module Swapper
    use Sim_parameters, only : wp, n_layers, theta, xi, eta_0, lambda_0
    use Layer_Class
    use S_Matrix_Class
    use Fields_Class
    use Math
    use data_global
    implicit none
    
    type(Layer), allocatable :: layers(:)
    type(S_Matrix), allocatable :: S_Matrices(:)
    type(Fields) :: inc_field
    
    abstract interface
        subroutine fun_temp(f, larray, inc_field, theta_in, xi_in, parameters)
            use Sim_parameters, only : wp
            use Fields_Class
            use Layer_Class
            real(wp), intent(in) :: f
            real(wp), intent(in), optional ::  theta_in, xi_in
            type(Layer), allocatable, intent(inout) :: larray(:)
            type(Fields), intent(inout) :: inc_field
            real(wp), intent(in), optional, dimension(:) :: parameters
        end subroutine fun_temp
    end interface
    
    contains
    ! swap freq and obtain the ref and tx coeff
    subroutine freq_swap(fun,freq_start,freq_end,n_points,file_name,output,savefig_flag,save_data_flag)
        use Plot_Pgplot
        use Sim_parameters, only : pec_flag
        implicit none
        real(wp), intent(in) :: freq_start, freq_end
        integer, intent(in) :: n_points
        procedure (fun_temp), pointer :: fun
        real(wp) :: freq_step, freq_cur
        complex(wp), dimension(2,2,2) :: tx_ref
        character(len=*), intent(in), optional :: file_name
        integer :: i, j
        logical :: file_exists
        type(Fields), allocatable :: fields_layer(:)
        real, allocatable :: freq_array(:), data_array(:,:), angle_array(:,:) 
        real(wp), optional, allocatable, intent(inout) :: output(:)
        integer, intent(in), optional :: savefig_flag, save_data_flag
        
        if (.not. allocated(freq_array)) then
            allocate(freq_array(n_points))
        end if 
        if (.not. allocated(data_array)) then
            allocate(data_array(2,n_points))
        end if
        if (.not. allocated(angle_array)) then
            allocate(angle_array(2,n_points))
        end if
        
        if (present(output) .and. (.NOT. allocated(output))) then
            allocate(output(2))
        end if
        
        if (n_points .EQ. 1) then
            freq_step = 0.0_wp
        else
            freq_step = (freq_end - freq_start) / real(n_points-1,wp)
        end if
        
        ! openfile
        if (present(file_name)) then
            inquire(FILE = file_name // '_tx_ref' // '.dat', EXIST=file_exists)
            if ( file_exists ) then
                open(11, file = file_name // '_tx_ref' // '.dat', status = 'old')
            else
                open(11, file = file_name // '_tx_ref' // '.dat', status = 'new')
            end if
        end if
        
        do j = 1, n_points
            !load parmeters
            freq_cur = real(j-1,wp) * freq_step + freq_start
            
            call fun(freq_cur,layers,inc_field)
            call compute_S_matrix(layers,S_Matrices)
            
            if (pec_flag .EQ. 1) then
                ! obtain reflection coeff of PEC backed
                tx_ref = trans_ref_coeff_pec(S_Matrices,layers(size(layers))%P_n)
            else
                 ! obtain reflection coeff of free space
                tx_ref = trans_ref_coeff_freespace(S_Matrices)
            end if
            
            ! write to file
            if (present(file_name)) then
                write(11,*) freq_cur, ABS(tx_ref(1,1,1))**2.0_wp, ABS(tx_ref(1,1,2))**2.0_wp, ABS(tx_ref(1,2,1))**2.0_wp, ABS(tx_ref(1,2,2))**2.0_wp, ABS(tx_ref(2,1,1))**2.0_wp, ABS(tx_ref(2,1,2))**2.0_wp, ABS(tx_ref(2,2,1))**2.0_wp, ABS(tx_ref(2,2,2))**2.0_wp
            end if
                
            freq_array(j) = real(freq_cur/1.0E12_wp)
            
            ! save data for plotting
            angle_array(1,j) = real(ellipse_angle(tx_ref(1,1,1),tx_ref(1,2,1))) / PI * 180.0
            data_array(1,j) = real((abs(tx_ref(1,1,1))**2.0 +  abs(tx_ref(1,2,1))**2.0)**0.5)
            !data_array(1,j) = (abs(tx_ref(2,1,1))**2.0 +  abs(tx_ref(2,2,1))**2.0)**0.5
            !data_array(2,j) = ellipticity(tx_ref(2,1,1),tx_ref(2,2,1))
            print*, data_array(1,j)
            ! print *, 'Freq: ', freq_cur/1.0E12_wp, ' THz'
            ! print *, 'Tx_Faraday_rot_angle: ', angle_array(j) , ' Degree'
            
            !print *, "tx_coeff_ee"
            !print *, ABS(tx_ref(1,1,1))
            !print *, "|tx_coeff_hh|^2"
            !print *, ABS(tx_ref(1,2,2))**2.0_wp
            !print *, "tx_coeff_he"
            !print *, (tx_ref(1,2,1))
            !print *, "|ref_coeff_ee|^2"
            !print *, ABS(tx_ref(2,1,1))**2.0_wp
            !print *, "|ref_coeff_hh|^2"
            !print *, ABS(tx_ref(2,2,2))**2.0_wp
            !print *, "|ref_coeff_he|^2"
            !print *, ABS(tx_ref(2,2,1))**2.0_wp
            
        end do
        
        ! phase unwrap
        !call phase_unwrap_1d(angle_array,1,90.0)
        !angle_array(2,:) = data_array(1,:) * (( (angle_array(1,:))**2.0 + (data_array(2,:) / PI * 180.0)**2.0  ) **0.5)
       
        ! save data for global plotting
        if (present(save_data_flag) .and. (save_data_flag .eq. 1)) then
            if (counter .eq. 1) then
                data_1 = freq_array
            end if
        
            ! save the data to global array
            data_2(counter+1,:) = -1.0*angle_array(1,:)
            data_3(counter+1,:) = data_array(1,:)
            counter = counter + 1
        end if
        
        !print *, 'FOM_BEST: ', maxval(angle_array(2,:))
        
        !print *, 'Rot_Angle: ', angle_array(1,2001)
        !print *, 'FOM: ', angle_array(2,2001)
        !print *, 'Ref: ', data_array(1,2001)
        !print *, 'Ell: ', tan(data_array(2,2001))
        
        ! call plotting subroutine
        if (present(savefig_flag) .and. (savefig_flag .EQ. 1) ) then 
            call plot_1d(freq_array, data_array,  x_label = 'Freq(THz)', y_label = 'Amp (A.U.) ', title = 'Ref Coeff Plot', dev='ref_coeff.ps/PS', color = (/1,2,2/), style = (/1,2,3/),legend = (/'Reflectance','Ellipticity'/))
            call plot_1d(freq_array, angle_array,  x_label = 'Freq(THz)', y_label = 'Angle (degrees) ', title = 'Kerr Rotation Angle & FOM', dev='Kerr_rot.ps/PS', style_flag = 1, color_flag = 1, legend = (/'Kerr Rot Angle','FOM'/))
            print*, 'Saved to Kerr_rot.ps'
            print*, 'Saved to ref_coeff.ps'
        else if (present(savefig_flag) .and. (savefig_flag .EQ. 2)) then
            call plot_1d(freq_array, angle_array,  x_label = 'Freq(THz)', y_label = 'Angle (degrees) ', title = 'Kerr Rotation Angle & FOM', dev='/WZ', style_flag = 1, color_flag = 1)
        end if
            
        if (present(file_name)) then
            close(11)
        end if
    
        deallocate(freq_array)
        deallocate(data_array)
        deallocate(angle_array)
    end subroutine freq_swap
    
    subroutine theta_swap(fun,freq,theta_start,theta_end,n_points)
        use Plot_Pgplot
        use Sim_parameters, only : pec_flag
        implicit none
        real(wp), intent(in) :: theta_start, theta_end, freq
        integer, intent(in) :: n_points
        procedure (fun_temp), pointer :: fun
        real(wp) :: theta_cur, theta_step
        complex(wp), dimension(2,2,2) :: tx_ref
        integer :: i, j
        logical :: file_exists
        type(Fields), allocatable :: fields_layer(:)
        real, allocatable :: theta_array(:), angle_array(:,:), data_array(:,:)
        
        allocate(theta_array(n_points))
        allocate(angle_array(2,n_points))
        allocate(data_array(3,n_points))
        
        if (n_points .EQ. 1) then
            theta_step = 0.0_wp
        else
            theta_step = (theta_end - theta_start) / real(n_points-1,wp)
        end if
    
        do j = 1, n_points
            !load parmeters
            theta_cur = real(j-1,wp) * theta_step + theta_start
            
            ! assume xi is 0
            call fun(freq, layers, inc_field,theta_in = theta_cur,xi_in = 0.0_wp )
                        
            call compute_S_matrix(layers,S_Matrices)
            
            if (pec_flag .EQ. 1) then
                ! obtain reflection coeff of PEC backed
                tx_ref = trans_ref_coeff_pec(S_Matrices,layers(size(layers))%P_n)
            else
                 ! obtain reflection coeff of free space
                tx_ref = trans_ref_coeff_freespace(S_Matrices)
            end if
            
            theta_array(j) = real(theta_cur) 
            
            !angle_array(1,j) = real(ellipse_angle(tx_ref(2,1,1),tx_ref(2,2,1))) / PI * 180.0
            data_array(1,j) = abs(tx_ref(1,1,1))**2.0 
            data_array(2,j) = abs(tx_ref(2,1,1))**2.0
            data_array(3,j) = 1-(abs(tx_ref(1,1,1))**2.0 + abs(tx_ref(2,1,1))**2.0)**0.5
            !angle_array(2,j) = data_array(1,j) * (angle_array(1,j) + data_array(2,j))
            !print *, 'ex: ', abs(tx_ref(2,1,1))**2.0
            !print *, 'ey: ', abs(tx_ref(2,2,1)* cos(theta_cur/pi*180))**2.0 
            !print *, 'Freq: ', freq_cur/1.0E12_wp, 'THz'
            !print *, 'Tx_Faraday_rot_angle: ', angle_array(j) / PI * 180.0, ' Degree'
            
            !print *, "tx_coeff_ee"
            !print *, ABS(tx_ref(1,1,1))
            !print *, "|tx_coeff_hh|^2"
            !print *, ABS(tx_ref(1,2,2))**2.0_wp
            !print *, "tx_coeff_he"
            !print *, (tx_ref(1,2,1))
            !print *, "|ref_coeff_ee|^2"
            !print *, ABS(tx_ref(2,1,1))**2.0_wp
            !print *, "|ref_coeff_hh|^2"
            !print *, ABS(tx_ref(2,2,2))**2.0_wp
            !print *, "|ref_coeff_he|^2"
            !print *, ABS(tx_ref(2,2,1))**2.0_wp
        end do
        
        ! call plotting subroutine
        ! call plot_1d(theta_array, angle_array, x_label = 'theta (degrees)',y_label = 'Angle (Degrees)',  dev = 'theta_plot.ps/PS',color_flag = 1, style_flag = 1,legend=(/'Kerr Rot Angle','FOM'/))
        call plot_1d(theta_array, data_array,  x_label = '\(0685) (degrees)', y_label = 'Amplitude (A.U.) ', title = '', dev='theta_ref_coeff.ps/CPS', color = (/1,2,3/), style = (/1,2,3/),legend = (/'Transmittance','Reflectance','Absorption'/))
        print*, 'theta plot saved'
        deallocate(theta_array)
        deallocate(data_array)
        deallocate(angle_array)
    end subroutine theta_swap
    
    ! obtain fields based on freq, fun stands for the function pointer of model's configuration
    subroutine fields_computation(fun,freq, n_points, file_name, savefig_flag)
        use Plot_Pgplot
        use Sim_parameters, only : pec_flag
        implicit none
        procedure (fun_temp), pointer :: fun
        real(wp), intent(in) :: freq
        integer, intent(in) :: n_points
        character(len=*), intent(in), optional :: file_name
        real(wp) :: totol_z=0.0_wp, current_z=0.0_wp, layer_z=0.0_wp, dz=0.0_wp, dz_round = 0.0_wp
        integer :: i, current_layer = 1
        complex(wp), dimension(2) :: fields_forward, fields_backward, E_field, H_field
        complex(wp), dimension(2,2) :: Pn_forward, Pn_backward
        logical :: file_exists
        ! z plotting range for the first and last layer
        real(wp) :: z_extension
        real, allocatable :: z_array(:), field_array(:,:)
        integer, intent(in), optional :: savefig_flag
        ! layer's field amplitude
        Type(Fields), allocatable :: fields_amp_layer(:)
        
        ! plotting data array
        allocate(z_array(n_points))
        allocate(field_array(2,n_points))
        
        ! load model paramters
        call fun(freq, layers, inc_field)
        call compute_S_matrix(layers,S_Matrices)
        
        ! compute each layer's fields amplitude
        call cal_fields_amp(S_Matrices, inc_field, fields_amp_layer)
        
        if (present(file_name)) then      
            ! output data into a file
            inquire(FILE=file_name // '_fields' // '.dat', EXIST = file_exists)
            if (file_exists) then
                open(10, file = file_name // '_fields' // '.dat', status = 'old')
            else
                open(10, file = file_name // '_fields' // '.dat', status = 'new')
            end if
        end if
        
        ! calculate the total thickness
        do i = 1,n_layers
            totol_z = totol_z + layers(i)%d
        end do
        
        ! extend the z by quarter total thickness
        if (totol_z .EQ. 0.0_wp)  then
            z_extension = lambda_0
        else
            z_extension = totol_z / 4.0_wp
        end if
        
        ! if pec backed, then 
        if (pec_flag .EQ. 1) then
            ! do not plot pec field
            totol_z = totol_z + z_extension
        else 
            totol_z = totol_z + 2.0_wp * z_extension
        end if
        dz = totol_z / real((n_points-1),wp)
        dz_round = dz/1E8
        
        layer_z = -1.0_wp * z_extension 
        current_z = -1.0_wp * z_extension 
        
        do i = 1, n_points
            ! update z-zn, dz_round to eliminate round-up error
            if ( (layer_z - layers(current_layer)%d > dz_round) .and. (current_layer .NE. n_layers)) then
                layer_z = layer_z - layers(current_layer)%d 
                current_layer = current_layer + 1
            end if
            
                        
            ! calculate forward backward field using phase
            ! current layer start from 1, which omits the layer of inc field
            Pn_forward(1,1) = exp(-1.0_wp*(0.0_wp,1.0_wp)*layers(current_layer)%kz_e*(layer_z))
            Pn_forward(1,2) = (0.0_wp,0.0_wp)
            Pn_forward(2,1) = (0.0_wp,0.0_wp)
            Pn_forward(2,2) = exp(-1.0_wp*(0.0_wp,1.0_wp)*layers(current_layer)%kz_h*(layer_z))
            
            Pn_backward(1,1) = exp((0.0_wp,1.0_wp)*layers(current_layer)%kz_e*(layer_z-layers(current_layer)%d))
            Pn_backward(2,1) = (0.0_wp,0.0_wp)
            Pn_backward(1,2) = (0.0_wp,0.0_wp)
            Pn_backward(2,2) = exp((0.0_wp,1.0_wp)*layers(current_layer)%kz_h*(layer_z-layers(current_layer)%d))
            
            fields_forward = MATMUL(Pn_forward, fields_amp_layer(current_layer)%Field_1)
            fields_backward = MATMUL(Pn_backward, fields_amp_layer(current_layer)%Field_2)
            E_field = fields_forward + fields_backward
            
            ! not sure which sign of eta should I choose
            H_field = MATMUL(layers(current_layer)%Y_n, fields_forward) - MATMUL(layers(current_layer)%Y_n, fields_backward)
            
            if (present(file_name)) then 
                write(10,*) REAL(current_z), REAL(Abs(E_field(1))), REAL(Abs(E_field(2))), REAL(Abs(H_field(1)*eta_0)), REAL(Abs(H_field(2)*eta_0))
            end if
            
            ! save the E field
            z_array(i) = REAL(current_z)
            field_array(1,i) = (real(E_field(1)))
            field_array(2,i) = (real(E_field(2)))
            
            current_z = current_z + dz
            layer_z = layer_z + dz
        end do

        ! call plotting subroutine
        if (present(savefig_flag) .and. (savefig_flag .eq. 1)) then
            call plot_1d(z_array,field_array, x_label = 'z(m)', y_label = 'E_Field(V/m)', title = 'E Fields Plot', style_flag = 1,color_flag = 1, dev='e_fields.ps/CPS')
            print*, 'E_field fig saved'
        else if (present(savefig_flag) .and. (savefig_flag .eq. 2)) then
            call plot_1d(z_array,field_array, x_label = 'z(m)', y_label = 'E_Field(V/m)', title = 'E Fields Plot', style_flag = 1,dev='/WZ')
        end if
        
        deallocate(z_array)
        deallocate(field_array)
        deallocate(fields_amp_layer)
        
        if (present(file_name)) then
            close(10)
        end if
    end subroutine fields_computation

    
    ! calculate fields' amplitude on each layer
    subroutine cal_fields_amp(S_matrices, inc_fields, fields_amp_layer)
        use Sim_parameters, only : pec_flag 
        implicit none
        type(S_Matrix), intent(in), allocatable :: S_matrices(:)
        type(Fields), intent(in) :: inc_fields
        type(Fields), allocatable, intent(inout) :: fields_amp_layer(:)
        integer :: n
        type(S_Matrix) :: S_matrices_Cascaded_1_n, S_matrices_Cascaded_n_N
        integer :: n_S_matrices
        complex(wp), dimension(2,2,2) :: tx_ref
        complex(wp), dimension(2) :: field_forward, field_backward
        n_S_matrices = SIZEOF(S_matrices)/SIZEOF(S_matrices(1))
        
        if (.not. allocated(fields_amp_layer)) then
            allocate(fields_amp_layer(n_layers))
        else
            deallocate(fields_amp_layer)
            allocate(fields_amp_layer(n_layers))
        end if
        
        if (pec_flag .EQ. 1) then
            ! pec back, update the reflected wave on the last layer
            tx_ref = trans_ref_coeff_pec(S_matrices,layers(size(layers))%P_n)
            ! calculate the amplitude of last layer
            fields_amp_layer(n_S_matrices+1)%Field_1 = MATMUL(tx_ref(1,:,:),inc_fields%Field_1)
            fields_amp_layer(n_S_matrices+1)%Field_2 = MATMUL(-layers(size(layers))%P_n, fields_amp_layer(n_S_matrices+1)%Field_1)
        else
            ! calculate last layer, free space and no inc field from the last layer
            tx_ref = trans_ref_coeff_freespace(S_matrices)
            ! calculate the amplitude of last layer
            fields_amp_layer(n_S_matrices+1)%Field_1 = MATMUL(tx_ref(1,:,:),inc_fields%Field_1)
            fields_amp_layer(n_S_matrices+1)%Field_2 = inc_fields%Field_2
        end if
        
        
        fields_amp_layer(1)%Field_1 = inc_fields%Field_1
        fields_amp_layer(1)%Field_2 = MATMUL(tx_ref(2,:,:),inc_fields%Field_1)

        !print*, 1
        !print*, 'v_forward: ', abs(fields_amp_layer(1)%Field_1)
        !print*, 'v_backward: ', abs(fields_amp_layer(1)%Field_2)
        
        do n = 1, n_S_matrices-1
            S_matrices_Cascaded_1_n = S_Matrices_Cascade(S_matrices,1,n)
            S_matrices_Cascaded_n_N = S_Matrices_Cascade(S_matrices,n+1,n_S_matrices)
            field_forward = MATMUL((unit_matrix - MATMUL( S_matrices_Cascaded_1_n%beta_n , S_matrices_Cascaded_n_N%alpha_n ))**-1 , MATMUL( S_matrices_Cascaded_1_n%delta_n,inc_fields%Field_1) + MATMUL( MATMUL(S_matrices_Cascaded_1_n%beta_n,S_matrices_Cascaded_n_N%gamma_n),fields_amp_layer(n_S_matrices+1)%Field_2) )
            field_backward = MATMUL((unit_matrix - MATMUL( S_matrices_Cascaded_n_N%alpha_n , S_matrices_Cascaded_1_n%beta_n ))**-1 , MATMUL( MATMUL ( S_matrices_Cascaded_n_N%alpha_n , S_matrices_Cascaded_1_n%delta_n), inc_fields%Field_1) + MATMUL( S_matrices_Cascaded_n_N%gamma_n, fields_amp_layer(n_S_matrices+1)%Field_2 ))
            fields_amp_layer( n + 1 ) = Fields(field_forward,field_backward)
            
            !print*, n + 1
            !print*, 'v_forward: ', abs(fields_amp_layer(n+1)%Field_1)
            !print*, 'v_backward: ', abs(fields_amp_layer(n+1)%Field_2)
        end do
        
        !print*, n+1
        !print*, 'v_forward: ', abs(fields_amp_layer(n_S_matrices+1)%Field_1)
        !print*, 'v_backward: ', abs(fields_amp_layer(n_S_matrices+1)%Field_2)
                
    end subroutine cal_fields_amp
    
        ! tx and ref coeff of last layer with free space 
    pure function trans_ref_coeff_freespace(S_matrices) result(trans_ref_coeff)
        implicit none    
        type(S_Matrix), intent(in), allocatable:: S_matrices(:)
        complex(wp), dimension(2,2,2) :: trans_ref_coeff
        type(S_Matrix) :: S_matrices_Cascaded
        integer :: n_S_matrices
        n_S_matrices = SIZEOF(S_matrices)/SIZEOF(S_matrices(1))
        ! cascade S Matrices
        S_matrices_Cascaded = S_Matrices_Cascade(S_matrices, 1, n_S_matrices)
        trans_ref_coeff(1,:,:) = S_matrices_Cascaded%delta_n
        trans_ref_coeff(2,:,:) = S_matrices_Cascaded%alpha_n
    end function trans_ref_coeff_freespace
    
    ! tx and ref coeff of last layer with PEC backed
    pure function trans_ref_coeff_pec(S_matrices, P_n) result(trans_ref_coeff)
        implicit none
        ! the S matrices should contain the last layer, namely S_1,N
        type(S_Matrix), intent(in), allocatable :: S_matrices(:)
        ! P_n is the last layer's(backed by pec) phase matrix
        complex(wp), dimension(2,2), intent(in) :: P_n
        complex(wp), dimension(2,2,2) :: trans_ref_coeff
        type(S_Matrix) :: S_matrices_Cascaded
        integer :: n_S_matrices
        
        n_S_matrices = SIZEOF(S_matrices)/SIZEOF(S_matrices(1))
        ! cascade S Matrices
        S_matrices_Cascaded = S_Matrices_Cascade(S_matrices, 1, n_S_matrices)
        trans_ref_coeff(1,:,:) = MATMUL((unit_matrix + MATMUL(S_matrices_Cascaded%beta_n , P_n))**-1 , S_matrices_Cascaded%delta_n)
        trans_ref_coeff(2,:,:) = S_matrices_Cascaded%alpha_n - MATMUL( MATMUL(S_matrices_Cascaded%gamma_n , P_n) , trans_ref_coeff(1,:,:))
    end function trans_ref_coeff_pec
  
end module Swapper