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
    implicit none
    
    type(Layer), allocatable :: layers(:)
    type(S_Matrix), allocatable :: S_Matrices(:)
    type(Fields) :: inc_field
    
    abstract interface
        subroutine fun_temp(f, larray, inc_field, theta, xi, parameters)
            use Sim_parameters, only : wp
            use Fields_Class
            use Layer_Class
            real(wp), intent(in) :: f
            real(wp), intent(in), optional ::  theta, xi
            type(Layer), allocatable, intent(inout) :: larray(:)
            type(Fields), intent(inout) :: inc_field
            real(wp), intent(in), optional, dimension(:) :: parameters
        end subroutine fun_temp
    end interface
    
    contains
    ! swap freq and obtain the ref and tx coeff
    subroutine freq_swap(fun,freq_start,freq_end,n_points,file_name)
        use Plot_Pgplot
        implicit none
        real(wp), intent(in) :: freq_start, freq_end
        integer, intent(in) :: n_points
        procedure (fun_temp), pointer :: fun
        real(wp) :: freq_step, freq_cur
        complex(wp), dimension(2,2,2) :: tx_ref
        character(len=*), intent(in) :: file_name
        integer :: i, j
        logical :: file_exists
        type(Fields), allocatable :: fields_layer(:)
        real, allocatable :: freq_array(:), coeff_array_1(:), coeff_array_2(:), coeff_array_3(:),angle_array(:)
        
        allocate(freq_array(n_points))
        allocate(coeff_array_1(n_points))
        allocate(coeff_array_2(n_points))
        allocate(coeff_array_3(n_points))
        allocate(angle_array(n_points))
        
        if (n_points .EQ. 1) then
            freq_step = 0.0_wp
        else
            freq_step = (freq_end - freq_start) / real(n_points-1,wp)
        end if
        
        ! openfile
        inquire(FILE = file_name // '_tx_ref' // '.dat', EXIST=file_exists)
        if ( file_exists ) then
            open(1, file = file_name // '_tx_ref' // '.dat', status = 'old')
        else
            open(1, file = file_name // '_tx_ref' // '.dat', status = 'new')
        end if
        
    
        do j = 1, n_points
            !load parmeters
            freq_cur = real(j-1,wp) * freq_step + freq_start
            
            call fun(freq_cur,layers,inc_field)
            call compute_S_matrix
            
            ! obtain reflection coeff
            tx_ref = trans_ref_coeff_freespace(S_Matrices)
            
            
            ! write to file
            write(1,*) freq_cur, ABS(tx_ref(1,1,1))**2.0_wp, ABS(tx_ref(1,1,2))**2.0_wp, ABS(tx_ref(1,2,1))**2.0_wp, ABS(tx_ref(1,2,2))**2.0_wp, ABS(tx_ref(2,1,1))**2.0_wp, ABS(tx_ref(2,1,2))**2.0_wp, ABS(tx_ref(2,2,1))**2.0_wp, ABS(tx_ref(2,2,2))**2.0_wp
            
            freq_array(j) = real(freq_cur/1.0E12_wp)
            
            !coeff_array_1(j) = real(ABS(tx_ref(2,1,1)))
            !coeff_array_2(j) = real(ABS(tx_ref(2,2,2)))
            !(j) = real(ABS(tx_ref(2,2,1)))
            
            angle_array(j) = ellipse_angle(tx_ref(1,1,1),tx_ref(1,2,1)) / PI * 180.0_wp * (abs(tx_ref(1,1,1))**2.0 + abs(tx_ref(1,2,1))**2.0) 
            coeff_array_1(j) = ellipse_angle(tx_ref(2,1,1),tx_ref(2,2,1)) / PI * 180.0_wp 
            coeff_array_2(j) = abs(tx_ref(2,1,1))**2.0 
            coeff_array_3(j) = abs(tx_ref(2,2,1))**2.0
            
            
            
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
    
        close(1) 
        
        ! call plotting subroutine
        call plot_1d(freq_array,coeff_array_1,  x_label = 'Freq(THz)', y_label = 'Angle (degrees) ', title = 'Faraday Rotation Angle', dev='/PS')
        print *,'Max Amp: ', maxval(coeff_array_2)
        print *,'Max Angle: ', maxval(coeff_array_1)
        ! call plot_1d(freq_array, coeff_array_1, 'Freq(THz)','Rx', 'Rx')
        ! call plot_1d(freq_array, coeff_array_2, 'Freq(THz)','Rx', 'Rx')
       ! call plot_1d(freq_array, coeff_array_1, coeff_array_2, coeff_array_3, 'Freq(THz)','OM Angle (Degrees)', '', dev = '/WZ')
        deallocate(freq_array)
        deallocate(coeff_array_1)
        deallocate(coeff_array_2)
        deallocate(coeff_array_3)
        deallocate(angle_array)
        
    end subroutine freq_swap
    
    subroutine theta_swap(fun,freq,theta_start,theta_end,n_points,file_name)
        use Plot_Pgplot
        implicit none
        real(wp), intent(in) :: theta_start, theta_end, freq
        integer, intent(in) :: n_points
        procedure (fun_temp), pointer :: fun
        real(wp) :: theta_cur, theta_step
        complex(wp), dimension(2,2,2) :: tx_ref
        character(len=*), intent(in) :: file_name
        integer :: i, j
        logical :: file_exists
        type(Fields), allocatable :: fields_layer(:)
        real, allocatable :: theta_array(:), coeff_array_1(:), coeff_array_2(:), coeff_array_3(:)
        
        allocate(theta_array(n_points))
        allocate(coeff_array_1(n_points))
        allocate(coeff_array_2(n_points))
        allocate(coeff_array_3(n_points))
        
        if (n_points .EQ. 1) then
            theta_step = 0.0_wp
        else
            theta_step = (theta_end - theta_start) / real(n_points-1,wp)
        end if
    
        do j = 1, n_points
            !load parmeters
            theta_cur = real(j-1,wp) * theta_step + theta_start
            
            ! assume xi is 0
            call fun(freq, layers, inc_field,theta_cur, 0.0_wp )
                        
            call compute_S_matrix
            
            ! obtain reflection coeff
            tx_ref = trans_ref_coeff_freespace(S_Matrices)
            
            theta_array(j) = real(theta_cur)
            !angle_array(j) = real(ellipse_angle(tx_ref(1,1,1),tx_ref(1,2,1)))
            coeff_array_1(j) = real(ABS(tx_ref(2,1,1))**2.0_wp)
            coeff_array_2(j) = real(ABS(tx_ref(1,1,1))**2.0_wp)
            coeff_array_3(j) = real( 1.0 - abs(tx_ref(1,1,1))**2.0_wp - abs(tx_ref(2,1,1))**2.0_wp)
            
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
        call plot_1d(theta_array, coeff_array_1, x_label = 'theta (degrees)',y_label = '',  dev = '/PS')
        
        deallocate(theta_array)
        deallocate(coeff_array_1)
        deallocate(coeff_array_2)
        deallocate(coeff_array_3)
    end subroutine theta_swap
    
    ! obtain fields based on freq, fun stands for the function pointer of model's configuration
    subroutine fields_computation(fun,freq_in, n_points, file_name)
        use Plot_Pgplot
        implicit none
        procedure (fun_temp), pointer :: fun
        real(wp), intent(in) :: freq_in
        integer, intent(in) :: n_points
        character(len=*), intent(in) :: file_name
        real(wp) :: totol_z=0.0_wp, current_z=0.0_wp, layer_z=0.0_wp, dz=0.0_wp, dz_round = 0.0_wp
        integer :: i, current_layer = 1
        complex(wp), dimension(2) :: fields_forward, fields_backward, E_field, H_field
        complex(wp), dimension(2,2) :: Pn_forward, Pn_backward
        logical :: file_exists
        ! z plotting range for the first and last layer
        real(wp) :: z_extension
        
        real, allocatable :: z_array(:), field_array(:)
        
        ! layer's field amplitude
        Type(Fields), allocatable :: fields_amp_layer(:)
        
        ! plotting data array
        allocate(z_array(n_points))
        allocate(field_array(n_points))
        
        ! load model paramters
        call fun(freq_in, layers, inc_field)
        call compute_S_matrix
        
        ! 1st layer amplitdue omitted 
        allocate(fields_amp_layer(n_layers))
        
        ! compute each layer's fields amplitude
        call cal_fields_amp(S_Matrices, inc_field, fields_amp_layer)
        
        ! output data into a file
        inquire(FILE=file_name // '_fields' // '.dat', EXIST = file_exists)
        if (file_exists) then
            open(10, file = file_name // '_fields' // '.dat', status = 'old')
        else
            open(10, file = file_name // '_fields' // '.dat', status = 'new')
        end if
        
        ! calculate the total thickness
        do i = 1,n_layers
            totol_z = totol_z + layers(i)%d
        end do
        
        if (totol_z .EQ. 0.0_wp)  then
            z_extension = lambda_0
        else
            z_extension = totol_z / 4.0_wp
        end if
        totol_z = totol_z + 2.0_wp * z_extension
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
            write(10,*) REAL(current_z), REAL(Abs(E_field(1))), REAL(Abs(E_field(2))), REAL(Abs(H_field(1)*eta_0)), REAL(Abs(H_field(2)*eta_0))
            
            ! plot the E_e field
            z_array(i) = REAL(current_z)
            field_array(i) = (real(E_field(1)))

            current_z = current_z + dz
            layer_z = layer_z + dz
        end do
        
        ! call plotting subroutine
        call plot_1d(z_array,field_array, x_label = 'z(m)', y_label = 'E_Field_e(v/m)', title = 'Fields Plot')
    end subroutine fields_computation
    
    subroutine compute_S_matrix()
        integer i
        ! allocate the S_Matrix
        if (.NOT. allocated(S_Matrices)) then
            allocate(S_Matrices(n_layers-1))
        end if
        ! compute S matrix
        do i = 1, n_layers-1
            ! assemble S Matrix from layer obj
            S_Matrices(i) = S_Matrix(layers(i),layers(i+1))
        end do
    end subroutine compute_S_matrix
    
    ! calculate fields' amplitude on each layer
    subroutine cal_fields_amp(S_matrices_in, inc_fields, fields_amp_layer_in)
        implicit none
        type(S_Matrix), intent(in), allocatable :: S_matrices_in(:)
        type(Fields), intent(in) :: inc_fields
        type(Fields), allocatable, intent(inout) :: fields_amp_layer_in(:)
        integer :: n
        type(S_Matrix) :: S_matrices_Cascaded_1_n, S_matrices_Cascaded_n_N
        integer :: n_S_matrices
        complex(wp), dimension(2,2,2) :: tx_ref
        complex(wp), dimension(2) :: field_forward, field_backward
        n_S_matrices = SIZEOF(S_matrices_in)/SIZEOF(S_matrices_in(1))
        
        ! calculate last layer, assume free space and no inc field from the last layer
        tx_ref = trans_ref_coeff_freespace(S_matrices_in)
        
        fields_amp_layer_in(1)%Field_1 = inc_fields%Field_1
        fields_amp_layer_in(1)%Field_2 = MATMUL(tx_ref(2,:,:),inc_fields%Field_1)

        !print*, 1
        !print*, 'v_forward: ', abs(fields_amp_layer_in(1)%Field_1)
        !print*, 'v_backward: ', abs(fields_amp_layer_in(1)%Field_2)
        
        do n = 1, n_S_matrices-1
            S_matrices_Cascaded_1_n = S_Matrices_Cascade(S_matrices_in,1,n)
            S_matrices_Cascaded_n_N = S_Matrices_Cascade(S_matrices_in,n+1,n_S_matrices)
            field_forward = MATMUL((unit_matrix - MATMUL( S_matrices_Cascaded_1_n%beta_n , S_matrices_Cascaded_n_N%alpha_n ))**-1 , MATMUL( S_matrices_Cascaded_1_n%delta_n,inc_fields%Field_1) )
            field_backward = MATMUL((unit_matrix - MATMUL( S_matrices_Cascaded_n_N%alpha_n , S_matrices_Cascaded_1_n%beta_n ))**-1 , MATMUL( MATMUL ( S_matrices_Cascaded_n_N%alpha_n , S_matrices_Cascaded_1_n%delta_n), inc_fields%Field_1) )
            fields_amp_layer_in( n + 1 ) = Fields(field_forward,field_backward)
            !print*, n + 1
            !print*, 'v_forward: ', abs(fields_amp_layer_in(n+1)%Field_1)
            !print*, 'v_backward: ', abs(fields_amp_layer_in(n+1)%Field_2)
        end do
        
        fields_amp_layer_in(n_S_matrices+1)%Field_1 = MATMUL(tx_ref(1,:,:),inc_fields%Field_1)
        fields_amp_layer_in(n_S_matrices+1)%Field_2 = inc_fields%Field_2
        !print*, n+1
        !print*, 'v_forward: ', abs(fields_amp_layer_in(n_S_matrices+1)%Field_1)
        !print*, 'v_backward: ', abs(fields_amp_layer_in(n_S_matrices+1)%Field_2)
        
        
    end subroutine cal_fields_amp
    
        ! tx and ref coeff of last layer with free space 
    pure function trans_ref_coeff_freespace(S_matrices_in) result(trans_ref_coeff)
        implicit none    
        type(S_Matrix), intent(in), allocatable:: S_matrices_in(:)
        complex(wp), dimension(2,2,2) :: trans_ref_coeff
        type(S_Matrix) :: S_matrices_Cascaded
        integer :: n_S_matrices
        n_S_matrices = SIZEOF(S_matrices_in)/SIZEOF(S_matrices_in(1))
        ! cascade S Matrices
        S_matrices_Cascaded = S_Matrices_Cascade(S_matrices_in, 1, n_S_matrices)
        trans_ref_coeff(1,:,:) = S_matrices_Cascaded%delta_n
        trans_ref_coeff(2,:,:) = S_matrices_Cascaded%alpha_n
    end function trans_ref_coeff_freespace
    
    ! tx and ref coeff of last layer with PEC backed
    pure function trans_ref_coeff_freespace_pec(S_matrices_in, P_n_in) result(trans_ref_coeff)
        implicit none
        type(S_Matrix), intent(in), allocatable :: S_matrices_in(:)
        complex(wp), dimension(2,2), intent(in) :: P_n_in
        complex(wp), dimension(2,2,2) :: trans_ref_coeff
        type(S_Matrix) :: S_matrices_Cascaded
        integer :: n_S_matrices
        
        n_S_matrices = SIZEOF(S_matrices_in)/SIZEOF(S_matrices_in(1))
        ! cascade S Matrices
        S_matrices_Cascaded = S_Matrices_Cascade(S_matrices_in, 1, n_S_matrices)
        trans_ref_coeff(1,:,:) = MATMUL((unit_matrix + MATMUL(S_matrices_Cascaded%beta_n , P_n_in))**-1 , S_matrices_Cascaded%delta_n)
        trans_ref_coeff(2,:,:) = S_matrices_Cascaded%alpha_n - MATMUL( MATMUL(S_matrices_Cascaded%gamma_n , P_n_in) , trans_ref_coeff(1,:,:))
    end function trans_ref_coeff_freespace_pec
  
end module Swapper