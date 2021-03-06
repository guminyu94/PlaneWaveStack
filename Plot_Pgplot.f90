!****************************************************************************
!
!   Module: Plot_Pgplot
!
!   PURPOSE: A wrapper for pgplot 1d plotting 
!
!****************************************************************************
Module Plot_Pgplot
    type dots 
        real, allocatable :: dots_data(:)
    end type dots
    
    contains
    subroutine plot_1d(x,y,x_label,y_label,title,dev,dots_x,dots_y, style, color, xrange, yrange, lw, style_flag, color_flag, font_size, font_style,legend, dots_style,data_points_style,dots_legend,dots_color)
        ! data
        real, allocatable, intent(in), dimension(:) :: x
        real, allocatable, intent(in), dimension(:,:) :: y
        integer :: n_p, n_l, n_d_p, n_d_l
        ! label
        character(len=*), intent(in), optional  :: x_label, y_label, title
        ! plot dev name
        character(len=*), intent(in), optional  :: dev
        real :: ymin, ymax, xmin, xmax
        ! mark point
        type(dots), allocatable, optional, intent(in) :: dots_x(:), dots_y(:)
        ! line style, color and width
        integer, allocatable, dimension(:) :: color_array, style_array, dots_style_array, dots_color_array
        integer, intent(in),optional, dimension(:)  :: style, color, dots_style, data_points_style, dots_color
        integer, intent(in), optional :: lw
        ! x y axis range
        real, intent(in), optional :: xrange(2),yrange(2)
        ! flag to determine whether style or line needs to be auto changed in differnet lines
        integer, intent(in), optional :: style_flag, color_flag
        ! label font size
        real, intent(in), optional :: font_size
        integer, intent(in), optional :: font_style
        ! legend
        character(len=*), intent(in), optional :: legend(:),dots_legend(:)
        ! legend length, 0 - 1
        real :: legend_length, legend_gap, legend_x_offset, legend_y_offset, legend_box_gap, legend_box_width
        
        n_p = size(y,2)
        n_l = size(y,1)
        
        if (present(dev)) then
            IER = pgbeg(0,dev,1,1)
        else
        ! initialize the pgplot, and select windows as plotting device   
            IER = pgbeg(0,'/WZ',1,1)
        end if
        
        ! adjust x and y lim automatically or based on input
        
        if (present(xrange)) then
            xmin = xrange(1)
            xmax = xrange(2)
        else
            xmin = minval(x)
            xmax = maxval(x)
        end if
        
        ! compute the range of y based on all y input
        if (present(yrange)) then
            ymin = yrange(1)
            ymax = yrange(2)
        else
            ymin = minval(minval(y,2),1)
            ymax = maxval(maxval(y,2),1)
        end if
        
        if (present(font_style)) then
            call pgscf(font_style)
        else
            ! font style
            call pgscf(2);
        end if
        
        if (present(font_size)) then 
            call PGSCH(font_size)
        else
            ! font size
            call PGSCH(1.2)
        end if
        
        ! line width of frame and label
        call pgslw(3)
        
        ! start a plot, xlim, ylim, box style, and x y scale
        CALL pgenv(xmin,xmax,ymin,ymax,0,1)
        
        
        if (present(x_label) .AND. present(y_label) .AND. present(title)) then
            ! labeling
            CALL pglab(x_label, y_label, title)
        else if ( present(x_label) .AND. present(y_label) ) then
            CALL pglab(x_label, y_label, '')
        end if
        
        ! determine style and color
        allocate(style_array(n_l))
        allocate(color_array(n_l))
        
        if (.NOT. present(style)) then
            if (present(style_flag) .and. (style_flag .EQ. 1)) then
                do i = 1,n_l
                    style_array(i) = i
                end do
            else 
                do i = 1,n_l
                    style_array(i) = 1
                end do
            end if
        else 
            do i = 1,n_l
                    style_array(i) = style(i)
            end do
        end if
        
        if (.NOT. present(color)) then
            if (present(color_flag) .and. (color_flag .EQ. 1)) then
                do i = 1,n_l
                    color_array(i) = i
                end do
            else 
                do i = 1,n_l
                    color_array(i) = 1
                end do
            end if
        else 
            do i = 1,n_l
                    color_array(i) = color(i)
            end do
        end if
        
        ! parameters for legend pos and distance
        legend_length = 0.1
        legend_gap = 0.05
        legend_x_offset = 0.56
        legend_y_offset = 0.32
        legend_yaxis_length = 38.0
        legend_box_width = 3.6
        legend_box_gap = 0.02
        
        do i = 1, n_l
            ! set line width
            if (present(lw)) then
                call pgslw(lw)
            else
                call pgslw(3)
            end if
            
            call PGSLS(style_array(i))
            call PGSCI(color_array(i))
            
            ! plot line, num of data points, data
            call PGLINE(n_p,x,y(i,:))        
            
            ! plot data points
            if (present(data_points_style)) then
                call pgslw(30)
                call PGPT(n_p, x, y(i,:), data_points_style(i))
            end if
            
            ! legend
            if (present(legend)) then
                ! set line width
                if (present(lw)) then
                    call pgslw(lw)
                else
                    call pgslw(3)
                end if
                call PGLINE(2,(/abs(xmax-xmin) * legend_x_offset+xmin, abs(xmax-xmin) * legend_x_offset + abs(xmax-xmin) * legend_length + xmin/), (/abs(ymax-ymin) * legend_y_offset - real(i-1)* abs(ymax-ymin) * legend_gap + ymin, abs(ymax-ymin) * legend_y_offset - real(i-1)* abs(ymax-ymin) * legend_gap + ymin/))
                !call pgarro(abs(xmax-xmin)*7.5/10.0 , abs(ymax-ymin)/2.0 * 8.5/10.0 - real(i-1)* abs(ymax-ymin) * 1.0/10.0 * legend_gap, abs(xmax-xmin)*7.5/10.0 + abs(xmax-xmin)*1.0/10.0*legend_length, abs(ymax-ymin)/2.0 * 8.5/10.0 - real(i-1)* abs(ymax-ymin) * 1.0/10.0 * legend_gap);
                
                ! legend line width
                call pgslw(1)
                ! legend font size
                call PGSCH(0.8)
                ! legend font style
                call pgscf(3)
                call PGMTXT('B',- (legend_yaxis_length * legend_y_offset) + real(i-1) * legend_gap * legend_yaxis_length , legend_x_offset + legend_length + 0.01,0.0,trim(legend(i)))
            
            end if
            
        end do
        
        ! plot data points
        if (present(dots_x) .and. present(dots_y)) then
            n_d_l = size(dots_y)
            allocate(dots_style_array(n_d_l))
            allocate(dots_color_array(n_d_l))
            do i = 1, n_d_l
                n_d_p = size(dots_y(i)%dots_data)
                if (present(dots_style)) then
                    dots_style_array(i) = dots_style(i)
                else
                    dots_style_array(i) = i
                end if
                if (present(dots_color)) then
                    dots_color_array(i) = dots_color(i)
                else
                    dots_color_array(i) = i
                end if
                
                ! plot data points
                call pgslw(50)
                call PGSCI(color_array(i))
                call PGPT(n_d_p, dots_x(i)%dots_data, dots_y(i)%dots_data, dots_style_array(i))
                ! dots legend
                if (present(dots_legend)) then
                    ! set line width
                    call PGPT(1,(/abs(xmax-xmin) * legend_x_offset+xmin+abs(xmax-xmin) * legend_length/2.0/),(/abs(ymax-ymin) * legend_y_offset - real(i-1)* abs(ymax-ymin) * legend_gap + ymin - real(n_l)* abs(ymax-ymin) * legend_gap /),dots_style_array(i))
                    
                    ! legend line width
                    call pgslw(1)
                    ! legend font size
                    call PGSCH(0.8)
                    ! legend font style
                    call pgscf(3)
                    call PGMTXT('B',- (legend_yaxis_length * legend_y_offset) + real(i-1) * legend_gap * legend_yaxis_length + real(n_l) * legend_gap * legend_yaxis_length, legend_x_offset + legend_length + 0.01,0.0,trim(dots_legend(i)))
            
                end if
            
            end do
            
        end if
        
        ! add legend bounding box
        if (present(legend)) then
            ! box line width
            call pgslw(2)
            ! box line color, filling, style
            call PGSLS(1)
            call PGSFS(2)
            call PGSCI(1)

            ! draw box
            call PGRECT(legend_x_offset * abs(xmax-xmin) + xmin - legend_box_gap * abs(xmax-xmin), abs(xmax-xmin) * legend_x_offset + xmin + legend_box_gap * abs(xmax-xmin) + abs(xmax-xmin) * legend_length * legend_box_width, abs(ymax-ymin) * legend_y_offset + ymin + abs(ymax-ymin) * legend_box_gap * 1.5, abs(ymax-ymin) * legend_y_offset - (size(legend)+size(dots_legend)-1.0) * abs(ymax-ymin) * legend_gap + ymin - abs(ymax-ymin) * legend_box_gap )
        end if
            
        call pgend
        deallocate(style_array)
        deallocate(color_array)
        
    end subroutine plot_1d
    
end module Plot_Pgplot