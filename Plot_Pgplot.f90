Module Plot_Pgplot
    
    contains
    subroutine plot_1d(x,y,x_label,y_label,title,dev,dots_x,dots_y, style, color, xrange, yrange, lw, style_flag, color_flag)
        real, allocatable, intent(in), dimension(:) :: x
        real, allocatable, intent(in), dimension(:,:) :: y
        character(len=*), intent(in), optional  :: x_label, y_label, title
        integer :: n_p, n_l
        ! plot dev name
        character(len=*), intent(in), optional  :: dev
        real :: ymin, ymax, xmin, xmax
        real, allocatable, optional, intent(in) :: dots_x(:,:), dots_y(:,:)
        integer, allocatable, dimension(:) :: color_array, style_array
        integer, intent(in),optional, dimension(:)  :: style, color
        real, intent(in), optional :: xrange(2)
        real, intent(in), optional :: yrange(2)
        integer, intent(in), optional :: lw
        ! flag to determine whether style or line needs to be auto changed in differnet lines
        integer, intent(in), optional :: style_flag, color_flag
        
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
            
        ! set line width
        if (present(lw)) then
            call pgslw(lw)
        else
            call pgslw(3)
        end if
        
        do i = 1, n_l
            call PGSLS(style_array(i))
            call PGSCI(color_array(i))
            
            ! plot line, num of data points, data
            call PGLINE(n_p,x,y(i,:))           
        end do
        
        call pgend
        deallocate(style_array)
        deallocate(color_array)
        
    end subroutine plot_1d
    
end module Plot_Pgplot