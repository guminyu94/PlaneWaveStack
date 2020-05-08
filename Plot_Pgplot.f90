Module Plot_Pgplot
    contains
    subroutine plot_1d(x,y,y2,y3,x_label,y_label,title,dev)
        real, allocatable, intent(in) :: x(:), y(:) 
        character(len=*), intent(in), optional  :: x_label, y_label, title
        integer :: n_p
        ! can plot 3 line
        real, allocatable, intent(in), optional :: y2(:), y3(:)
        character(len=*), intent(in), optional  :: dev
        real :: ymin, ymax
        
        n_p = size(x)
        
        if (present(dev)) then
            IER = pgbeg(0,dev,1,1)
        else
        ! initialize the pgplot, and select windows as plotting device   
            IER = pgbeg(0,'/WZ',1,1)
        end if
        
        ! adjust x and y lim automatically
        ! start a plot, xlim, ylim, box style, and x y scale    
        
        ! compute the range of y based on all y input
        ymin = minval(y)
        ymax = maxval(y)
        if (present(y2)) then
            ymin = min(ymin,minval(y2))
            ymax = max(ymax,maxval(y2))
        end if
        
        if (present(y3)) then
            ymin = min(ymin,minval(y3))
            ymax = max(ymax,maxval(y3))
        end if
        
        CALL pgenv(minval(x),maxval(x),ymin,ymax,maxval(y),0,1)
        
        if (present(x_label) .AND. present(y_label) .AND. present(title)) then
            ! labeling
            CALL pglab(x_label, y_label, title)
        end if
        
        ! draw point mark
        ! CALL PGPT(5,XS,YS,9)

        ! plot line, num of data points, data
        ! CALL PGSLS(1)
        CALL PGLINE(n_p,x,y)
        
        if (present(y2)) then
            !CALL PGSLS(3)
            CALL PGSCI(2)
            CALL PGLINE(n_p,x,y2)
            ! CALL PGSLS(1)
        end if
        
        if (present(y3)) then
            CALL PGSLS(2)
            CALL PGSCI(1)
            CALL PGLINE(n_p,x,y3)
        end if
    
        CALL pgend
    end subroutine plot_1d
    
    ! plot dot mark
    subroutine plot_dot(x,y,style)
        real, allocatable, intent(in) :: x(:), y(:)
        integer, intent(in), optional :: style
        ! draw point mark
        if (present(style)) then
            CALL PGPT(size(x),x,y,style)
        else 
            CALL PGPT(size(x),x,y,9)
        end if
    end subroutine plot_dot
    
end module Plot_Pgplot