      subroutine kernel(n, t, t0, dt, softr, val)

        implicit none

        integer, intent(in) :: n
        double precision, dimension(n), intent(in) :: t, softr
        double precision, intent(in) :: t0, dt
        double precision, intent(out) :: val

        integer :: i
        double precision :: norm, k

        val = 0.d0
        norm = 0.d0
        do i=1,n

          if (t(i) - t0 .ge. -dt .and. t(i) - t0 .le. dt) then

            if (t(i) .ge. t0) then
              k = ((t(i) - t0) / dt - 1)
              k = k * k
            else
              k = ((t(i) - t0) / dt + 1)
              k = -k * k
            endif

            norm = norm + k * k
            k = softr(i) * k
            val = val + k

          endif

        enddo

        val = val / norm
        val = val * val

      end subroutine

      subroutine discontinuities(n, t, chi, dt, Q, thresh, ind)

        implicit none

        integer, intent(in) :: n
        double precision, dimension(n), intent(in) :: t, chi
        double precision, intent(in) :: dt, Q, thresh
        integer, intent(out) :: ind

        integer :: i
        double precision :: tmid, maxv, val
        double precision, dimension(n) :: softr

        softr(:) = dsqrt(Q / (Q + chi * chi)) * chi

        ind = -1
        maxv = 0.d0
        do i=1,n - 1

          tmid = 0.5 * (t(i) + t(i + 1))

          val = 0.d0
          call kernel(n, t, tmid, dt, softr, val)

          if (val .ge. thresh .and. val .ge. maxv) then
            ind = i - 1
            maxv = val
          endif

        enddo

      end subroutine
