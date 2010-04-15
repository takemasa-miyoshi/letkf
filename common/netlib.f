      subroutine daxpy(n,da,dx,incx,dy,incy)
c
c     constant times a vector plus a vector.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double precision dx(*),dy(*),da
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if (da .eq. 0.0d0) return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dy(iy) + da*dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dy(i) + da*dx(i)
   30 continue
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,4
        dy(i) = dy(i) + da*dx(i)
        dy(i + 1) = dy(i + 1) + da*dx(i + 1)
        dy(i + 2) = dy(i + 2) + da*dx(i + 2)
        dy(i + 3) = dy(i + 3) + da*dx(i + 3)
   50 continue
      return
      end
      subroutine  dcopy(n,dx,incx,dy,incy)
c
c     copies a vector, x, to a vector, y.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double precision dx(*),dy(*)
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,7)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dx(i)
   30 continue
      if( n .lt. 7 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,7
        dy(i) = dx(i)
        dy(i + 1) = dx(i + 1)
        dy(i + 2) = dx(i + 2)
        dy(i + 3) = dx(i + 3)
        dy(i + 4) = dx(i + 4)
        dy(i + 5) = dx(i + 5)
        dy(i + 6) = dx(i + 6)
   50 continue
      return
      end
      double precision function ddot(n,dx,incx,dy,incy)
c
c     forms the dot product of two vectors.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double precision dx(*),dy(*),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
c
      ddot = 0.0d0
      dtemp = 0.0d0
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = dtemp + dx(ix)*dy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      ddot = dtemp
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dtemp + dx(i)*dy(i)
   30 continue
      if( n .lt. 5 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) +
     *   dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
   50 continue
   60 ddot = dtemp
      return
      end
      subroutine dspdi(ap,n,kpvt,det,inert,work,job)
      integer n,job
      double precision ap(1),work(1)
      double precision det(2)
      integer kpvt(1),inert(3)
c
c     dspdi computes the determinant, inertia and inverse
c     of a double precision symmetric matrix using the factors from
c     dspfa, where the matrix is stored in packed form.
c
c     on entry
c
c        ap      double precision (n*(n+1)/2)
c                the output from dspfa.
c
c        n       integer
c                the order of the matrix a.
c
c        kpvt    integer(n)
c                the pivot vector from dspfa.
c
c        work    double precision(n)
c                work vector.  contents ignored.
c
c        job     integer
c                job has the decimal expansion  abc  where
c                   if  c .ne. 0, the inverse is computed,
c                   if  b .ne. 0, the determinant is computed,
c                   if  a .ne. 0, the inertia is computed.
c
c                for example, job = 111  gives all three.
c
c     on return
c
c        variables not requested by job are not used.
c
c        ap     contains the upper triangle of the inverse of
c               the original matrix, stored in packed form.
c               the columns of the upper triangle are stored
c               sequentially in a one-dimensional array.
c
c        det    double precision(2)
c               determinant of original matrix.
c               determinant = det(1) * 10.0**det(2)
c               with 1.0 .le. dabs(det(1)) .lt. 10.0
c               or det(1) = 0.0.
c
c        inert  integer(3)
c               the inertia of the original matrix.
c               inert(1)  =  number of positive eigenvalues.
c               inert(2)  =  number of negative eigenvalues.
c               inert(3)  =  number of zero eigenvalues.
c
c     error condition
c
c        a division by zero will occur if the inverse is requested
c        and  dspco  has set rcond .eq. 0.0
c        or  dspfa  has set  info .ne. 0 .
c
c     linpack. this version dated 08/14/78 .
c     james bunch, univ. calif. san diego, argonne nat. lab.
c
c     subroutines and functions
c
c     blas daxpy,dcopy,ddot,dswap
c     fortran dabs,iabs,mod
c
c     internal variables.
c
      double precision akkp1,ddot,temp
      double precision ten,d,t,ak,akp1
      integer ij,ik,ikp1,iks,j,jb,jk,jkp1
      integer k,kk,kkp1,km1,ks,ksj,kskp1,kstep
      logical noinv,nodet,noert
c
      noinv = mod(job,10) .eq. 0
      nodet = mod(job,100)/10 .eq. 0
      noert = mod(job,1000)/100 .eq. 0
c
      if (nodet .and. noert) go to 140
         if (noert) go to 10
            inert(1) = 0
            inert(2) = 0
            inert(3) = 0
   10    continue
         if (nodet) go to 20
            det(1) = 1.0d0
            det(2) = 0.0d0
            ten = 10.0d0
   20    continue
         t = 0.0d0
         ik = 0
         do 130 k = 1, n
            kk = ik + k
            d = ap(kk)
c
c           check if 1 by 1
c
            if (kpvt(k) .gt. 0) go to 50
c
c              2 by 2 block
c              use det (d  s)  =  (d/t * c - t) * t  ,  t = dabs(s)
c                      (s  c)
c              to avoid underflow/overflow troubles.
c              take two passes through scaling.  use  t  for flag.
c
               if (t .ne. 0.0d0) go to 30
                  ikp1 = ik + k
                  kkp1 = ikp1 + k
                  t = dabs(ap(kkp1))
                  d = (d/t)*ap(kkp1+1) - t
               go to 40
   30          continue
                  d = t
                  t = 0.0d0
   40          continue
   50       continue
c
            if (noert) go to 60
               if (d .gt. 0.0d0) inert(1) = inert(1) + 1
               if (d .lt. 0.0d0) inert(2) = inert(2) + 1
               if (d .eq. 0.0d0) inert(3) = inert(3) + 1
   60       continue
c
            if (nodet) go to 120
               det(1) = d*det(1)
               if (det(1) .eq. 0.0d0) go to 110
   70             if (dabs(det(1)) .ge. 1.0d0) go to 80
                     det(1) = ten*det(1)
                     det(2) = det(2) - 1.0d0
                  go to 70
   80             continue
   90             if (dabs(det(1)) .lt. ten) go to 100
                     det(1) = det(1)/ten
                     det(2) = det(2) + 1.0d0
                  go to 90
  100             continue
  110          continue
  120       continue
            ik = ik + k
  130    continue
  140 continue
c
c     compute inverse(a)
c
      if (noinv) go to 270
         k = 1
         ik = 0
  150    if (k .gt. n) go to 260
            km1 = k - 1
            kk = ik + k
            ikp1 = ik + k
            kkp1 = ikp1 + k
            if (kpvt(k) .lt. 0) go to 180
c
c              1 by 1
c
               ap(kk) = 1.0d0/ap(kk)
               if (km1 .lt. 1) go to 170
                  call dcopy(km1,ap(ik+1),1,work,1)
                  ij = 0
                  do 160 j = 1, km1
                     jk = ik + j
                     ap(jk) = ddot(j,ap(ij+1),1,work,1)
                     call daxpy(j-1,work(j),ap(ij+1),1,ap(ik+1),1)
                     ij = ij + j
  160             continue
                  ap(kk) = ap(kk) + ddot(km1,work,1,ap(ik+1),1)
  170          continue
               kstep = 1
            go to 220
  180       continue
c
c              2 by 2
c
               t = dabs(ap(kkp1))
               ak = ap(kk)/t
               akp1 = ap(kkp1+1)/t
               akkp1 = ap(kkp1)/t
               d = t*(ak*akp1 - 1.0d0)
               ap(kk) = akp1/d
               ap(kkp1+1) = ak/d
               ap(kkp1) = -akkp1/d
               if (km1 .lt. 1) go to 210
                  call dcopy(km1,ap(ikp1+1),1,work,1)
                  ij = 0
                  do 190 j = 1, km1
                     jkp1 = ikp1 + j
                     ap(jkp1) = ddot(j,ap(ij+1),1,work,1)
                     call daxpy(j-1,work(j),ap(ij+1),1,ap(ikp1+1),1)
                     ij = ij + j
  190             continue
                  ap(kkp1+1) = ap(kkp1+1)
     *                         + ddot(km1,work,1,ap(ikp1+1),1)
                  ap(kkp1) = ap(kkp1)
     *                       + ddot(km1,ap(ik+1),1,ap(ikp1+1),1)
                  call dcopy(km1,ap(ik+1),1,work,1)
                  ij = 0
                  do 200 j = 1, km1
                     jk = ik + j
                     ap(jk) = ddot(j,ap(ij+1),1,work,1)
                     call daxpy(j-1,work(j),ap(ij+1),1,ap(ik+1),1)
                     ij = ij + j
  200             continue
                  ap(kk) = ap(kk) + ddot(km1,work,1,ap(ik+1),1)
  210          continue
               kstep = 2
  220       continue
c
c           swap
c
            ks = iabs(kpvt(k))
            if (ks .eq. k) go to 250
               iks = (ks*(ks - 1))/2
               call dswap(ks,ap(iks+1),1,ap(ik+1),1)
               ksj = ik + ks
               do 230 jb = ks, k
                  j = k + ks - jb
                  jk = ik + j
                  temp = ap(jk)
                  ap(jk) = ap(ksj)
                  ap(ksj) = temp
                  ksj = ksj - (j - 1)
  230          continue
               if (kstep .eq. 1) go to 240
                  kskp1 = ikp1 + ks
                  temp = ap(kskp1)
                  ap(kskp1) = ap(kkp1)
                  ap(kkp1) = temp
  240          continue
  250       continue
            ik = ik + k
            if (kstep .eq. 2) ik = ik + k + 1
            k = k + kstep
         go to 150
  260    continue
  270 continue
      return
      end
      subroutine dspfa(ap,n,kpvt,info)
      integer n,kpvt(1),info
      double precision ap(1)
c
c     dspfa factors a double precision symmetric matrix stored in
c     packed form by elimination with symmetric pivoting.
c
c     to solve  a*x = b , follow dspfa by dspsl.
c     to compute  inverse(a)*c , follow dspfa by dspsl.
c     to compute  determinant(a) , follow dspfa by dspdi.
c     to compute  inertia(a) , follow dspfa by dspdi.
c     to compute  inverse(a) , follow dspfa by dspdi.
c
c     on entry
c
c        ap      double precision (n*(n+1)/2)
c                the packed form of a symmetric matrix  a .  the
c                columns of the upper triangle are stored sequentially
c                in a one-dimensional array of length  n*(n+1)/2 .
c                see comments below for details.
c
c        n       integer
c                the order of the matrix  a .
c
c     output
c
c        ap      a block diagonal matrix and the multipliers which
c                were used to obtain it stored in packed form.
c                the factorization can be written  a = u*d*trans(u)
c                where  u  is a product of permutation and unit
c                upper triangular matrices , trans(u) is the
c                transpose of  u , and  d  is block diagonal
c                with 1 by 1 and 2 by 2 blocks.
c
c        kpvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if the k-th pivot block is singular. this is
c                     not an error condition for this subroutine,
c                     but it does indicate that dspsl or dspdi may
c                     divide by zero if called.
c
c     packed storage
c
c          the following program segment will pack the upper
c          triangle of a symmetric matrix.
c
c                k = 0
c                do 20 j = 1, n
c                   do 10 i = 1, j
c                      k = k + 1
c                      ap(k)  = a(i,j)
c             10    continue
c             20 continue
c
c     linpack. this version dated 08/14/78 .
c     james bunch, univ. calif. san diego, argonne nat. lab.
c
c     subroutines and functions
c
c     blas daxpy,dswap,idamax
c     fortran dabs,dmax1,dsqrt
c
c     internal variables
c
      double precision ak,akm1,bk,bkm1,denom,mulk,mulkm1,t
      double precision absakk,alpha,colmax,rowmax
      integer idamax,ij,ijj,ik,ikm1,im,imax,imaxp1,imim,imj,imk
      integer j,jj,jk,jkm1,jmax,jmim,k,kk,km1,km1k,km1km1,km2,kstep
      logical swap
c
c
c     initialize
c
c     alpha is used in choosing pivot block size.
      alpha = (1.0d0 + dsqrt(17.0d0))/8.0d0
c
      info = 0
c
c     main loop on k, which goes from n to 1.
c
      k = n
      ik = (n*(n - 1))/2
   10 continue
c
c        leave the loop if k=0 or k=1.
c
c     ...exit
         if (k .eq. 0) go to 200
         if (k .gt. 1) go to 20
            kpvt(1) = 1
            if (ap(1) .eq. 0.0d0) info = 1
c     ......exit
            go to 200
   20    continue
c
c        this section of code determines the kind of
c        elimination to be performed.  when it is completed,
c        kstep will be set to the size of the pivot block, and
c        swap will be set to .true. if an interchange is
c        required.
c
         km1 = k - 1
         kk = ik + k
         absakk = dabs(ap(kk))
c
c        determine the largest off-diagonal element in
c        column k.
c
         imax = idamax(k-1,ap(ik+1),1)
         imk = ik + imax
         colmax = dabs(ap(imk))
         if (absakk .lt. alpha*colmax) go to 30
            kstep = 1
            swap = .false.
         go to 90
   30    continue
c
c           determine the largest off-diagonal element in
c           row imax.
c
            rowmax = 0.0d0
            imaxp1 = imax + 1
            im = imax*(imax - 1)/2
            imj = im + 2*imax
            do 40 j = imaxp1, k
               rowmax = dmax1(rowmax,dabs(ap(imj)))
               imj = imj + j
   40       continue
            if (imax .eq. 1) go to 50
               jmax = idamax(imax-1,ap(im+1),1)
               jmim = jmax + im
               rowmax = dmax1(rowmax,dabs(ap(jmim)))
   50       continue
            imim = imax + im
            if (dabs(ap(imim)) .lt. alpha*rowmax) go to 60
               kstep = 1
               swap = .true.
            go to 80
   60       continue
            if (absakk .lt. alpha*colmax*(colmax/rowmax)) go to 70
               kstep = 1
               swap = .false.
            go to 80
   70       continue
               kstep = 2
               swap = imax .ne. km1
   80       continue
   90    continue
         if (dmax1(absakk,colmax) .ne. 0.0d0) go to 100
c
c           column k is zero.  set info and iterate the loop.
c
            kpvt(k) = k
            info = k
         go to 190
  100    continue
         if (kstep .eq. 2) go to 140
c
c           1 x 1 pivot block.
c
            if (.not.swap) go to 120
c
c              perform an interchange.
c
               call dswap(imax,ap(im+1),1,ap(ik+1),1)
               imj = ik + imax
               do 110 jj = imax, k
                  j = k + imax - jj
                  jk = ik + j
                  t = ap(jk)
                  ap(jk) = ap(imj)
                  ap(imj) = t
                  imj = imj - (j - 1)
  110          continue
  120       continue
c
c           perform the elimination.
c
            ij = ik - (k - 1)
            do 130 jj = 1, km1
               j = k - jj
               jk = ik + j
               mulk = -ap(jk)/ap(kk)
               t = mulk
               call daxpy(j,t,ap(ik+1),1,ap(ij+1),1)
               ijj = ij + j
               ap(jk) = mulk
               ij = ij - (j - 1)
  130       continue
c
c           set the pivot array.
c
            kpvt(k) = k
            if (swap) kpvt(k) = imax
         go to 190
  140    continue
c
c           2 x 2 pivot block.
c
            km1k = ik + k - 1
            ikm1 = ik - (k - 1)
            if (.not.swap) go to 160
c
c              perform an interchange.
c
               call dswap(imax,ap(im+1),1,ap(ikm1+1),1)
               imj = ikm1 + imax
               do 150 jj = imax, km1
                  j = km1 + imax - jj
                  jkm1 = ikm1 + j
                  t = ap(jkm1)
                  ap(jkm1) = ap(imj)
                  ap(imj) = t
                  imj = imj - (j - 1)
  150          continue
               t = ap(km1k)
               ap(km1k) = ap(imk)
               ap(imk) = t
  160       continue
c
c           perform the elimination.
c
            km2 = k - 2
            if (km2 .eq. 0) go to 180
               ak = ap(kk)/ap(km1k)
               km1km1 = ikm1 + k - 1
               akm1 = ap(km1km1)/ap(km1k)
               denom = 1.0d0 - ak*akm1
               ij = ik - (k - 1) - (k - 2)
               do 170 jj = 1, km2
                  j = km1 - jj
                  jk = ik + j
                  bk = ap(jk)/ap(km1k)
                  jkm1 = ikm1 + j
                  bkm1 = ap(jkm1)/ap(km1k)
                  mulk = (akm1*bk - bkm1)/denom
                  mulkm1 = (ak*bkm1 - bk)/denom
                  t = mulk
                  call daxpy(j,t,ap(ik+1),1,ap(ij+1),1)
                  t = mulkm1
                  call daxpy(j,t,ap(ikm1+1),1,ap(ij+1),1)
                  ap(jk) = mulk
                  ap(jkm1) = mulkm1
                  ijj = ij + j
                  ij = ij - (j - 1)
  170          continue
  180       continue
c
c           set the pivot array.
c
            kpvt(k) = 1 - k
            if (swap) kpvt(k) = -imax
            kpvt(k-1) = kpvt(k)
  190    continue
         ik = ik - (k - 1)
         if (kstep .eq. 2) ik = ik - (k - 2)
         k = k - kstep
      go to 10
  200 continue
      return
      end
      subroutine  dswap (n,dx,incx,dy,incy)
c
c     interchanges two vectors.
c     uses unrolled loops for increments equal one.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double precision dx(*),dy(*),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c       code for unequal increments or equal increments not equal
c         to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = dx(ix)
        dx(ix) = dy(iy)
        dy(iy) = dtemp
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c       code for both increments equal to 1
c
c
c       clean-up loop
c
   20 m = mod(n,3)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp
   30 continue
      if( n .lt. 3 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,3
        dtemp = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp
        dtemp = dx(i + 1)
        dx(i + 1) = dy(i + 1)
        dy(i + 1) = dtemp
        dtemp = dx(i + 2)
        dx(i + 2) = dy(i + 2)
        dy(i + 2) = dtemp
   50 continue
      return
      end
      integer function idamax(n,dx,incx)
c
c     finds the index of element having max. absolute value.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double precision dx(*),dmax
      integer i,incx,ix,n
c
      idamax = 0
      if( n.lt.1 .or. incx.le.0 ) return
      idamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      dmax = dabs(dx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(dabs(dx(ix)).le.dmax) go to 5
         idamax = i
         dmax = dabs(dx(ix))
    5    ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 dmax = dabs(dx(1))
      do 30 i = 2,n
         if(dabs(dx(i)).le.dmax) go to 30
         idamax = i
         dmax = dabs(dx(i))
   30 continue
      return
      end
      double precision function pythag(a,b)
      double precision a,b
c
c     finds dsqrt(a**2+b**2) without overflow or destructive underflow
c
      double precision p,r,s,t,u
      p = dmax1(dabs(a),dabs(b))
      if (p .eq. 0.0d0) go to 20
      r = (dmin1(dabs(a),dabs(b))/p)**2
   10 continue
         t = 4.0d0 + r
         if (t .eq. 4.0d0) go to 20
         s = r/t
         u = 1.0d0 + 2.0d0*s
         p = u*p
         r = (s/u)**2 * r
      go to 10
   20 pythag = p
      return
      end
      subroutine rs(nm,n,a,w,matz,z,fv1,fv2,ierr)
c
      integer n,nm,ierr,matz
      double precision a(nm,n),w(n),z(nm,n),fv1(n),fv2(n)
c
c     this subroutine calls the recommended sequence of
c     subroutines from the eigensystem subroutine package (eispack)
c     to find the eigenvalues and eigenvectors (if desired)
c     of a real symmetric matrix.
c
c     on input
c
c        nm  must be set to the row dimension of the two-dimensional
c        array parameters as declared in the calling program
c        dimension statement.
c
c        n  is the order of the matrix  a.
c
c        a  contains the real symmetric matrix.
c
c        matz  is an integer variable set equal to zero if
c        only eigenvalues are desired.  otherwise it is set to
c        any non-zero integer for both eigenvalues and eigenvectors.
c
c     on output
c
c        w  contains the eigenvalues in ascending order.
c
c        z  contains the eigenvectors if matz is not zero.
c
c        ierr  is an integer output variable set equal to an error
c           completion code described in the documentation for tqlrat
c           and tql2.  the normal completion code is zero.
c
c        fv1  and  fv2  are temporary storage arrays.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      if (n .le. nm) go to 10
      ierr = 10 * n
      go to 50
c
   10 if (matz .ne. 0) go to 20
c     .......... find eigenvalues only ..........
      call  tred1(nm,n,a,w,fv1,fv2)
*  tqlrat encounters catastrophic underflow on the Vax
*     call  tqlrat(n,w,fv2,ierr)
      call  tql1(n,w,fv1,ierr)
      go to 50
c     .......... find both eigenvalues and eigenvectors ..........
   20 call  tred2(nm,n,a,w,fv1,z)
      call  tql2(nm,n,w,fv1,z,ierr)
   50 return
      end
      subroutine tql1(n,d,e,ierr)
c
      integer i,j,l,m,n,ii,l1,l2,mml,ierr
      double precision d(n),e(n)
      double precision c,c2,c3,dl1,el1,f,g,h,p,r,s,s2,tst1,tst2,pythag
c
c     this subroutine is a translation of the algol procedure tql1,
c     num. math. 11, 293-306(1968) by bowdler, martin, reinsch, and
c     wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 227-240(1971).
c
c     this subroutine finds the eigenvalues of a symmetric
c     tridiagonal matrix by the ql method.
c
c     on input
c
c        n is the order of the matrix.
c
c        d contains the diagonal elements of the input matrix.
c
c        e contains the subdiagonal elements of the input matrix
c          in its last n-1 positions.  e(1) is arbitrary.
c
c      on output
c
c        d contains the eigenvalues in ascending order.  if an
c          error exit is made, the eigenvalues are correct and
c          ordered for indices 1,2,...ierr-1, but may not be
c          the smallest eigenvalues.
c
c        e has been destroyed.
c
c        ierr is set to
c          zero       for normal return,
c          j          if the j-th eigenvalue has not been
c                     determined after 30 iterations.
c
c     calls pythag for  dsqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      ierr = 0
      if (n .eq. 1) go to 1001
c
      do 100 i = 2, n
  100 e(i-1) = e(i)
c
      f = 0.0d0
      tst1 = 0.0d0
      e(n) = 0.0d0
c
      do 290 l = 1, n
         j = 0
         h = dabs(d(l)) + dabs(e(l))
         if (tst1 .lt. h) tst1 = h
c     .......... look for small sub-diagonal element ..........
         do 110 m = l, n
            tst2 = tst1 + dabs(e(m))
            if (tst2 .eq. tst1) go to 120
c     .......... e(n) is always zero, so there is no exit
c                through the bottom of the loop ..........
  110    continue
c
  120    if (m .eq. l) go to 210
  130    if (j .eq. 30) go to 1000
         j = j + 1
c     .......... form shift ..........
         l1 = l + 1
         l2 = l1 + 1
         g = d(l)
         p = (d(l1) - g) / (2.0d0 * e(l))
         r = pythag(p,1.0d0)
         d(l) = e(l) / (p + dsign(r,p))
         d(l1) = e(l) * (p + dsign(r,p))
         dl1 = d(l1)
         h = g - d(l)
         if (l2 .gt. n) go to 145
c
         do 140 i = l2, n
  140    d(i) = d(i) - h
c
  145    f = f + h
c     .......... ql transformation ..........
         p = d(m)
         c = 1.0d0
         c2 = c
         el1 = e(l1)
         s = 0.0d0
         mml = m - l
c     .......... for i=m-1 step -1 until l do -- ..........
         do 200 ii = 1, mml
            c3 = c2
            c2 = c
            s2 = s
            i = m - ii
            g = c * e(i)
            h = c * p
            r = pythag(p,e(i))
            e(i+1) = s * r
            s = e(i) / r
            c = p / r
            p = c * d(i) - s * g
            d(i+1) = h + s * (c * g + s * d(i))
  200    continue
c
         p = -s * s2 * c3 * el1 * e(l) / dl1
         e(l) = s * p
         d(l) = c * p
         tst2 = tst1 + dabs(e(l))
         if (tst2 .gt. tst1) go to 130
  210    p = d(l) + f
c     .......... order eigenvalues ..........
         if (l .eq. 1) go to 250
c     .......... for i=l step -1 until 2 do -- ..........
         do 230 ii = 2, l
            i = l + 2 - ii
            if (p .ge. d(i-1)) go to 270
            d(i) = d(i-1)
  230    continue
c
  250    i = 1
  270    d(i) = p
  290 continue
c
      go to 1001
c     .......... set error -- no convergence to an
c                eigenvalue after 30 iterations ..........
 1000 ierr = l
 1001 return
      end
      subroutine tql2(nm,n,d,e,z,ierr)
c
      integer i,j,k,l,m,n,ii,l1,l2,nm,mml,ierr
      double precision d(n),e(n),z(nm,n)
      double precision c,c2,c3,dl1,el1,f,g,h,p,r,s,s2,tst1,tst2,pythag
c
c     this subroutine is a translation of the algol procedure tql2,
c     num. math. 11, 293-306(1968) by bowdler, martin, reinsch, and
c     wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 227-240(1971).
c
c     this subroutine finds the eigenvalues and eigenvectors
c     of a symmetric tridiagonal matrix by the ql method.
c     the eigenvectors of a full symmetric matrix can also
c     be found if  tred2  has been used to reduce this
c     full matrix to tridiagonal form.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        d contains the diagonal elements of the input matrix.
c
c        e contains the subdiagonal elements of the input matrix
c          in its last n-1 positions.  e(1) is arbitrary.
c
c        z contains the transformation matrix produced in the
c          reduction by  tred2, if performed.  if the eigenvectors
c          of the tridiagonal matrix are desired, z must contain
c          the identity matrix.
c
c      on output
c
c        d contains the eigenvalues in ascending order.  if an
c          error exit is made, the eigenvalues are correct but
c          unordered for indices 1,2,...,ierr-1.
c
c        e has been destroyed.
c
c        z contains orthonormal eigenvectors of the symmetric
c          tridiagonal (or full) matrix.  if an error exit is made,
c          z contains the eigenvectors associated with the stored
c          eigenvalues.
c
c        ierr is set to
c          zero       for normal return,
c          j          if the j-th eigenvalue has not been
c                     determined after 30 iterations.
c
c     calls pythag for  dsqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      ierr = 0
      if (n .eq. 1) go to 1001
c
      do 100 i = 2, n
  100 e(i-1) = e(i)
c
      f = 0.0d0
      tst1 = 0.0d0
      e(n) = 0.0d0
c
      do 240 l = 1, n
         j = 0
         h = dabs(d(l)) + dabs(e(l))
         if (tst1 .lt. h) tst1 = h
c     .......... look for small sub-diagonal element ..........
         do 110 m = l, n
            tst2 = tst1 + dabs(e(m))
            if (tst2 .eq. tst1) go to 120
c     .......... e(n) is always zero, so there is no exit
c                through the bottom of the loop ..........
  110    continue
c
  120    if (m .eq. l) go to 220
  130    if (j .eq. 30) go to 1000
         j = j + 1
c     .......... form shift ..........
         l1 = l + 1
         l2 = l1 + 1
         g = d(l)
         p = (d(l1) - g) / (2.0d0 * e(l))
         r = pythag(p,1.0d0)
         d(l) = e(l) / (p + dsign(r,p))
         d(l1) = e(l) * (p + dsign(r,p))
         dl1 = d(l1)
         h = g - d(l)
         if (l2 .gt. n) go to 145
c
         do 140 i = l2, n
  140    d(i) = d(i) - h
c
  145    f = f + h
c     .......... ql transformation ..........
         p = d(m)
         c = 1.0d0
         c2 = c
         el1 = e(l1)
         s = 0.0d0
         mml = m - l
c     .......... for i=m-1 step -1 until l do -- ..........
         do 200 ii = 1, mml
            c3 = c2
            c2 = c
            s2 = s
            i = m - ii
            g = c * e(i)
            h = c * p
            r = pythag(p,e(i))
            e(i+1) = s * r
            s = e(i) / r
            c = p / r
            p = c * d(i) - s * g
            d(i+1) = h + s * (c * g + s * d(i))
c     .......... form vector ..........
            do 180 k = 1, n
               h = z(k,i+1)
               z(k,i+1) = s * z(k,i) + c * h
               z(k,i) = c * z(k,i) - s * h
  180       continue
c
  200    continue
c
         p = -s * s2 * c3 * el1 * e(l) / dl1
         e(l) = s * p
         d(l) = c * p
         tst2 = tst1 + dabs(e(l))
         if (tst2 .gt. tst1) go to 130
  220    d(l) = d(l) + f
  240 continue
c     .......... order eigenvalues and eigenvectors ..........
      do 300 ii = 2, n
         i = ii - 1
         k = i
         p = d(i)
c
         do 260 j = ii, n
            if (d(j) .ge. p) go to 260
            k = j
            p = d(j)
  260    continue
c
         if (k .eq. i) go to 300
         d(k) = d(i)
         d(i) = p
c
         do 280 j = 1, n
            p = z(j,i)
            z(j,i) = z(j,k)
            z(j,k) = p
  280    continue
c
  300 continue
c
      go to 1001
c     .......... set error -- no convergence to an
c                eigenvalue after 30 iterations ..........
 1000 ierr = l
 1001 return
      end
      subroutine tred1(nm,n,a,d,e,e2)
c
      integer i,j,k,l,n,ii,nm,jp1
      double precision a(nm,n),d(n),e(n),e2(n)
      double precision f,g,h,scale
c
c     this subroutine is a translation of the algol procedure tred1,
c     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
c
c     this subroutine reduces a real symmetric matrix
c     to a symmetric tridiagonal matrix using
c     orthogonal similarity transformations.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        a contains the real symmetric input matrix.  only the
c          lower triangle of the matrix need be supplied.
c
c     on output
c
c        a contains information about the orthogonal trans-
c          formations used in the reduction in its strict lower
c          triangle.  the full upper triangle of a is unaltered.
c
c        d contains the diagonal elements of the tridiagonal matrix.
c
c        e contains the subdiagonal elements of the tridiagonal
c          matrix in its last n-1 positions.  e(1) is set to zero.
c
c        e2 contains the squares of the corresponding elements of e.
c          e2 may coincide with e if the squares are not needed.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      do 100 i = 1, n
         d(i) = a(n,i)
         a(n,i) = a(i,i)
  100 continue
c     .......... for i=n step -1 until 1 do -- ..........
      do 300 ii = 1, n
         i = n + 1 - ii
         l = i - 1
         h = 0.0d0
         scale = 0.0d0
         if (l .lt. 1) go to 130
c     .......... scale row (algol tol then not needed) ..........
         do 120 k = 1, l
  120    scale = scale + dabs(d(k))
c
         if (scale .ne. 0.0d0) go to 140
c
         do 125 j = 1, l
            d(j) = a(l,j)
            a(l,j) = a(i,j)
            a(i,j) = 0.0d0
  125    continue
c
  130    e(i) = 0.0d0
         e2(i) = 0.0d0
         go to 300
c
  140    do 150 k = 1, l
            d(k) = d(k) / scale
            h = h + d(k) * d(k)
  150    continue
c
         e2(i) = scale * scale * h
         f = d(l)
         g = -dsign(dsqrt(h),f)
         e(i) = scale * g
         h = h - f * g
         d(l) = f - g
         if (l .eq. 1) go to 285
c     .......... form a*u ..........
         do 170 j = 1, l
  170    e(j) = 0.0d0
c
         do 240 j = 1, l
            f = d(j)
            g = e(j) + a(j,j) * f
            jp1 = j + 1
            if (l .lt. jp1) go to 220
c
            do 200 k = jp1, l
               g = g + a(k,j) * d(k)
               e(k) = e(k) + a(k,j) * f
  200       continue
c
  220       e(j) = g
  240    continue
c     .......... form p ..........
         f = 0.0d0
c
         do 245 j = 1, l
            e(j) = e(j) / h
            f = f + e(j) * d(j)
  245    continue
c
         h = f / (h + h)
c     .......... form q ..........
         do 250 j = 1, l
  250    e(j) = e(j) - h * d(j)
c     .......... form reduced a ..........
         do 280 j = 1, l
            f = d(j)
            g = e(j)
c
            do 260 k = j, l
  260       a(k,j) = a(k,j) - f * e(k) - g * d(k)
c
  280    continue
c
  285    do 290 j = 1, l
            f = d(j)
            d(j) = a(l,j)
            a(l,j) = a(i,j)
            a(i,j) = f * scale
  290    continue
c
  300 continue
c
      return
      end
      subroutine tred2(nm,n,a,d,e,z)
c
      integer i,j,k,l,n,ii,nm,jp1
      double precision a(nm,n),d(n),e(n),z(nm,n)
      double precision f,g,h,hh,scale
c
c     this subroutine is a translation of the algol procedure tred2,
c     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
c
c     this subroutine reduces a real symmetric matrix to a
c     symmetric tridiagonal matrix using and accumulating
c     orthogonal similarity transformations.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        a contains the real symmetric input matrix.  only the
c          lower triangle of the matrix need be supplied.
c
c     on output
c
c        d contains the diagonal elements of the tridiagonal matrix.
c
c        e contains the subdiagonal elements of the tridiagonal
c          matrix in its last n-1 positions.  e(1) is set to zero.
c
c        z contains the orthogonal transformation matrix
c          produced in the reduction.
c
c        a and z may coincide.  if distinct, a is unaltered.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      do 100 i = 1, n
c
         do 80 j = i, n
   80    z(j,i) = a(j,i)
c
         d(i) = a(n,i)
  100 continue
c
      if (n .eq. 1) go to 510
c     .......... for i=n step -1 until 2 do -- ..........
      do 300 ii = 2, n
         i = n + 2 - ii
         l = i - 1
         h = 0.0d0
         scale = 0.0d0
         if (l .lt. 2) go to 130
c     .......... scale row (algol tol then not needed) ..........
         do 120 k = 1, l
  120    scale = scale + dabs(d(k))
c
         if (scale .ne. 0.0d0) go to 140
  130    e(i) = d(l)
c
         do 135 j = 1, l
            d(j) = z(l,j)
            z(i,j) = 0.0d0
            z(j,i) = 0.0d0
  135    continue
c
         go to 290
c
  140    do 150 k = 1, l
            d(k) = d(k) / scale
            h = h + d(k) * d(k)
  150    continue
c
         f = d(l)
         g = -dsign(dsqrt(h),f)
         e(i) = scale * g
         h = h - f * g
         d(l) = f - g
c     .......... form a*u ..........
         do 170 j = 1, l
  170    e(j) = 0.0d0
c
         do 240 j = 1, l
            f = d(j)
            z(j,i) = f
            g = e(j) + z(j,j) * f
            jp1 = j + 1
            if (l .lt. jp1) go to 220
c
            do 200 k = jp1, l
               g = g + z(k,j) * d(k)
               e(k) = e(k) + z(k,j) * f
  200       continue
c
  220       e(j) = g
  240    continue
c     .......... form p ..........
         f = 0.0d0
c
         do 245 j = 1, l
            e(j) = e(j) / h
            f = f + e(j) * d(j)
  245    continue
c
         hh = f / (h + h)
c     .......... form q ..........
         do 250 j = 1, l
  250    e(j) = e(j) - hh * d(j)
c     .......... form reduced a ..........
         do 280 j = 1, l
            f = d(j)
            g = e(j)
c
            do 260 k = j, l
  260       z(k,j) = z(k,j) - f * e(k) - g * d(k)
c
            d(j) = z(l,j)
            z(i,j) = 0.0d0
  280    continue
c
  290    d(i) = h
  300 continue
c     .......... accumulation of transformation matrices ..........
      do 500 i = 2, n
         l = i - 1
         z(n,l) = z(l,l)
         z(l,l) = 1.0d0
         h = d(i)
         if (h .eq. 0.0d0) go to 380
c
         do 330 k = 1, l
  330    d(k) = z(k,i) / h
c
         do 360 j = 1, l
            g = 0.0d0
c
            do 340 k = 1, l
  340       g = g + z(k,i) * z(k,j)
c
            do 360 k = 1, l
               z(k,j) = z(k,j) - g * d(k)
  360    continue
c
  380    do 400 k = 1, l
  400    z(k,i) = 0.0d0
c
  500 continue
c
  510 do 520 i = 1, n
         d(i) = z(n,i)
         z(n,i) = 0.0d0
  520 continue
c
      z(n,n) = 1.0d0
      e(1) = 0.0d0
      return
      end
