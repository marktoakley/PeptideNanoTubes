! Analyse hydrogen-bonding pattern in cyclic peptide nanotubes

module hb
   implicit none
! peptides array contains coordinates of atoms
!   residue number 
!   atom number: 1 O, 2 C, 3 N, 4 H
!   axis: 1 x, 2 y, 3 z
   double precision :: peptides(100,4,3)=0d0
! number of rings and number of residues per ring
   integer rings,residues
! counters for different types of hydrogen bonds
   integer :: num_hb=0,num_anti=0,num_par=0,max_hb=0
contains

! Get the indices either side of residue ires. These are normally ires+1 and
! ires-1 unless ires is one of the residues that connects to close the ring.
   subroutine getneighbours(ires,i1,i2)
      implicit none
      integer, intent(in):: ires ! Residue index
      integer, intent(out):: i1,i2 ! Indices of neighbours
      i1=ires-1
      i2=ires+1
      if (mod(i1,residues).eq.0) then
         i1=i1+residues
      endif
      if (mod(ires,residues).eq.0) then
         i2=i2-residues
      endif
      return
   end subroutine getneighbours

! Get distance between two atoms
   function get_r(ires,iat,jres,jat)
      implicit none
      integer ires,jres ! Residue indices
      integer iat,jat ! Atom indices
      double precision get_r
      integer i
      get_r=0d0
      do i=1,3
         get_r=get_r+(peptides(ires,iat,i)-peptides(jres,jat,i))**2
      enddo
      get_r=dsqrt(get_r)
      return
      end

! Count intermolecular hydrogen bonds in a structure containing two or more cyclic peptides.
   subroutine count_hb()
      implicit none
      integer i,j,k ! Loop indices
      double precision ivec(3),jvec(3) ! Vectors used in parallel/antiparallel calculation
      integer i1,i2,j1,j2 ! Residue indices used in parallel/antiparallel calculation
      double precision e_dssp ! DSSP energy
      double precision dotp ! A dot product
      double precision dist ! Length of hydrogen bond

! Atom indices in peptides array
! 1 O, 2 C, 3 N, 4 H

      max_hb=residues*(rings-1)
! Loop through acceptor O atoms
      do i=1,rings*residues
! Loop through donor H atoms
         do j=1,rings*residues
! Check both atoms not in same ring
            if (((i-1)/residues).ne.((j-1)/residues)) then
               dist=get_r(i,1,j,4)
! DSSP energy calculation
               e_dssp=0.084*332*(1/get_r(i,1,j,3)+1/get_r(i,2,j,4)-1/get_r(i,1,j,4)-1/get_r(i,2,j,3))
               if (e_dssp.lt.-0.5) then! True if this is a hydrogen bond
                  num_hb=num_hb+1
! Check for parallel/antiparallel
                  call getneighbours(i,i1,i2)
                  call getneighbours(j,j1,j2)
                  dotp=0d0
                  do k=1,3
                     ivec(k)=peptides(i1,2,k)-peptides(i2,2,k)
                     jvec(k)=peptides(j1,2,k)-peptides(j2,2,k)
                     dotp=dotp+ivec(k)*jvec(k)
                  enddo
                  if (dotp.lt.0) then
                     num_anti=num_anti+1
                  else if (dotp.gt.0) then
                     num_par=num_par+1
                  endif
               endif
            endif
         enddo
      enddo
      return
   end subroutine count_hb

   function centroid(ring_no)
! Calculate the position of the centroid of a ring
! This calculation is based on the peptide N and C atoms
      integer ring_no
      integer i,j
      double precision centroid(3)
      centroid=0d0
      do j =1,3
         do i=residues*(ring_no-1)+1,residues*ring_no
            centroid(j)=centroid(j)+peptides(i,2,j) ! C atom
            centroid(j)=centroid(j)+peptides(i,3,j) ! N atom
         enddo
         centroid(j)=centroid(j)/(2*residues)
      enddo
      return
   end function centroid

   function centroid_distance(r1,r2)
! Calculate the distance between the centroids of two peptide rings
      integer r1,r2
      integer i
      double precision cent1(3), cent2(3)
      double precision centroid_distance
      !double precision centroid(3)
      centroid_distance=0d0
      cent1=centroid(r1)
      cent2=centroid(r2)
      do i=1,3
         centroid_distance=centroid_distance+(cent1(i)-cent2(i))**2
      enddo
      centroid_distance=dsqrt(centroid_distance)
      return
   end function centroid_distance

end module hb
