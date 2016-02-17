! Analyse hydrogen-bonding pattern in cyclic peptide nanotubes

module hb
implicit none
! Atom indices in peptides array
! 1 O, 2 C, 3 N, 4 H
double precision :: peptides(100,4,3)=0d0 !Dynamically allocate this
integer residues,rings
integer :: num_hb=0,num_anti=0,num_par=0
end module hb

subroutine getneighbours(ires,residues,i1,i2)
implicit none
integer ires
integer i1,i2
integer residues
i1=ires-1
i2=ires+1
if (mod(i1,residues).eq.0) then
   i1=i1+residues
endif
if (mod(ires,residues).eq.0) then
   i2=i2-residues
endif
return
end

function get_r(ires,iat,jres,jat)
use hb
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

subroutine count_hb()
use hb
implicit none
integer i,j,k
double precision ivec(3),jvec(3)
double precision get_r
integer i1,i2,j1,j2
double precision e_dssp
double precision dotp
double precision dist

! Atom indices in peptides array
! 1 O, 2 C, 3 N, 4 H

do i=1,rings*residues
   do j=1,rings*residues
      if (((i-1)/residues).ne.((j-1)/residues)) then
         dist=get_r(i,1,j,4)
         e_dssp=0.084*332*(1/get_r(i,1,j,3)+1/get_r(i,2,j,4)-1/get_r(i,1,j,4)-1/get_r(i,2,j,3))
!         write(*,*) dist,e_dssp
         if (e_dssp.lt.-0.5) then
            num_hb=num_hb+1
            call getneighbours(i,residues,i1,i2)
            call getneighbours(j,residues,j1,j2)
            dotp=0d0
            do k=1,3
               ivec(k)=peptides(i1,2,k)-peptides(i2,2,k)
               jvec(k)=peptides(j1,2,k)-peptides(j2,2,k)
               dotp=dotp+ivec(k)*jvec(k)
            enddo
!            write(*,*) dotp
            if (dotp.lt.0) then
               num_anti=num_anti+1
            else if (dotp.gt.0) then
               num_par=num_par+1
            endif
         endif
            !write(*,*) i1,i,i2
      endif
   enddo
enddo
return
end
