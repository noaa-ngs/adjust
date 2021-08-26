C %P%

      integer function gfirst(i,nx)
********************************************************************************
*  return head of list pointer
********************************************************************************
      CHARACTER*99 SCCSID
      dimension nx(*)
      common/kathy/n,n1,n2,iemx,ifree
C      SCCSID='$Id: newhog.for 44197 2010-07-13 12:26:03Z bruce.tran $	20$Date: 2007/11/20 15:24:49 $ NGS'
      gfirst=nx(i)
      return 
      end
*
      integer function gnext(iptr,nx)
********************************************************************************
* return the item at this node and update the pointer
********************************************************************************
      dimension nx(*)
      integer gnext1
      common/kathy/n,n1,n2,iemx,ifree
      gnext=gnext1(iptr,nx(1),nx(n1),nx(n1+iemx))
      return
      end
*
      integer function gnext1(iptr,ihead,nbr,link)
********************************************************************************
* return the item at this node and update the pointer
********************************************************************************
      dimension ihead(*),nbr(*),link(*)
      common/kathy/n,n1,n2,iemx,ifree
      gnext1=nbr(iptr)
      iptr=link(iptr)
      return 
      end
*
      subroutine dumpt(nx,lunit)
      dimension nx(*)
      common/kathy/n,n1,n2,iemx,ifree
      call dumpt1(lunit,nx(1),nx(n1),nx(n1+iemx))
      return
      end
*
      subroutine dumpt1(lunit,ihead,nbr,link)
      dimension ihead(*),nbr(*),link(*)
      common/kathy/n,n1,n2,iemx,ifree
      do 100 i=1,n
        write(lunit,6010) i,ihead(i)
  100 continue
 6010 format(2i6)

      do 200 i=1,ifree
        write(lunit,6020) i,nbr(i),link(i)
  200 continue
 6020 format(3i6)
      return
      end
