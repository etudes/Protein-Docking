c
c ===========================================================================
c
      subroutine gkuser (usernm)
c
      implicit none
c
c --- GKUSER (...) => returns user's login name
c
c --- G J Kleywegt @ 920403
c
      character usernm*(*)
c
code ...
c
c Disabled MRH11 April 2011
cmrh      call getlog (usernm)
      usernm = "Anonymous"
c
      return
      end
