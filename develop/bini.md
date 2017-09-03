Computation of the Newton Polygon {#bini}
============================================================
\author Dario Andrea Bini,\n University of Pisa, Italy
\note E-mail: bini@dm.unipi.it
\version 1.4
\date June 1996 

\par License
\verbatim
All the software  contained in this library  is protected by copyright. 
Permission  to use, copy, modify, and  distribute this software for any 
purpose without fee is hereby granted, provided that this entire notice 
is included  in all copies  of any software which is or includes a copy 
or modification  of this software  and in all copies  of the supporting 
documentation for such software. 
  
THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT ANY EXPRESS OR IMPLIED 
WARRANTY. IN NO EVENT, NEITHER  THE AUTHORS, NOR THE PUBLISHER, NOR ANY 
MEMBER  OF THE EDITORIAL BOARD OF  THE JOURNAL  "NUMERICAL ALGORITHMS", 
NOR ITS EDITOR-IN-CHIEF, BE  LIABLE FOR ANY ERROR  IN THE SOFTWARE, ANY 
MISUSE  OF IT  OR ANY DAMAGE ARISING OUT OF ITS USE. THE ENTIRE RISK OF 
USING THE SOFTWARE LIES WITH THE PARTY DOING SO.
  
ANY USE  OF THE SOFTWARE  CONSTITUTES  ACCEPTANCE  OF THE TERMS  OF THE 
ABOVE STATEMENT.
\endverbatim

\par Subroutines
\verbatim
SUBROUTINES AND FUNCTIONS

 CNVEX  :  Computes the convex hull
 CMERGE :  Used by CNVEX 
 LEFT   :  Used by CMERGE
 RIGHT  :  Used by CMERGE
 CTEST  :  Convexity test, Used by CMERGE
\endverbatim

\see 
<a href="https://link.springer.com/article/10.1007/BF02207694">Numerical computation of polynomial zeros by means of Aberth's method</a>\n
Numerical Algorithms, 13 (1996), PP. 179-200

\remarks 
This version, which is compatible with Lahey's free ELF90 compiler
by Alan Miller, CSIRO Mathematical & Information Sciences,
Private Bag 10, Clayton South MDC, Victoria, Australia 3169.
Alan.Miller @ mel.dms.csiro.au http://www.mel.dms.csiro.au/~alan
Latest revision of ELF90 version - 5 May 1997 \n
Work performed under the support of the ESPRIT BRA project 6846 POSSO 



  
  
