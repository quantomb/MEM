	  A  �   k820309    <
          11.1        �]L                                                                                                           
       dqmc_util.F90 DQMC_UTIL                                                    
                        �                                 
       #         @                                     	              #DAXPY%KIND    #N    #A    #X    #INCX    #Y 	   #INCY 
                                                    KIND                                                                                                     
              @                                                     
     p          1     1                                                                                @                                 	                    
     p          1     1                                                             
            #         @                                     	              #DSCAL%KIND    #N    #DA    #DX    #INCX                                                     KIND                                                                                                     
              @                                                     
     p          1     1                                                                                                                                                                                                                                 KIND                                                  
                 
                                 0.0D0                                                 
                 
                       �?        1.0D0                                                 
                 
                        @        2.0D0                                                 
                 
                       �?        0.5D0                                                                                                    0                                                                                                   6                                                                                                   5                                                 
                                                   C(a30, i12)                                                                                                                    C(a30, f19.6)                                                                                                                    C(a30, f12.6,' +- ',f12.6)                                                                                                                    C(i3,i3)                                                                 	                            
                       C(76('='))                                                                 	                            
                       C(76('-'))                %         @                                                     
       #DQMC_MATDIFF%SQRT !   #N "   #A #   #B $                                              !     SQRT           
@ @                               "                    
@ @                              #                    
      p        5 � p        r "   p          5 � p        r "     5 � p        r "       5 � p        r "     5 � p        r "                              
@ @                              $                    
      p        5 � p        r "   p          5 � p        r "     5 � p        r "       5 � p        r "     5 � p        r "                     %         @                                 %                   
       #DQMC_MATNORM%ABS &   #DQMC_MATNORM%SQRT '   #N (   #A )                                              &     ABS                                            '     SQRT           
                                  (                    
                                 )                    
      p        5 � p        r (   p          5 � p        r (     5 � p        r (       5 � p        r (     5 � p        r (                     #         @                                  *                   #N +   #A ,             
                                  +                    
D                                ,                    
       p        5 � p        r +   p          5 � p        r +     5 � p        r +       5 � p        r +     5 � p        r +                     #         @                                  -                   #N .   #AT /   #A 0             
                                  .                    
D                                /                    
       p        5 � p        r .   p          5 � p        r .     5 � p        r .       5 � p        r .     5 � p        r .                              
                                 0                    
      p        5 � p        r .   p          5 � p        r .     5 � p        r .       5 � p        r .     5 � p        r .                     #         @                                  1                   #N 2   #A 3   #D 4             
@ @                               2                    
D @                              3                    
       p        5 � p        r 2   p          5 � p        r 2     5 � p        r 2       5 � p        r 2     5 � p        r 2                              
@ @                              4                    
    p          5 � p        r 2       5 � p        r 2                     #         @                                  5                   #N 6   #A 7   #D 8             
@ @                               6                    
D @                              7                    
 	      p        5 � p        r 6   p          5 � p        r 6     5 � p        r 6       5 � p        r 6     5 � p        r 6                              
@ @                              8                    
 
   p          5 � p        r 6       5 � p        r 6                     #         @                                  9                   #N :   #A ;   #D <             
@ @                               :                    
D @                              ;                    
       p        5 � p        r :   p          5 � p        r :     5 � p        r :       5 � p        r :     5 � p        r :                              
                                 <                    
    p          5 � p        r :       5 � p        r :                     #         @                                  =                   #N >   #A ?   #D @             
@ @                               >                    
D @                              ?                    
       p        5 � p        r >   p          5 � p        r >     5 � p        r >       5 � p        r >     5 � p        r >                              
                                 @                    
    p          5 � p        r >       5 � p        r >                     #         @                                  A                  #DQMC_SIGNJACKKNIFE%SUM B   #DQMC_SIGNJACKKNIFE%ABS C   #DQMC_SIGNJACKKNIFE%SQRT D   #N E   #AVG F   #ERR G   #X H   #Y I   #SGN J   #SUM_SGN K                                              B     SUM                                            C     ABS                                            D     SQRT           
                                  E                     D @                              F     
                 D @                              G     
                 
  @                              H                   
              &                                                     
D @                              I                   
               &                                                     
                                 J                   
              &                                                     
                                 K     
      #         @                                  L                  #DQMC_JACKKNIFE%SUM M   #DQMC_JACKKNIFE%ABS N   #DQMC_JACKKNIFE%SQRT O   #N P   #AVG Q   #ERR R   #X S   #Y T   #SGN U   #SUM_SGN V                                              M     SUM                                            N     ABS                                            O     SQRT           
                                  P                     D @                              Q     
                 D @                              R     
                
  @                              S                    
    p          5 � p        r P       5 � p        r P                              
D @                              T                    
     p          5 � p        r P       5 � p        r P                              
D                                U                    
     p          5 � p        r P       5 � p        r P                               
D                                V     
       #         @                                  W                  #DQMC_GETERR%SUM X   #DQMC_GETERR%ABS Y   #DQMC_GETERR%SQRT Z   #N [   #ERR \   #AVG ]   #LIST ^                                              X     SUM                                            Y     ABS                                            Z     SQRT           
                                  [                     D @                              \     
                 
  @                              ]     
               
D @                              ^                    
     p          5 � p        r [       5 � p        r [                     #         @                                  _                  #DQMC_GETERR1%SUM `   #DQMC_GETERR1%ABS a   #DQMC_GETERR1%SQRT b   #N c   #DATA d   #AVG e   #ERR f                                              `     SUM                                            a     ABS                                            b     SQRT           
                                  c                    
  @                              d                    
    p          5 � p        r c       5 � p        r c                               D @                              e     
                 D                                f     
       #         @                                  g                  #DQMC_GETERR2%ABS h   #DQMC_GETERR2%SQRT i   #N j   #SM k   #SSQ l   #AVG m   #ERR n                                              h     ABS                                            i     SQRT           
                                  j                     
                                 k     
                
                                 l     
                D @                              m     
                 D                                n     
       #         @                                  o                   #MESSAGE p   #NO q             
                                 p                    1           
                                  q           #         @                                  r                   #MESSAGE s   #NO t             
                                 s                    1           
                                  t           #         @                                  u                   #N v   #VAR w   #SEED x             
                                  v                    D                                w                    
     p          5 � p        r v       5 � p        r v                               
                                 x                        p          p            p                          #         @                                  y                   #N z   #VAR {   #SEED |             
                                  z                    D                                {                    
     p          5 � p        r z       5 � p        r z                               
                                 |                        p          p            p                          #         @                                  }                   #M ~   #N    #A �   #OPT �             
                                  ~                     
                                                      
      �                           �                    
      p        5 � p        r ~   p          & p        5 � p        r ~     & p        5 � p        r        5 � p        r ~     5 � p        r                                
                                  �           #         @                                  �                   #N �   #M �   #TITLE �   #LABEL �   #AVG �   #ERR �   #OPT �             
                                  �                     
                                  �                     
                                 �                    1 ,          
                                 �                                  &                                           1           
                                 �                   
              &                   &                                                     
                                 �                   
              &                   &                                                     
                                  �           #         @                                  �                   #N �   #M �   #TITLE �   #LABEL �   #AVG �   #ERR �   #OPT �             
                                  �                     
                                  �                     
                                 �                    1 ,          
                                 �                                   &                                           1           
                                 �                   
 !             &                   &                                                     
                                 �                   
 "             &                   &                                                     
                                  �              �          fn#fn    �   @   J   LAPACK_MOD       @   J   BLAS_MOD    @  �       DAXPY+BLAS_MOD $   �  =      DAXPY%KIND+BLAS_MOD !     @   a   DAXPY%N+BLAS_MOD !   E  @   a   DAXPY%A+BLAS_MOD !   �  �   a   DAXPY%X+BLAS_MOD $   	  @   a   DAXPY%INCX+BLAS_MOD !   I  �   a   DAXPY%Y+BLAS_MOD $   �  @   a   DAXPY%INCY+BLAS_MOD      y       DSCAL+BLAS_MOD $   �  =      DSCAL%KIND+BLAS_MOD !   �  @   a   DSCAL%N+BLAS_MOD "     @   a   DSCAL%DA+BLAS_MOD "   C  �   a   DSCAL%DX+BLAS_MOD $   �  @   a   DSCAL%INCX+BLAS_MOD      p       WP    w  =       KIND    �  u       ZERO    )  u       ONE    �  u       TWO      u       HALF    �  q       STDERR    �  q       STDOUT    j	  q       STDIN    �	  �       FMT_STRINT    f
  �       FMT_STRDBL    �
  �       FMT_VALERR    �  �       FMT_INTPAR      �       FMT_DBLINE    �  �       FMT_SGLINE    )  |       DQMC_MATDIFF "   �  =      DQMC_MATDIFF%SQRT    �  @   a   DQMC_MATDIFF%N    "  $  a   DQMC_MATDIFF%A    F  $  a   DQMC_MATDIFF%B    j  �       DQMC_MATNORM !   �  <      DQMC_MATNORM%ABS "   1  =      DQMC_MATNORM%SQRT    n  @   a   DQMC_MATNORM%N    �  $  a   DQMC_MATNORM%A    �  V       DQMC_EYE    (  @   a   DQMC_EYE%N    h  $  a   DQMC_EYE%A    �  ^       DQMC_TRANS    �  @   a   DQMC_TRANS%N    *  $  a   DQMC_TRANS%AT    N  $  a   DQMC_TRANS%A    r  ]       DQMC_SCALECOL     �  @   a   DQMC_SCALECOL%N       $  a   DQMC_SCALECOL%A     3  �   a   DQMC_SCALECOL%D    �  ]       DQMC_SCALEROW     D  @   a   DQMC_SCALEROW%N     �  $  a   DQMC_SCALEROW%A     �  �   a   DQMC_SCALEROW%D !   \  ]       DQMC_SCALECOLINV #   �  @   a   DQMC_SCALECOLINV%N #   �  $  a   DQMC_SCALECOLINV%A #     �   a   DQMC_SCALECOLINV%D !   �  ]       DQMC_SCALEROWINV #   .  @   a   DQMC_SCALEROWINV%N #   n  $  a   DQMC_SCALEROWINV%A #   �   �   a   DQMC_SCALEROWINV%D #   F!  �       DQMC_SIGNJACKKNIFE '    "  <      DQMC_SIGNJACKKNIFE%SUM '   \"  <      DQMC_SIGNJACKKNIFE%ABS (   �"  =      DQMC_SIGNJACKKNIFE%SQRT %   �"  @   a   DQMC_SIGNJACKKNIFE%N '   #  @   a   DQMC_SIGNJACKKNIFE%AVG '   U#  @   a   DQMC_SIGNJACKKNIFE%ERR %   �#  �   a   DQMC_SIGNJACKKNIFE%X %   !$  �   a   DQMC_SIGNJACKKNIFE%Y '   �$  �   a   DQMC_SIGNJACKKNIFE%SGN +   9%  @   a   DQMC_SIGNJACKKNIFE%SUM_SGN    y%  �       DQMC_JACKKNIFE #   G&  <      DQMC_JACKKNIFE%SUM #   �&  <      DQMC_JACKKNIFE%ABS $   �&  =      DQMC_JACKKNIFE%SQRT !   �&  @   a   DQMC_JACKKNIFE%N #   <'  @   a   DQMC_JACKKNIFE%AVG #   |'  @   a   DQMC_JACKKNIFE%ERR !   �'  �   a   DQMC_JACKKNIFE%X !   p(  �   a   DQMC_JACKKNIFE%Y #   $)  �   a   DQMC_JACKKNIFE%SGN '   �)  @   a   DQMC_JACKKNIFE%SUM_SGN    *  �       DQMC_GETERR     �*  <      DQMC_GETERR%SUM     �*  <      DQMC_GETERR%ABS !   ;+  =      DQMC_GETERR%SQRT    x+  @   a   DQMC_GETERR%N     �+  @   a   DQMC_GETERR%ERR     �+  @   a   DQMC_GETERR%AVG !   8,  �   a   DQMC_GETERR%LIST    �,  �       DQMC_GETERR1 !   �-  <      DQMC_GETERR1%SUM !   �-  <      DQMC_GETERR1%ABS "   .  =      DQMC_GETERR1%SQRT    O.  @   a   DQMC_GETERR1%N "   �.  �   a   DQMC_GETERR1%DATA !   C/  @   a   DQMC_GETERR1%AVG !   �/  @   a   DQMC_GETERR1%ERR    �/  �       DQMC_GETERR2 !   b0  <      DQMC_GETERR2%ABS "   �0  =      DQMC_GETERR2%SQRT    �0  @   a   DQMC_GETERR2%N     1  @   a   DQMC_GETERR2%SM !   [1  @   a   DQMC_GETERR2%SSQ !   �1  @   a   DQMC_GETERR2%AVG !   �1  @   a   DQMC_GETERR2%ERR    2  ]       DQMC_ERROR #   x2  L   a   DQMC_ERROR%MESSAGE    �2  @   a   DQMC_ERROR%NO    3  ]       DQMC_WARNING %   a3  L   a   DQMC_WARNING%MESSAGE     �3  @   a   DQMC_WARNING%NO    �3  b       RAN0    O4  @   a   RAN0%N    �4  �   a   RAN0%VAR    C5  �   a   RAN0%SEED    �5  b       RAN1    96  @   a   RAN1%N    y6  �   a   RAN1%VAR    -7  �   a   RAN1%SEED    �7  f       DUMPA    '8  @   a   DUMPA%M    g8  @   a   DUMPA%N    �8  D  a   DUMPA%A    �9  @   a   DUMPA%OPT %   +:  �       DQMC_PRINT_REALARRAY '   �:  @   a   DQMC_PRINT_REALARRAY%N '   �:  @   a   DQMC_PRINT_REALARRAY%M +   2;  L   a   DQMC_PRINT_REALARRAY%TITLE +   ~;  �   a   DQMC_PRINT_REALARRAY%LABEL )   <  �   a   DQMC_PRINT_REALARRAY%AVG )   �<  �   a   DQMC_PRINT_REALARRAY%ERR )   V=  @   a   DQMC_PRINT_REALARRAY%OPT (   �=  �       DQMC_PRINT_COMPLEXARRAY *   >  @   a   DQMC_PRINT_COMPLEXARRAY%N *   ]>  @   a   DQMC_PRINT_COMPLEXARRAY%M .   �>  L   a   DQMC_PRINT_COMPLEXARRAY%TITLE .   �>  �   a   DQMC_PRINT_COMPLEXARRAY%LABEL ,   y?  �   a   DQMC_PRINT_COMPLEXARRAY%AVG ,   @  �   a   DQMC_PRINT_COMPLEXARRAY%ERR ,   �@  @   a   DQMC_PRINT_COMPLEXARRAY%OPT 