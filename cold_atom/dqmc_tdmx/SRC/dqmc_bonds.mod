	  `)  M   k820309    <
          11.1        �]L                                                                                                           
       dqmc_bonds.F90 DQMC_BONDS                                                    
                                                          
                         @                                '�                   #NSITES    #NATOM    #NCELL    #NDIM    #SC    #AC 	   #SCC 
   #POS    #CARTPOS    #XAT    #PHASE    #TRANSLATION    #NCLASS    #MYCLASS    #CLASS_SIZE    #CLASS_LABEL    #OLABEL    #INITIALIZED    #CONSTRUCTED    #ANALYZED                �                                                               �                                                              �                                                              �                                                              �                                    	                          p          p          p            p          p                                       �                              	     	       8                 
  p          p          p            p          p                                       �                              
     	       �                 
  p          p          p            p          p                                      �                                         �                 
            &                   &                                                       �                                         (             	   
            &                   &                                                       �                                         �             
   
            &                   &                                                       �                                         �                
            &                                                       �                                         0                
            &                   &                                                        �                                    �                        �                                          �                            &                   &                                                       �                                          �                            &                                                       �                                         @                
            &                   &                                           .            �                                          �                            &                                                                �                                    �                         �                                    �                         �                                    �                                                                                                               3%         @                                                         #MOVE_TO_RECORD%INDEX    #STRING    #IUNIT                                                     INDEX                                                                 1                                                                     @                                           
                                                            TWp          n                         �              C#NDIM       n                          �              C#MU         n                          �              C#PRIM       n                          �              C#SUPER      n                          �              C#ORB        n                          �              C#HAMILT     n                          �              C#K-POINT    n                          �              C#SYMM       n                          �              C#PHASE      n                          �              C#BONDS      n                          �              C#PAIR       n                          �              C#DILUT      h  p          p          p            p                                                                                                                                                                                                                                                                                                                                                                   
               10           @@                                                      @                                                          p          p            p                                                                       !                                                      11%         @                               "                         #LATTICE_T    #HOPTOWHO%SUM #   #HOPTOWHO%NINT $   #IAT %   #DELTA &   #JAT '   #LATTICE (                                               #     SUM                                             $     NINT           
                                  %                     
                                 &                   
 ,   p          p            p                                    
                                  '                     
                                  (     �             #LATTICE_T                      @                           )     '�             
      #NTOTBOND *   #BOND_LABEL +   #BOND_ORIGIN ,   #BOND_TARGET -   #XXBOND .   #NCLASS_B /   #MYCLASS_B 0   #CLASS_SIZE_B 1   #INITIALIZED 2   #ANALYZED 3                �                               *                               �                              +                                         &                                                       �                              ,            P                             &                                                       �                              -            �                             &                                                       �                             .            �                 
            &                   &                                                        �                               /     @                        �                              0            H                            &                   &                                                       �                              1            �                            &                                                        �                               2     �      	                   �                               3     �      
                        @                           4     '`                   #NWAVE 5   #NBOND 6   #NBONDV 7   #BOND_ORIGIN 8   #BOND_END 9   #BOND_MAP :   #PAIR_MAP ;   #BOND_NUMBER <   #BOND_WGT =   #WAVE_LABEL >   #MYCLASS_P ?   #NCLASS_P @   #CLASS_SIZE_P A   #INITIALIZED B                �                               5                                �                               6                              �                              7                                         &                                                       �                              8            P                             &                   &                                                       �                              9            �                             &                   &                                                       �                              :                                        &                                                       �                              ;            X                            &                                                       �                              <            �                            &                   &                                                       �                             =                          	   
            &                   &                                           .            �                              >            `             
               &                                                               �                              ?            �                            &                   &                                                        �                               @                             �                              A                                        &                                                        �                               B     X            #         @                                  C                  #READ_BONDS%SUM D   #BONDS E                                              D     SUM           D                                 E     �              #BONDS_T )   #         @                                  F                 #LATTICE_T    #CONSTRUCT_PAIRS%MOD G   #CONSTRUCT_PAIRS%SQRT H   #CONSTRUCT_PAIRS%SUM I   #BONDS J   #PAIRS K   #LATTICE L                                              G     MOD                                            H     SQRT                                            I     SUM                                            J     �              #BONDS_T )             D                                 K     `              #PAIRING 4              @                               L     �              #LATTICE_T       �   "      fn#fn     �   @   J   DQMC_GEOM_PARAM      @   J   DQMC_LATT $   B  F      LATTICE_T+DQMC_LATT +   �  H   a   LATTICE_T%NSITES+DQMC_LATT *   �  H   a   LATTICE_T%NATOM+DQMC_LATT *     H   a   LATTICE_T%NCELL+DQMC_LATT )   `  H   a   LATTICE_T%NDIM+DQMC_LATT '   �  �   a   LATTICE_T%SC+DQMC_LATT '   d  �   a   LATTICE_T%AC+DQMC_LATT (      �   a   LATTICE_T%SCC+DQMC_LATT (   �  �   a   LATTICE_T%POS+DQMC_LATT ,   �  �   a   LATTICE_T%CARTPOS+DQMC_LATT (   4  �   a   LATTICE_T%XAT+DQMC_LATT *   �  �   a   LATTICE_T%PHASE+DQMC_LATT 0   t  �   a   LATTICE_T%TRANSLATION+DQMC_LATT +    	  H   a   LATTICE_T%NCLASS+DQMC_LATT ,   h	  �   a   LATTICE_T%MYCLASS+DQMC_LATT /   
  �   a   LATTICE_T%CLASS_SIZE+DQMC_LATT 0   �
  �   a   LATTICE_T%CLASS_LABEL+DQMC_LATT +   T  �   a   LATTICE_T%OLABEL+DQMC_LATT 0   �  H   a   LATTICE_T%INITIALIZED+DQMC_LATT 0   8  H   a   LATTICE_T%CONSTRUCTED+DQMC_LATT -   �  H   a   LATTICE_T%ANALYZED+DQMC_LATT %   �  q       RDIM+DQMC_GEOM_PARAM /   9  �       MOVE_TO_RECORD+DQMC_GEOM_PARAM 5   �  >      MOVE_TO_RECORD%INDEX+DQMC_GEOM_PARAM 6   �  L   a   MOVE_TO_RECORD%STRING+DQMC_GEOM_PARAM 5   D  @   a   MOVE_TO_RECORD%IUNIT+DQMC_GEOM_PARAM -   �  �      INPUT_FIELDS+DQMC_GEOM_PARAM (   l  r       BONDS_F+DQMC_GEOM_PARAM (   �  @       INPUNIT+DQMC_GEOM_PARAM ,     �       FOUND_FIELD+DQMC_GEOM_PARAM (   �  r       PAIRS_F+DQMC_GEOM_PARAM #   $  �       HOPTOWHO+DQMC_LATT '   �  <      HOPTOWHO%SUM+DQMC_LATT (     =      HOPTOWHO%NINT+DQMC_LATT '   K  @   a   HOPTOWHO%IAT+DQMC_LATT )   �  �   a   HOPTOWHO%DELTA+DQMC_LATT '     @   a   HOPTOWHO%JAT+DQMC_LATT +   _  W   a   HOPTOWHO%LATTICE+DQMC_LATT    �  �       BONDS_T !   �  H   a   BONDS_T%NTOTBOND #   �  �   a   BONDS_T%BOND_LABEL $   |  �   a   BONDS_T%BOND_ORIGIN $     �   a   BONDS_T%BOND_TARGET    �  �   a   BONDS_T%XXBOND !   P  H   a   BONDS_T%NCLASS_B "   �  �   a   BONDS_T%MYCLASS_B %   D  �   a   BONDS_T%CLASS_SIZE_B $   �  H   a   BONDS_T%INITIALIZED !      H   a   BONDS_T%ANALYZED    h        PAIRING    �  H   a   PAIRING%NWAVE    �  H   a   PAIRING%NBOND      �   a   PAIRING%NBONDV $   �  �   a   PAIRING%BOND_ORIGIN !   T   �   a   PAIRING%BOND_END !    !  �   a   PAIRING%BOND_MAP !   �!  �   a   PAIRING%PAIR_MAP $   ("  �   a   PAIRING%BOND_NUMBER !   �"  �   a   PAIRING%BOND_WGT #   �#  �   a   PAIRING%WAVE_LABEL "   $  �   a   PAIRING%MYCLASS_P !   �$  H   a   PAIRING%NCLASS_P %   %  �   a   PAIRING%CLASS_SIZE_P $   �%  H   a   PAIRING%INITIALIZED    �%  g       READ_BONDS    S&  <      READ_BONDS%SUM !   �&  U   a   READ_BONDS%BONDS     �&  �       CONSTRUCT_PAIRS $   �'  <      CONSTRUCT_PAIRS%MOD %   �'  =      CONSTRUCT_PAIRS%SQRT $   #(  <      CONSTRUCT_PAIRS%SUM &   _(  U   a   CONSTRUCT_PAIRS%BONDS &   �(  U   a   CONSTRUCT_PAIRS%PAIRS (   	)  W   a   CONSTRUCT_PAIRS%LATTICE 