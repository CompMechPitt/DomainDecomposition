
    F   k820309    P          16.0        ­tMW                                                                                                           
       ../Sources/precond.f90 PRECOND_CLASS                                                     
                                                           
                                                           
                         @              @                'X                    #N    #DATA    #DESC                                                                                                                                         
            &                                                                                           4       P              #DATA_DESC                      @                               '4                    #NX_G 	   #NY_G 
   #NX_L    #NY_L    #NODE_NUM    #ELEM_NUM    #ME    #NP    #UP_PID    #UP_OFF    #DOWN_PID    #DOWN_OFF    #MPI_COMM                                                	                                                               
                                                                                                                                                                                                                                                                                                                                                                                                                                                               	                                                       $       
                                                       (                                                              ,                                                              0                               @               @                'p                    #M    #N    #DATA    #DESC                                                                                                                                                                                                       
            &                   &                                                                                           4       h              #DATA_DESC    #         @                                                      #DESC    #V              
                                       4              #DATA_DESC                                                    X               #VECTOR    #         @                                                      #X    #Y               
                                       X              #VECTOR              
                                       X               #VECTOR    #         @                                  !                    #X "             
                                 "     X               #VECTOR    #         @                                  #                    #V $             
                                 $     X               #VECTOR                      @               À           %     'È                    #DESC &   #TYPE '   #INV_DIAG (   #CHOL )                                              &     4                      #DATA_DESC                                                '                                                              (     X                     #VECTOR                                             )            h                 
            &                   &                                                                                        *                                                       0                                             +                                                      1                                             ,                                                      2#         @                                  -     	               #UPLO .   #N /   #A 0   #LDA 1   #INFO 2             
                                 .                                     
                                 /                  B  
                               0                    
       p        5 O p        p          5 O p          1     5 O p          1                             
                                 1                     
                                 2           #         @                                  3     	               #UPLO 4   #N 5   #NRHS 6   #A 7   #LDA 8   #B 9   #LDB :   #INFO ;             
                                 4                                     
                                 5                     
                                 6                  B  
                                7                    
      p        5 O p        p          5 O p          1     5 O p          1                             
                                 8                  B  
                               9                    
       p        5 O p        p          5 O p          1     5 O p          1                             
                                 :                     
                                 ;           #         @                                   <                    #A =   #TYPE >   #PREC ?             
                                  =     p              #MATRIX              
                                  >                     D @                               ?     È               #PRECOND %   #         @                                   @                    #M A   #R B   #Z C             
                                  A     È              #PRECOND %             
  @                               B     X              #VECTOR              
D @                               C     X               #VECTOR    #         @                                   D                    #PREC E             
D @                               E     È               #PRECOND %          -      fn#fn     Í   @   J   DATA_DESC_CLASS      @   J   MATRIX_CLASS    M  @   J   VECTOR_CLASS $     k       VECTOR+VECTOR_CLASS &   ø  H   a   VECTOR%N+VECTOR_CLASS )   @     a   VECTOR%DATA+VECTOR_CLASS )   Ô  _   a   VECTOR%DESC+VECTOR_CLASS *   3  æ       DATA_DESC+DATA_DESC_CLASS /     H   a   DATA_DESC%NX_G+DATA_DESC_CLASS /   a  H   a   DATA_DESC%NY_G+DATA_DESC_CLASS /   ©  H   a   DATA_DESC%NX_L+DATA_DESC_CLASS /   ñ  H   a   DATA_DESC%NY_L+DATA_DESC_CLASS 3   9  H   a   DATA_DESC%NODE_NUM+DATA_DESC_CLASS 3     H   a   DATA_DESC%ELEM_NUM+DATA_DESC_CLASS -   É  H   a   DATA_DESC%ME+DATA_DESC_CLASS -     H   a   DATA_DESC%NP+DATA_DESC_CLASS 1   Y  H   a   DATA_DESC%UP_PID+DATA_DESC_CLASS 1   ¡  H   a   DATA_DESC%UP_OFF+DATA_DESC_CLASS 3   é  H   a   DATA_DESC%DOWN_PID+DATA_DESC_CLASS 3   1  H   a   DATA_DESC%DOWN_OFF+DATA_DESC_CLASS 3   y  H   a   DATA_DESC%MPI_COMM+DATA_DESC_CLASS $   Á  r       MATRIX+MATRIX_CLASS &   3  H   a   MATRIX%M+MATRIX_CLASS &   {  H   a   MATRIX%N+MATRIX_CLASS )   Ã  ¬   a   MATRIX%DATA+MATRIX_CLASS )   o	  _   a   MATRIX%DESC+MATRIX_CLASS +   Î	  Y       VECTOR_CREATE+VECTOR_CLASS 0   '
  W   a   VECTOR_CREATE%DESC+VECTOR_CLASS -   ~
  T   a   VECTOR_CREATE%V+VECTOR_CLASS )   Ò
  V       VECTOR_COMM+VECTOR_CLASS +   (  T   a   VECTOR_COMM%X+VECTOR_CLASS +   |  T   a   VECTOR_COMM%Y+VECTOR_CLASS +   Ð  O       VECTOR_WEIGHT+VECTOR_CLASS -     T   a   VECTOR_WEIGHT%X+VECTOR_CLASS ,   s  O       VECTOR_DESTROY+VECTOR_CLASS .   Â  T   a   VECTOR_DESTROY%V+VECTOR_CLASS      |       PRECOND      _   a   PRECOND%DESC    ñ  H   a   PRECOND%TYPE !   9  \   a   PRECOND%INV_DIAG      ¬   a   PRECOND%CHOL    A  q       IDENTITY    ²  q       DIAGONAL    #  q       NEUMANN      s       DPOTRF      P   a   DPOTRF%UPLO    W  @   a   DPOTRF%N      Ü   a   DPOTRF%A    s  @   a   DPOTRF%LDA    ³  @   a   DPOTRF%INFO    ó         DPOTRS      P   a   DPOTRS%UPLO    Ð  @   a   DPOTRS%N      @   a   DPOTRS%NRHS    P  Ü   a   DPOTRS%A    ,  @   a   DPOTRS%LDA    l  Ü   a   DPOTRS%B    H  @   a   DPOTRS%LDB      @   a   DPOTRS%INFO    È  c       PRECOND_CREATE !   +  T   a   PRECOND_CREATE%A $     @   a   PRECOND_CREATE%TYPE $   ¿  U   a   PRECOND_CREATE%PREC      ]       PRECOND_APPLY     q  U   a   PRECOND_APPLY%M     Æ  T   a   PRECOND_APPLY%R       T   a   PRECOND_APPLY%Z     n  R       PRECOND_DESTROY %   À  U   a   PRECOND_DESTROY%PREC 