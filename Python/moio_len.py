
mo_distributons = {
              "A<B"  = "1"
              "a<b"  = "2"
              "I<J"  = "3"
              "A<=B" = "4"
              "a<=b" = "5"
              "I<=J" = "6"
              "i<=j" = "7"
              "AI"   = "8"
              "ai"   = "9"
              "Ai"   = "10"
              "aI"   = "11"
              "Ab"   = "12"
              "Ij"   = "13"
              "ia"   = "14"
              "IA"   = "15"
              "ia"   = "16"
              "Ia"   = "17"
              "AB"   = "18"
              "ab"   = "19"
              "IJ"   = "20"
              "ij"   = "21"
}

mo_iolists = {
             "W<AB|IJ>:AI,BJ" = 5
             "W<ab|ij>:ai,bj" = 6

             "W<IJ||KL>:I<J,K<L" = 11
             "W<ij||kl>:i<j,k<l" = 12
             "W<iJ||kL>:ij,kl"   = 13
               
             "W<IJ||KA>:I<J,KA" = 7
             "W<ij||ka>:i<j,ka" = 8
             "W<iJ||Ak>:Ij,Ak"  = 9
             "W<iJ||Ka>:Ij,Ka"  = 10 

             "W<AB||IJ>:A<B,I<J" = 14
	     "W<ab||ij>:a<b,i<j" = 15
             "W<Ab||Ij>:ab,Ij"   = 16
             "W<Ab||Ij>:bj,AI" = 17
             "W<Ab||Ij>:AI,bj" = 18
             "W<AB||IJ>:AI,BJ" = 19
             "W<ab||ij>:ai,bj" = 20
             "W<Ab||Ij>:Aj,bI" = 21
             "W<Ab||Ij>:bI,Aj" = 22

             "W<IA||JB>:BI,AJ" = 23 
             "W<ia||jb>:bi,aj" = 24
             "W<iA||jB>:Bi,Aj" = 25 
             "W<aI||bJ>:bI,aJ" = 26

             "W<AB||CI>:A<B,CI" = 27
             "W<ab||ci>:a<b,ci" = 28
             "W<Ab||Ic>:Ab,Ic"  = 29
             "W<Ab||Ci>:Ab,Ci"  = 30
            
             "T2(IJ,AB):AJ,BI"  = 34
             "T2(ij,ab):aj,bi"  = 35
             "T2(Ij,Ab):bj,AI"  = 36
             "T2(Ij,Ab):AI,bj"  = 37
             "T2(Ij,Ab):bI,Aj"  = 38
             "T2(Ij,Ab):Aj,bI"  = 39

             "T2(IJ,AB):A<B,I<J"  = 44
             "T2(ij,ab):i<j,a<b"  = 45
             "T2(Ij,Ab):Ab,Ij"    = 46

             "T1(AI):AI"            = 901 
             "T1(ai):Ai"            = 902 
             "F(MI):MI"             = 913
             "F(mi):mi"             = 914
             "F(EA):EA"             = 923
             "F(ea):ea"             = 924 
             "F(AI):AI"             = 933
             "F(ai):ai"             = 934 


             "HB(MI):MI"             = 911
             "HB(MI):MI"             = 912
             "HB(AE):AE"             = 921
             "HB(ae):ae"             = 922
             "HB(IA):IA"             = 931
             "HB(Ia):ia"             = 932

             "HB(MN,IJ):M<N,I<J"      = 51
             "HB(mn,ij):m<n,i<j"      = 52
             "HB(Mn,Ij):M<n,I<j"      = 53

             "HB(AI,BC):AI,B<C"      = 27
             "HB(ai,bc):ai,b<c"      = 28
             "HB(Ai,Bc):Ai,Bc"       =  29
             "HB(Ai,Bc):Ai,cB"       =  30

             "HB(IJ,KA):I<J,KA"     = 7
             "HB(ij,ka):i<j,ka"     = 8
             "HB(Ij,Ka):Ij,Ka"      = 9
             "HB(Ij,Ak):IjkA"       = 10

             "HB(MB,EJ):MBEJ"       = 54
             "HB(mb,ej):mbej"       = 55
             "HB(Mb,Ej):MbEj"       = 56
             "HB(Mb,Ej):MbEj"       = 57
             "HB(Mb,Ej):MbEj"       = 58
             "HB(Mb,Ej):MbEj"       = 59 

             "HB(AB,CI):A<B,CI"      = 127
             "HB(ab,ci):a<b,CI"      = 128
             "HB(Ab,Ci):Ab,CI"      =  129
             "HB(Ab,Ic):Ab,CI"      =  130


             "HB(IA,JK):IA,J<K"    = 107
             "HB(ia,jk):ia,j<k"    = 108
             "HB(Ia,Jk):Ia,Jk"     = 109
             "HB(Ia,Kj):Ia,jK"     = 110

             "HB(AB,CD):A<B,C<D"   = 231
             "HB(ab,cd):a<b,c<d"   = 232
             "HB(Ab,Cd):Ab,cd"    = 233

             "L1AI:AI"          = 1901 
             "L1ai:Ai"          = 1902 
             "L2(IJ,AB):AJ,BI"  = 134
             "L2(ij,ab):aj,bi"  = 135
             "L2(Ij,Ab):bj,AI"  = 136
             "L2(Ij,Ab):AI,bj"  = 137
             "L2(Ij,Ab):bI,Aj"  = 138
             "L2(Ij,Ab):Aj,bI"  = 139

             "L2(IJ,AB):A<B,I<J"  = 144
             "L2(ij,ab):i<j,a<b"  = 145
             "L2(Ij,Ab):Ab,Ij"    = 146

              "R2(AB,IJ):A<B,I<J" = 444
              "R2(AB,IJ):a<b,i<j" = 445
              "R2(Ab,Ij):Ab,Ij"   = 446
  
             "R1AI:AI"            = 4901
             "E1ai:Ai"            = 4902
}



