*-----------------------------------------------------------
* Title      :
* Written by :
* Date       :
* Description:
*-----------------------------------------------------------
    ORG    $1000
START:          

    LEA filename, A1
    MOVE.B #51, D0
    TRAP #15
    
    MOVE.L D3,D2
    MOVE.B #55,D0
    TRAP #15
    
    LEA buffer,A1
    MOVE.W #datasize,D2
    MOVE.B #53,D0
    TRAP #15
    MOVEA.L A1,A4
    MOVE.B D0,D7
     
loop1 CMP #$F0,D6
      BEQ done
      JSR convert
      MOVE.B (A1)+,D6
      BRA loop1

done MOVEA.L A4,A1
    CLR D1
    
loop2 MOVE.B #16,D2
    MOVE.B (A1)+,D1
    MOVE.B #15,D0
    TRAP #15
    MOVE.B (A1),D1
    CMP #$0D,D1
    BNE loop2
        
    MOVE.B #56, D0
    TRAP #15
        
    SIMHALT

* Put variables and constants here
filename dc.b 'test1.txt', 0
datasize equ 2000
buffer ds.b 2000

convert MOVE.B (A1),D5
       CMP #$40,D5
       BGT let
       SUBI.L #$30,D5
       BRA exit
let    SUBI.L #$37,D5
exit   MOVE.B D5,(A1)
       RTS

    END    START        ; last line of source




*~Font name~Courier New~
*~Font size~10~
*~Tab type~1~
*~Tab size~4~
