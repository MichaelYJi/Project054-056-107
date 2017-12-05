*-----------------------------------------------------------
* Title      :
* Written by :
* Date       :
* Description:
*-----------------------------------------------------------
    ORG    $1000
START:          

    LEA filename,A1
    MOVE.B #51,D0
    TRAP #15
    
    MOVE.L D3,D2
    MOVE.B #55,D0
    TRAP #15
    
    LEA buffer,A1
    MOVE.W #datasize,D2
    MOVE.B #53,D0
    TRAP #15
    MOVEA.L A1,A3
    MOVE.B D0,D7
     
loop1 MOVE.B (A1),D3
      CMP #$FF,D3
      BEQ done
      JSR convert
      ADDQ #1,A1
      BRA loop1

done CLR D3
    MOVE.B (A3)+,D3
    
loop2
    CLR D5
    JSR checkop
    CMP #$CC,D5 *CC represents invalid opcode
    BEQ clr
    JSR printdata 
    JSR print
clr CMP #$DD,D3 *DD represents new line
    BNE continue
    MOVE.B (A3)+,D3
    BRA clr
continue CMP #$FF,D3 *FF represents end of stream
    BNE loop2
        
    SIMHALT

* Put variables and constants here
filename dc.b 'test1.txt', 0
datasize equ 2000
buffer ds.b 2000
data dc.b 'DATA ',0
source dc.b 0
dest dc.b 0
sourcem dc.b 0
destm dc.b 0

pmove_b dc.b 'MOVE.B ',0
pmove_w dc.b 'MOVE.W ',0
pmove_l dc.b 'MOVE.L ',0
pmovea_w dc.b 'MOVEA.W ',0
pmovea_l dc.b 'MOVEA.L ',0
pmoveq dc.b 'MOVEQ ',0
pnop dc.b 'NOP ',0
padd dc.b 'ADD',0
padda dc.b 'ADDA',0
paddi dc.b 'ADDI',0
prts dc.b 'RTS',0

pas dc.b 'AS',0
pls dc.b 'LS',0
pro dc.b 'RO',0
pright dc.b 'R',0
pleft dc.b 'L',0
space dc.b ' ',0

dregister dc.b 'D',0
aregister dc.b 'A',0
hashtag dc.b '#',0
dollarsign dc.b '$',0
comma dc.b ',',0
dotb dc.b '.B ',0
dotw dc.b '.W ',0
dotl dc.b '.L ',0
leftp dc.b '(',0
rightp dc.b ')',0
plus dc.b   '+',0
minus dc.b '-',0

error dc.b 'error',0

convert CMP #$FF,D3
       BEQ exit
       CMP #$30,D3
       BLT fake
       CMP #$46,D3
       BGT fake
       CMP #$39,D3
       BLE not
       CMP #$41,D3
       BGE not
       BRA fake
fake   MOVE.B #$DD,D3
       BRA exit
not    CMP #$40,D3
       BGT let
       SUBI.L #$30,D3
       BRA exit
let    SUBI.L #$37,D3
exit   MOVE.B D3,(A1)
       RTS
       
newline MOVE #$0D,D1
        MOVE #6,D0
        TRAP #15
        MOVE #$0A,D1
        MOVE #6,D0
        TRAP #15
        RTS

print        
loop    CMP #$DD,D3
        BEQ done1
        CMP #$FF,D3
        BEQ done2
        MOVE.B D3,D1
        MOVE.B #16,D2
        MOVE.B #15,D0
        TRAP #15
        MOVE.B (A3)+,D3
        BRA loop
done1   JSR newline
done2   RTS

        
printdata LEA data,A1
        MOVE.B #14,D0
        TRAP #15
        RTS 
        
checkop MOVEA.W A3,A4

        MOVE.B -(A4),D4
        
        CMP #$00,D4     *0000
        BEQ ADDI_
        
        CMP #$01,D4     *0001
        BEQ MOVE_B
        
        CMP #$02,D4     *0010
        BEQ MOVE_L
        
        CMP #$03,D4     *0011
        BEQ MOVE_W
        
        CMP #$04,D4     *0100
        BEQ code0100
        
        CMP #$07,D4     *0111
        BEQ MOVEQ_
        
        CMP #$0D,D4     *1101
        BEQ ADD_
        
        CMP #$0E,D4     *1110
        BEQ rl_check
        
        *First level jump table goes here
        RTS
        
code0100 

        *NOP check
        MOVEA.W A3,A4
        ADDQ #2,A4
        MOVE.B (A4),D4
        CMP #$01,D4
        BEQ NOP_   
        CMP #$5,D4
        BEQ RTS_
        RTS

ADDI_   JSR IsolateRegisters
        
        *Isolate size
        MOVEA.W A3,A4
        *ADDQ #1,A4
        MOVE.B (A4)+,D5
        CMP #$6,D5
        BNE addiexit
        MOVE.B (A4),D5
        MOVE.B #28,D7
        LSL.L D7,D5
        ADDQ #2,D7
        LSR.L D7,D5
        
        *Check size
        CMP #$0,D5
        BNE addin1
        LEA dotb,A5
        MOVE.B #2,D4
        BRA prnaddi
addin1  CMP #$1,D5
        BNE addin2
        LEA dotw,A5
        MOVE.B #4,D4
        BRA prnaddi
addin2  CMP #$2,D5
        BNE addiexit
        LEA dotl,A5
        MOVE.B #8,D4
        
        *Print
prnaddi LEA paddi,A1
        MOVE.B #14,D0
        TRAP #15
        LEA (A5),A1
        TRAP #15
        LEA hashtag,A1
        TRAP #15
        LEA dollarsign,A1
        TRAP #15
        ADDQ #2,A4
        CLR D5
        MOVE.B #16,D2
        MOVE.B #15,D0
addilp  ADDQ #1,D5
        MOVE.B (A4)+,D1
        TRAP #15
        CMP D4,D5
        BNE addilp
        LEA comma,A1
        MOVE.B #14,D0
        TRAP #15
        MOVE.B (source),D5
        MOVE.B (sourcem),D6
        JSR EA

        MOVE.B #$CC,D5
        ADDQ #4,A3
        ADDA.L D4,A3
        MOVE.B (A3),D3
        JSR newline
addiexit RTS
        
MOVE_B  JSR IsolateRegisters
                
        *Check valid opcode
        MOVE.B #$1,D7    
        MOVE.B (destm),D6
        CMP D6,D7
        BEQ mbexit
        
        *Print
        LEA pmove_b,A1
        MOVE.B #14,D0
        TRAP #15
        CLR D5
        CLR D6
        MOVE.B (source),D5
        MOVE.B (sourcem),D6
        JSR EA
        LEA comma,A1
        MOVE.B #14,D0
        TRAP #15
        MOVE.B (dest),D5
        MOVE.B (destm),D6
        JSR EA
        
        MOVE.B #$CC,D5
        ADDQ #4,A3
        MOVE.B (A3),D3
        JSR newline
mbexit  RTS

MOVE_W  JSR IsolateRegisters
        
        *Check movea or move
        MOVE.B #$1,D7
        MOVE.B (destm),D6
        CMP D6,D7
        BNE printnaw
        LEA pmovea_w,A1
        MOVE.B #14,D0
        TRAP #15
        BRA printmw
        
        *Print
printnaw LEA pmove_w,A1
        MOVE.B #14,D0
        TRAP #15
printmw CLR D5
        CLR D6
        MOVE.B (source),D5
        MOVE.B (sourcem),D6
        JSR EA
        LEA comma,A1
        MOVE.B #14,D0
        TRAP #15
        MOVE.B (dest),D5
        MOVE.B (destm),D6
        JSR EA
        
        MOVE.B #$CC,D5
        ADDQ #4,A3
        MOVE.B (A3),D3
        JSR newline
mwexit  RTS

MOVE_L  JSR IsolateRegisters 
        
        *Check movea or move
        MOVE.B #$1,D7
        MOVE.B (destm),D6
        CMP D6,D7
        BNE printnal
        LEA pmovea_l,A1
        MOVE.B #14,D0
        TRAP #15
        BRA printml
        
        *Print
printnal LEA pmove_l,A1
        MOVE.B #14,D0
        TRAP #15
printml CLR D5
        CLR D6
        MOVE.B (source),D5
        MOVE.B (sourcem),D6
        JSR EA
        LEA comma,A1
        MOVE.B #14,D0
        TRAP #15
        MOVE.B (dest),D5
        MOVE.B (destm),D6
        JSR EA
        
        MOVE.B #$CC,D5
        ADDQ #4,A3
        MOVE.B (A3),D3
        JSR newline
mlexit  RTS

MOVEQ_  CLR D5
        CLR D6
        CLR D7

        *Check 8th bit for 0
        ADDQ #1,A4
        MOVE.B (A4),D5
        MOVE.L #32,D7
        ASL.L D7,D5
        ASR.L D7,D5
        CMP #$1,D5
        BEQ mqexit
        
        *Get data register
        MOVE.B (A4)+,D5
        MOVE.B #28,D7
        LSL.L D7,D5
        ADDQ #1,D7
        LSR.L D7,D5
        MOVE.B D5,dest
        
        *Get and print value
        LEA pmoveq,A1
        MOVE.B #14,D0
        TRAP #15
        LEA hashtag,A1
        TRAP #15
        LEA dollarsign,A1
        TRAP #15
        MOVE.B (A4)+,D1
        CMP #$0,D1
        BEQ mqpnext
        MOVE.B #16,D2
        MOVE.B #15,D0
        TRAP #15
mqpnext MOVE.B (A4),D1
        MOVE.B #16,D2
        MOVE.B #15,D0
        TRAP #15
        LEA comma,A1
        MOVE.B #14,D0
        TRAP #15
        *Print data register
        CLR D5
        MOVE.B #$0,D6
        MOVE.B (dest),D5
        JSR EA
        
        CLR D5
        MOVE.B #$CC,D5
        ADDQ #4,A3
        MOVE.B (A3),D3
        JSR newline
mqexit  RTS

ADD_    JSR IsolateRegisters
        MOVE.B (destm),D5
        
        *EA + Dn to Dn
        CMP #$0,D5
        BNE addn1
        MOVE.B (sourcem),D6
        CMP #$1,D6
        BEQ addexit
        LEA dotb,A5
        MOVE.B #$0,destm
        BRA printadd
addn1   CMP #$1,D5
        BNE addn2
        LEA dotw,A5
        MOVE.B #$0,destm
        BRA printadd
addn2   CMP #$2,D5
        BNE addn3
        LEA dotl,A5
        MOVE.B #$0,destm
        BRA printadd
        
        *Dn + EA to EA
addn3   CMP #$4,D5
        BNE addn4
        LEA dotb,A5
        MOVE.B #$0,destm
        BRA printadd
addn4   CMP #$5,D5
        BNE addn5
        LEA dotw,A5
        MOVE.B #$0,destm
        BRA printadd
addn5   CMP #$6,D5
        BNE addn6
        LEA dotl,A5
        MOVE.B #$0,destm
        BRA printadd
        
        *ADDA
addn6   CMP #$3,D5
        BNE addn7
        LEA dotw,A5
        MOVE.B #$1,destm
        BRA printadda
addn7   LEA dotl,A5
        MOVE.B #$1,destm
        BRA printadda
        
printadda LEA padda,A1
        BRA addrest
printadd LEA padd,A1
addrest MOVE.B #14,D0
        TRAP #15
        LEA (A5),A1
        MOVE.B #14,D0
        TRAP #15
        CLR D5
        CLR D6
        MOVE.B (source),D5
        MOVE.B (sourcem),D6
        JSR EA
        LEA comma,A1
        MOVE.B #14,D0
        TRAP #15
        MOVE.B (dest),D5
        MOVE.B (destm),D6
        JSR EA

        MOVE.B #$CC,D5
        ADDQ #4,A3
        MOVE.B (A3),D3
        JSR newline
addexit RTS

NOP_    MOVEA.W A3,A4
        MOVE.B -(A4),D4

        MOVE.B (A4)+,D4
        CMP #$04,D4
        BEQ n1NOP
        BRA notNOP
n1NOP   MOVE.B (A4)+,D4
        CMP #$0E,D4
        BEQ n2NOP
        BRA notNOP
n2NOP   MOVE.B (A4)+,D4
        CMP #$07,D4
        BEQ n3NOP
        BRA notNOP
n3NOP   MOVE.B (A4),D4
        CMP #$01,D4
        BNE notNOP

        LEA pnop,A1
        MOVE.B #14,D0
        TRAP #15
        
        MOVE.B #$CC,D5
        ADDQ #4,A3
        MOVE.B (A3),D3
        JSR newline
notNOP  RTS

RTS_    MOVEA.W A3,A4
        MOVE.B -(A4),D4

        MOVE.B (A4)+,D4
        CMP #$04,D4
        BEQ r1RTS
        BRA notRTS
r1RTS   MOVE.B (A4)+,D4
        CMP #$0E,D4
        BEQ r2RTS
        BRA notRTS
r2RTS   MOVE.B (A4)+,D4
        CMP #$07,D4
        BEQ r3RTS
        BRA notRTS
r3RTS   MOVE.B (A4),D4
        CMP #$05,D4
        BNE notRTS

        LEA prts,A1
        MOVE.B #14,D0
        TRAP #15
        
        MOVE.B #$CC,D5
        ADDQ #4,A3
        MOVE.B (A3),D3
        JSR newline
notRTS  RTS

*============================ASL, ASR, LSL, LSR, ROL, ROR================================

rl_check    *for registers
    CLR.L D4
	ADDQ #$2,A4	*gets digits 7->4 into D4
	MOVE.B (A4),D4
	
	LSR.L #2,D4	*isolates the 7->6 spot, if it == 11 branch
	CMP #$3,D4
    BEQ rl_check2   *branch to check2 for ea
    
    CLR.L D4
	MOVE.B (A4),D4	*7->4 back into D4
	MOVE.B #31,D7
	LSL.L D7,D4	*isolates 4
	LSR.L D7,D4
	LSL.L #1,D4	*preps to get added
	CLR.L D2
	MOVE.B D4,D2	*stores 4 in D4

	ADDQ #$1,A4	*SHIFT gets digits 3->0 into D4
	MOVE.B (A4),D4
	LSR.L #3,D4	*isolate 3
	
	ADD.B D2,D4	*gets 4->3 in D4
	
	CMP #$2,D4  *if 4->3 == 2, bad
    BEQ not_rl_shift 
    
    BRA rl_shift
    
rl_check2   *ea check
    CLR.L D4
    SUBQ #1,A4  *shift 11->8
    MOVE.B (A4),D4  *get 11->8
    
    LSR.L #1,D4 *shift right to isolate 11->9
    
	CMP #$2,D4  *cant be 010
    BEQ not_rl_shift
    
    CMP #$4,D4  cant be 100
    BEQ not_rl_shift
    
    CMP #$5,D4
    BEQ not_rl_shift 
    
    CMP #$6,D4
    BEQ not_rl_shift 
    
    CMP #$7,D4
    BEQ not_rl_shift
    
    CLR.L D4
	ADDQ #$1,A4	*gets digits 7->4 into D4
	MOVE.B (A4),D4

    MOVE.B #30,D7
	LSL.L D7,D4	*isolate 5->4
	LSR.L D7,D4

	LSL.L #1,D4	*prep to add 3 onto it
	CLR.L D2
	MOVE.L D4,D2	*store 5->4 on D2 to read

	CLR.L D4
	ADDQ #1,A4	*gets 3->0
	MOVE.B (A4),D4
	
	LSR.L #3,D4	*isolate 3

	ADD.L D4,D2	*add shifted spot 5->4 to spot 3, D2 = 5->3
	
	CMP #$0,D2
	BEQ not_rl_shift
	
	CMP #$1,D2
	BEQ not_rl_shift
	
	CMP #$5,D2
	BEQ not_rl_shift
	
	CMP #$6,D2
	BEQ not_rl_shift
	
	CMP #$2,D2
	BEQ rl_shift
	
	CMP #$3,D2
	BEQ rl_shift
	
	CMP #$4,D2
	BEQ rl_shift
	
	*D4 equals 111 CONTINUE
	
    CLR.L D4
    ADDQ #1,A4  * spot -1 first bit after 3->0
    MOVE.B (A4),D4
    CMP #$FF,D4
	BEQ not_rl_shift
    CMP #$DD,D4
    BEQ not_rl_shift
    
    SUBQ #1,A4  *gets range 3->0
	
	CLR.L D4
	*ADDQ #3,A4  *get range 3->0
	MOVE.B (A4),D4  *load 3->0 onto D4
	
	CLR.L D7
	MOVE.B #29,D7
	LSL.L D7,D4     *isolate 2->0 on D4
	LSR.L D7,D4

    *SUBQ #3,A4  *reset counter back to 15->12
    
    CMP #$0,D4
    BEQ rl_shift
    
    CMP #$1,D4
    BEQ rl_shift
    
    BRA not_rl_shift

rl_shift	*is it memory or register shift?
    SUBQ #3,A4  *reset counter back to 15->12
	CLR.L D4
	ADDQ #$2,A4	*gets digits 7->4 into D4
	MOVE.B (A4),D4

	LSR.L #2,D4	*isolates the 7->6 spot
	CMP #$3,D4
	BNE rl_shift_b	*branch to part b(register)if necessary

*continue with ea shift

	CLR.L D4
	SUBQ #1,A4	*gets 11->8 number
	MOVE.B (A4),D4

	LSR.L #1,D4	*isolate 11->9

	
	CMP #$0,D4	*if 11->9 == 000, jump to print AS
	BEQ print_as_a

	CMP #$1,D4	*if 11->9 == 001, jump to print LS
	BEQ print_ls_a

	CMP #$2,D4	*if 11->9 == 010, bad op 
	BEQ not_rl_shift 

	CMP #$3,D4	*if 11->9 == 011. jump to print RO
	BEQ print_ro_a

	*could probably just BRA not_rl_shift
	
	CMP #$4,D4	*if 11->9 == 100, bad op 
	BEQ not_rl_shift

	CMP #$5,D4	*if 11->9 == 101, bad op 
	BEQ not_rl_shift 

	CMP #$6,D4	*if 11->9 == 110, bad op 
	BEQ not_rl_shift 

	CMP #$7,D4	*if 11->9 == 111, bad op 
	BEQ not_rl_shift

	BRA not_rl_shift

print_as_a  *print AS

	LEA pas,A1
	MOVE.B #14,D0
	TRAP #15
	BRA rl_shift_a1

print_ls_a  *print LS
	LEA pls,A1
	MOVE.B #14,D0
	TRAP #15
	BRA rl_shift_a1

print_ro_a  *print RO
	LEA pro,A1
	MOVE.B #14,D0
	TRAP #15
	BRA rl_shift_a1

rl_shift_a1	*ea
	CLR.L D4
	MOVE.B (A4),D4	*gets 11->8 into D4
	MOVE.B #31,D7
	LSL.L D7,D4	*isolate 8
	LSR.L D7,D4
	
	CMP #$0,D4	*if spot 8 == 0,go to print R
	BEQ print_r_a
	
	CMP #$1,D4	*if spot 8 == 1, go to print L
	BEQ print_l_a

	BRA not_rl_shift

print_r_a   *print R
	LEA pright,A1
	MOVE.B #14,D0
	TRAP #15
	BRA rl_shift_a2

print_l_a   *print L
	LEA pleft,A1
	MOVE.B #14,D0
	TRAP #15
	BRA rl_shift_a2

rl_shift_a2	*register shift
	CLR.L D4
	ADDQ #$1,A4	*gets digits 7->4 into D4
	MOVE.B (A4),D4

    MOVE.B #30,D7
	LSL.L D7,D4	*isolate 5->4
	LSR.L D7,D4

	LSL.L #1,D4	*prep to add 3 onto it
	CLR.L D2
	MOVE.L D4,D2	*store 5->4 on D2 to read

	CLR.L D4
	ADDQ #1,A4	*gets 3->0
	MOVE.B (A4),D4
	
	LSR.L #3,D4	*isolate 3

	ADD.L D4,D2	*add shifted spot 5->4 to spot 3, D2 = 5->3

	CLR.L D4

	MOVE.B (A4),D4	*load 3->0 on D4
	MOVE.B #29,D7
	LSL.L D7,D4	*isolate 2->0
	LSR.L D7,D4
	CLR.L D1	*load 2->0 into D1
	MOVE.B D4,D1

	LEA space,A1
	MOVE.B #14,D0
	TRAP #15

	CMP #$2,D2	*if == 010, print (An)
	BEQ print_pAnp_rl_a

	CMP #$3,D2	*if == 011, print (An)+
	BEQ print_pAnpplus_rl_a

	CMP #$4,D2	*if == 100, print -(An)
	BEQ print_minuspAnp_rl_a

	CMP #$7,D2	*if == 111, print ea
	BEQ print_ea_rl_a
	
	BRA not_rl_shift

print_pAnp_rl_a     *print (An)
	LEA leftp,A1
	MOVE.B #14,D0
	TRAP #15

	LEA aregister,A1
	MOVE.B #14,D0
	TRAP #15
	
	MOVE.B #3,D0
	TRAP #15
    
	LEA rightp,A1
	MOVE.B #14,D0
	TRAP #15

	MOVE.B #$CC,D5      *DONE
	ADDQ #4,A3
	MOVE.B (A3),D3
	BRA newline

print_pAnpplus_rl_a *print (An)+
	LEA leftp,A1
	MOVE.B #14,D0
	TRAP #15

	LEA aregister,A1
	MOVE.B #14,D0
	TRAP #15
	
	MOVE.B #3,D0
	TRAP #15

	LEA rightp,A1
	MOVE.B #14,D0
	TRAP #15

	LEA plus,A1
	MOVE.B #14,D0
	TRAP #15

	MOVE.B #$CC,D5      *DONE
	ADDQ #4,A3
	MOVE.B (A3),D3
	BRA newline

print_minuspAnp_rl_a    *print -(An)

	LEA minus,A1
	MOVE.B #14,D0
	TRAP #15

	LEA leftp,A1
	MOVE.B #14,D0
	TRAP #15

	LEA aregister,A1
	MOVE.B #14,D0
	TRAP #15
	
	MOVE.B #3,D0
	TRAP #15

	LEA rightp,A1
	MOVE.B #14,D0
	TRAP #15

	MOVE.B #$CC,D5      *DONE
	ADDQ #4,A3
	MOVE.B (A3),D3
	BRA newline
	
print_ea_rl_a   *print ea 
	CMP #$2,D1
	BEQ not_rl_shift

	CMP #$3,D1
	BEQ not_rl_shift

	CMP #$4,D1
	BEQ not_rl_shift

	CMP #$5,D1
	BEQ not_rl_shift

	CMP #$6,D1
	BEQ not_rl_shift

	CMP #$7,D1
	BEQ not_rl_shift

	LEA dollarsign,A1
	MOVE.B #14,D0
	TRAP #15

    CLR.L D4
    ADDQ #1,A4
    MOVE.B (A4),D4
    CMP #$FF,D4
	BEQ not_rl_shift
    CMP #$DD,D4
    BEQ not_rl_shift
    
    BRA print_rl_ea
    
print_rl_ea
    CLR.L D4
    MOVE.B (A4),D4
	CLR.L D2
	MOVE.B #16,D2
	CLR.L D1
	MOVE.L D4,D1	*print extra ea from D4
	MOVE.B #15,D0
	TRAP #15
	
	ADDQ #1,A4
	MOVE.B (A4),D4
	
	ADDQ #1,A3
	
	CMP #$FF,D4
	BEQ rl_shift_done
	CMP #$DD,D4
	BNE print_rl_ea

	MOVE.B #$CC,D5  *DONE
	ADDQ #4,A3
	MOVE.B (A3),D3
	BRA newline
	
rl_shift_done
	MOVE.B #$CC,D5  *DONE
	ADDQ #4,A3
	MOVE.B (A3),D3
	BRA newline

rl_shift_b  *=====================PART B FOR MEMORY LOCATIONS==================
    CLR.L D4
	MOVE.B (A4),D4	*7->4 back into D4
	MOVE.B #31,D7
	LSL.L D7,D4	*isolates 4
	LSR.L D7,D4
	LSL.L #1,D4	*preps to get added
	CLR.L D2
	MOVE.B D4,D2	*stores 4 in D4

	ADDQ #$1,A4	*SHIFT gets digits 3->0 into D4
	MOVE.B (A4),D4
	LSR.L #3,D4	*isolate 3
	
	ADD.B D2,D4	*gets 4->3
	
	CMP #$0,D4
	BEQ print_as_b

	CMP #$1,D4
	BEQ print_ls_b

	CMP #$3,D4
	BEQ print_ro_b

	BRA not_rl_shift
	
print_as_b      *print AS
	LEA pas,A1
	MOVE.B #14,D0
	TRAP #15
	BRA rl_shift_b1

print_ls_b      *print LS
	LEA pls,A1
	MOVE.B #14,D0
	TRAP #15
	BRA rl_shift_b1

print_ro_b  *print RO
	LEA pro,A1
	MOVE.B #14,D0
	TRAP #15
	BRA rl_shift_b1

rl_shift_b1
	SUBQ #2,A4      *SHIFT
	MOVE.B (A4),D4	*gets 11->8 into D4
	MOVE.B #31,D7
	LSL.L D7,D4	    *isolate 8
	LSR.L D7,D4

	CMP #$0,D4	    *if spot 8 == 0, print R
	BEQ print_r_b

	CMP #$1,D4	    *if spot 8 == 1, print L
	BEQ print_l_b

	BRA not_rl_shift

print_r_b   *print R
	LEA pright,A1
	MOVE.B #14,D0
	TRAP #15
	BRA rl_shift_b2

print_l_b   *print L
	LEA pleft,A1
	MOVE.B #14,D0
	TRAP #15
	BRA rl_shift_b2

rl_shift_b2
	ADDQ #1,A4      *SHIFT
	MOVE.B (A4),D4	*gets 7->4 into D4
	LSR.L #2,D4	    *isolate 7->6

	CMP #$0,D4
	BEQ print_byte_rl_b

	CMP #$1,D4
	BEQ print_word_rl_b

	CMP #$2,D4
	BEQ print_long_rl_b
	
	BRA not_rl_shift

print_byte_rl_b *print .B
	LEA dotb,A1
	MOVE.B #14,D0
	TRAP #15
	BRA rl_shift_b3

print_word_rl_b *print .W
	LEA dotw,A1
	MOVE.B #14,D0
	TRAP #15
	BRA rl_shift_b3

print_long_rl_b *print .L
	LEA dotl,A1
	MOVE.B #14,D0
	TRAP #15
	BRA rl_shift_b3

rl_shift_b3
    CLR.L D4
	MOVE.B (A4),D4	*gets 7->4 into D4
	LSR.L #1,D4	    *isolate 5
	MOVE.B #31,D7
	LSL.L D7,D4
	LSR.L D7,D4

	CMP #$0,D4
	BEQ print_number_rl_b

	CMP #$1,D4
	BEQ print_d_rl_b

	BRA not_rl_shift

print_number_rl_b   *prints #number
	SUBQ #1,A4	    *SHIFT gets 11->8 into D4
	MOVE.B (A4),D4

	LSR.L #1,D4	    *isolates 11->9
	CLR.L D1
	MOVE.L D4,D1
	
	LEA hashtag,A1  *print #
	MOVE.B #14,D0
	TRAP #15

	MOVE.B #3,D0    *print number in D1
	TRAP #15
	
	LEA comma,A1    *print ,
	MOVE.B #14,D0
	TRAP #15

	LEA dregister,A1	*print D
	MOVE.B #14,D0
	TRAP #15

	ADDQ #2,A4	    *SHIFT get 3->0 into D4
	MOVE.B (A4),D4

	CLR.L D1
	MOVE.B #29,D7
	LSL.L D7,D4	    *isolate 2->0
	LSR.L D7,D4
	MOVE.L D4,D1

	MOVE.B #3,D0	*print number in D1, dregister
	TRAP #15

	MOVE.B #$CC,D5  *DONE
	ADDQ #4,A3
	MOVE.B (A3),D3
	BRA newline

print_d_rl_b
	SUBQ #1,A4	    *SHIFT gets 11->8 into D4
	MOVE.B (A4),D4

	LSR.L #1,D4	    *isolates 11->9
	CLR.L D1
	MOVE.L D4,D1
	
	LEA dregister,A1	*print D
	MOVE.B #14,D0
	Trap #15

	MOVE.B #3,D0	*print number in D1
	TRAP #15

	LEA comma,A1	*print ,
	MOVE.B #14,D0
	Trap #15

	LEA dregister,A1	*print D
	MOVE.B #14,D0
	TRAP #15

	ADDQ #2,A4  	*get 3->0 into D4
	MOVE.B (A4),D4

	CLR.L D1
	MOVE.B #29,D7
	LSL.L D7,D4	    *isolate 2->0
	LSR.L D7,D4
	MOVE.L D4,D1

	MOVE.B #3,D0	*print number in D1, dregister
	TRAP #15

	MOVE.B #$CC,D5  *DONE
	ADDQ #4,A3
	MOVE.B (A3),D3
	BRA  newline

not_rl_shift
	RTS
*===================================================================================

*Pre-conditions:    
*D6 needs to contain addressing mode
*D5 needs to contain EA/register number          
EA      CMP #$0,D6
        BEQ EA_000
        
        CMP #$1,D6
        BEQ EA_001
        
        RTS

EA_000  LEA dregister,A1
        MOVE.B #14,D0
        TRAP #15
        MOVE.B D5,D1
        MOVE.B #16,D2
        MOVE.B #15,D0
        TRAP #15
        RTS
        
EA_001  LEA aregister,A1
        MOVE.B #14,D0
        TRAP #15
        MOVE.B D5,D1
        MOVE.B #16,D2
        MOVE.B #15,D0
        TRAP #15
        RTS

EA_010  LEA leftp,A1
        MOVE.B #14,D0
        LEA aregister,A1
        TRAP #15
        MOVE.B D5,D1
        MOVE.B #16,D2
        MOVE.B #15,D0
        TRAP #15
        LEA rightp,A1
        MOVE.B #14,D0
        LEA aregister,A1
        RTS

EA_011  NOP
        RTS 

EA_100  NOP
        RTS    

EA_101  NOP
        RTS  
        
EA_110  NOP
        RTS

EA_111  NOP
        RTS 

IsolateRegisters
        CLR D5
        CLR D6
        CLR D7
        
        *Isolate and get source mode
        ADDQ #2,A4
        MOVE.B (A4),D5
        ADDQ #1,A4
        MOVE.B (A4),D6
        MOVE.B #30,D7
        LSL.L D7,D5
        LSR.L D7,D5
        LSR.B #$3,D6
        ADD.B D6,D5
        MOVE.B D5,sourcem
        
        *Isolate and get source register
        MOVEA.W A3,A4
        ADDQ #2,A4
        MOVE.B (A4),D5
        MOVE.B #29,D7
        LSL.L D7,D5
        LSR.L D7,D5
        MOVE.B D5,source
        
        *Isolate and get destination mode
        MOVEA.W A3,A4
        MOVE.B (A4),D5
        ADDQ #1,A4
        MOVE.B (A4),D6
        MOVE.B #31,D7
        LSL.L D7,D5
        SUBI.B #2,D7
        LSR.L D7,D5
        LSR.B #$2,D6
        ADD.B D5,D6
        MOVE.B D6,destm
        
        *Isolate and get destination register
        MOVEA.W A3,A4
        MOVE.B (A4),D5
        LSR.B #$1,D5
        MOVE.B D5,dest

        RTS 
      
    END    START        ; last line of source


























*~Font name~Courier New~
*~Font size~10~
*~Tab type~1~
*~Tab size~4~
