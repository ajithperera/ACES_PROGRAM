	.text
	.file	"oed__ovl_prepare_ctr.f"
	.section	.rodata.cst8,"aM",@progbits,8
	.p2align	3, 0x0
.LCPI0_0:
	.quad	0x3fe0000000000000
.LCPI0_1:
	.quad	0x3fe8000000000000
	.text
	.globl	oed__ovl_prepare_ctr_
	.p2align	4
	.type	oed__ovl_prepare_ctr_,@function
oed__ovl_prepare_ctr_:
.Lfunc_begin0:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
	pushq	%r15
	pushq	%r14
	pushq	%r13
	pushq	%r12
	pushq	%rbx
	subq	$88, %rsp
	.cfi_offset %rbx, -56
	.cfi_offset %r12, -48
	.cfi_offset %r13, -40
	.cfi_offset %r14, -32
	.cfi_offset %r15, -24
	movq	%r9, -80(%rbp)
	movq	%r8, %r15
	movq	%rcx, -64(%rbp)
	movq	%rdx, -72(%rbp)
	movq	%rdi, -104(%rbp)
	movq	96(%rbp), %rax
	movq	%rax, -88(%rbp)
	movq	80(%rbp), %r14
	movq	72(%rbp), %r13
	movq	56(%rbp), %rdi
	movq	48(%rbp), %rbx
	movq	40(%rbp), %rax
	movq	%rax, -112(%rbp)
	movq	32(%rbp), %rax
	movq	%rax, -120(%rbp)
.Ltmp0:
	movq	(%rsi), %rdx
	movq	24(%rbp), %rax
	movq	%rax, -96(%rbp)
	movq	16(%rbp), %r12
.Ltmp1:
	testq	%rdx, %rdx
	movq	%rdi, -56(%rbp)
	jle	.LBB0_2
	movq	88(%rbp), %rdi
	movq	64(%rbp), %rsi
	shlq	$3, %rdx
.Ltmp2:
	callq	memcpy@PLT
	movq	-56(%rbp), %rdi
.LBB0_2:
.Ltmp3:
	cvtsi2sdq	(%r15), %xmm0
	mulsd	.LCPI0_0(%rip), %xmm0
	addsd	.LCPI0_1(%rip), %xmm0
	movsd	%xmm0, -48(%rbp)
.Ltmp4:
	testb	$1, (%rbx)
	movq	-72(%rbp), %rax
	movq	(%rax), %rbx
	jne	.LBB0_3
	testq	%rbx, %rbx
.Ltmp5:
	jle	.LBB0_12
	xorl	%r15d, %r15d
	.p2align	4
.LBB0_11:
.Ltmp6:
	movsd	(%r12,%r15,8), %xmm0
	movsd	-48(%rbp), %xmm1
	callq	pow@PLT
	movsd	%xmm0, (%r13,%r15,8)
.Ltmp7:
	incq	%r15
	cmpq	%r15, %rbx
	jne	.LBB0_11
.LBB0_12:
	movq	-64(%rbp), %rax
.Ltmp8:
	movq	(%rax), %rbx
	testq	%rbx, %rbx
	movq	-96(%rbp), %r12
	jle	.LBB0_6
	movq	-80(%rbp), %rax
.Ltmp9:
	xorps	%xmm0, %xmm0
	cvtsi2sdq	(%rax), %xmm0
	mulsd	.LCPI0_0(%rip), %xmm0
	addsd	.LCPI0_1(%rip), %xmm0
	movsd	%xmm0, -48(%rbp)
	xorl	%r15d, %r15d
	.p2align	4
.LBB0_14:
.Ltmp10:
	movsd	(%r12,%r15,8), %xmm0
	movsd	-48(%rbp), %xmm1
	callq	pow@PLT
	movsd	%xmm0, (%r14,%r15,8)
.Ltmp11:
	incq	%r15
	cmpq	%r15, %rbx
	jne	.LBB0_14
	jmp	.LBB0_6
.LBB0_3:
	testq	%rbx, %rbx
	jle	.LBB0_7
	xorl	%r15d, %r15d
	.p2align	4
.LBB0_5:
.Ltmp12:
	movsd	(%r12,%r15,8), %xmm0
	movsd	-48(%rbp), %xmm1
	callq	pow@PLT
	movsd	%xmm0, (%r13,%r15,8)
.Ltmp13:
	movsd	%xmm0, (%r14,%r15,8)
.Ltmp14:
	incq	%r15
	cmpq	%r15, %rbx
	jne	.LBB0_5
.LBB0_6:
	movq	-72(%rbp), %rax
.Ltmp15:
	movq	(%rax), %rbx
	movq	-56(%rbp), %rdi
.LBB0_7:
	movq	-112(%rbp), %rax
.Ltmp16:
	movsd	(%rax), %xmm0
	movq	-120(%rbp), %rax
	mulsd	(%rax), %xmm0
	movq	-64(%rbp), %rax
.Ltmp17:
	movq	(%rax), %rax
	cmpq	%rax, %rbx
	jle	.LBB0_8
.Ltmp18:
	testq	%rax, %rax
	jle	.LBB0_20
	xorl	%ecx, %ecx
	.p2align	4
.LBB0_19:
	movsd	(%r14,%rcx,8), %xmm1
.Ltmp19:
	mulsd	%xmm0, %xmm1
	movsd	%xmm1, (%r14,%rcx,8)
.Ltmp20:
	incq	%rcx
	cmpq	%rcx, %rax
	jne	.LBB0_19
	jmp	.LBB0_20
.LBB0_8:
.Ltmp21:
	testq	%rbx, %rbx
	jle	.LBB0_20
	xorl	%eax, %eax
	.p2align	4
.LBB0_10:
	movsd	(%r13,%rax,8), %xmm1
.Ltmp22:
	mulsd	%xmm0, %xmm1
	movsd	%xmm1, (%r13,%rax,8)
.Ltmp23:
	incq	%rax
	cmpq	%rax, %rbx
	jne	.LBB0_10
.LBB0_20:
.Ltmp24:
	testb	$1, (%rdi)
	je	.LBB0_22
	movq	-104(%rbp), %rax
.Ltmp25:
	movq	(%rax), %rdx
	testq	%rdx, %rdx
	jle	.LBB0_22
	shlq	$3, %rdx
	movq	-88(%rbp), %rdi
.Ltmp26:
	xorl	%esi, %esi
	addq	$88, %rsp
	popq	%rbx
	popq	%r12
	popq	%r13
	popq	%r14
	popq	%r15
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	jmp	memset@PLT
.LBB0_22:
	.cfi_def_cfa %rbp, 16
.Ltmp27:
	addq	$88, %rsp
	popq	%rbx
	popq	%r12
	popq	%r13
	popq	%r14
	popq	%r15
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	retq
.Ltmp28:
.Lfunc_end0:
	.size	oed__ovl_prepare_ctr_, .Lfunc_end0-oed__ovl_prepare_ctr_
	.cfi_endproc

	.section	.trace,"a",@progbits
.Lmodule_begin0:
	.byte	10
	.short	2
	.byte	0
	.long	.Lmodule_end0-.Lmodule_begin0
	.quad	oed__ovl_prepare_ctr_
	.long	1
	.long	.Lfunc_end0-oed__ovl_prepare_ctr_
	.short	0
	.short	22
	.ascii	"oed__ovl_prepare_ctr.f"
	.p2align	3, 0x0
	.byte	12
	.byte	0
	.short	20
	.quad	oed__ovl_prepare_ctr_
	.ascii	"oed__ovl_prepare_ctr"
	.byte	4
	.byte	26
	.byte	9
	.long	(.Ltmp1-oed__ovl_prepare_ctr_)-1
	.byte	4
	.byte	127
	.byte	9
	.long	(.Ltmp2-.Ltmp1)-1
	.byte	4
	.byte	1
	.byte	9
	.long	(.Ltmp3-.Ltmp2)-1
	.byte	4
	.byte	9
	.byte	9
	.long	(.Ltmp4-.Ltmp3)-1
	.byte	4
	.byte	1
	.byte	9
	.long	(.Ltmp5-.Ltmp4)-1
	.byte	4
	.byte	6
	.byte	9
	.long	(.Ltmp6-.Ltmp5)-1
	.byte	4
	.byte	1
	.byte	9
	.long	(.Ltmp7-.Ltmp6)-1
	.byte	4
	.byte	1
	.byte	9
	.long	(.Ltmp8-.Ltmp7)-1
	.byte	4
	.byte	2
	.byte	9
	.long	(.Ltmp9-.Ltmp8)-1
	.byte	4
	.byte	-1
	.byte	9
	.long	(.Ltmp10-.Ltmp9)-1
	.byte	4
	.byte	2
	.byte	9
	.long	(.Ltmp11-.Ltmp10)-1
	.byte	4
	.byte	1
	.byte	9
	.long	(.Ltmp12-.Ltmp11)-1
	.byte	4
	.byte	-10
	.byte	9
	.long	(.Ltmp13-.Ltmp12)-1
	.byte	4
	.byte	1
	.byte	9
	.long	(.Ltmp14-.Ltmp13)-1
	.byte	4
	.byte	1
	.byte	9
	.long	(.Ltmp15-.Ltmp14)-1
	.byte	4
	.byte	18
	.byte	9
	.long	(.Ltmp16-.Ltmp15)-1
	.byte	4
	.byte	-2
	.byte	9
	.long	(.Ltmp17-.Ltmp16)-1
	.byte	4
	.byte	2
	.byte	9
	.long	(.Ltmp18-.Ltmp17)-1
	.byte	4
	.byte	5
	.byte	9
	.long	(.Ltmp19-.Ltmp18)-1
	.byte	4
	.byte	1
	.byte	9
	.long	(.Ltmp20-.Ltmp19)-1
	.byte	4
	.byte	1
	.byte	9
	.long	(.Ltmp21-.Ltmp20)-1
	.byte	4
	.byte	-6
	.byte	9
	.long	(.Ltmp22-.Ltmp21)-1
	.byte	4
	.byte	1
	.byte	9
	.long	(.Ltmp23-.Ltmp22)-1
	.byte	4
	.byte	1
	.byte	9
	.long	(.Ltmp24-.Ltmp23)-1
	.byte	4
	.byte	11
	.byte	9
	.long	(.Ltmp25-.Ltmp24)-1
	.byte	4
	.byte	1
	.byte	9
	.long	(.Ltmp26-.Ltmp25)-1
	.byte	4
	.byte	1
	.byte	9
	.long	(.Ltmp27-.Ltmp26)-1
	.byte	4
	.byte	9
	.byte	9
	.long	(.Lfunc_end0-.Ltmp27)-1
.Lmodule_end0:
	.section	".note.GNU-stack","",@progbits
