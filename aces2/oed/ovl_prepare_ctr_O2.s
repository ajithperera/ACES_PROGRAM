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
	subq	$72, %rsp
	.cfi_offset %rbx, -56
	.cfi_offset %r12, -48
	.cfi_offset %r13, -40
	.cfi_offset %r14, -32
	.cfi_offset %r15, -24
	movq	%r9, -88(%rbp)
	movq	96(%rbp), %r15
	movq	80(%rbp), %r13
	movq	72(%rbp), %r12
	movq	56(%rbp), %rax
	movq	%rax, -112(%rbp)
	movq	48(%rbp), %r10
	movq	40(%rbp), %rax
	movq	%rax, -96(%rbp)
	movq	32(%rbp), %rax
	movq	%rax, -104(%rbp)
.Ltmp0:
	movq	(%rdx), %r9
	movq	(%rcx), %rax
	movq	%rax, -56(%rbp)
	movq	(%rsi), %rdx
	movq	(%rdi), %rax
	movq	%rax, -72(%rbp)
	movq	24(%rbp), %r14
	movq	16(%rbp), %rbx
.Ltmp1:
	testq	%rdx, %rdx
	movq	%r9, -48(%rbp)
	jle	.LBB0_1
	movq	88(%rbp), %rdi
	movq	64(%rbp), %rsi
	cmpq	$13, %rdx
	jb	.LBB0_28
.Ltmp2:
	shlq	$3, %rdx
	movq	%r15, -64(%rbp)
	movq	%r8, %r15
	movq	%r10, -80(%rbp)
	callq	_intel_fast_memcpy@PLT
	movq	-80(%rbp), %r10
	movq	-48(%rbp), %r9
	movq	%r15, %r8
	movq	-64(%rbp), %r15
	jmp	.LBB0_1
.LBB0_28:
	xorl	%eax, %eax
	.p2align	4
.LBB0_29:
	movsd	(%rsi,%rax,8), %xmm0
	movsd	%xmm0, (%rdi,%rax,8)
.Ltmp3:
	incq	%rax
	cmpq	%rax, %rdx
	jne	.LBB0_29
.LBB0_1:
.Ltmp4:
	cvtsi2sdq	(%r8), %xmm8
	mulsd	.LCPI0_0(%rip), %xmm8
	addsd	.LCPI0_1(%rip), %xmm8
.Ltmp5:
	testb	$1, (%r10)
	jne	.LBB0_2
	testq	%r9, %r9
.Ltmp6:
	jle	.LBB0_5
	movq	%r9, %rdi
	movabsq	$9223372036854775806, %rax
	andq	%rax, %rdi
	je	.LBB0_10
	movq	%r15, -64(%rbp)
.Ltmp7:
	movapd	%xmm8, %xmm9
	unpcklpd	%xmm8, %xmm9
	xorl	%esi, %esi
	movq	__svml_pow2@GOTPCREL(%rip), %r15
	.p2align	4
.LBB0_36:
	movups	(%rbx,%rsi,8), %xmm0
	movapd	%xmm9, %xmm1
	callq	*%r15
	movups	%xmm0, (%r12,%rsi,8)
.Ltmp8:
	addq	$2, %rsi
	cmpq	%rdi, %rsi
	jb	.LBB0_36
.Ltmp9:
	cmpq	%rdi, -48(%rbp)
	movq	-64(%rbp), %r15
	jne	.LBB0_38
	jmp	.LBB0_5
.LBB0_2:
	testq	%r9, %r9
.Ltmp10:
	jle	.LBB0_11
	movq	-48(%rbp), %rdi
	movabsq	$9223372036854775806, %rax
	andq	%rax, %rdi
	je	.LBB0_4
.Ltmp11:
	movapd	%xmm8, %xmm9
	unpcklpd	%xmm8, %xmm9
	xorl	%r14d, %r14d
	movq	__svml_pow2@GOTPCREL(%rip), %rsi
	.p2align	4
.LBB0_31:
	movups	(%rbx,%r14,8), %xmm0
	movapd	%xmm9, %xmm1
	callq	*%rsi
	movups	%xmm0, (%r12,%r14,8)
.Ltmp12:
	movups	%xmm0, (%r13,%r14,8)
.Ltmp13:
	addq	$2, %r14
	cmpq	%rdi, %r14
	jb	.LBB0_31
.Ltmp14:
	cmpq	%rdi, -48(%rbp)
	jne	.LBB0_33
	jmp	.LBB0_11
.LBB0_10:
	xorl	%edi, %edi
.LBB0_38:
	movq	__svml_pow1@GOTPCREL(%rip), %rsi
	.p2align	4
.LBB0_39:
.Ltmp15:
	movsd	(%rbx,%rdi,8), %xmm0
.Ltmp16:
	movapd	%xmm8, %xmm1
	callq	*%rsi
.Ltmp17:
	movsd	%xmm0, (%r12,%rdi,8)
.Ltmp18:
	incq	%rdi
	cmpq	%rdi, -48(%rbp)
	jne	.LBB0_39
.LBB0_5:
.Ltmp19:
	cmpq	$0, -56(%rbp)
	jle	.LBB0_11
	movq	-88(%rbp), %rax
.Ltmp20:
	xorps	%xmm8, %xmm8
	cvtsi2sdq	(%rax), %xmm8
	mulsd	.LCPI0_0(%rip), %xmm8
	addsd	.LCPI0_1(%rip), %xmm8
	movq	-56(%rbp), %rdi
.Ltmp21:
	movabsq	$9223372036854775806, %rax
	andq	%rax, %rdi
	je	.LBB0_7
.Ltmp22:
	movapd	%xmm8, %xmm9
	unpcklpd	%xmm8, %xmm9
	xorl	%esi, %esi
	movq	__svml_pow2@GOTPCREL(%rip), %rbx
	.p2align	4
.LBB0_41:
	movups	(%r14,%rsi,8), %xmm0
	movapd	%xmm9, %xmm1
	callq	*%rbx
	movups	%xmm0, (%r13,%rsi,8)
.Ltmp23:
	addq	$2, %rsi
	cmpq	%rdi, %rsi
	jb	.LBB0_41
.Ltmp24:
	cmpq	%rdi, -56(%rbp)
	jne	.LBB0_43
	jmp	.LBB0_11
.LBB0_4:
	xorl	%edi, %edi
.LBB0_33:
	movq	__svml_pow1@GOTPCREL(%rip), %r14
	.p2align	4
.LBB0_34:
.Ltmp25:
	movsd	(%rbx,%rdi,8), %xmm0
.Ltmp26:
	movapd	%xmm8, %xmm1
	callq	*%r14
.Ltmp27:
	movsd	%xmm0, (%r12,%rdi,8)
.Ltmp28:
	movsd	%xmm0, (%r13,%rdi,8)
.Ltmp29:
	incq	%rdi
	cmpq	%rdi, -48(%rbp)
	jne	.LBB0_34
	jmp	.LBB0_11
.LBB0_7:
	xorl	%edi, %edi
.LBB0_43:
	movq	__svml_pow1@GOTPCREL(%rip), %rsi
	.p2align	4
.LBB0_44:
.Ltmp30:
	movsd	(%r14,%rdi,8), %xmm0
.Ltmp31:
	movapd	%xmm8, %xmm1
	callq	*%rsi
.Ltmp32:
	movsd	%xmm0, (%r13,%rdi,8)
.Ltmp33:
	incq	%rdi
	cmpq	%rdi, -56(%rbp)
	jne	.LBB0_44
.LBB0_11:
	movq	-96(%rbp), %rax
.Ltmp34:
	movsd	(%rax), %xmm0
	movq	-104(%rbp), %rax
	mulsd	(%rax), %xmm0
	movq	-48(%rbp), %rdi
	movq	-56(%rbp), %r8
.Ltmp35:
	cmpq	%r8, %rdi
	jle	.LBB0_12
.Ltmp36:
	testq	%r8, %r8
	movq	-72(%rbp), %rdx
	movabsq	$9223372036854775806, %rsi
	jle	.LBB0_22
	cmpq	$8, %r8
	jb	.LBB0_21
.Ltmp37:
	movq	%r8, %rax
	shrq	$3, %rax
.Ltmp38:
	movapd	%xmm0, %xmm1
	unpcklpd	%xmm0, %xmm1
	leaq	48(%r13), %rcx
	.p2align	4
.LBB0_20:
	movupd	-48(%rcx), %xmm2
	movupd	-32(%rcx), %xmm3
	movupd	-16(%rcx), %xmm4
	movupd	(%rcx), %xmm5
	mulpd	%xmm1, %xmm2
	movupd	%xmm2, -48(%rcx)
	mulpd	%xmm1, %xmm3
	movupd	%xmm3, -32(%rcx)
	mulpd	%xmm1, %xmm4
	movupd	%xmm4, -16(%rcx)
	mulpd	%xmm1, %xmm5
	movupd	%xmm5, (%rcx)
.Ltmp39:
	addq	$64, %rcx
	decq	%rax
	jne	.LBB0_20
.LBB0_21:
	addq	$-6, %rsi
	andq	%r8, %rsi
	andl	$7, %r8d
	jmpq	*.LJTI0_1(,%r8,8)
.LBB0_52:
	movsd	48(%r13,%rsi,8), %xmm1
.Ltmp40:
	mulsd	%xmm0, %xmm1
	movsd	%xmm1, 48(%r13,%rsi,8)
.LBB0_53:
	movsd	40(%r13,%rsi,8), %xmm1
	mulsd	%xmm0, %xmm1
	movsd	%xmm1, 40(%r13,%rsi,8)
.LBB0_54:
	movsd	32(%r13,%rsi,8), %xmm1
	mulsd	%xmm0, %xmm1
	movsd	%xmm1, 32(%r13,%rsi,8)
.LBB0_55:
	movsd	24(%r13,%rsi,8), %xmm1
	mulsd	%xmm0, %xmm1
	movsd	%xmm1, 24(%r13,%rsi,8)
.LBB0_56:
	movsd	16(%r13,%rsi,8), %xmm1
	mulsd	%xmm0, %xmm1
	movsd	%xmm1, 16(%r13,%rsi,8)
.LBB0_57:
	movsd	8(%r13,%rsi,8), %xmm1
	mulsd	%xmm0, %xmm1
	movsd	%xmm1, 8(%r13,%rsi,8)
.LBB0_58:
	mulsd	(%r13,%rsi,8), %xmm0
	movsd	%xmm0, (%r13,%rsi,8)
	jmp	.LBB0_22
.LBB0_12:
	testq	%rdi, %rdi
	movq	-72(%rbp), %rdx
	movabsq	$9223372036854775806, %rsi
.Ltmp41:
	jle	.LBB0_22
	cmpq	$8, %rdi
	jb	.LBB0_16
.Ltmp42:
	movq	%rdi, %rax
	shrq	$3, %rax
.Ltmp43:
	movapd	%xmm0, %xmm1
	unpcklpd	%xmm0, %xmm1
	leaq	48(%r12), %rcx
	.p2align	4
.LBB0_15:
	movupd	-48(%rcx), %xmm2
	movupd	-32(%rcx), %xmm3
	movupd	-16(%rcx), %xmm4
	movupd	(%rcx), %xmm5
	mulpd	%xmm1, %xmm2
	movupd	%xmm2, -48(%rcx)
	mulpd	%xmm1, %xmm3
	movupd	%xmm3, -32(%rcx)
	mulpd	%xmm1, %xmm4
	movupd	%xmm4, -16(%rcx)
	mulpd	%xmm1, %xmm5
	movupd	%xmm5, (%rcx)
.Ltmp44:
	addq	$64, %rcx
	decq	%rax
	jne	.LBB0_15
.LBB0_16:
	addq	$-6, %rsi
	andq	%rdi, %rsi
	andl	$7, %edi
	jmpq	*.LJTI0_0(,%rdi,8)
.LBB0_45:
	movsd	48(%r12,%rsi,8), %xmm1
.Ltmp45:
	mulsd	%xmm0, %xmm1
	movsd	%xmm1, 48(%r12,%rsi,8)
.LBB0_46:
	movsd	40(%r12,%rsi,8), %xmm1
	mulsd	%xmm0, %xmm1
	movsd	%xmm1, 40(%r12,%rsi,8)
.LBB0_47:
	movsd	32(%r12,%rsi,8), %xmm1
	mulsd	%xmm0, %xmm1
	movsd	%xmm1, 32(%r12,%rsi,8)
.LBB0_48:
	movsd	24(%r12,%rsi,8), %xmm1
	mulsd	%xmm0, %xmm1
	movsd	%xmm1, 24(%r12,%rsi,8)
.LBB0_49:
	movsd	16(%r12,%rsi,8), %xmm1
	mulsd	%xmm0, %xmm1
	movsd	%xmm1, 16(%r12,%rsi,8)
.LBB0_50:
	movsd	8(%r12,%rsi,8), %xmm1
	mulsd	%xmm0, %xmm1
	movsd	%xmm1, 8(%r12,%rsi,8)
.LBB0_51:
	mulsd	(%r12,%rsi,8), %xmm0
	movsd	%xmm0, (%r12,%rsi,8)
.LBB0_22:
	movq	-112(%rbp), %rax
.Ltmp46:
	testb	$1, (%rax)
	je	.LBB0_61
	testq	%rdx, %rdx
	jle	.LBB0_61
	cmpq	$12, %rdx
	jbe	.LBB0_59
.Ltmp47:
	shlq	$3, %rdx
	movq	%r15, %rdi
	xorl	%esi, %esi
	addq	$72, %rsp
	popq	%rbx
	popq	%r12
	popq	%r13
	popq	%r14
	popq	%r15
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	jmp	_intel_fast_memset@PLT
.LBB0_59:
	.cfi_def_cfa %rbp, 16
	xorl	%eax, %eax
	.p2align	4
.LBB0_60:
	movq	$0, (%r15,%rax,8)
.Ltmp48:
	incq	%rax
	cmpq	%rax, %rdx
	jne	.LBB0_60
.LBB0_61:
.Ltmp49:
	addq	$72, %rsp
	popq	%rbx
	popq	%r12
	popq	%r13
	popq	%r14
	popq	%r15
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	retq
.Ltmp50:
.Lfunc_end0:
	.size	oed__ovl_prepare_ctr_, .Lfunc_end0-oed__ovl_prepare_ctr_
	.cfi_endproc
	.section	.rodata,"a",@progbits
	.p2align	3, 0x0
.LJTI0_0:
	.quad	.LBB0_22
	.quad	.LBB0_51
	.quad	.LBB0_50
	.quad	.LBB0_49
	.quad	.LBB0_48
	.quad	.LBB0_47
	.quad	.LBB0_46
	.quad	.LBB0_45
.LJTI0_1:
	.quad	.LBB0_22
	.quad	.LBB0_58
	.quad	.LBB0_57
	.quad	.LBB0_56
	.quad	.LBB0_55
	.quad	.LBB0_54
	.quad	.LBB0_53
	.quad	.LBB0_52

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
	.byte	1
	.byte	9
	.long	(.Ltmp4-.Ltmp3)-1
	.byte	4
	.byte	8
	.byte	9
	.long	(.Ltmp5-.Ltmp4)-1
	.byte	4
	.byte	1
	.byte	9
	.long	(.Ltmp6-.Ltmp5)-1
	.byte	4
	.byte	6
	.byte	9
	.long	(.Ltmp7-.Ltmp6)-1
	.byte	4
	.byte	1
	.byte	9
	.long	(.Ltmp8-.Ltmp7)-1
	.byte	4
	.byte	1
	.byte	9
	.long	(.Ltmp9-.Ltmp8)-1
	.byte	4
	.byte	-2
	.byte	9
	.long	(.Ltmp10-.Ltmp9)-1
	.byte	4
	.byte	-5
	.byte	9
	.long	(.Ltmp11-.Ltmp10)-1
	.byte	4
	.byte	1
	.byte	9
	.long	(.Ltmp12-.Ltmp11)-1
	.byte	4
	.byte	1
	.byte	9
	.long	(.Ltmp13-.Ltmp12)-1
	.byte	4
	.byte	1
	.byte	9
	.long	(.Ltmp14-.Ltmp13)-1
	.byte	4
	.byte	-3
	.byte	9
	.long	(.Ltmp15-.Ltmp14)-1
	.byte	4
	.byte	6
	.byte	9
	.long	(.Ltmp16-.Ltmp15)-1
	.byte	4
	.byte	-18
	.byte	9
	.long	(.Ltmp17-.Ltmp16)-1
	.byte	4
	.byte	18
	.byte	9
	.long	(.Ltmp18-.Ltmp17)-1
	.byte	4
	.byte	1
	.byte	9
	.long	(.Ltmp19-.Ltmp18)-1
	.byte	4
	.byte	2
	.byte	9
	.long	(.Ltmp20-.Ltmp19)-1
	.byte	4
	.byte	-1
	.byte	9
	.long	(.Ltmp21-.Ltmp20)-1
	.byte	4
	.byte	1
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
	.byte	-2
	.byte	9
	.long	(.Ltmp25-.Ltmp24)-1
	.byte	4
	.byte	-8
	.byte	9
	.long	(.Ltmp26-.Ltmp25)-1
	.byte	4
	.byte	-13
	.byte	9
	.long	(.Ltmp27-.Ltmp26)-1
	.byte	4
	.byte	13
	.byte	9
	.long	(.Ltmp28-.Ltmp27)-1
	.byte	4
	.byte	1
	.byte	9
	.long	(.Ltmp29-.Ltmp28)-1
	.byte	4
	.byte	1
	.byte	9
	.long	(.Ltmp30-.Ltmp29)-1
	.byte	4
	.byte	7
	.byte	9
	.long	(.Ltmp31-.Ltmp30)-1
	.byte	4
	.byte	-22
	.byte	9
	.long	(.Ltmp32-.Ltmp31)-1
	.byte	4
	.byte	22
	.byte	9
	.long	(.Ltmp33-.Ltmp32)-1
	.byte	4
	.byte	1
	.byte	9
	.long	(.Ltmp34-.Ltmp33)-1
	.byte	4
	.byte	8
	.byte	9
	.long	(.Ltmp35-.Ltmp34)-1
	.byte	4
	.byte	2
	.byte	9
	.long	(.Ltmp36-.Ltmp35)-1
	.byte	4
	.byte	5
	.byte	9
	.long	(.Ltmp37-.Ltmp36)-1
	.byte	4
	.byte	-38
	.byte	9
	.long	(.Ltmp38-.Ltmp37)-1
	.byte	4
	.byte	39
	.byte	9
	.long	(.Ltmp39-.Ltmp38)-1
	.byte	4
	.byte	1
	.byte	9
	.long	(.Ltmp40-.Ltmp39)-1
	.byte	4
	.byte	-1
	.byte	9
	.long	(.Ltmp41-.Ltmp40)-1
	.byte	4
	.byte	-5
	.byte	9
	.long	(.Ltmp42-.Ltmp41)-1
	.byte	4
	.byte	-34
	.byte	9
	.long	(.Ltmp43-.Ltmp42)-1
	.byte	4
	.byte	35
	.byte	9
	.long	(.Ltmp44-.Ltmp43)-1
	.byte	4
	.byte	1
	.byte	9
	.long	(.Ltmp45-.Ltmp44)-1
	.byte	4
	.byte	-1
	.byte	9
	.long	(.Ltmp46-.Ltmp45)-1
	.byte	4
	.byte	12
	.byte	9
	.long	(.Ltmp47-.Ltmp46)-1
	.byte	4
	.byte	2
	.byte	9
	.long	(.Ltmp48-.Ltmp47)-1
	.byte	4
	.byte	1
	.byte	9
	.long	(.Ltmp49-.Ltmp48)-1
	.byte	4
	.byte	8
	.byte	9
	.long	(.Lfunc_end0-.Ltmp49)-1
.Lmodule_end0:
	.section	".note.GNU-stack","",@progbits
