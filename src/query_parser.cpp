/*******************************************************************************
 * Copyright (c) 2017 Wojciech Migda
 * All rights reserved
 * Distributed under the terms of the MIT License
 *******************************************************************************
 *
 * Filename: query_parser.cpp
 *
 * Description:
 *      description
 *
 * Authors:
 *          Wojciech Migda (wm)
 *
 *******************************************************************************
 * History:
 * --------
 * Date         Who  Ticket     Description
 * ----------   ---  ---------  ------------------------------------------------
 * 2017-02-01   wm              Initial version
 *
 ******************************************************************************/

//#define USE_QI_PARSER
#if defined USE_QI_PARSER
#include "boost/spirit/include/qi_action.hpp"
#include "boost/spirit/include/qi_operator.hpp"

#include "boost/spirit/include/qi_parse.hpp"
#include "boost/spirit/include/qi_uint.hpp"
#include "boost/spirit/include/qi_eps.hpp"
#include "boost/spirit/include/qi_char.hpp"
#endif


#include "query_parser.hpp"

#include "likely.h"

#include <cstddef>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <string>

struct query_view
{
    std::uint32_t * vec_p;
    std::uint32_t size;
    std::uint32_t capacity;
};

extern "C"
bool parse_query_(const char * first, std::size_t sz, query_view * optr);

#if defined USE_QI_PARSER
bool parse_query_(const char * first, std::size_t sz, query_view * optr)
{
    using namespace boost::spirit;

    const char * last = first + sz;

    bool ok = qi::phrase_parse(
         first, last,
         (
             qi::uint_[([&optr](std::uint32_t x)
                 {
                     if (UNLIKELY(optr->size == optr->capacity))
                     {
                         optr->capacity *= 2;
                         optr->vec_p = static_cast<std::uint32_t *>(realloc(optr->vec_p, optr->capacity));
                     }
                     else
                     {
                         optr->vec_p[optr->size++] = x;
                     }
                 })]
         ) % ',',
         qi::ascii::blank
     );

    return ok;
}
#endif

std::vector<std::uint32_t> parse_query(std::string const & s)
{
    constexpr auto INIT_CAP = 250u;

    struct query_view view = {static_cast<std::uint32_t *>(malloc(INIT_CAP * sizeof (std::uint32_t))), 0, INIT_CAP};

    parse_query_(s.c_str(), s.size(), &view);

    std::vector<std::uint32_t> ret(view.vec_p, view.vec_p + view.size);

    free(view.vec_p);

    return ret;
}

#ifndef USE_QI_PARSER
// version with qi::ascii::blank
asm(R"(
    .section    .text.unlikely,"ax",@progbits
parse_query_.LCOLDB0:
    .text
parse_query_.LHOTB0:
    .p2align 4,,15
    .globl  parse_query_
    .type   parse_query_, @function
parse_query_:
parse_query_.LFB13551:
    .cfi_startproc
    pushq   %r14
    .cfi_def_cfa_offset 16
    .cfi_offset 14, -16
    pushq   %r13
    .cfi_def_cfa_offset 24
    .cfi_offset 13, -24
    pushq   %r12
    .cfi_def_cfa_offset 32
    .cfi_offset 12, -32
    pushq   %rbp
    .cfi_def_cfa_offset 40
    .cfi_offset 6, -40
    leaq    (%rdi,%rsi), %rbp
    pushq   %rbx
    .cfi_def_cfa_offset 48
    .cfi_offset 3, -48
    subq    $16, %rsp
    .cfi_def_cfa_offset 64
    movq    %rdx, 8(%rsp)
    cmpq    %rdi, %rbp
    je  parse_query_.L2
    movq    %rdi, %rbx
    jmp parse_query_.L4
    .p2align 4,,10
    .p2align 3
parse_query_.L100:
    addq    $1, %rbx
    cmpq    %rbx, %rbp
    je  parse_query_.L2
parse_query_.L4:
    movsbl  (%rbx), %eax
    cmpl    $9, %eax
    sete    %dl
    cmpl    $32, %eax
    sete    %al
    orb %dl, %al
    jne parse_query_.L100
    cmpq    %rbx, %rbp
    je  parse_query_.L91
    xorl    %ecx, %ecx
    .p2align 4,,10
    .p2align 3
parse_query_.L8:
    movsbl  (%rbx), %edx
    cmpb    $48, %dl
    jne parse_query_.L6
    addq    $1, %rbx
    addq    $1, %rcx
    cmpq    %rbx, %rbp
    jne parse_query_.L8
    testq   %rcx, %rcx
    je  parse_query_.L91
    movq    %rbp, %rbx
    xorl    %edx, %edx
parse_query_.L9:
    movq    8(%rsp), %r12
    movl    8(%r12), %eax
    cmpl    12(%r12), %eax
    je  parse_query_.L101
    movq    (%r12), %rcx
    leal    1(%rax), %esi
    movl    %esi, 8(%r12)
    movl    %edx, (%rcx,%rax,4)
parse_query_.L19:
    leaq    -2(%rbp), %r13
    leaq    -1(%rbp), %r12
parse_query_.L47:
    cmpq    %rbx, %rbp
    je  parse_query_.L51
    .p2align 4,,10
    .p2align 3
parse_query_.L23:
    movsbl  (%rbx), %eax
    cmpl    $9, %eax
    je  parse_query_.L69
    cmpl    $32, %eax
    je  parse_query_.L69
    cmpq    %rbx, %rbp
    je  parse_query_.L51
    cmpb    $44, %al
    jne parse_query_.L51
parse_query_.L93:
    addq    $1, %rbx
    cmpq    %rbx, %rbp
    je  parse_query_.L51
    movsbl  (%rbx), %eax
    cmpl    $9, %eax
    je  parse_query_.L93
    cmpl    $32, %eax
    je  parse_query_.L93
    cmpq    %rbx, %rbp
    je  parse_query_.L51
    xorl    %edx, %edx
    .p2align 4,,10
    .p2align 3
parse_query_.L30:
    movsbl  (%rbx), %eax
    cmpb    $48, %al
    jne parse_query_.L28
    addq    $1, %rbx
    addq    $1, %rdx
    cmpq    %rbx, %rbp
    jne parse_query_.L30
    testq   %rdx, %rdx
    je  parse_query_.L51
    movq    %rbp, %rbx
    xorl    %eax, %eax
parse_query_.L31:
    movq    8(%rsp), %r14
    movl    8(%r14), %edx
    cmpl    12(%r14), %edx
    je  parse_query_.L102
    movq    (%r14), %rcx
    leal    1(%rdx), %esi
    movl    %esi, 8(%r14)
    movl    %eax, (%rcx,%rdx,4)
    jmp parse_query_.L47
    .p2align 4,,10
    .p2align 3
parse_query_.L2:
    xorl    %eax, %eax
parse_query_.L91:
    addq    $16, %rsp
    .cfi_remember_state
    .cfi_def_cfa_offset 48
    popq    %rbx
    .cfi_def_cfa_offset 40
    popq    %rbp
    .cfi_def_cfa_offset 32
    popq    %r12
    .cfi_def_cfa_offset 24
    popq    %r13
    .cfi_def_cfa_offset 16
    popq    %r14
    .cfi_def_cfa_offset 8
    ret
    .p2align 4,,10
    .p2align 3
parse_query_.L69:
    .cfi_restore_state
    addq    $1, %rbx
    cmpq    %rbx, %rbp
    jne parse_query_.L23
parse_query_.L51:
    addq    $16, %rsp
    .cfi_remember_state
    .cfi_def_cfa_offset 48
    movl    $1, %eax
    popq    %rbx
    .cfi_def_cfa_offset 40
    popq    %rbp
    .cfi_def_cfa_offset 32
    popq    %r12
    .cfi_def_cfa_offset 24
    popq    %r13
    .cfi_def_cfa_offset 16
    popq    %r14
    .cfi_def_cfa_offset 8
    ret
    .p2align 4,,10
    .p2align 3
parse_query_.L28:
    .cfi_restore_state
    leal    -48(%rax), %ecx
    cmpb    $9, %cl
    jbe parse_query_.L48
    testq   %rdx, %rdx
    je  parse_query_.L51
    xorl    %eax, %eax
    jmp parse_query_.L31
parse_query_.L6:
    leal    -48(%rdx), %esi
    cmpb    $9, %sil
    jbe parse_query_.L103
    testq   %rcx, %rcx
    je  parse_query_.L91
    xorl    %edx, %edx
    jmp parse_query_.L9
parse_query_.L48:
    leaq    1(%rbx), %rdx
    subl    $48, %eax
    cmpq    %rdx, %rbp
    je  parse_query_.L67
    movsbl  1(%rbx), %ecx
    leal    -48(%rcx), %esi
    cmpb    $9, %sil
    ja  parse_query_.L62
    xorl    %esi, %esi
    jmp parse_query_.L33
    .p2align 4,,10
    .p2align 3
parse_query_.L104:
    leal    (%rax,%rax,4), %eax
    leal    -48(%rcx,%rax,2), %eax
parse_query_.L35:
    leaq    1(%rdx), %r8
    cmpq    %r12, %rdx
    je  parse_query_.L66
    movsbl  1(%rdx), %ecx
    leal    -48(%rcx), %edi
    cmpb    $9, %dil
    ja  parse_query_.L66
    movq    %rdx, %rdi
    subq    %rbx, %rdi
    cmpq    $7, %rdi
    ja  parse_query_.L36
    leal    (%rax,%rax,4), %eax
    leal    -48(%rcx,%rax,2), %eax
parse_query_.L37:
    leaq    2(%rdx), %r8
    cmpq    %r13, %rdx
    je  parse_query_.L66
    movsbl  2(%rdx), %ecx
    leal    -48(%rcx), %edi
    cmpb    $9, %dil
    ja  parse_query_.L66
    leaq    2(%rsi), %rdi
    cmpq    $7, %rdi
    ja  parse_query_.L38
    leal    (%rax,%rax,4), %eax
    leal    -48(%rcx,%rax,2), %eax
parse_query_.L39:
    addq    $3, %rdx
    addq    $3, %rsi
    cmpq    %rdx, %rbp
    je  parse_query_.L67
    movsbl  (%rdx), %ecx
    leal    -48(%rcx), %edi
    cmpb    $9, %dil
    ja  parse_query_.L62
parse_query_.L33:
    cmpq    $7, %rsi
    jbe parse_query_.L104
    cmpl    $429496729, %eax
    ja  parse_query_.L51
    leal    (%rax,%rax,4), %eax
    leal    (%rax,%rax), %edi
    movsbl  %cl, %eax
    subl    $48, %eax
    movl    %eax, %ecx
    notl    %ecx
    cmpl    %ecx, %edi
    ja  parse_query_.L51
    addl    %edi, %eax
    jmp parse_query_.L35
    .p2align 4,,10
    .p2align 3
parse_query_.L36:
    cmpl    $429496729, %eax
    ja  parse_query_.L51
    leal    (%rax,%rax,4), %eax
    leal    (%rax,%rax), %edi
    movsbl  %cl, %eax
    subl    $48, %eax
    movl    %eax, %ecx
    notl    %ecx
    cmpl    %ecx, %edi
    ja  parse_query_.L51
    addl    %edi, %eax
    jmp parse_query_.L37
    .p2align 4,,10
    .p2align 3
parse_query_.L38:
    cmpl    $429496729, %eax
    ja  parse_query_.L51
    leal    (%rax,%rax,4), %eax
    leal    (%rax,%rax), %edi
    movsbl  %cl, %eax
    subl    $48, %eax
    movl    %eax, %ecx
    notl    %ecx
    cmpl    %ecx, %edi
    ja  parse_query_.L51
    addl    %edi, %eax
    jmp parse_query_.L39
parse_query_.L103:
    leaq    1(%rbx), %rcx
    subl    $48, %edx
    cmpq    %rcx, %rbp
    je  parse_query_.L59
    movsbl  1(%rbx), %esi
    leal    -48(%rsi), %edi
    cmpb    $9, %dil
    ja  parse_query_.L54
    leaq    -2(%rbp), %r9
    xorl    %edi, %edi
    leaq    -1(%rbp), %r8
    jmp parse_query_.L11
    .p2align 4,,10
    .p2align 3
parse_query_.L105:
    leal    (%rdx,%rdx,4), %edx
    leal    -48(%rsi,%rdx,2), %edx
parse_query_.L13:
    leaq    1(%rcx), %r11
    cmpq    %r8, %rcx
    je  parse_query_.L58
    movsbl  1(%rcx), %esi
    leal    -48(%rsi), %r10d
    cmpb    $9, %r10b
    ja  parse_query_.L58
    movq    %rcx, %r10
    subq    %rbx, %r10
    cmpq    $7, %r10
    ja  parse_query_.L14
    leal    (%rdx,%rdx,4), %edx
    leal    -48(%rsi,%rdx,2), %edx
parse_query_.L15:
    leaq    2(%rcx), %r11
    cmpq    %r9, %rcx
    je  parse_query_.L58
    movsbl  2(%rcx), %esi
    leal    -48(%rsi), %r10d
    cmpb    $9, %r10b
    ja  parse_query_.L58
    leaq    2(%rdi), %r10
    cmpq    $7, %r10
    ja  parse_query_.L16
    leal    (%rdx,%rdx,4), %edx
    leal    -48(%rsi,%rdx,2), %edx
parse_query_.L17:
    addq    $3, %rcx
    addq    $3, %rdi
    cmpq    %rcx, %rbp
    je  parse_query_.L59
    movsbl  (%rcx), %esi
    leal    -48(%rsi), %r10d
    cmpb    $9, %r10b
    ja  parse_query_.L54
parse_query_.L11:
    cmpq    $7, %rdi
    jbe parse_query_.L105
    cmpl    $429496729, %edx
    ja  parse_query_.L91
    leal    (%rdx,%rdx,4), %edx
    leal    (%rdx,%rdx), %r10d
    movsbl  %sil, %edx
    subl    $48, %edx
    movl    %edx, %esi
    notl    %esi
    cmpl    %esi, %r10d
    ja  parse_query_.L91
    addl    %r10d, %edx
    jmp parse_query_.L13
parse_query_.L14:
    cmpl    $429496729, %edx
    ja  parse_query_.L91
    leal    (%rdx,%rdx,4), %edx
    leal    (%rdx,%rdx), %r10d
    movsbl  %sil, %edx
    subl    $48, %edx
    movl    %edx, %esi
    notl    %esi
    cmpl    %esi, %r10d
    ja  parse_query_.L91
    addl    %r10d, %edx
    jmp parse_query_.L15
parse_query_.L16:
    cmpl    $429496729, %edx
    ja  parse_query_.L91
    leal    (%rdx,%rdx,4), %edx
    leal    (%rdx,%rdx), %r10d
    movsbl  %sil, %edx
    subl    $48, %edx
    movl    %edx, %esi
    notl    %esi
    cmpl    %esi, %r10d
    ja  parse_query_.L91
    addl    %r10d, %edx
    jmp parse_query_.L17
parse_query_.L66:
    movq    %r8, %rbx
    jmp parse_query_.L31
parse_query_.L58:
    movq    %r11, %rbx
    jmp parse_query_.L9
parse_query_.L62:
    movq    %rdx, %rbx
    jmp parse_query_.L31
parse_query_.L67:
    movq    %rbp, %rbx
    jmp parse_query_.L31
parse_query_.L102:
    movq    (%r14), %rdi
    leal    (%rdx,%rdx), %esi
    movl    %esi, 12(%r14)
    call    realloc
    movq    %rax, (%r14)
    jmp parse_query_.L47
parse_query_.L101:
    movq    (%r12), %rdi
    leal    (%rax,%rax), %esi
    movl    %esi, 12(%r12)
    call    realloc
    movq    %rax, (%r12)
    jmp parse_query_.L19
parse_query_.L54:
    movq    %rcx, %rbx
    jmp parse_query_.L9
parse_query_.L59:
    movq    %rbp, %rbx
    jmp parse_query_.L9
    .cfi_endproc
parse_query_.LFE13551:
    .size   parse_query_, .-parse_query_
    .section    .text.unlikely
parse_query_.LCOLDE0:
    .text
parse_query_.LHOTE0:
    .section    .text.unlikely
parse_query_.LCOLDB1:
    .section    .text.startup,"ax",@progbits
parse_query_.LHOTB1:
    .p2align 4,,15
    .type   _GLOBAL__sub_I_parse_query, @function
_GLOBAL__sub_I_parse_query:
parse_query_.LFB14939:
    .cfi_startproc
    subq    $8, %rsp
    .cfi_def_cfa_offset 16
    movl    $_ZStL8__ioinit, %edi
    call    _ZNSt8ios_base4InitC1Ev
    movl    $__dso_handle, %edx
    movl    $_ZStL8__ioinit, %esi
    movl    $_ZNSt8ios_base4InitD1Ev, %edi
    addq    $8, %rsp
    .cfi_def_cfa_offset 8
    jmp __cxa_atexit
    .cfi_endproc
parse_query_.LFE14939:
    .size   _GLOBAL__sub_I_parse_query, .-_GLOBAL__sub_I_parse_query
    .section    .text.unlikely
parse_query_.LCOLDE1:
    .section    .text.startup
parse_query_.LHOTE1:
    .section    .init_array,"aw"
    .align 8
    .quad   _GLOBAL__sub_I_parse_query
    .local  _ZStL8__ioinit
    .comm   _ZStL8__ioinit,1,1
    .hidden __dso_handle
)");
#endif
