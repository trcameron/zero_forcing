%token _BANG_t
%token _BANG_EQUALS_t
%token _AMPER_t
%token _AMPERAMPER_t
%token _LPAREN_t
%token _RPAREN_t
%token _PLUS_EQUALS_t
%token _COLON_t
%token _SEMIC_t
%token _LANGLE_t
%token _LANGLE_EQUALS_t
%token _EQUALS_t
%token _RANGLE_t
%token _RANGLE_EQUALS_t
%token _QUESTION_EQUALS_t
%token _LBRACKET_t
%token _RBRACKET_t
%token ACTIONS_t
%token BIND_t
%token BREAK_t
%token CASE_t
%token CLASS_t
%token CONTINUE_t
%token DEFAULT_t
%token ELSE_t
%token EXISTING_t
%token FOR_t
%token IF_t
%token IGNORE_t
%token IN_t
%token INCLUDE_t
%token LOCAL_t
%token MODULE_t
%token ON_t
%token PIECEMEAL_t
%token QUIETLY_t
%token RETURN_t
%token RULE_t
%token SWITCH_t
%token TOGETHER_t
%token UPDATED_t
%token WHILE_t
%token _LBRACE_t
%token _BAR_t
%token _BARBAR_t
%token _RBRACE_t
/*
 * Copyright 1993, 2000 Christopher Seiwald.
 *
 * This file is part of Jam - see jam.c for Copyright information.
 */

/*  This file is ALSO:
 *  Copyright 2001-2004 David Abrahams.
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or http://www.boost.org/LICENSE_1_0.txt)
 */

/*
 * jamgram.yy - jam grammar
 *
 * 04/13/94 (seiwald) - added shorthand L0 for null list pointer
 * 06/01/94 (seiwald) - new 'actions existing' does existing sources
 * 08/23/94 (seiwald) - Support for '+=' (append to variable)
 * 08/31/94 (seiwald) - Allow ?= as alias for "default =".
 * 09/15/94 (seiwald) - if conditionals take only single arguments, so
 *			that 'if foo == bar' gives syntax error (use =).
 * 02/11/95 (seiwald) - when scanning arguments to rules, only treat
 *			punctuation keywords as keywords.  All arg lists
 *			are terminated with punctuation keywords.
 *
 * 09/11/00 (seiwald) - Support for function calls:
 *
 *		Rules now return lists (LIST *), rather than void.
 *
 *		New "[ rule ]" syntax evals rule into a LIST.
 *
 *		Lists are now generated by compile_list() and
 *		compile_append(), and any other rule that indirectly
 *		makes a list, rather than being built directly here,
 *		so that lists values can contain rule evaluations.
 *
 *		New 'return' rule sets the return value, though
 *		other statements also may have return values.
 *
 *		'run' production split from 'block' production so
 *		that empty blocks can be handled separately.
 */

%token ARG STRING

%left _BARBAR_t _BAR_t
%left _AMPERAMPER_t _AMPER_t
%left _EQUALS_t _BANG_EQUALS_t IN_t
%left _LANGLE_t _LANGLE_EQUALS_t _RANGLE_t _RANGLE_EQUALS_t
%left _BANG_t

%{
#include "jam.h"

#include "lists.h"
#include "parse.h"
#include "scan.h"
#include "compile.h"
#include "object.h"
#include "rules.h"

# define YYINITDEPTH 5000 /* for C++ parsing */
# define YYMAXDEPTH 10000	/* for OSF and other less endowed yaccs */

# define F0 -1
# define P0 (PARSE *)0
# define S0 (OBJECT *)0

# define pappend( l,r )    	parse_make( PARSE_APPEND,l,r,P0,S0,S0,0 )
# define peval( c,l,r )	parse_make( PARSE_EVAL,l,r,P0,S0,S0,c )
# define pfor( s,l,r,x )    	parse_make( PARSE_FOREACH,l,r,P0,s,S0,x )
# define pif( l,r,t )	  	parse_make( PARSE_IF,l,r,t,S0,S0,0 )
# define pincl( l )       	parse_make( PARSE_INCLUDE,l,P0,P0,S0,S0,0 )
# define plist( s )	  	parse_make( PARSE_LIST,P0,P0,P0,s,S0,0 )
# define plocal( l,r,t )  	parse_make( PARSE_LOCAL,l,r,t,S0,S0,0 )
# define pmodule( l,r )	  	parse_make( PARSE_MODULE,l,r,P0,S0,S0,0 )
# define pclass( l,r )	  	parse_make( PARSE_CLASS,l,r,P0,S0,S0,0 )
# define pnull()	  	parse_make( PARSE_NULL,P0,P0,P0,S0,S0,0 )
# define pon( l,r )	  	parse_make( PARSE_ON,l,r,P0,S0,S0,0 )
# define prule( s,p )     	parse_make( PARSE_RULE,p,P0,P0,s,S0,0 )
# define prules( l,r )	  	parse_make( PARSE_RULES,l,r,P0,S0,S0,0 )
# define pset( l,r,a )          parse_make( PARSE_SET,l,r,P0,S0,S0,a )
# define pset1( l,r,t,a )	parse_make( PARSE_SETTINGS,l,r,t,S0,S0,a )
# define psetc( s,p,a,l )     	parse_make( PARSE_SETCOMP,p,a,P0,s,S0,l )
# define psete( s,l,s1,f ) 	parse_make( PARSE_SETEXEC,l,P0,P0,s,s1,f )
# define pswitch( l,r )   	parse_make( PARSE_SWITCH,l,r,P0,S0,S0,0 )
# define pwhile( l,r )   	parse_make( PARSE_WHILE,l,r,P0,S0,S0,0 )
# define preturn( l )       parse_make( PARSE_RETURN,l,P0,P0,S0,S0,0 )
# define pbreak()           parse_make( PARSE_BREAK,P0,P0,P0,S0,S0,0 )
# define pcontinue()        parse_make( PARSE_CONTINUE,P0,P0,P0,S0,S0,0 )

# define pnode( l,r )    	parse_make( F0,l,r,P0,S0,S0,0 )
# define psnode( s,l )     	parse_make( F0,l,P0,P0,s,S0,0 )

%}

%%

run	: /* empty */
		/* do nothing */
	| rules
		{ parse_save( $1.parse ); }
	;

/*
 * block - zero or more rules
 * rules - one or more rules
 * rule - any one of jam's rules
 * right-recursive so rules execute in order.
 */

block	: null
                { $$.parse = $1.parse; }
	| rules
		{ $$.parse = $1.parse; }
	;

rules	: rule
		{ $$.parse = $1.parse; }
	| rule rules
		{ $$.parse = prules( $1.parse, $2.parse ); }
	| LOCAL_t { yymode( SCAN_ASSIGN ); } list assign_list_opt _SEMIC_t { yymode( SCAN_NORMAL ); } block
		{ $$.parse = plocal( $3.parse, $4.parse, $7.parse ); }
	;

null    : /* empty */
        { $$.parse = pnull(); }
        ;

assign_list_opt : _EQUALS_t { yymode( SCAN_PUNCT ); } list
                { $$.parse = $3.parse; $$.number = ASSIGN_SET; }
        | null
		{ $$.parse = $1.parse; $$.number = ASSIGN_APPEND; }
        ;

arglist_opt : _LPAREN_t lol _RPAREN_t
                { $$.parse = $2.parse; }
        |
                { $$.parse = P0; }
        ;

local_opt : LOCAL_t
                { $$.number = 1; }
          | /* empty */
                { $$.number = 0; }
          ;

else_opt : ELSE_t rule
                { $$.parse = $2.parse; }
          | /* empty */
                { $$.parse = pnull(); }

rule	: _LBRACE_t block _RBRACE_t
		{ $$.parse = $2.parse; }
	| INCLUDE_t { yymode( SCAN_PUNCT ); } list _SEMIC_t
		{ $$.parse = pincl( $3.parse ); yymode( SCAN_NORMAL ); }
	| ARG { yymode( SCAN_PUNCT ); } lol _SEMIC_t
		{ $$.parse = prule( $1.string, $3.parse ); yymode( SCAN_NORMAL ); }
	| arg assign { yymode( SCAN_PUNCT ); } list _SEMIC_t
		{ $$.parse = pset( $1.parse, $4.parse, $2.number ); yymode( SCAN_NORMAL ); }
	| arg ON_t { yymode( SCAN_ASSIGN ); } list assign { yymode( SCAN_PUNCT ); } list _SEMIC_t
		{ $$.parse = pset1( $1.parse, $4.parse, $7.parse, $5.number ); yymode( SCAN_NORMAL ); }
	| RETURN_t { yymode( SCAN_PUNCT ); } list _SEMIC_t
		{ $$.parse = preturn( $3.parse ); yymode( SCAN_NORMAL ); }
    | BREAK_t _SEMIC_t
        { $$.parse = pbreak(); }
    | CONTINUE_t _SEMIC_t
        { $$.parse = pcontinue(); }
	| FOR_t local_opt ARG IN_t { yymode( SCAN_PUNCT ); } list _LBRACE_t { yymode( SCAN_NORMAL ); } block _RBRACE_t
		{ $$.parse = pfor( $3.string, $6.parse, $9.parse, $2.number ); }
	| SWITCH_t { yymode( SCAN_PUNCT ); } list _LBRACE_t { yymode( SCAN_NORMAL ); } cases _RBRACE_t
		{ $$.parse = pswitch( $3.parse, $6.parse ); }
	| IF_t { yymode( SCAN_CONDB ); } expr _LBRACE_t { yymode( SCAN_NORMAL ); } block _RBRACE_t else_opt
		{ $$.parse = pif( $3.parse, $6.parse, $8.parse ); }
	| MODULE_t { yymode( SCAN_PUNCT ); } list _LBRACE_t { yymode( SCAN_NORMAL ); } block _RBRACE_t
		{ $$.parse = pmodule( $3.parse, $6.parse ); }
	| CLASS_t { yymode( SCAN_PUNCT ); } lol _LBRACE_t { yymode( SCAN_NORMAL ); } block _RBRACE_t
		{ $$.parse = pclass( $3.parse, $6.parse ); }
	| WHILE_t { yymode( SCAN_CONDB ); } expr { yymode( SCAN_NORMAL ); } _LBRACE_t block _RBRACE_t
		{ $$.parse = pwhile( $3.parse, $6.parse ); }
     | local_opt RULE_t { yymode( SCAN_PUNCT ); } ARG { yymode( SCAN_PARAMS ); } arglist_opt { yymode( SCAN_NORMAL ); } rule
		{ $$.parse = psetc( $4.string, $8.parse, $6.parse, $1.number ); }
	| ON_t arg rule
		{ $$.parse = pon( $2.parse, $3.parse ); }
	| ACTIONS_t eflags ARG bindlist _LBRACE_t
		{ yymode( SCAN_STRING ); }
	  STRING
		{ yymode( SCAN_NORMAL ); }
	  _RBRACE_t
		{ $$.parse = psete( $3.string,$4.parse,$7.string,$2.number ); }
	;

/*
 * assign - = or +=
 */

assign	: _EQUALS_t
		{ $$.number = ASSIGN_SET; }
	| _PLUS_EQUALS_t
		{ $$.number = ASSIGN_APPEND; }
	| _QUESTION_EQUALS_t
		{ $$.number = ASSIGN_DEFAULT; }
	| DEFAULT_t _EQUALS_t
		{ $$.number = ASSIGN_DEFAULT; }
	;

/*
 * expr - an expression for if
 */
expr	: arg
		{ $$.parse = peval( EXPR_EXISTS, $1.parse, pnull() ); yymode( SCAN_COND ); }
	| expr _EQUALS_t { yymode( SCAN_CONDB ); } expr
		{ $$.parse = peval( EXPR_EQUALS, $1.parse, $4.parse ); }
	| expr _BANG_EQUALS_t { yymode( SCAN_CONDB ); } expr
		{ $$.parse = peval( EXPR_NOTEQ, $1.parse, $4.parse ); }
	| expr _LANGLE_t { yymode( SCAN_CONDB ); } expr
		{ $$.parse = peval( EXPR_LESS, $1.parse, $4.parse ); }
	| expr _LANGLE_EQUALS_t { yymode( SCAN_CONDB ); } expr
		{ $$.parse = peval( EXPR_LESSEQ, $1.parse, $4.parse ); }
	| expr _RANGLE_t { yymode( SCAN_CONDB ); } expr
		{ $$.parse = peval( EXPR_MORE, $1.parse, $4.parse ); }
	| expr _RANGLE_EQUALS_t { yymode( SCAN_CONDB ); } expr
		{ $$.parse = peval( EXPR_MOREEQ, $1.parse, $4.parse ); }
	| expr _AMPER_t { yymode( SCAN_CONDB ); } expr
		{ $$.parse = peval( EXPR_AND, $1.parse, $4.parse ); }
	| expr _AMPERAMPER_t { yymode( SCAN_CONDB ); } expr
		{ $$.parse = peval( EXPR_AND, $1.parse, $4.parse ); }
	| expr _BAR_t { yymode( SCAN_CONDB ); } expr
		{ $$.parse = peval( EXPR_OR, $1.parse, $4.parse ); }
	| expr _BARBAR_t { yymode( SCAN_CONDB ); } expr
		{ $$.parse = peval( EXPR_OR, $1.parse, $4.parse ); }
	| arg IN_t { yymode( SCAN_PUNCT ); } list
		{ $$.parse = peval( EXPR_IN, $1.parse, $4.parse ); yymode( SCAN_COND ); }
	| _BANG_t { yymode( SCAN_CONDB ); } expr
		{ $$.parse = peval( EXPR_NOT, $3.parse, pnull() ); }
	| _LPAREN_t { yymode( SCAN_CONDB ); } expr _RPAREN_t
		{ $$.parse = $3.parse; }
	;


/*
 * cases - action elements inside a 'switch'
 * case - a single action element inside a 'switch'
 * right-recursive rule so cases can be examined in order.
 */

cases	: /* empty */
		{ $$.parse = P0; }
	| case cases
		{ $$.parse = pnode( $1.parse, $2.parse ); }
	;

case	: CASE_t { yymode( SCAN_CASE ); } ARG _COLON_t { yymode( SCAN_NORMAL ); } block
		{ $$.parse = psnode( $3.string, $6.parse ); }
	;

/*
 * lol - list of lists
 * right-recursive rule so that lists can be added in order.
 */

lol	: list
		{ $$.parse = pnode( P0, $1.parse ); }
	| list _COLON_t lol
		{ $$.parse = pnode( $3.parse, $1.parse ); }
	;

/*
 * list - zero or more args in a LIST
 * listp - list (in puncutation only mode)
 * arg - one ARG or function call
 */

list	: listp
		{ $$.parse = $1.parse; }
	;

listp	: /* empty */
		{ $$.parse = pnull(); }
	| listp arg
        	{ $$.parse = pappend( $1.parse, $2.parse ); }
	;

arg	: ARG
		{ $$.parse = plist( $1.string ); }
	| _LBRACKET_t { $$.number = yymode( SCAN_CALL ); } func _RBRACKET_t
		{ $$.parse = $3.parse; yymode( $2.number ); }
	;

/*
 * func - a function call (inside [])
 * This needs to be split cleanly out of 'rule'
 */

func	: ARG { yymode( SCAN_PUNCT ); } lol
		{ $$.parse = prule( $1.string, $3.parse ); }
	| ON_t arg ARG { yymode( SCAN_PUNCT ); } lol
		{ $$.parse = pon( $2.parse, prule( $3.string, $5.parse ) ); }
	| ON_t arg RETURN_t { yymode( SCAN_PUNCT ); } list
		{ $$.parse = pon( $2.parse, $5.parse ); }
	;


/*
 * eflags - zero or more modifiers to 'executes'
 * eflag - a single modifier to 'executes'
 */

eflags	: /* empty */
		{ $$.number = 0; }
	| eflags eflag
		{ $$.number = $1.number | $2.number; }
	;

eflag	: UPDATED_t
		{ $$.number = EXEC_UPDATED; }
	| TOGETHER_t
		{ $$.number = EXEC_TOGETHER; }
	| IGNORE_t
		{ $$.number = EXEC_IGNORE; }
	| QUIETLY_t
		{ $$.number = EXEC_QUIETLY; }
	| PIECEMEAL_t
		{ $$.number = EXEC_PIECEMEAL; }
	| EXISTING_t
		{ $$.number = EXEC_EXISTING; }
	;


/*
 * bindlist - list of variable to bind for an action
 */

bindlist : /* empty */
		{ $$.parse = pnull(); }
	| BIND_t { yymode( SCAN_PUNCT ); } list
		{ $$.parse = $3.parse; }
	;

