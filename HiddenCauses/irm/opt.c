#define OPT_PACKAGE "opt"
#define OPT_VERSION "3.19"
#include "opt.h"
/* OPT v3.19: options parsing tool */
/*
 * This SOFTWARE has been authored by an employee of the University of
 * California, operator of the Los Alamos National Laboratory under
 * Contract No. W-7405-ENG-36 with the U.S. Department of Energy. The
 * U.S. Government has rights to use, reproduce, and distribute this
 * SOFTWARE. Neither the Government nor the University makes any
 * warranty, express or implied, or assumes any liability or
 * responsibility for the use of this SOFTWARE. If SOFTWARE is modified
 * to produce derivative works, such modified SOFTWARE should be clearly
 * marked, so as not to confuse it with the version available from LANL.
 * 
 * Additionally, this program is free software; you can distribute it
 * and/or modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of the
 * License, or any later version.  Accordingly, this program is
 * distributed in the hope that it will be useful, but WITHOUT A
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 * for more details (http://www.gnu.org/copyleft/gpl.txt).
 * 
 * 
 * The OPT package is available from
 * http://nis-www.lanl.gov/~jt/Software
 */
/******************************************************************

   This SOFTWARE has been authored by an employee of the University of
   California, operator of the Los Alamos National Laboratory under
   Contract No. W-7405-ENG-36 with the U.S. Department of Energy. The
   U.S. Government has rights to use, reproduce, and distribute this
   SOFTWARE. Neither the Government nor the University makes any
   warranty, express or implied, or assumes any liability or
   responsibility for the use of this SOFTWARE. If SOFTWARE is
   modified to produce derivative works, such modified SOFTWARE should
   be clearly marked, so as not to confuse it with the version
   available from LANL.

   Additionally, this program is free software; you can distribute it
   and/or modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2 of
   the License, or any later version.  Accordingly, this program is
   distributed in the hope that it will be useful, but WITHOUT A
   WARRANTY; without even the implied warranty of MERCHANTABILITY or
   FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
   License for more details (http://www.gnu.org/copyleft/gpl.txt).

   You can find the latest version of OPT on the web
   http://nis-www.lanl.gov/~jt/Software
   Or you can email the author: jt@lanl.gov

*********************************************************************/

/* $Id: opt.c,v 1.1.1.1 2005/10/29 19:20:46 ckemp Exp $ */
#define SINGLEFILE 1

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

/** You may want to edit the #define's below **/

#ifndef HAVE_STRDUP
#define HAVE_STRDUP 1
#endif

#ifndef HAVE_UNISTD_H
#define HAVE_UNISTD_H 0
#endif

#ifndef HAVE_LONGJMP
#define HAVE_LONGJMP 0
#endif

#ifndef HAVE_SIGNAL
#define HAVE_SIGNAL 0
#endif

#ifndef RETSIGTYPE
#define RETSIGTYPE void
#endif


#ifdef HAVE_READLINE
#define READLINE_LIBRARY
#include <readline/readline.h>
#include <readline/history.h>
#endif

/** end opt_sf.h **/





/* OPT v3.19: options parsing tool */
/*
 * This SOFTWARE has been authored by an employee of the University of
 * California, operator of the Los Alamos National Laboratory under
 * Contract No. W-7405-ENG-36 with the U.S. Department of Energy. The
 * U.S. Government has rights to use, reproduce, and distribute this
 * SOFTWARE. Neither the Government nor the University makes any
 * warranty, express or implied, or assumes any liability or
 * responsibility for the use of this SOFTWARE. If SOFTWARE is modified
 * to produce derivative works, such modified SOFTWARE should be clearly
 * marked, so as not to confuse it with the version available from LANL.
 * 
 * Additionally, this program is free software; you can distribute it
 * and/or modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of the
 * License, or any later version.  Accordingly, this program is
 * distributed in the hope that it will be useful, but WITHOUT A
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 * for more details (http://www.gnu.org/copyleft/gpl.txt).
 * 
 * 
 * The OPT package is available from
 * http://nis-www.lanl.gov/~jt/Software
 */
/* $Id: opt.c,v 1.1.1.1 2005/10/29 19:20:46 ckemp Exp $ */
/* ag.h */


#ifndef _AG_H
#define _AG_H

#ifndef	TRUE
#define	TRUE	1
#define	FALSE	0
#endif

/*************************************
 * ARGVECTOR structure
 * Basically, an (argc,argv) construct
 * with indices to which character of 
 * which word is the current position
 */

typedef struct {
    char    **v;            /* argument vector */
    int     c;              /* argument count  */
    int     iw;             /* current word    */
    int     ic;             /* current character */
}   ARGVECTOR;

/* 
 * Function prototypes 
 */


extern  ARGVECTOR *ag_new(int,char **);
extern  void       ag_free(ARGVECTOR *);
extern  void       ag_fprint(FILE *, ARGVECTOR *);
extern  int        ag_enstring(char *,ARGVECTOR *,int);
extern  void       ag_reset(ARGVECTOR *);
extern  int        ag_w_number(ARGVECTOR *);
extern  void       ag_w_advance(ARGVECTOR *);
extern  int        ag_eow(ARGVECTOR *);
extern  int        ag_end(ARGVECTOR *);
extern  char       ag_c(ARGVECTOR *);
extern  char       ag_c_next(ARGVECTOR *);
extern  char       ag_cnn_next(ARGVECTOR *);
extern  char       ag_c_advance(ARGVECTOR *);
extern  char       ag_backspace(ARGVECTOR *);
extern  void       ag_backword(ARGVECTOR *);
extern  char      *ag_s(ARGVECTOR *);
extern  char      *ag_s_next(ARGVECTOR *);
extern  char      *ag_s_advance(ARGVECTOR *);
extern  int        ag_clear(ARGVECTOR *);
extern  int        ag_argc(ARGVECTOR *);
extern  char     **ag_argv(ARGVECTOR *);
extern  ARGVECTOR *ag_copy(ARGVECTOR *);
extern  void       ag_freeall(ARGVECTOR *);
extern  ARGVECTOR *ag_prepend_argv0(ARGVECTOR *, char *);
extern  char      *argnext(ARGVECTOR *);


#endif /* _AG_H */








/* OPT v3.19: options parsing tool */
/*
 * This SOFTWARE has been authored by an employee of the University of
 * California, operator of the Los Alamos National Laboratory under
 * Contract No. W-7405-ENG-36 with the U.S. Department of Energy. The
 * U.S. Government has rights to use, reproduce, and distribute this
 * SOFTWARE. Neither the Government nor the University makes any
 * warranty, express or implied, or assumes any liability or
 * responsibility for the use of this SOFTWARE. If SOFTWARE is modified
 * to produce derivative works, such modified SOFTWARE should be clearly
 * marked, so as not to confuse it with the version available from LANL.
 * 
 * Additionally, this program is free software; you can distribute it
 * and/or modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of the
 * License, or any later version.  Accordingly, this program is
 * distributed in the hope that it will be useful, but WITHOUT A
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 * for more details (http://www.gnu.org/copyleft/gpl.txt).
 * 
 * 
 * The OPT package is available from
 * http://nis-www.lanl.gov/~jt/Software
 */
/* $Id: opt.c,v 1.1.1.1 2005/10/29 19:20:46 ckemp Exp $ */
/* opt_p.h */
/*
 *	Private header file for OPT package.
 */

#ifndef _OPT_P_H
#define _OPT_P_H

/***
 * Enable linking with C++
 ***/
#ifdef __cplusplus
extern "C"
{
#endif

#ifndef	OPT_TRUE
#define	OPT_TRUE	1
#define	OPT_FALSE	0
#endif

#ifndef DEBUG
#define DEBUG 0
#endif

/* 
 * These are delimiter characters 
 */
#define DELIM		'-'		/* option delimiter character */
#define	ALTDELIM	'/'		/* alternate delimiter character */
#define	OPTFROMFILE	'@'		/* denotes options are in a file */
#define	OPTTOFILE	'%'		/* says to put options in a file */
#define	DOCUMENT	'-'		/* write document to file */
#define	INTERACT	'$'		/* Flags interactive menu */
#define IGNOREEOL       ';'		/* Ignore until the end of line */

/*
 * These are not available on command line 
 * But may be invoked from a file or the menu
 */
#define	RUN             '='		/* Says to just run to completion */
#define	QUITCH		'.'		/* Quit character */
/*
 * These are not available on command line or from a file
 * But may be invoked from the menu
 */
#define	BANG		'!'		/* Shell escape character */
#define	ADDITIONAL_OPTS	'+'		/* Additional options */

/*
 * These are not available from a file, but may be used on
 * the command line or from the menu.
 */
#define	HELPCH		'?'		/* Help character */
#define REASONCH        '~'             /* Displays deprecation reason */

#define ISEMPTYSTRING(s) (!((s) && *(s)))

/******************
 * Global variables
 */
extern	int	opt_nreg;
extern	int	opt_menuflag;
extern  int     opt_exit_number;

/*********************
 * Function Prototypes
 */

extern	void	opt_get_help(char);
extern  char   *optgetTitle(void);
extern  int     opt_char_number(char);

extern  void	opt_tofile(FILE *);
extern	char	*optstrval(int);
extern	char	*optstrtyp(int);
extern  int     opt_parse_positional(ARGVECTOR *);
extern  int     opt_parse_delim(ARGVECTOR *);
extern  int     opt_parse_longdelim(ARGVECTOR *);

extern void short_usage(void);
extern void long_usage(void);

  
#define OPT_NUM_VALID   0x1
#define OPT_NUM_FLOAT   0x2
#define OPT_NUM_NEG     0x4

extern double        opt_atof(char *);
extern long          opt_atoi(char *);
extern unsigned long opt_atou(char *);
extern int           opt_isvalidnumber(char *);
extern void          opt_setstring(char **, char *);
extern char         *opt_justify(char* s, int width, int indent, 
                                 int initial, char* prefix); 

		/* --------------- */
		/* Process Options */
		/* --------------- */

extern	char	  *opt_program_name;
extern  OPT_PFI    opt_additional_usage_fcn; 

extern	void       opt_wr_title(void);
extern  int        opt_fromfname(char *);
extern	int	   opt_tofname(char *);
extern  int        opt_getline(char *,FILE *);
extern	char	  *opt_mstring(int);
extern  ARGVECTOR *opt_process(int,char **);
extern  char      *short_progname(char *);
extern  int        opt_lineprocess(char *);
extern	void	   opt_menu(void);
extern  void       opt_menu_getline(char *, char *, int);
extern  void       opt_freestrings(void);

extern	void	opt_help(char *);
extern	void	opt_quit(void);
extern	void	opt_usage(void);
extern  void    opt_readline_init(char *);

extern  void    opt_write_xml_file(FILE *fp);

#define opt_mess_1(fmt,var1)                printf(fmt,var1)
#define opt_mess_2(fmt,var1,var2)           printf(fmt,var1,var2)
#define opt_mess_3(fmt,var1,var2,var3)      printf(fmt,var1,var2,var3)
#define opt_mess_4(fmt,var1,var2,var3,var4) printf(fmt,var1,var2,var3,var4)

#if HAVE_SNPRINTF /* if snprintf exists, then use it */
#define opt_snprintf_1(buf,n,fmt,v1)          snprintf(buf,n,fmt,v1)
#define opt_snprintf_2(buf,n,fmt,v1,v2)       snprintf(buf,n,fmt,v1,v2)
#define opt_snprintf_3(buf,n,fmt,v1,v2,v3)    snprintf(buf,n,fmt,v1,v2,v3)
#define opt_snprintf_4(buf,n,fmt,v1,v2,v3,v4) snprintf(buf,n,fmt,v1,v2,v3,v4)
#else /* if no snprintf, then use sprintf and take your chances! */
#define opt_snprintf_1(buf,n,fmt,v1)           sprintf(buf,fmt,v1)
#define opt_snprintf_2(buf,n,fmt,v1,v2)        sprintf(buf,fmt,v1,v2)
#define opt_snprintf_3(buf,n,fmt,v1,v2,v3)     sprintf(buf,fmt,v1,v2,v3)
#define opt_snprintf_4(buf,n,fmt,v1,v2,v3,v4)  sprintf(buf,fmt,v1,v2,v3,v4)
#endif

#define OPT_ERRMAXSTRLEN 1024  /* shouldn't be fixed length, but it is! */
  
#define opt_warn_1(fmt,var1) do { \
    char gstr[OPT_ERRMAXSTRLEN]; \
        opt_snprintf_1(gstr,OPT_ERRMAXSTRLEN,fmt,var1); \
        opt_warning(gstr); } while(0)
#define opt_warn_2(fmt,var1,var2) do { \
    char gstr[OPT_ERRMAXSTRLEN]; \
        opt_snprintf_2(gstr,OPT_ERRMAXSTRLEN,fmt,var1,var2); \
	opt_warning(gstr); } while(0)
#define opt_warn_3(fmt,var1,var2,var3) do { \
    char gstr[OPT_ERRMAXSTRLEN]; \
        opt_snprintf_3(gstr,OPT_ERRMAXSTRLEN,fmt,var1,var2,var3); \
	opt_warning(gstr); } while(0)

#define opt_fatal_1(fmt,var1) do { \
    char gstr[OPT_ERRMAXSTRLEN]; \
        opt_snprintf_1(gstr,OPT_ERRMAXSTRLEN,fmt,var1); \
	opt_fatal(gstr); } while(0)
#define opt_fatal_2(fmt,var1,var2) do { \
    char gstr[OPT_ERRMAXSTRLEN]; \
        opt_snprintf_2(gstr,OPT_ERRMAXSTRLEN,fmt,var1,var2); \
	opt_fatal(gstr); } while(0)
#define opt_fatal_3(fmt,var1,var2,var3) do { \
    char gstr[OPT_ERRMAXSTRLEN]; \
        opt_snprintf_3(gstr,OPT_ERRMAXSTRLEN,fmt,var1,var2,var3); \
	opt_fatal(gstr); } while(0)

#define OPT_FREE(s) do { if (s) { free(s); s=NULL; } } while (0)


/* ...well, assume that you do have strdup, since if you do and you
 * say you don't you can get into trouble.  if you really don't, you
 * should '#define HAVE_STRDUP 0' (or -DHAVE_STRDUP=0 on command line)
 */
#ifndef HAVE_STRDUP
#define HAVE_STRDUP 1
#endif

#if !HAVE_STRDUP
char *strdup(char *s);
#endif

/***
 * Enable linking with C++
 ***/
#ifdef __cplusplus
}
#endif
#endif /* _OPT_P_H */






/* OPT v3.19: options parsing tool */
/*
 * This SOFTWARE has been authored by an employee of the University of
 * California, operator of the Los Alamos National Laboratory under
 * Contract No. W-7405-ENG-36 with the U.S. Department of Energy. The
 * U.S. Government has rights to use, reproduce, and distribute this
 * SOFTWARE. Neither the Government nor the University makes any
 * warranty, express or implied, or assumes any liability or
 * responsibility for the use of this SOFTWARE. If SOFTWARE is modified
 * to produce derivative works, such modified SOFTWARE should be clearly
 * marked, so as not to confuse it with the version available from LANL.
 * 
 * Additionally, this program is free software; you can distribute it
 * and/or modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of the
 * License, or any later version.  Accordingly, this program is
 * distributed in the hope that it will be useful, but WITHOUT A
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 * for more details (http://www.gnu.org/copyleft/gpl.txt).
 * 
 * 
 * The OPT package is available from
 * http://nis-www.lanl.gov/~jt/Software
 */
/* $Id: opt.c,v 1.1.1.1 2005/10/29 19:20:46 ckemp Exp $ */
/**************************
 *  some useful utilities *
 **************************/
/*
 *  *** NO LONGER AVAILABLE!!
 *  *** 1) Variable argument lists: for error messages
 *  ***   opt_errmess(), gives message and then aborts run
 *     
 *  2) Long Jump: for aborting run() without losing the entire job
 *     optBeginRun()
 *     optAbortRun()
 *     2a) Signal: trapping ^C so that it aborts run()
 *  4) Readline: interface to GNU readline
 *     
 */
/*  Can enable longjmp with: #define HAVE_LONGJMP 1
 *  Only if longjmp is enabled can signal be enabled, using
 *                           #define HAVE_SIGNAL 1
 *  Finally, readline can be enabled with
 *                           #define HAVE_READLINE 1
 *  If this version of opt came with a configure script, then those
 *  values should be set automatically by the system at configure time.
 */
#ifndef SINGLEFILE
#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ag.h"
#ifdef __cplusplus
}
#endif
#include "opt.h"
#include "opt_p.h"
#endif /* SINGLEFILE */


/* If these are not defined on the compile command line, then assume
 * the worst, that you don't have anything!
 */
#ifndef HAVE_LONGJMP
#define HAVE_LONGJMP 0
#endif
#ifndef HAVE_SIGNAL
#define HAVE_SIGNAL 0
#endif
#ifndef RETSIGTYPE
#define RETSIGTYPE void
#endif
#ifndef HAVE_READLINE
#define HAVE_READLINE 0
#endif

/* Provide message/warning/fatal functions */
/* Ordinary messages go to stdout, but warn/fatal messages to stderr;
 * (This convention adopted at suggestion of Jason V. Morgan) */
void
opt_message(char *s)
{
    if (s) fputs(s,stdout);
}
void
opt_warning(char *s)
{
    fputs("OPT Warning: ",stderr);
    if (s) fputs(s,stderr);
    fputs("\n",stderr);
}
void
opt_fatal(char *s)
{
    fputs("OPT Fatal error: ",stderr);
    if (s) fputs(s,stderr);
    exit(opt_exit_number);
}

void
opt_setstring(char **s, char *t)
{
  if (*s) free(*s);
  *s = strdup(t);
}

/* A function to justify text to a given width, useful for      */
/* printing long description strings. Returns a pointer to      */
/* a malloc'd string buffer which must be free'd by the caller. */
char*
opt_justify(char* s, int width, int indent, int initial, char* prefix) 
{
    int i, j, n, p, len, approxLines, bufSize;
    size_t k;
    char* buf;

    /* If s is null or empty, then return an empty string */
    if (ISEMPTYSTRING(s)) {
      buf = (char*) malloc(1);
      buf[0] = '\0';
      return buf;
    }

    /* Allocate a nice big buffer for the justified string */
    /* (This is very generous estimate...) */
    len = strlen(s);
    approxLines = len / width + 1;
    bufSize = initial + width + 1 +         /* The first line, inc. newline */
	(1 + indent + width) * approxLines * 3; /* Justfied lines, assuming bad line breaks */

    buf = (char*) malloc(bufSize + 1);
    if (!buf)
	opt_fatal("unable to allocate memory for string processing");

    j = 0;
    
    /* Insert initial indent and prefix */
    for (n = 0; n < initial; ++n)
	buf[j++] = ' ';
    if (prefix) {
	strcpy(buf + j, prefix);
	j += strlen(prefix);
    }

    /* Scan through input, breaking up into words */
    p = 0;   /* Position with justified column */
    for (i = 0; i < len; ++i) {
	/* Skip any leading whitespace */
	k = strspn(s + i, " \t\r\n");
	i += k;
	if (i >= len)
	    break;
	/* Search for next whitespace char */
	k = strcspn(s + i, " \t\r\n");
	/* Do we need a line break? */
	if (p + k + 1 >= width && p != 0) {
	    buf[j++] = '\n';
	    for (n = 0; n < indent; ++n)
		buf[j++] = ' ';
	    if (prefix) {
		strcpy(buf + j, prefix);
		j += strlen(prefix);
	    }
	    p = 0;
	}
	/* Add the word to the output buffer */
	/* If not at the beginning of the line, add a space */
	if (p != 0) {
	    buf[j++] = ' ';
	    ++p;
	}
	/* Copy over the word */
	strncpy(buf + j, s + i, k);
	p += k;
	i += k;
	j += k;
    }
    /* Terminate output buffer */
    buf[j] = '\0';

    /* Return the justified buffer - remember to free this! */
    return buf;
}


#if HAVE_LONGJMP
#include <setjmp.h>
static jmp_buf	opt_jumpstart;
static int	opt_jump_set=FALSE;

/* prototypes */
static RETSIGTYPE set_signal(void);
static RETSIGTYPE unset_signal(void);

/**********
 * SIGNAL *
 **********/
#if HAVE_SIGNAL
#include <signal.h>

static void
catch_signal(int sig
#if 0                           /* someday make this portable! */
             ,int code, struct sigcontext *scp, char *addr
#endif                         
             )
{
    opt_message("\nOPT Interrupted:\n");
    longjmp(opt_jumpstart,1);
    perror("Did not long-jump");
    exit(0);
}

static RETSIGTYPE (*old_catcher)(int) = NULL;
static RETSIGTYPE
set_signal(void)
{
    old_catcher = signal(SIGINT,catch_signal);
}
static RETSIGTYPE
unset_signal(void)
{
    signal(SIGINT,old_catcher);
}
#else

static RETSIGTYPE   set_signal(void) {}
static RETSIGTYPE unset_signal(void) {}

#endif /* end SIGNAL */


/* optBeginRun, optAbortRun: these functions are for running
 * the run hook defined by the programmer in optRun().
 */

int
optBeginRun(OPT_PFI run)
{
    int value;

    opt_jump_set=OPT_TRUE;
    set_signal();

    value = setjmp(opt_jumpstart);
    if(value!=0)
        opt_message("OPT: Run aborted...try again\n");
    else
        value = (*run)();
	
    unset_signal();
    opt_jump_set=OPT_FALSE;

    return value;
}

void
optAbortRun(void)
{
    if(opt_jump_set)
        longjmp(opt_jumpstart,1);
    else
        /** should this be some opt_exit() ?? **/
        opt_quit();
}

int optBeginMain(OPT_PFI_ARG run_args,int argc, char **argv)
{
    int value;

    opt_jump_set=OPT_TRUE;
    set_signal();

    value = setjmp(opt_jumpstart);
    if(value!=0)
        opt_message("OPT: Run aborted...try again\n");
    else
        value = (*run_args)(argc,argv);
	
    unset_signal();
    opt_jump_set=OPT_FALSE;

    return value;
}    

#else  /* If LONGJMP is unavailable, or disabled, then these functions
        * are rather trivial.
        */
int  optBeginRun(OPT_PFI run) { return (*run)(); }
void optAbortRun() { exit(opt_exit_number); }
int  optBeginMain(OPT_PFI_ARG run_args,int argc, char **argv)
{
    return (*run_args)(argc,argv);
}

#endif /* end LONGJMP */

void optAbortMain() { optAbortRun(); }


/************
 * Readline *
 ************/

#if HAVE_READLINE
#ifdef __cplusplus
extern "C" {
#include <readline/readline.h>
#include <readline/history.h>
}
#else
#include <readline/readline.h>
#include <readline/history.h>
#endif
#endif

void
opt_readline_init(char *s)
{
#if HAVE_READLINE
    rl_readline_name = s;
#endif
    return;
}

void 
opt_menu_getline(char *prompt, char *respond, int maxlen)
{
#if HAVE_READLINE
    char *line;
    line = readline(prompt);
    if (line && *line) {
	add_history(line);
	strncpy(respond,line,maxlen);
	free(line);
    } else
	respond[0] = '\0';
	/* double terminate the line */
    respond[strlen(respond)+1] = '\0';
#else
    opt_message(prompt);
    opt_getline(respond,stdin);
#endif
}


#if !HAVE_STRDUP
/* Provided here in case it is not available on system */     
char *strdup(char *s)
{
    char *sdup;
    int i,n;
    n = strlen(s);
    sdup = (char *)malloc(n+1);
    for (i=0; i<n; ++i)
        sdup[i] = s[i];
    sdup[n]='\0';
    return sdup;
}
#endif
    






/* OPT v3.19: options parsing tool */
/*
 * This SOFTWARE has been authored by an employee of the University of
 * California, operator of the Los Alamos National Laboratory under
 * Contract No. W-7405-ENG-36 with the U.S. Department of Energy. The
 * U.S. Government has rights to use, reproduce, and distribute this
 * SOFTWARE. Neither the Government nor the University makes any
 * warranty, express or implied, or assumes any liability or
 * responsibility for the use of this SOFTWARE. If SOFTWARE is modified
 * to produce derivative works, such modified SOFTWARE should be clearly
 * marked, so as not to confuse it with the version available from LANL.
 * 
 * Additionally, this program is free software; you can distribute it
 * and/or modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of the
 * License, or any later version.  Accordingly, this program is
 * distributed in the hope that it will be useful, but WITHOUT A
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 * for more details (http://www.gnu.org/copyleft/gpl.txt).
 * 
 * 
 * The OPT package is available from
 * http://nis-www.lanl.gov/~jt/Software
 */
/* $Id: opt.c,v 1.1.1.1 2005/10/29 19:20:46 ckemp Exp $ */
/* opt_reg.c */

/************************ register options ********************************
 *
 *    options are stored in an array,
 *    each element of which is a structure
 *
 *    The members of the structure contain information about
 *        1) the variable which is altered
 *        2) the type of variable
 *        3) a character by which the variable is invoked
 *        4) a string by which the variable can be invoked
 *        5) the positional value of the option if applicable
 *        6) a brief description of the variable's role in the program
 *        7) a longer help message
 *        8) a flag indicating that an option has been invoked
 *        9) a "hook", a function that is called when variable is invoked
 *        10) a pointer to an array structure if the option is for an array
 *            
 *    The structure elements are not assigned directly, but are
 *    "registered" with a call to the function optreg().
 *        
 *    In this file are the routines which access the structure of
 *    registered options.
 *
 *          I) Register options 
 *         II) Get number corresponding to name of option
 *        III) Print options out
 *             a) to a file which can be used as an input command line
 *             b) to the screen as part of a usage message
 *             c) to the screen as part of an interactive menu
 *         IV) Process single delimited and positional options
 */

#ifndef	SINGLEFILE
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "opt.h"
#include "ag.h"
#include "opt_p.h"
#endif

#define OPT_DEPRECATED_OK 1

/**
 * If you would like to remove these warnings then re-compile opt with
 * gcc -DOPT_DEPRECATED_OK
 */

static void
opt_CSTRING_deprecated()
{
#ifndef OPT_DEPRECATED_OK
  opt_warning("OPT_CSTRING is deprecated...you should use OPT_STRING instead");
  opt_warning("If this warning message is a problem, then recompile the source");
  opt_warning("with the symbol OPT_DEPRECATED_OK defined.");
#endif
}

static void
opt_UNDELIM_deprecated()
{
#ifndef OPT_DEPRECATED_OK
  opt_warning("OPT_UNDELIM and OPT_UNDELIMC are deprecated...");
  opt_warning("You should use positional options instead.");
  opt_warning("If this warning message is a problem, then recompile the source");
  opt_warning("with the symbol OPT_DEPRECATED_OK defined.");
#endif
}



/* Default is to indicate whether a flag is true or false by using + or -
 * But, by #define'ing FLAG_ZERO_ONE, we will use 1 or 0 instead
 */

#ifdef FLAG_ZERO_ONE
#define FLAG_TRUE '1'
#define FLAG_FALSE '0'        
#else
#define FLAG_TRUE '+'
#define FLAG_FALSE '-'        
#endif

int opt_exit_number=0;          /* if we have to exit(), then exit
                                 * with this number.
                                 */

        /* -------------------------- */
        /* List of Registered Options */
        /* -------------------------- */


typedef struct {
  void **value;      /* pointer to array of data */
  int *size;         /* pointer to size of data (number of elements in array) */
  char delim;        /* char delimits elements of array: eg "," in "1,2,3" */
  opt_TYPE base_type;  /* what the array of data is... */
}   OptArray_t;


typedef struct {
  void      *value;         /* pointer to value of option */
  opt_TYPE   type;          /* will be cast according to type */
  char       name;          /* name by which option is invoked */
  char      *longname;      /* long version of the name */
  opt_MODE   mode;          /* delimited, positional or flexible? */
  char      *descript;      /* a brief description */
  char      *help;          /* longer help message */
  int        invoked;       /* counts number of times option is invoked */
  OPT_HOOK   hook;          /* code to evaluate if option is invoked */
  OptArray_t  *array;       /* set to null if not an array option */
} Option_t;

static	Option_t *optlist=NULL;             /* array of options */
int              opt_nreg=0;                /* number of registered opts */

/* Prototypes for static (local) functions */
static void opt_fileheader(FILE *);
static int opt_action(int, ARGVECTOR *);
static int opt_number(void *);
static int opt_longname_number(const char *);
static char *optarray_getstrval(int, void *, opt_TYPE, char);
static char *opt_getstrval( void *, opt_TYPE );
static int optsizeof(opt_TYPE);
static void opt_setvalue(void *, opt_TYPE, char *);
static void optarray_action(OptArray_t *, char *);

static char opt_array_delim=',';    /* the delimiter for arrays. (eads)*/

    /* ---------- */
    /* OPT macros */
    /* ---------- */

#define OPT_isvalidnum(n) ((n)>=0 && (n)<opt_nreg)

#define	OPT_isflagtype(o) \
        ((o)==OPT_TOGGLE || (o)==OPT_BOOL || \
         (o)==OPT_NEGTOGGLE || (o)==OPT_NEGBOOL)


/* OPT_GETVALUE: returns value of what 'void *' points to */
/* OPT_SETVALUE: sets the value that the 'void *' points to */

#define OPT_GETVALUE(typ,v)      ((typ)(*((typ *)(v))))
#define OPT_SETVALUE(typ,v,val)  do { typ *xptr;       \
                                      xptr = (typ *)v; \
                                      *xptr = (typ)val; } while (0)

/* OPT_GET_ITH_VALUE: gets value of i'th option */
/* OPT_SET_ITH_VALUE: sets value of i'th option */

#define	OPT_GET_ITH_VALUE(typ,i)     OPT_GETVALUE(typ,optlist[i].value)
#define OPT_SET_ITH_VALUE(typ,i,val) OPT_SETVALUE(typ,optlist[i].value,val)


void
opt_free()
{
  int i;
  for (i=0; i<opt_nreg; ++i) {
    OPT_FREE(optlist[i].longname);
    OPT_FREE(optlist[i].help);
    OPT_FREE(optlist[i].descript);
    if (optlist[i].array && optlist[i].type == OPT_STRING) {
      char *s = *((char **)(optlist[i].value));
      OPT_FREE(s);
      OPT_FREE(optlist[i].value);
    }
    OPT_FREE(optlist[i].array);
  }
  OPT_FREE(optlist);
  opt_nreg=0;
  opt_freestrings();
}

static char *
optarray_getstrval( int size, void *val, opt_TYPE otype, char delim )
{
  int n;
  int len=0;
  char *str_val = NULL; /* string with all values in array 
                          * will be re-allocated for each item
                          */
  /* Convert the array 'val' into a string 'str_val' */
  for ( n = 0; n < size; n++ ) {
    /* this next line may look a little arcane...
     * 1. Recall val is in real life a pointer to an array, that is type **
     *    but we don't know the type, so we pass it through the argument list
     *    as a void *
     * 2. *(void **)val recasts val to a void ** and then dereferences it
     *    so that it is the pointer that is the array, but its type is
     *    `void *' instead of `type *'
     * 3. (char *)(*(void **)val) recasts this generic `void *' to be
     *    a `char *' because we know what the size of a char is (it's 1)
     *    and we can do arithmetic with it
     * 4. n*optsizeof(otype) is the offset of the n'th element of the array
     * 5. nth_val points to the n'th value in the array and if we were to
     *    recast it as *(type *)nth_val it would be the explicit value.
     *    Of course we don't know type, but we do know otype, and so
     *    we let opt_getstrval() take care of that
     */
    char *nth_val = (char *)(*(void **)val) + n*optsizeof( otype );
    char *str_nth_val = opt_getstrval( nth_val, otype );
    char str_delim[2] = { delim, '\0' }; 

    if ( n == size-1 ) str_delim[0]='\0';  /* no delimiter at the end */
      
    /* Make room in str_val for the next item */
    len += strlen(str_nth_val) + strlen(str_delim);
    if (str_val) {
      str_val = (char *)realloc(str_val,sizeof(char)*(len+1));
    } else {
      str_val = (char *)malloc(sizeof(char)*(len+1));
      str_val[0]='\0';
    }
    strcat(str_val,str_nth_val);
    strcat(str_val,str_delim);
  }
  return str_val;
}

        /* ------------------------------------------- */
        /* Routine to register each option as an array */
        /* Added by: Damian Eads <dre9227@cs.rit.edu>  */
        /* ------------------------------------------- */

int
optrega_array(  int *size, void *val, opt_TYPE otype, char name,
                char *longname, char *descript )
{
#if 0  
  static char **str_val_array=NULL;
  static int    str_val_array_n=0;

    ++str_val_array_n;
    str_val_array = (char **)realloc(str_val_array,
                                     str_val_array_n*sizeof(char *));
    str_val_array[str_val_array_n-1] = optarray_getstrval(*size,val,otype,opt_array_delim);
  
    /* Register the option as an OPT_STRING */
    optrega( &str_val_array[str_val_array_n-1],
             OPT_STRING, name, longname, descript );
#else
    char **s = (char **)malloc(sizeof(char *));
    *s = optarray_getstrval(*size,val,otype,opt_array_delim);
    optrega( s, OPT_STRING, name, longname, descript );
#endif     

    /* Warn if inappropriate base_type is specified */
    switch (otype) {
      /* These are all okay */
    case OPT_INT:
    case OPT_UINT:
    case OPT_SHORT:
    case OPT_USHORT:
    case OPT_LONG:
    case OPT_ULONG:
    case OPT_FLOAT:
    case OPT_DOUBLE:
      break;

      /* These _should_ be okay too */
    case OPT_STRING:
      break;

      /* These might be problematic */
    case OPT_CSTRING:    
    case OPT_UNDELIMC:
      opt_CSTRING_deprecated(); /* no break */
    case OPT_CHAR:
    case OPT_UCHAR:
    case OPT_INTLEVEL:
    case OPT_TOGGLE:
    case OPT_BOOL:
    case OPT_NEGTOGGLE:
    case OPT_NEGBOOL:
    case OPT_UNDELIM:
      opt_UNDELIM_deprecated(); /* no break */
    case OPT_NUL:
      opt_warn_1("Undefined behavior for %s array",optstrtyp(opt_nreg-1));
      break;

    default:
        opt_fatal("opt_getstrval: Undefined o-type");
        break;
    }

    /* Allocate an OptArray_t structure */
    optlist[ opt_nreg-1 ].array = (OptArray_t *)malloc(sizeof(OptArray_t));
                                   
    /* Then add array-specific information to the registered option */
    optlist[ opt_nreg-1 ].array->value     = (void **)val;
    optlist[ opt_nreg-1 ].array->size      = size;
    optlist[ opt_nreg-1 ].array->base_type = otype;
    optlist[ opt_nreg-1 ].array->delim     = opt_array_delim;
    
    return opt_nreg-1;
}

int optreg_array( int *size, void *val, opt_TYPE otype, char name,
		  char *descript ) {
    return optrega_array(size,val,otype,name,NULL,descript);
}
int optregc_array( int *size, void *val, opt_TYPE otype, char name ) {
    return optrega_array(size,val,otype,name,NULL,NULL);
}
int optregcb_array( int *size, void *val, opt_TYPE otype, char name,
		    char *descript ) {
    return optrega_array(size,val,otype,name,NULL,descript);
}
int optregs_array( int *size, void *val, opt_TYPE otype, char *longname ) {
    return optrega_array(size,val,otype,'\0',longname,longname);
}
int optregsb_array( int *size, void *val, opt_TYPE otype, char *longname,
		    char *descript) {
    return optrega_array(size,val,otype,'\0',longname,descript);
}
int optregp_array( int *size, void *val, opt_TYPE otype, char* longname, char* descript) {
    int n = optrega_array(size,val,otype,'\0',longname,descript);
    optmode_n(n, OPT_POSITIONAL);
    return n;
}
int optregf_array( int *size, void *val, opt_TYPE otype, char name, char* longname, char* descript) {
    int n = optrega_array(size,val,otype,name,longname,descript);
    optmode_n(n, OPT_FLEXIBLE);
    return n;
}


        /* ------------------------------- */
        /* Routine to register each option */
        /* ------------------------------- */

int optreg(void *val, opt_TYPE otype, char name, char *descript) {
    return optrega(val,otype,name,NULL,descript);
}
int optregc(void *val, opt_TYPE otype, char name) {
    return optrega(val,otype,name,NULL,NULL);
}
int optregcb(void *val, opt_TYPE otype, char name, char *descript) {
    return optrega(val,otype,name,NULL,descript);
}
int optregs(void *val, opt_TYPE otype, char *longname) {
    return optrega(val,otype,'\0',longname,longname);
}
int optregsb(void *val, opt_TYPE otype, char *longname, char *descript) {
    return optrega(val,otype,'\0',longname,descript);
}
int optregp(void *val, opt_TYPE otype, char* longname, char* descript) {
  int n = optrega(val, otype, '\0', longname, descript);
  optmode_n(n, OPT_POSITIONAL);
  return n;
}
int optregf(void *val, opt_TYPE otype, char name, char* longname, char* descript) {
  int n = optrega(val, otype, name, longname, descript);
  optmode_n(n, OPT_FLEXIBLE);
  return n;
}
int optexec(char *longname, OPT_HOOK fcn, char *descript)
{
    int n;
    n = optrega(NULL,OPT_NUL,'\0',longname,descript);
    opthook_n(n,fcn);
    return n;
}
     
int
optrega(void *val, opt_TYPE otype, char name, const char *longname, 
	const char *descript)
{
    opt_MODE mode = OPT_DELIMITED;

    /* write a warning if that option name has been used before */
    if( name != '\0' && opt_char_number(name) != -1 ) {
        opt_warn_1("optrega: Duplicate option name \'%c\'",name);
    }
    if (!ISEMPTYSTRING(longname) && opt_longname_number(longname) != -1) {
        opt_warn_1("optrega: Duplicate option longname \"%s\"",longname);
    }

    /* val points to the memory where the value of the option
     * is maintained.  This cannot be NULL.
     */
    if( val == NULL && otype != OPT_NUL ) {
        if (name)     { opt_warn_1("Name: -%c", name); }
        if (longname) { opt_warn_1("Name: --%s",longname); }
        opt_fatal("optrega: invalid first argument=NULL\n");
    }
    
    /* This is probably a mistake to print out a warning message,
       since it will be printed every time the program is run.
       Perhaps it should only be run if the user types 'run --testopt'
       */
#if 0    
    if( val != NULL  && opt_val_number(val) != -1) {
        opt_warn_1("optrega: This variable [%s] has already been registered",
                   thename);
    }
#endif    
    
    /* Programmer may optionally signal that an option is to be a
     * positional option by setting the name to '\0' and the longname
     * to NULL or an empty string.
     */
    if(name=='\0' && ISEMPTYSTRING(longname)) {
      mode = OPT_POSITIONAL;
    }

    /* The old opt concept of undelimited strings has been */
    /* replaced with that of positional options. The old   */
    /* OPT_UNDELIM type for instance, corresponds to a new */
    /* option of type OPT_STRING and mode OPT_POSITIONAL.  */
    if (otype == OPT_UNDELIM) {
      opt_UNDELIM_deprecated();
      otype = OPT_STRING;
      mode = OPT_POSITIONAL;
    }
    else if (otype == OPT_UNDELIMC) {
      opt_UNDELIM_deprecated();
      otype = OPT_CSTRING;
      mode = OPT_POSITIONAL;
    }

    /* Make sure that the longname doesn't have any whitespace */
    if (!ISEMPTYSTRING(longname)) {
     	const char *t = longname;
        while (*t != '\0') {
          if (isspace(*t)) {
            opt_fatal_1("invalid longname [%s] contains whitespace\n",longname);
          }
          ++t;
        }
    }

    /* reallocate another struct for the new opt */
    ++opt_nreg;
    optlist = (Option_t *)realloc(optlist,opt_nreg*sizeof(Option_t));

    /* Having checked for various warnings, now register the options */
    
    optlist[opt_nreg-1].value         = val;
    optlist[opt_nreg-1].type          = otype;
    optlist[opt_nreg-1].name          = name;
    optlist[opt_nreg-1].longname      = longname ? strdup(longname) : NULL;
    optlist[opt_nreg-1].mode          = mode;
    optlist[opt_nreg-1].descript      = descript ? strdup(descript) : NULL;

    optlist[opt_nreg-1].invoked       = 0;
    optlist[opt_nreg-1].help          = NULL;
    optlist[opt_nreg-1].hook          = NULL;
    optlist[opt_nreg-1].array         = NULL;

    return opt_nreg-1;

}

void optchar_n(int n, char name)
{
    if (OPT_isvalidnum(n)) optlist[n].name = name;
}    
void optmode_n(int n, opt_MODE mode)
{
  if (OPT_isvalidnum(n)) {
    if (!(mode == OPT_POSITIONAL || mode == OPT_DELIMITED || mode == OPT_FLEXIBLE))
      opt_warning("optmode() was passed an invalid mode");
    optlist[n].mode = mode;
  }
}
void optdescript_n(int n, char *descript)
{
    if (OPT_isvalidnum(n)) opt_setstring(&optlist[n].descript,descript);
}
void opthelp_n(int n, char *help)
{
    if (OPT_isvalidnum(n)) opt_setstring(&optlist[n].help,help);
}
void optlongname_n(int n, char *longname)
{
    if (OPT_isvalidnum(n)) opt_setstring(&optlist[n].longname,longname);
}
void opthook_n(int n, OPT_HOOK hook)
{
    if (OPT_isvalidnum(n)) optlist[n].hook = hook;
}
int optinvoked_n(int n)
{
    if (OPT_isvalidnum(n))
        return optlist[n].invoked;
    return 0;
}

int optinvoked(void *v) { return optinvoked_n( opt_number(v) ); }

void optchar(void *v,  char name) { optchar_n( opt_number(v),   name ); }
void opthelp(void *v,  char *help) { opthelp_n( opt_number(v),  help ); }
void opthook(void *v,  OPT_HOOK hook) { opthook_n( opt_number(v), hook ); }      
void optlongname(void *v, char *longname)
{
    optlongname_n( opt_number(v), longname );
}
void optmode(void* v, opt_MODE mode) 
{
  optmode_n( opt_number(v), mode );
}
void optdescript(void *v, char *descript)
{
    optdescript_n( opt_number(v), descript );
}
void optarraydelim_n(int n, char delim)
{
    if (OPT_isvalidnum(n) && optlist[n].array) {
      if (optlist[n].array->delim != delim) {
        optlist[n].array->delim = delim;
        if (*(optlist[n].array->size) > 0) {
          OPT_SET_ITH_VALUE(char *,n,\
                            optarray_getstrval(*(optlist[n].array->size),
                                               optlist[n].array->value,
                                               optlist[n].array->base_type,
                                               optlist[n].array->delim));
        }
        
      }
    }
}  
void optarraydelim(void *v, char delim)
{
    int n;
    if (v == NULL) {
      /* then reset default array delim and all opt delims*/
      opt_array_delim = delim;
      for (n=0; n<opt_nreg; ++n) {
        optarraydelim_n(n,delim);
      }
    } else {
      /* just reset for the specific option */
      n = opt_number(v);
      optarraydelim_n(n,delim);
    }
}
    
void optExitNumber(int n) {
    opt_exit_number = n;
}


void
opt_get_help(char c)
{
    int n;
    n = opt_char_number(c);
    if (OPT_isvalidnum(n) && optlist[n].help != NULL)
        opt_mess_2("%c: %s\n",c,optlist[n].help);
    else
        opt_mess_1("Help unavailable for \'%c\'\n",c);
}
	
	

/* opt_XXX_number(): Get number corresponding to option name; return
 * a value of -1 if the option does not exist.
 *  opt_char_number(c) takes as input the single-character name;
 *  opt_longname_number(longname) takes as input the long name
 *  opt_number(void *) takes as input pointer to variable
 */

int
opt_char_number(char c)
{
    int i;                         /* see which registered option */
    for(i=0; i<opt_nreg; ++i) {
        if( c == optlist[i].name )
            return i;
    }
    return -1;
}
    
static int
opt_number(void *v)
{
    int i;                         /* see which registered option */
    for(i=0; i<opt_nreg; ++i) {
      if (optlist[i].array && v == optlist[i].array->value) return i;
      if (!optlist[i].array && v == optlist[i].value ) return i;
    }
    return -1;                     /* to signify not an option */
}
static int
opt_longname_number(const char *s)
{
    int i;                         /* see which registered option */
    if ( ISEMPTYSTRING(s) ) return -1;
    for(i=0; i<opt_nreg; ++i) {
        if ( !ISEMPTYSTRING(optlist[i].longname) ) {
            if (0==strcmp(s,optlist[i].longname)) {
                return i;
            }
        }
    }
    return -1;                     /* to signify not an option */
}

	/* -------------------------------------------- */
	/* Print value of registered option to a string */
	/* -------------------------------------------- */

/* optstrval:
 * returns a string corresponding to the ith option
 * Thus if x = 12.6, this returns the string "12.6"
 * Calls opt_getstrval for all types except INTLEVEL.
 * The return value is a pointer to a static string which is
 * overwritten with each call
 * FLAG values are returned as "+" or "-"
 * Previously: INTLEVEL values are returned as "- -xxxx..." depending on value
 * Now, more straightforward, treat INTLEVEL just like INT
 */

char *
optstrval(int i)
{
  return opt_getstrval(optlist[i].value,optlist[i].type);
}/*optstrval*/

	/* -------------------------------------------- */
	/* Print value of registered option to a string */
        /* without specifying an index. (eads)          */
	/* -------------------------------------------- */

static char *
opt_getstrval( void *val, opt_TYPE typ ) {
    int maybe;
    static char stval_buf[80];

    switch( typ ) {
    case OPT_INT:
        opt_snprintf_1(stval_buf,80,"%d", OPT_GETVALUE(int,val) );
        break;
    case OPT_UINT:
        opt_snprintf_1(stval_buf,80,"%u", OPT_GETVALUE(unsigned int,val) );
        break;
    case OPT_SHORT:
        opt_snprintf_1(stval_buf,80,"%d", OPT_GETVALUE(short,val) );
        break;
    case OPT_USHORT:
        opt_snprintf_1(stval_buf,80,"%u", OPT_GETVALUE(unsigned short,val) );
        break;
    case OPT_LONG:
        opt_snprintf_1(stval_buf,80,"%ld", OPT_GETVALUE(long,val) );
        break;
    case OPT_ULONG:
        opt_snprintf_1(stval_buf,80,"%lu", OPT_GETVALUE(unsigned long,val) );
        break;
    case OPT_CHAR:
        opt_snprintf_1(stval_buf,80,"%c", OPT_GETVALUE(char,val));
        break;
    case OPT_UCHAR:
        opt_snprintf_1(stval_buf,80,"%c", OPT_GETVALUE(unsigned char,val));
        break;
    case OPT_FLOAT:
        opt_snprintf_1(stval_buf,80,"%g", OPT_GETVALUE(float,val) );
        break;
    case OPT_DOUBLE:
        opt_snprintf_1(stval_buf,80,"%g", OPT_GETVALUE(double,val) );
        break;
    case OPT_INTLEVEL:
        if (1) {
          int v = OPT_GETVALUE(int,val);
          if (v < 0) {
            opt_fatal_1("Invalid (negative) value [v=%d] for OPT_INTLEVEL",v);
          }
          if (v < 10) {
            opt_snprintf_1(stval_buf,80,"%1d",v);
          }
          if (v >= 10) {
            /* Could set it to "- -xxxxx" but need to know x */
            /* and it only works if 'x' is a one-character name */
            opt_warn_1("Cannot set string for OPT_INTLEVEL v=%d > 9",v);
          }
        }
        break;
    case OPT_TOGGLE:
    case OPT_BOOL:
        maybe = OPT_GETVALUE(int,val);
        opt_snprintf_1(stval_buf,80,"%c", (maybe==OPT_TRUE ? FLAG_TRUE : FLAG_FALSE) );
        break;
    case OPT_NEGTOGGLE:
    case OPT_NEGBOOL:
        maybe = OPT_GETVALUE(int,val);
        opt_snprintf_1(stval_buf,80,"%c", (maybe==OPT_TRUE ? FLAG_FALSE : FLAG_TRUE) );
        break;
    case OPT_CSTRING:    
        opt_CSTRING_deprecated();
        opt_snprintf_1(stval_buf,80,"\"%s\"",(char *)val);
        break;
    case OPT_STRING:
        if( val == NULL )	/* this should never happen */
            sprintf(stval_buf,"\"\"");
        else 
            if( OPT_GETVALUE(char *,val) == NULL )
                sprintf(stval_buf,"\"\"");
            else
                opt_snprintf_1(stval_buf,80,"\"%s\"",OPT_GETVALUE(char *,val));
        break;
    case OPT_NUL:
        stval_buf[0]='\0';
        break;
    default:
        opt_fatal("opt_getstrval: Undefined o-type");
        break;
    }
    /* ---- Return value is static buffer ---- */
    return stval_buf;
}/* opt_getstrval */

static int
optsizeof( opt_TYPE typ )
{
    int retval=0;
    
    switch( typ ) {
    case OPT_INT:
        retval = sizeof(int);
        break;
    case OPT_UINT:
        retval = sizeof(unsigned int);
        break;
    case OPT_SHORT:
        retval = sizeof(short);
        break;
    case OPT_USHORT:
        retval = sizeof(unsigned short);
        break;
    case OPT_LONG:
        retval = sizeof(long);
        break;
    case OPT_ULONG:
        retval = sizeof(unsigned long);
        break;
    case OPT_CHAR:
        retval = sizeof(char);
        break;
    case OPT_UCHAR:
        retval = sizeof(unsigned char);
        break;
    case OPT_FLOAT:
        retval = sizeof(float);
        break;
    case OPT_DOUBLE:
        retval = sizeof(double);
        break;
    case OPT_INTLEVEL:
        retval = sizeof(int);
        break;
    case OPT_TOGGLE:
    case OPT_BOOL:
        retval = sizeof(int);
        break;
    case OPT_NEGTOGGLE:
    case OPT_NEGBOOL:
        retval = sizeof(int);
        break;
    case OPT_CSTRING:
        opt_CSTRING_deprecated();
        retval = sizeof(char*);
        break;
    case OPT_STRING:
        retval = sizeof(char*);
        break;
    case OPT_NUL:
      /* FIXME */
        retval = sizeof(int(*)());
        break;
    default:
        opt_fatal("optsizeof: Undefined o-type");
        break;
    }
    return retval;
}/*optsizeof*/

/* optstrtyp:
 * returns the type of the ith registered option as a string
 * Thus if x = 12.6 is a float , this returns the string "<FLOAT>"
 * The return value is a pointer to a static string which is
 * overwritten with each call
 */
char *
optstrtyp(int i)
{
    static char	sttyp_buf[80];
    opt_TYPE o;

    o = (optlist[i].array) ? optlist[i].array->base_type : optlist[i].type;
    
    switch( o ) {
        /* Integers */
    case OPT_INT:         strcpy(sttyp_buf,"<INT"); break;
    case OPT_UINT:        strcpy(sttyp_buf,"<UINT"); break;
    case OPT_SHORT:       strcpy(sttyp_buf,"<SHORT"); break;
    case OPT_LONG:        strcpy(sttyp_buf,"<LONG"); break;
    case OPT_USHORT:      strcpy(sttyp_buf,"<USHORT"); break;
    case OPT_ULONG:       strcpy(sttyp_buf,"<ULONG"); break;
        /* Floating point */
    case OPT_FLOAT:       strcpy(sttyp_buf,"<FLOAT"); break;
    case OPT_DOUBLE:      strcpy(sttyp_buf,"<DOUBLE"); break;
        /* Char */
    case OPT_CHAR:        strcpy(sttyp_buf,"<CHAR"); break;
    case OPT_UCHAR:       strcpy(sttyp_buf,"<UCHAR"); break;
        /* Flag */
    case OPT_INTLEVEL:    strcpy(sttyp_buf,"<INTLEVEL"); break;
    case OPT_TOGGLE:      strcpy(sttyp_buf,"<TOGGLE"); break;
    case OPT_BOOL:        strcpy(sttyp_buf,"<BOOL"); break;
    case OPT_NEGTOGGLE:   strcpy(sttyp_buf,"<NEGTOGGLE"); break;
    case OPT_NEGBOOL:     strcpy(sttyp_buf,"<NEGBOOL"); break;
        /* String */
    case OPT_STRING:      strcpy(sttyp_buf,"<STRING"); break;
    case OPT_CSTRING:     strcpy(sttyp_buf,"<STRING"); break;
        /* Nul */
    case OPT_NUL:         strcpy(sttyp_buf,"<NUL"); break;
        
    default:
        opt_fatal("usage: undefined o-type");
    }
    /* If it's an array, then add " ARRAY" to the name of the type */
    if ( optlist[i].array ) {
      strcat(sttyp_buf," ARRAY");
    }
    /* finish with closing bracket */
    strcat(sttyp_buf, ">");
    return sttyp_buf;
}

		/* ------------------------------ */
		/* Put registered options to file */
		/* ------------------------------ */

void
opt_tofile(FILE *fp)
{
    int	i;
    char *format;
    /* Short Option Names */
    char *fnumstr="-%c %-33s ";	/* numbers and strings  */
    char *fflgchr="-%c%-33s  ";	/* flags and characters */

    opt_fileheader(fp);
    /* Short Option Names First */
    for(i=0; i<opt_nreg; ++i) {
        char *s;
        if (optlist[i].name == '\0')
            continue;
        if (optlist[i].type == OPT_NUL)
            continue;
	if (optlist[i].mode == OPT_POSITIONAL)
	  continue;

        if( OPT_isflagtype( optlist[i].type )
            || optlist[i].type == OPT_INTLEVEL
            || optlist[i].type == OPT_CHAR )
            format = fflgchr;
        else
            format = fnumstr;

        fprintf(fp,format,optlist[i].name,optstrval(i));
        s = optlist[i].descript;
        if (strlen(s) > 37) {
          s = opt_justify(optlist[i].descript,55,20,20,"; ");
          fprintf(fp,"\n%s\n",s);
          free(s);
        } else {
          fprintf(fp,"; %s\n",s);
        }
          
        /* If an option has both a short and a long option name
         * then write it as a comment '    ; -b == --option-b
         */
        if (!ISEMPTYSTRING(optlist[i].longname)) {
          if (strlen(optlist[i].descript) > 37) {
            fprintf(fp,"%-20s%c -%c == --%s\n","",IGNOREEOL,
                    optlist[i].name,optlist[i].longname);
          } else {
            fprintf(fp,"%-37s%c -%c == --%s\n","",IGNOREEOL,
                    optlist[i].name,optlist[i].longname);
          }
        }
    }
    /* Now Long Option Names (thanks to Jason V. Morgan) */
    for(i=0; i<opt_nreg; ++i) {
        int len;
        char *s;
        /* do not output the long name if a short name exists */
        if (optlist[i].name != '\0')
           continue;
        if (!optlist[i].longname) 
           continue;
	if (optlist[i].mode == OPT_POSITIONAL)
	  continue;
        /* The following arguably applies to short options too ?? */
        /* don't want anything with a hook (do we?) */
        if (optlist[i].type == OPT_NUL)
            continue;
        if (optlist[i].hook != NULL)
            continue;

        len = 34 - strlen(optlist[i].longname);
        if (len<1) len=1;
        fprintf(fp, "--%s=%-*s", optlist[i].longname, len, optstrval(i));

        s = optlist[i].descript;
        if (strlen(s) > 35) {
          s = opt_justify(optlist[i].descript,55,20,20,"; ");
          fprintf(fp,"\n%s\n",s);
          free(s);
        } else {
          fprintf(fp,"; %s\n",s);
        }
    }
}
static void
opt_fileheader(FILE *fp)
{
    extern char *opt_program_name;
    int	i,n;
    char buf[80];

    opt_snprintf_1(buf,80,"Options file created by: [%.30s]",opt_program_name);
    n = strlen(buf);

    fputc(IGNOREEOL,fp);
    fprintf(fp,"%s",buf);
    fputc(IGNOREEOL,fp); fputc('\n',fp);

    fputc(IGNOREEOL,fp);
    for(i=0;i<n;++i) fprintf(fp,"-");
	fputc(IGNOREEOL,fp); fputc('\n',fp);
}
		/* ----------------------------------------- */
		/* write options out to the interactive menu */
		/* ----------------------------------------- */

/* Produces a malloc'd string that can be used to display    */
/* the essential info about an option in a line or two. Note */
/* the string must be free'd by the caller!                  */

char *
opt_mstring(int i)
{
    char* mstring;
    char* justified_descript;
    char	flgstr[10];
    char	*s;
    char	*format = "-%c %s%*s%s";
    char        *lformat = "--%s=%s%*s%s";

    /* If we have an unnamed option then don't bother returning it */
    if (optlist[i].name == '\0' && ISEMPTYSTRING(optlist[i].longname))
      return NULL;

    /* If the option is purely positional then don't return it either */
    if (optlist[i].mode == OPT_POSITIONAL)
      return NULL;

    s = optstrval(i);
    if( optlist[i].type == OPT_INTLEVEL ) {
        int f;
        f = OPT_GET_ITH_VALUE(int,i);
        if(f==0)
            strcpy(flgstr,"OFF");
        else if (f==1)
            strcpy(flgstr,"ON");
        else
            opt_snprintf_1(flgstr,10,"ON:%d",f);
        s=flgstr;
    } else if( optlist[i].type == OPT_CHAR ) {
        switch (s[0]) {
        case '\0':	sprintf(flgstr,"\'\\0\'"); break;
        case '\n':	sprintf(flgstr,"\'\\n\'"); break;
        case '\t':	sprintf(flgstr,"\'\\t\'"); break;
        case ' ':	sprintf(flgstr,"\'  \'");  break;
        default:	sprintf(flgstr,"%c",s[0]); break;
        }
        s = flgstr;
    } else if( OPT_isflagtype( optlist[i].type ) ) {
        if( *s == '+' || *s == '1')
            strcpy(flgstr,"TRUE");
        if( *s == '-' || *s == '0')
            strcpy(flgstr,"FALSE");
        s=flgstr;
    }

    /* Note that opt_justify() always returns a valid string, even */
    /* if descript is NULL or empty.                               */
    justified_descript = opt_justify(optlist[i].descript, 40, 35, 0, 0);

    mstring = (char*) malloc(strlen(justified_descript) + 200);
    if (!mstring)
	opt_fatal("unable to allocate memory for menu text strings");

    if (optlist[i].name)
        sprintf(mstring,format,optlist[i].name,s,
		35-strlen(s)-3,"",justified_descript);
    else if (optlist[i].longname)
        sprintf(mstring,lformat,optlist[i].longname,s,
		35-strlen(optlist[i].longname)-strlen(s)-3,"",justified_descript);

    free (justified_descript);
    return mstring;
}

/* opt_usage: Write a usage statement describing registered options.
 * This is what the program writes in response to '--help' on the
 * command line.
 */

#define OPTUSAGEDEBUG 0
#define OPTDESCRIPTWIDTH  43
#define OPTDESCRIPTOFFSET 38
#define OPTMAXLONGNAMELEN 14
void
opt_usage(void)
{
    static char *dformat =        " -%c                   %-14s  %s\n";
    static char *lformat =                 " -%c, --%-14s %-14s  %s\n";
    static char *pformat =                        " %-20s %-14s  %s\n";
    static char *uformat =  " %-10s                              %s\n";
    static char *aformat =                  "     --%-14s %-14s  %s\n";
    int  i;
    int initial;
    char* justified_descript;

    /* First we print the positional options */
    initial = 1;
    for (i = 0; i < opt_nreg; ++i) {
      if (optlist[i].mode != OPT_DELIMITED) {
	if (initial) {
	  opt_message("The arguments are:\n");
	  initial = 0;
	}
	/* Get a suitably formatted version of the description  */
	/* descript may be NULL, but opt_justify always returns */
	/* a valid (perhaps empty) string.                      */
	justified_descript = opt_justify(optlist[i].descript,
                                         OPTDESCRIPTWIDTH,
                                         OPTDESCRIPTOFFSET,0,0);
        if(OPTUSAGEDEBUG) fprintf(stderr,"%2d u ",i);
	if (!ISEMPTYSTRING(optlist[i].descript))
	  opt_mess_3(pformat,optlist[i].longname,optstrtyp(i),justified_descript);
	else
	  opt_mess_2(uformat,optstrtyp(i),justified_descript);
	free(justified_descript);
      }
    }

    /* Second we print the delimited options */
    initial = 1;
    for(i=0; i<opt_nreg; ++i) {
      if (optlist[i].mode != OPT_POSITIONAL && 
	  (optlist[i].name != '\0' || !ISEMPTYSTRING(optlist[i].longname))) {
	if (initial) {
	  opt_message("The options are:\n");
	  initial = 0;
	}
	/* Get a suitably formatted version of the description */
	/* If the option is positional, then refer to the argument description */
	if (optlist[i].mode == OPT_FLEXIBLE) {
	  justified_descript = opt_justify("(See argument description above)",
                                           OPTDESCRIPTWIDTH,
                                           OPTDESCRIPTOFFSET,0,0);
        } else if (optlist[i].longname &&
                   strlen(optlist[i].longname) > OPTMAXLONGNAMELEN) {
	  justified_descript = opt_justify(optlist[i].descript,
                                           OPTDESCRIPTWIDTH,
                                           OPTDESCRIPTOFFSET,
                                           OPTDESCRIPTOFFSET,0);
        } else {
	  justified_descript = opt_justify(optlist[i].descript,
                                           OPTDESCRIPTWIDTH,
                                           OPTDESCRIPTOFFSET,0,0);
        }
        
	if (optlist[i].name == '\0') {
	  if(OPTUSAGEDEBUG) fprintf(stderr,"a ");
          if (strlen(optlist[i].longname) > OPTMAXLONGNAMELEN) {
            opt_mess_3(aformat,optlist[i].longname,
                       optstrtyp(i),"");
            opt_mess_1("%s\n",justified_descript);
          } else {
            opt_mess_3(aformat,optlist[i].longname,
                       optstrtyp(i),justified_descript);
          }            
	} else {
	  if (!ISEMPTYSTRING(optlist[i].longname)) {
	    if(OPTUSAGEDEBUG) fprintf(stderr,"l ");
            if (strlen(optlist[i].longname) > OPTMAXLONGNAMELEN) {
              opt_mess_4(lformat,optlist[i].name,optlist[i].longname,
                         optstrtyp(i),"");
              opt_mess_1("%s\n",justified_descript);
            } else {
              opt_mess_4(lformat,optlist[i].name,optlist[i].longname,
                         optstrtyp(i),justified_descript);
            }            
	  } else {
	    if(OPTUSAGEDEBUG) fprintf(stderr,"d ");
	    opt_mess_3(dformat,optlist[i].name,
		       optstrtyp(i),justified_descript);
	  }
	}
	free(justified_descript);
      }
    }

    if (opt_additional_usage_fcn != NULL) {
        /* precede additional usage by -- to indicate that we
         * are finished with the options
         */
        opt_message(" --\n");
        /* return value is ignored */
        (*opt_additional_usage_fcn)();
    }
}

	/* ----------------------- */
	/* Get and process options */
	/* ----------------------- */

	/*
	 * The routines below make use of the ARGVECTOR structure
	 * defined in "ag.h" to process the command line options
	 */
		
	/* -------------------------------------- */
	/* get and process a positional option */
	/* -------------------------------------- */

int
opt_parse_positional(ARGVECTOR *ag)
{
    int     i;
 
    /* Find the first positional option that has not been invoked and
     * set that to the argument string.  Note that once a positional
     * option has been invoked, it cannot be un-invoked.
     */

    for (i=0; i<opt_nreg; ++i) {
        if (optlist[i].invoked || optlist[i].mode == OPT_DELIMITED) 
	  continue;
	return opt_action(i,ag);
    }
    /* If we get here, we didn't find an uninvoked positional option */
    return 0;
}

		/* ---------------------------------- */
		/* get and process a delimited option */
		/* ---------------------------------- */

int
opt_parse_delim(ARGVECTOR *ag)
{
    int      i;
    char     c;
 
    c = ag_c_advance(ag);   /* first character gives name of option */
    i = opt_char_number(c);      /* this is number of option w/ name 'c' */
 
    if(! OPT_isvalidnum(i)) { /* if invalid option name */
        opt_warn_1("%c not a registered option name\n",c);
        return 0;
    }
    return opt_action(i,ag);
}


int
opt_parse_longdelim(ARGVECTOR *ag)
{
    int      i;
    char     c;
    char     *sc,*scptr;

    /* this code parses input options of the form --var=value */
    
    sc = strdup(ag_s(ag));      /* long option name is a string */
    scptr = sc;
    while( !ag_eow(ag) ) {
        c=ag_c_advance(ag);
        if (c == '=') {
            *scptr = '\0';
            break;
        }
        ++scptr;
    }
    
    i  = opt_longname_number(sc);   /* this is number of option w/ name `sc' */
    if(! OPT_isvalidnum(i)) {
        opt_warn_1("%s not a registered option name\n",sc);
        return 0;
    }
    free((char *)sc); 
    return opt_action(i,ag);
}

static void
opt_setvalue(void *v, opt_TYPE o, char *s)
{
  switch(o) {
  case OPT_INT:
    OPT_SETVALUE(int,v, opt_atoi(s));
    break;
  case OPT_UINT:
    OPT_SETVALUE(unsigned int, v, opt_atou(s));
    break;
  case OPT_SHORT:
    OPT_SETVALUE(short, v, opt_atoi(s));
    break;
  case OPT_LONG:
    OPT_SETVALUE(long, v, opt_atoi(s));
    break;
  case OPT_USHORT:
    OPT_SETVALUE(unsigned short, v, opt_atou(s));
    break;
  case OPT_ULONG:
    OPT_SETVALUE(unsigned long, v, opt_atou(s));
    break;
  case OPT_FLOAT:
    OPT_SETVALUE(float, v, opt_atof(s));
    break;
  case OPT_DOUBLE:
    OPT_SETVALUE(double, v, opt_atof(s));
    break;
  case OPT_CHAR:
    OPT_SETVALUE(char, v, s[0]);
    break;
  case OPT_UCHAR:
    OPT_SETVALUE(unsigned char, v, s[0]);
    break;
        
  case OPT_STRING: /* warning: memory leak! */
    OPT_SETVALUE(char *, v,strdup(s));
    break;
  case OPT_CSTRING:
    opt_CSTRING_deprecated();
    strncpy( (char *)v , s, OPT_MAXSTRLEN );
    break;
        
  default:
    opt_fatal("Invalid type");
  }
}
    
static void
optarray_action(OptArray_t *arr, char *stmp)
{
  char *s = strdup(stmp);
  char *t = s;
  char *u = s;
  char cu = '\0';
  char *v;  /* why is this a char * and not a void *, even though it
             * is meant as a generic pointer?  see the long comment
             * above in the optarray_getstrval() function
             */

  /* Reset the array values */
  if (*(arr->value)) {
#if 0 /* Trying to free some memory, but I get "FREE DANGLING" */   
    if (arr->base_type == OPT_STRING) {
      int i;
      for (i=0; i < *(arr->size); ++i) {
        char *ith_val = (char *)(*(arr->value)) + i;
        free( ith_val );
      }
    }
#endif    
    free(*(arr->value));
  }
  *(arr->value)=NULL;
  *(arr->size)=0;

  while (*t) {
    for (u=t; !(*u == '\0' || *u == arr->delim); ++u) {
      /* read until the first delimiter */
      /* t is first character of word */
      /* u is delimiter, so u-1 is last character of t's word */
    }
    cu = *u;
    *u = '\0';  /* now t is null-terminated string */

    /* add new element to array */
    *(arr->size) += 1;
    *(arr->value) = realloc(*(arr->value),
                            *(arr->size)*optsizeof(arr->base_type));
    /* v points to the next element in the array */
    v = (char *)(*(arr->value)) + (*(arr->size)-1)*optsizeof(arr->base_type);
    /* set the value of v according to the string t */
    opt_setvalue(v,arr->base_type,t);
    /* now go to the next word */
    if (cu) t = u+1;
    else    t = u; /* unless that was the last word */
  }
  free(s);
}
                
static int
opt_action(int i, ARGVECTOR *ag)
{
    /* ------------- */
    /* act on option */
    /* ------------- */
    opt_TYPE  o;
    char      c;
    char      *s;
    int       yes,no,maybe,toggle;    /* flag values */

    if (!OPT_isvalidnum(i)) {
        opt_warn_1("opt_action: invalid number %d\n",i);
        return 0; /* failure */
    }

    switch( o=optlist[i].type ) {
    case OPT_INT:
    case OPT_UINT:
    case OPT_SHORT:
    case OPT_LONG:
    case OPT_USHORT:
    case OPT_ULONG:
    case OPT_FLOAT:
    case OPT_DOUBLE:
      opt_setvalue(optlist[i].value,o,argnext(ag));
      break;
    case OPT_CHAR:
      OPT_SET_ITH_VALUE(char,i,ag_c_advance(ag));
      break;
    case OPT_UCHAR:
      OPT_SET_ITH_VALUE(unsigned char,i,ag_c_advance(ag));
      break;

    case OPT_TOGGLE:
    case OPT_BOOL:
    case OPT_NEGTOGGLE:
    case OPT_NEGBOOL:
        /* define terms:
         * yes: value of flag when  ON
         * no:  value of flag when OFF
         * maybe:  present value of flag
         * toggle: negative of maybe
         */
      maybe  = OPT_GET_ITH_VALUE(int,i);
      toggle = (maybe==OPT_TRUE ? OPT_FALSE : OPT_TRUE);
      if( o == OPT_TOGGLE || o==OPT_BOOL )
        yes=OPT_TRUE;
      else
        yes=OPT_FALSE;
      no = (yes==OPT_TRUE ? OPT_FALSE : OPT_TRUE);
      c = ag_c(ag);           /* character following 'c' */
      switch(c) {
        case '+':
        case '1':
            OPT_SET_ITH_VALUE(int,i,yes);
            ag_c_advance(ag);       /* eat the '+' */
            break;
        case '-':
        case '0':
            OPT_SET_ITH_VALUE(int,i,no);
            ag_c_advance(ag);       /* eat the '-' */
            break;
        default:
            if( o==OPT_BOOL || o==OPT_NEGBOOL ) {
                OPT_SET_ITH_VALUE(int,i,yes);
            } else {    /* otherwise toggle value */
                OPT_SET_ITH_VALUE(int,i,toggle);
            }
            break;
      }
      break;

    case OPT_INTLEVEL:
      /* default is to increment */
        OPT_SET_ITH_VALUE(int,i,OPT_GET_ITH_VALUE(int,i) + 1);  
        c = ag_c(ag);
        switch(c) {
        case '+':
            /* we've already incremented */
            ag_c_advance(ag);
            break;
        case '-':
            OPT_SET_ITH_VALUE(int,i,0);
            ag_c_advance(ag);
            break;
        case '0': case '1': case '2': case '3': case '4':
        case '5': case '6': case '7': case '8': case '9':
            OPT_SET_ITH_VALUE(int,i,c-'0');
            ag_c_advance(ag);
        default:
            break;
        }
        break;
        
    case OPT_STRING:
        s = argnext(ag);
        if (optinvoked_n(i) && optlist[i].value != NULL) {
          /* if this option has been invoked before, then we can be
           * confident that it can be safely freed, before reassigning
           * it to a new value
           */
#if 0
          fprintf(stderr,"invoked opt[%d] for %d times:s=%s\n",i,optinvoked_n(i),s);
#endif          
          free(*((char **)(optlist[i].value)));
        }
        OPT_SET_ITH_VALUE(char *,i,strdup(s));
        /* If this type is really an array, then act accordingly */
        if (optlist[i].array) {
          optarray_action(optlist[i].array,s);
        }
        break;
    case OPT_CSTRING:
        opt_CSTRING_deprecated();
        strncpy( (char *)optlist[i].value , argnext(ag), OPT_MAXSTRLEN );
        break;
        
    case OPT_NUL:
        break;
        
    default:
        opt_fatal("opt_delim: don't know this type");
        break;
    }

    /* Indicate that the option was invoked */
    optlist[i].invoked += 1;

    /* Run the hook if its defined */
    if (optlist[i].hook) {
        int retval;
        retval = optlist[i].hook(optlist[i].value);

        if (retval == OPT_ERROR) {
            /* do nothing for now */;
        }
        if (retval == OPT_ABORT && !opt_menuflag)
            exit(opt_exit_number);
        if (retval == OPT_QUIT)
            opt_quit();
        if (retval == OPT_EXIT)
            exit(opt_exit_number);
    }
    
    return 1;                  /* success */
}

static char *optquotedstrtyp(int i)
{
  char *s,*ss;
  ss = s = optstrtyp(i);
  while (*s) {
    if (*s == '<') { *s = '\"'; ++s; }
    else if (*s == '>') { *s = '\"'; ++s; *s='\0'; }
    else ++s;
  }
  return ss;
}
  
static void opt_xml(int i, FILE *fp)
{
  fprintf(fp,"  <option type=%s",optquotedstrtyp(i));
  if (optlist[i].mode == OPT_POSITIONAL) {
    fprintf(fp," mode=\"POSITIONAL\"");
  } else if (optlist[i].mode == OPT_FLEXIBLE) {
    fprintf(fp," mode=\"FLEXIBLE\"");
  }    
  fprintf(fp,">\n");
  if (optlist[i].name)     fprintf(fp,"    <short>%c</short>\n",optlist[i].name);
  if (optlist[i].longname) fprintf(fp,"    <long>%s</long>\n",optlist[i].longname);
  if (optlist[i].descript) fprintf(fp,"    <description>%s</description>\n",optlist[i].descript);
  if (optlist[i].help)     fprintf(fp,"    <help>%s</help>\n",optlist[i].help);
  if (optlist[i].array) {
    fprintf(fp,"    <array-delimiter>%c</array-delimiter>\n",
            optlist[i].array->delim);
  }
  fprintf(fp,"  </option>\n");
}
void opt_write_xml_file(FILE *fp)
{
  int i;
  fprintf(fp,"<opt-interface program=\"%s\">\n",opt_program_name);
  for (i=0; i<opt_nreg; ++i) {
    opt_xml(i,fp);
  }
  fprintf(fp,"</opt-interface>\n");
}  














/* OPT v3.19: options parsing tool */
/*
 * This SOFTWARE has been authored by an employee of the University of
 * California, operator of the Los Alamos National Laboratory under
 * Contract No. W-7405-ENG-36 with the U.S. Department of Energy. The
 * U.S. Government has rights to use, reproduce, and distribute this
 * SOFTWARE. Neither the Government nor the University makes any
 * warranty, express or implied, or assumes any liability or
 * responsibility for the use of this SOFTWARE. If SOFTWARE is modified
 * to produce derivative works, such modified SOFTWARE should be clearly
 * marked, so as not to confuse it with the version available from LANL.
 * 
 * Additionally, this program is free software; you can distribute it
 * and/or modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of the
 * License, or any later version.  Accordingly, this program is
 * distributed in the hope that it will be useful, but WITHOUT A
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 * for more details (http://www.gnu.org/copyleft/gpl.txt).
 * 
 * 
 * The OPT package is available from
 * http://nis-www.lanl.gov/~jt/Software
 */
/* opt_regt.c */
/* Note: this file generated automatically from makeopth.pl */
/* Edits to this file may be lost! */
#include "opt.h"
int
optrega_CHAR(char *v, char c, char *n, char *B) {
      return optrega(v,OPT_CHAR,c,n,B);
}
int
optreg_CHAR(char *v, char c, char *B) {
      return optreg(v,OPT_CHAR,c,B);
}
int
optregc_CHAR(char *v, char c) {
      return optregc(v,OPT_CHAR,c);
}
int
optregs_CHAR(char *v, char *n) {
      return optregs(v,OPT_CHAR,n);
}
int
optregcb_CHAR(char *v, char c, char *B) {
      return optregcb(v,OPT_CHAR,c,B);
}
int
optregsb_CHAR(char *v, char *s, char *B) {
      return optregsb(v,OPT_CHAR,s,B);
}
int
optrega_CSTRING(char *v, char c, char *n, char *B) {
      return optrega(v,OPT_CSTRING,c,n,B);
}
int
optreg_CSTRING(char *v, char c, char *B) {
      return optreg(v,OPT_CSTRING,c,B);
}
int
optregc_CSTRING(char *v, char c) {
      return optregc(v,OPT_CSTRING,c);
}
int
optregs_CSTRING(char *v, char *n) {
      return optregs(v,OPT_CSTRING,n);
}
int
optregcb_CSTRING(char *v, char c, char *B) {
      return optregcb(v,OPT_CSTRING,c,B);
}
int
optregsb_CSTRING(char *v, char *s, char *B) {
      return optregsb(v,OPT_CSTRING,s,B);
}
int
optrega_UNDELIM(char **v, char c, char *n, char *B) {
      return optrega(v,OPT_UNDELIM,c,n,B);
}
int
optreg_UNDELIM(char **v, char c, char *B) {
      return optreg(v,OPT_UNDELIM,c,B);
}
int
optregc_UNDELIM(char **v, char c) {
      return optregc(v,OPT_UNDELIM,c);
}
int
optregs_UNDELIM(char **v, char *n) {
      return optregs(v,OPT_UNDELIM,n);
}
int
optregcb_UNDELIM(char **v, char c, char *B) {
      return optregcb(v,OPT_UNDELIM,c,B);
}
int
optregsb_UNDELIM(char **v, char *s, char *B) {
      return optregsb(v,OPT_UNDELIM,s,B);
}
int
optrega_VSTRING(char **v, char c, char *n, char *B) {
      return optrega(v,OPT_VSTRING,c,n,B);
}
int
optreg_VSTRING(char **v, char c, char *B) {
      return optreg(v,OPT_VSTRING,c,B);
}
int
optregc_VSTRING(char **v, char c) {
      return optregc(v,OPT_VSTRING,c);
}
int
optregs_VSTRING(char **v, char *n) {
      return optregs(v,OPT_VSTRING,n);
}
int
optregcb_VSTRING(char **v, char c, char *B) {
      return optregcb(v,OPT_VSTRING,c,B);
}
int
optregsb_VSTRING(char **v, char *s, char *B) {
      return optregsb(v,OPT_VSTRING,s,B);
}
int
optrega_DOUBLE(double *v, char c, char *n, char *B) {
      return optrega(v,OPT_DOUBLE,c,n,B);
}
int
optreg_DOUBLE(double *v, char c, char *B) {
      return optreg(v,OPT_DOUBLE,c,B);
}
int
optregc_DOUBLE(double *v, char c) {
      return optregc(v,OPT_DOUBLE,c);
}
int
optregs_DOUBLE(double *v, char *n) {
      return optregs(v,OPT_DOUBLE,n);
}
int
optregcb_DOUBLE(double *v, char c, char *B) {
      return optregcb(v,OPT_DOUBLE,c,B);
}
int
optregsb_DOUBLE(double *v, char *s, char *B) {
      return optregsb(v,OPT_DOUBLE,s,B);
}
int
optrega_NUL(void *v, char c, char *n, char *B) {
      return optrega(v,OPT_NUL,c,n,B);
}
int
optreg_NUL(void *v, char c, char *B) {
      return optreg(v,OPT_NUL,c,B);
}
int
optregc_NUL(void *v, char c) {
      return optregc(v,OPT_NUL,c);
}
int
optregs_NUL(void *v, char *n) {
      return optregs(v,OPT_NUL,n);
}
int
optregcb_NUL(void *v, char c, char *B) {
      return optregcb(v,OPT_NUL,c,B);
}
int
optregsb_NUL(void *v, char *s, char *B) {
      return optregsb(v,OPT_NUL,s,B);
}
int
optrega_BOOL(int *v, char c, char *n, char *B) {
      return optrega(v,OPT_BOOL,c,n,B);
}
int
optreg_BOOL(int *v, char c, char *B) {
      return optreg(v,OPT_BOOL,c,B);
}
int
optregc_BOOL(int *v, char c) {
      return optregc(v,OPT_BOOL,c);
}
int
optregs_BOOL(int *v, char *n) {
      return optregs(v,OPT_BOOL,n);
}
int
optregcb_BOOL(int *v, char c, char *B) {
      return optregcb(v,OPT_BOOL,c,B);
}
int
optregsb_BOOL(int *v, char *s, char *B) {
      return optregsb(v,OPT_BOOL,s,B);
}
int
optrega_UNDELIMV(char **v, char c, char *n, char *B) {
      return optrega(v,OPT_UNDELIMV,c,n,B);
}
int
optreg_UNDELIMV(char **v, char c, char *B) {
      return optreg(v,OPT_UNDELIMV,c,B);
}
int
optregc_UNDELIMV(char **v, char c) {
      return optregc(v,OPT_UNDELIMV,c);
}
int
optregs_UNDELIMV(char **v, char *n) {
      return optregs(v,OPT_UNDELIMV,n);
}
int
optregcb_UNDELIMV(char **v, char c, char *B) {
      return optregcb(v,OPT_UNDELIMV,c,B);
}
int
optregsb_UNDELIMV(char **v, char *s, char *B) {
      return optregsb(v,OPT_UNDELIMV,s,B);
}
int
optrega_LONG(long *v, char c, char *n, char *B) {
      return optrega(v,OPT_LONG,c,n,B);
}
int
optreg_LONG(long *v, char c, char *B) {
      return optreg(v,OPT_LONG,c,B);
}
int
optregc_LONG(long *v, char c) {
      return optregc(v,OPT_LONG,c);
}
int
optregs_LONG(long *v, char *n) {
      return optregs(v,OPT_LONG,n);
}
int
optregcb_LONG(long *v, char c, char *B) {
      return optregcb(v,OPT_LONG,c,B);
}
int
optregsb_LONG(long *v, char *s, char *B) {
      return optregsb(v,OPT_LONG,s,B);
}
int
optrega_UINT(unsigned int *v, char c, char *n, char *B) {
      return optrega(v,OPT_UINT,c,n,B);
}
int
optreg_UINT(unsigned int *v, char c, char *B) {
      return optreg(v,OPT_UINT,c,B);
}
int
optregc_UINT(unsigned int *v, char c) {
      return optregc(v,OPT_UINT,c);
}
int
optregs_UINT(unsigned int *v, char *n) {
      return optregs(v,OPT_UINT,n);
}
int
optregcb_UINT(unsigned int *v, char c, char *B) {
      return optregcb(v,OPT_UINT,c,B);
}
int
optregsb_UINT(unsigned int *v, char *s, char *B) {
      return optregsb(v,OPT_UINT,s,B);
}
int
optrega_UNSINT(unsigned int *v, char c, char *n, char *B) {
      return optrega(v,OPT_UNSINT,c,n,B);
}
int
optreg_UNSINT(unsigned int *v, char c, char *B) {
      return optreg(v,OPT_UNSINT,c,B);
}
int
optregc_UNSINT(unsigned int *v, char c) {
      return optregc(v,OPT_UNSINT,c);
}
int
optregs_UNSINT(unsigned int *v, char *n) {
      return optregs(v,OPT_UNSINT,n);
}
int
optregcb_UNSINT(unsigned int *v, char c, char *B) {
      return optregcb(v,OPT_UNSINT,c,B);
}
int
optregsb_UNSINT(unsigned int *v, char *s, char *B) {
      return optregsb(v,OPT_UNSINT,s,B);
}
int
optrega_STRING(char **v, char c, char *n, char *B) {
      return optrega(v,OPT_STRING,c,n,B);
}
int
optreg_STRING(char **v, char c, char *B) {
      return optreg(v,OPT_STRING,c,B);
}
int
optregc_STRING(char **v, char c) {
      return optregc(v,OPT_STRING,c);
}
int
optregs_STRING(char **v, char *n) {
      return optregs(v,OPT_STRING,n);
}
int
optregcb_STRING(char **v, char c, char *B) {
      return optregcb(v,OPT_STRING,c,B);
}
int
optregsb_STRING(char **v, char *s, char *B) {
      return optregsb(v,OPT_STRING,s,B);
}
int
optrega_UCHAR(unsigned char *v, char c, char *n, char *B) {
      return optrega(v,OPT_UCHAR,c,n,B);
}
int
optreg_UCHAR(unsigned char *v, char c, char *B) {
      return optreg(v,OPT_UCHAR,c,B);
}
int
optregc_UCHAR(unsigned char *v, char c) {
      return optregc(v,OPT_UCHAR,c);
}
int
optregs_UCHAR(unsigned char *v, char *n) {
      return optregs(v,OPT_UCHAR,n);
}
int
optregcb_UCHAR(unsigned char *v, char c, char *B) {
      return optregcb(v,OPT_UCHAR,c,B);
}
int
optregsb_UCHAR(unsigned char *v, char *s, char *B) {
      return optregsb(v,OPT_UCHAR,s,B);
}
int
optrega_SHORT(short *v, char c, char *n, char *B) {
      return optrega(v,OPT_SHORT,c,n,B);
}
int
optreg_SHORT(short *v, char c, char *B) {
      return optreg(v,OPT_SHORT,c,B);
}
int
optregc_SHORT(short *v, char c) {
      return optregc(v,OPT_SHORT,c);
}
int
optregs_SHORT(short *v, char *n) {
      return optregs(v,OPT_SHORT,n);
}
int
optregcb_SHORT(short *v, char c, char *B) {
      return optregcb(v,OPT_SHORT,c,B);
}
int
optregsb_SHORT(short *v, char *s, char *B) {
      return optregsb(v,OPT_SHORT,s,B);
}
int
optrega_USHORT(unsigned short *v, char c, char *n, char *B) {
      return optrega(v,OPT_USHORT,c,n,B);
}
int
optreg_USHORT(unsigned short *v, char c, char *B) {
      return optreg(v,OPT_USHORT,c,B);
}
int
optregc_USHORT(unsigned short *v, char c) {
      return optregc(v,OPT_USHORT,c);
}
int
optregs_USHORT(unsigned short *v, char *n) {
      return optregs(v,OPT_USHORT,n);
}
int
optregcb_USHORT(unsigned short *v, char c, char *B) {
      return optregcb(v,OPT_USHORT,c,B);
}
int
optregsb_USHORT(unsigned short *v, char *s, char *B) {
      return optregsb(v,OPT_USHORT,s,B);
}
int
optrega_INT(int *v, char c, char *n, char *B) {
      return optrega(v,OPT_INT,c,n,B);
}
int
optreg_INT(int *v, char c, char *B) {
      return optreg(v,OPT_INT,c,B);
}
int
optregc_INT(int *v, char c) {
      return optregc(v,OPT_INT,c);
}
int
optregs_INT(int *v, char *n) {
      return optregs(v,OPT_INT,n);
}
int
optregcb_INT(int *v, char c, char *B) {
      return optregcb(v,OPT_INT,c,B);
}
int
optregsb_INT(int *v, char *s, char *B) {
      return optregsb(v,OPT_INT,s,B);
}
int
optrega_NEGBOOL(int *v, char c, char *n, char *B) {
      return optrega(v,OPT_NEGBOOL,c,n,B);
}
int
optreg_NEGBOOL(int *v, char c, char *B) {
      return optreg(v,OPT_NEGBOOL,c,B);
}
int
optregc_NEGBOOL(int *v, char c) {
      return optregc(v,OPT_NEGBOOL,c);
}
int
optregs_NEGBOOL(int *v, char *n) {
      return optregs(v,OPT_NEGBOOL,n);
}
int
optregcb_NEGBOOL(int *v, char c, char *B) {
      return optregcb(v,OPT_NEGBOOL,c,B);
}
int
optregsb_NEGBOOL(int *v, char *s, char *B) {
      return optregsb(v,OPT_NEGBOOL,s,B);
}
int
optrega_UNDELIMC(char *v, char c, char *n, char *B) {
      return optrega(v,OPT_UNDELIMC,c,n,B);
}
int
optreg_UNDELIMC(char *v, char c, char *B) {
      return optreg(v,OPT_UNDELIMC,c,B);
}
int
optregc_UNDELIMC(char *v, char c) {
      return optregc(v,OPT_UNDELIMC,c);
}
int
optregs_UNDELIMC(char *v, char *n) {
      return optregs(v,OPT_UNDELIMC,n);
}
int
optregcb_UNDELIMC(char *v, char c, char *B) {
      return optregcb(v,OPT_UNDELIMC,c,B);
}
int
optregsb_UNDELIMC(char *v, char *s, char *B) {
      return optregsb(v,OPT_UNDELIMC,s,B);
}
int
optrega_ULONG(unsigned long *v, char c, char *n, char *B) {
      return optrega(v,OPT_ULONG,c,n,B);
}
int
optreg_ULONG(unsigned long *v, char c, char *B) {
      return optreg(v,OPT_ULONG,c,B);
}
int
optregc_ULONG(unsigned long *v, char c) {
      return optregc(v,OPT_ULONG,c);
}
int
optregs_ULONG(unsigned long *v, char *n) {
      return optregs(v,OPT_ULONG,n);
}
int
optregcb_ULONG(unsigned long *v, char c, char *B) {
      return optregcb(v,OPT_ULONG,c,B);
}
int
optregsb_ULONG(unsigned long *v, char *s, char *B) {
      return optregsb(v,OPT_ULONG,s,B);
}
int
optrega_FLOAT(float *v, char c, char *n, char *B) {
      return optrega(v,OPT_FLOAT,c,n,B);
}
int
optreg_FLOAT(float *v, char c, char *B) {
      return optreg(v,OPT_FLOAT,c,B);
}
int
optregc_FLOAT(float *v, char c) {
      return optregc(v,OPT_FLOAT,c);
}
int
optregs_FLOAT(float *v, char *n) {
      return optregs(v,OPT_FLOAT,n);
}
int
optregcb_FLOAT(float *v, char c, char *B) {
      return optregcb(v,OPT_FLOAT,c,B);
}
int
optregsb_FLOAT(float *v, char *s, char *B) {
      return optregsb(v,OPT_FLOAT,s,B);
}
int
optrega_NEGTOGGLE(int *v, char c, char *n, char *B) {
      return optrega(v,OPT_NEGTOGGLE,c,n,B);
}
int
optreg_NEGTOGGLE(int *v, char c, char *B) {
      return optreg(v,OPT_NEGTOGGLE,c,B);
}
int
optregc_NEGTOGGLE(int *v, char c) {
      return optregc(v,OPT_NEGTOGGLE,c);
}
int
optregs_NEGTOGGLE(int *v, char *n) {
      return optregs(v,OPT_NEGTOGGLE,n);
}
int
optregcb_NEGTOGGLE(int *v, char c, char *B) {
      return optregcb(v,OPT_NEGTOGGLE,c,B);
}
int
optregsb_NEGTOGGLE(int *v, char *s, char *B) {
      return optregsb(v,OPT_NEGTOGGLE,s,B);
}
int
optrega_INTLEVEL(int *v, char c, char *n, char *B) {
      return optrega(v,OPT_INTLEVEL,c,n,B);
}
int
optreg_INTLEVEL(int *v, char c, char *B) {
      return optreg(v,OPT_INTLEVEL,c,B);
}
int
optregc_INTLEVEL(int *v, char c) {
      return optregc(v,OPT_INTLEVEL,c);
}
int
optregs_INTLEVEL(int *v, char *n) {
      return optregs(v,OPT_INTLEVEL,n);
}
int
optregcb_INTLEVEL(int *v, char c, char *B) {
      return optregcb(v,OPT_INTLEVEL,c,B);
}
int
optregsb_INTLEVEL(int *v, char *s, char *B) {
      return optregsb(v,OPT_INTLEVEL,s,B);
}
int
optrega_TOGGLE(int *v, char c, char *n, char *B) {
      return optrega(v,OPT_TOGGLE,c,n,B);
}
int
optreg_TOGGLE(int *v, char c, char *B) {
      return optreg(v,OPT_TOGGLE,c,B);
}
int
optregc_TOGGLE(int *v, char c) {
      return optregc(v,OPT_TOGGLE,c);
}
int
optregs_TOGGLE(int *v, char *n) {
      return optregs(v,OPT_TOGGLE,n);
}
int
optregcb_TOGGLE(int *v, char c, char *B) {
      return optregcb(v,OPT_TOGGLE,c,B);
}
int
optregsb_TOGGLE(int *v, char *s, char *B) {
      return optregsb(v,OPT_TOGGLE,s,B);
}
/* OPT v3.19: options parsing tool */
/*
 * This SOFTWARE has been authored by an employee of the University of
 * California, operator of the Los Alamos National Laboratory under
 * Contract No. W-7405-ENG-36 with the U.S. Department of Energy. The
 * U.S. Government has rights to use, reproduce, and distribute this
 * SOFTWARE. Neither the Government nor the University makes any
 * warranty, express or implied, or assumes any liability or
 * responsibility for the use of this SOFTWARE. If SOFTWARE is modified
 * to produce derivative works, such modified SOFTWARE should be clearly
 * marked, so as not to confuse it with the version available from LANL.
 * 
 * Additionally, this program is free software; you can distribute it
 * and/or modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of the
 * License, or any later version.  Accordingly, this program is
 * distributed in the hope that it will be useful, but WITHOUT A
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 * for more details (http://www.gnu.org/copyleft/gpl.txt).
 * 
 * 
 * The OPT package is available from
 * http://nis-www.lanl.gov/~jt/Software
 */
/* $Id: opt.c,v 1.1.1.1 2005/10/29 19:20:46 ckemp Exp $ */
/* opt_menu.c */


#ifndef	SINGLEFILE
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "opt.h"
#include "ag.h"
#include "opt_p.h"
#endif /* SINGLEFILE */

/* --------------------- menu flag ------------------- */

int opt_menuflag=OPT_FALSE;

static char mgstring[160];	/* global menu string */

#define menu_wr_string(str)       opt_message(str)

static void write_the_menu(int iopfrom, int iopto);
static int auto_prefix_delim(char *r);


/* ----------------------------------------------------	*/
/*  opt_menu:	Get options from an interactive menu	*/
/* ----------------------------------------------------	*/

#define	MAXRESPONDLINE	280
#ifndef	MAXOPTSINMENU
#define	MAXOPTSINMENU	20
#endif

char *MENUPROMPT="-> ";

void
opt_menu(void)
{
    char respond[MAXRESPONDLINE+2];
    static int maxoptsinmenu=MAXOPTSINMENU;
    int	iopfrom,iopto;

    opt_menuflag=OPT_TRUE;	/* turn on MENUFLAG in case it is not set already */
    
    iopfrom = 0;
    iopto = ( opt_nreg < maxoptsinmenu ? opt_nreg : maxoptsinmenu );
    
    respond[0]='\0';

    opt_mess_1("%s\n",optgetTitle());
    write_the_menu(iopfrom,iopto);

    while( opt_menuflag ) {
        opt_menu_getline(MENUPROMPT,respond,MAXRESPONDLINE);
		
        switch(*respond) {
	case REASONCH:
	    opt_message( "Option Information: To be implemented.\n" );
	    opt_message( "\n   Name:        " );
	    opt_message( "\n   Type:        " );
	    opt_message( "\n   Reason:      " );
	    break;
        case ADDITIONAL_OPTS:
            if( respond[1] != '\0' && respond[1] != ' ' ) {
                maxoptsinmenu = atoi(respond+1);
                if(maxoptsinmenu < 1)
                    maxoptsinmenu = opt_nreg;
                opt_snprintf_1(mgstring,160,"Scroll %d options\n",maxoptsinmenu);
                menu_wr_string(mgstring);
                iopfrom = 0;
                iopto = 
                    ( opt_nreg < maxoptsinmenu ? opt_nreg : maxoptsinmenu );
            } else {
                iopfrom += maxoptsinmenu;
                if( iopfrom > opt_nreg)
                    iopfrom = 0;
                iopto = iopfrom + maxoptsinmenu;
                if( iopto > opt_nreg )
                    iopto = opt_nreg;
            }
            write_the_menu(iopfrom,iopto);
            break;
        case INTERACT:
            opt_menuflag=OPT_FALSE;
            break;
        case BANG:
            system( respond+1 );
            break;
        case '\0':
            write_the_menu(iopfrom,iopto);
            break;
        case QUITCH:
            /* Only quit if the QUITCH is followed by whitespace.  In
             * other words, if respond = '.m=5', don't quit.  However,
             * note that respond = '. m=5' will cause a quit.
             */
            if ( respond[1]=='\0' || respond[1]==' ' )
                opt_quit();
            else
                opt_mess_1("Invalid line: [%s]\n",respond);
            break;
        case DELIM:
#ifdef PERMIT_ALTDELIM            
        case ALTDELIM:
#endif            
            opt_lineprocess(respond);
            break;
        default:
            auto_prefix_delim(respond);
            opt_lineprocess(respond);
            break;
        }
    }
}/* opt_menu */


/**********
 * write_the_menu:
 *	write the menu including options from iopfrom to iopto.
 */
static void
write_the_menu(int iopfrom, int iopto)
{
    int iop;
    int fullmenu;

    fullmenu = ((iopfrom==0 && iopto==opt_nreg) ? OPT_TRUE : OPT_FALSE );

    if( !fullmenu ) {
        opt_snprintf_3(mgstring,160,"menu: %d->%d [%d]\n",iopfrom,iopto,opt_nreg);
        menu_wr_string(mgstring);
    }

    for(iop=iopfrom; iop<iopto; ++iop) {
        char *s;
        s = opt_mstring(iop);
        if (s!=NULL) {
            menu_wr_string(s);
            menu_wr_string("\n");
	    free(s);
        }
    }
    if (!fullmenu) {
        opt_snprintf_1(mgstring,160,"%c Additional options\n",ADDITIONAL_OPTS);
        menu_wr_string(mgstring);
    }
    opt_snprintf_1(mgstring,160,"(Type %c for Help)\n",HELPCH);
    menu_wr_string(mgstring);
}

/*	auto_prefix_delim:	
 *		this is a fru-fru piece of code that automatically
 *		sticks a DELIM character onto the front of a string
 *		in cases where it imagines that that is what the user
 *		really meant.  Thus
 *		-> m4
 *		gives the same effect as
 *		-> -m4
 *		But be warned that this only applies in limited cases.
 *		Eg.,
 *		-> m4 b3
 *		is not the same as
 *		-> -m4 -b3
 *
 *		a '-' will be prepended in the case that 
 *		the first character is a registered name
 */
static int
auto_prefix_delim(char *r)
{
    if( opt_char_number( *r ) != -1 ) {
        int len;
        len = strlen(r)+1;	/* +1 since double terminated */
        while(len>=0) {
            r[len+1]=r[len];
            --len;
        }
        r[0]=DELIM;
        return(1);
    }
    else
        return(0);
}/* auto_prefix_delim */


/* OPT v3.19: options parsing tool */
/*
 * This SOFTWARE has been authored by an employee of the University of
 * California, operator of the Los Alamos National Laboratory under
 * Contract No. W-7405-ENG-36 with the U.S. Department of Energy. The
 * U.S. Government has rights to use, reproduce, and distribute this
 * SOFTWARE. Neither the Government nor the University makes any
 * warranty, express or implied, or assumes any liability or
 * responsibility for the use of this SOFTWARE. If SOFTWARE is modified
 * to produce derivative works, such modified SOFTWARE should be clearly
 * marked, so as not to confuse it with the version available from LANL.
 * 
 * Additionally, this program is free software; you can distribute it
 * and/or modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of the
 * License, or any later version.  Accordingly, this program is
 * distributed in the hope that it will be useful, but WITHOUT A
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 * for more details (http://www.gnu.org/copyleft/gpl.txt).
 * 
 * 
 * The OPT package is available from
 * http://nis-www.lanl.gov/~jt/Software
 */

/* $Id: opt.c,v 1.1.1.1 2005/10/29 19:20:46 ckemp Exp $ */
/* opt_proc.c */
/*
 *	process options in the OPT package
 */

#ifndef	SINGLEFILE
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "opt.h"
#include "ag.h"
#include "opt_p.h"
#endif /* SINGLEFILE */

#ifndef OPT_VERSION
#define OPT_VERSION VERSION
#endif

/******************
 * Global variables 
 */
OPT_PFI  opt_additional_usage_fcn=NULL;
char	*opt_program_name = NULL;

/* Global to this file */

static int   opt_fileflag = OPT_FALSE;
static int   opt_menu_enabled = OPT_TRUE;
static char *opt_filename = NULL;
static char *opt_defaultfile = NULL;
static char *opt_titlestring = NULL;
static char *opt_usagestring = NULL;
static char *opt_envstring=NULL;
static char *opt_defstring=NULL;
static char *opt_pkg_version=NULL;

static OPT_PFI	   opt_quit_fcn=NULL;
static OPT_PFI	   opt_run_fcn=NULL;
static OPT_PFI_ARG opt_main_fcn=NULL;

/****************************
 * Static Function prototypes 
 */

static char *append_string(char *, char *);
static int  line2argv(int *, char ***, char *);
static int  l2a(char *, char **);
static int  break_word(int, char *);
static int  backup_file(char *fname);

/**********************
 * free various strings
 */

void
opt_freestrings()
{
  OPT_FREE( opt_filename );
  OPT_FREE( opt_defaultfile );
  OPT_FREE( opt_titlestring );
  OPT_FREE( opt_usagestring );
  OPT_FREE( opt_envstring );
  OPT_FREE( opt_defstring );
  OPT_FREE( opt_pkg_version );
}

/***************************
 * set title, usage, etc.  *
 ***************************/
void
optTitle(char *s)
{
  opt_setstring(&opt_titlestring,s);
}
char *
optgetTitle(void)
{
    if (!ISEMPTYSTRING(opt_titlestring))
        return opt_titlestring;
    else
        return opt_program_name;
}
void
optVersion(char *s)
{
  if (!ISEMPTYSTRING(s)) opt_setstring(&opt_pkg_version,s);
}
void
optUsage(char *s)
{
  opt_setstring(&opt_usagestring,s);
}
void
optEnvVarName(char *s)
{
  opt_setstring(&opt_envstring,s);
}
void
optDefaultString(char *s)
{
  opt_setstring(&opt_defstring,s);
}
void
optDefaultFile(char *s)
{
  opt_setstring(&opt_defaultfile,s);
}
void
optQuit(OPT_PFI fcn)
{
    opt_quit_fcn = fcn;
}
void
optRun(OPT_PFI fcn)
{
    opt_run_fcn = fcn;
}
void
optMain(OPT_PFI_ARG fcn)
{
    opt_main_fcn = fcn;
}
void
optAdditionalUsage(OPT_PFI fcn)
{
    opt_additional_usage_fcn = fcn;
}
void
optDisableMenu()
{
    opt_menu_enabled = OPT_FALSE;
}


static void
optWriteVersion()
{
  char *opt_pkg;
  opt_pkg = optgetTitle();
  
#ifdef VERSION
  if (ISEMPTYSTRING(opt_pkg_version)) {
    opt_pkg_version=VERSION;
  }
#endif

  if (!ISEMPTYSTRING(opt_pkg_version)) {
    opt_mess_2("%s: Version %s\n",opt_pkg,opt_pkg_version);
  } else {
    opt_mess_1("%s: no version information\n",opt_pkg);
  }
}



/* ------------------------------------------------------------------ */

/* getopts:
 * This routine is maintained only for backward compatibility.  It
 * takes argc, argv and returns the new value of argc; Now, it is
 * recommended that you use opt(), which takes &argc,&argv as arguments,
 * and then returns new argc,argv through the pointers on the arg list.
 */
int
getopts(int argc, char **argv)
{
    
	/*
	 * Before processing, set the global variables
	 * opt_program_name : name of routine that is running
	 * opt_filename     : default name of options file
	 */
    if (ISEMPTYSTRING(opt_program_name))
      opt_program_name = short_progname(argv[0]);
    opt_filename = append_string(opt_program_name,OPT_EXT);
    opt_readline_init(opt_program_name);
                                               
	/* Begin options processing */
	/* ------------------------ */

    /* First process default string */
    opt_lineprocess(opt_defstring);

    /* Second, check environment variable */
    if(opt_envstring != NULL) {
	char	*s;
	s = getenv(opt_envstring);
	opt_lineprocess(s);
    }

    /* Finally, parse the command line */
    opt_process(argc-1,argv+1);

    return(argc);
}
/*
 * opt: This is the routine that does all the work.  It is generally called
 * from main(), after all calls to optreg() have been made.
 * It returns nothing, but new argc/argv are returned via the arglist
 */
void
opt(int *argcptr, char ***argvptr)
{
    int	 argc;
    char **argv;
    ARGVECTOR *ag;
    
    argc = *argcptr;
    argv = *argvptr;

    /*
     * Before processing, set the global variables
     * opt_program_name : name of routine that is running
     * opt_filename     : default name of options file
     */
    if (ISEMPTYSTRING(opt_program_name))
        opt_program_name = short_progname(argv[0]);
    opt_filename = append_string(opt_program_name,OPT_EXT);
    opt_readline_init(opt_program_name);
                                               
    /* Begin options processing */
    /* ------------------------ */

    /* First process default string */
    opt_lineprocess(opt_defstring);

    /* Second, check the default file (eg ~/.programrc) */
    if (!ISEMPTYSTRING(opt_defaultfile)) {
        opt_fromfname(opt_defaultfile);
    }
    
    /* Third, check environment variable */
    if(opt_envstring != NULL)
    {
	char	*s;
	s = getenv(opt_envstring);
	opt_lineprocess(s);
    }

    /* Finally, parse the command line */
    if (argc == 0) return;      /* assuming there is one... */
    ag = opt_process(argc-1,argv+1);
    /* And, return the leftover argc/argv, but first prepend argv[0]
     * to this!
     */
    if (ag == NULL) {
        *argcptr = 1;
        *argvptr = argv;
    } else {
        /* make a new argc,argv */
        ARGVECTOR *agnew;
        agnew = ag_prepend_argv0(ag,argv[0]);
        *argcptr = ag_argc(agnew);
        *argvptr = ag_argv(agnew);
        ag_free(agnew);         /* don't ag_freeall(agnew) because
                                 * that will ruin the new argv strings
                                 */
    }
    ag_free(ag);

    /* Now that we are done processing, free what we no longer need */
    OPT_FREE(opt_filename);
}
static char *
append_string(char *s, char *t)
{
    /* input two strings "s" and "t":
     * concatenates them to get string "st"
     * which it allocates space for and returns
     */
    char *st;
    st = (char *)malloc( strlen(s)+strlen(t)+1 );
    strcpy(st,s);
    strcat(st,t);
    return(st);
}

/* Allow the user to explcitly set the program name */
/* Note that this routine strips off the leading directory */
/* so that you can safely use optProgName(argv[0]) */
void optProgName(char* progname) {
  opt_program_name = short_progname(progname);
}


char *
short_progname(char *pname)
{
    char *p;
    int	len;

    /*
     * inelegant routine
     * returns the shortened name of the program named [pname].
     * It (recursively) strips off leading directory or drive number
     * and trailing ".EXE" in the case of MS-DOS executable
     */

    p = pname;
    while( *p != '/' && *p!='\\' && *p!=':' && *p!='\0' )
    	++p;
    if( *p=='\0') {	/* remove .EXE if it exists */
    	len = strlen(pname);
    	if( len > 4 ) {
    	    if( ( 0==strcmp(pname+len-4,".EXE") ||
    	    	  0==strcmp(pname+len-4,".exe") ) )
    	    		pname[len-4]='\0';
	}
    	return( pname );
    }
    else
    	return( short_progname(++p) );
}
/* ----------------------------------------------------------------- */

	/* ------- */
	/* process */
	/* ------- */

/****************************************************************
 * opt_process(): 
 *	INPUT:
 *		argc, argv;
 *
 *	this is the central receiving facility for the options package.
 *	opt_process takes an argument list given by (argc,argv) and
 *	and processes it,
 *	Although the input may come from either of
 *		system command line
 *		options file
 *		interactive command line
 *	the behavior is slightly different if the global opt_menuflag is set
 */

ARGVECTOR *
opt_process(int argc, char **argv)
{
    char	c;
    int		iword;
    ARGVECTOR	*ag;
    
    /*
     * convert the argument vector (argc,argv) into
     * an ARGVECTOR structure, to ease manipulations
     */
    
    if( argc==0 )
        return (ARGVECTOR *)NULL;
    
    ag = ag_new(argc,argv);

    if(DEBUG) {
        ag_fprint(stderr,ag);
    }

    /* Loop through the options in the argument vector */
    
    while( !ag_end(ag) ) {	
        if( ag_eow(ag)	)		/* if necessary...      */
            ag_w_advance(ag);	        /* advance to next word */
		
        if(DEBUG)
            if( ag->ic != 0 )		/* this should never happen */
                opt_warning("opt_proc: not on first character");
		
        c=ag_c_advance(ag);		/* first character of word */
		
        /* ------ */
        /* switch */
        /* ------ */
		
        switch( c ) {
        case DELIM:
#ifdef PERMIT_ALTDELIM            
        case ALTDELIM:
#endif            
            if( ag_eow(ag) ) {
                /* unadorned "-" will be treated as a positional option */
                /* Try to process it as an positional option */
	        ag_backspace(ag);
                if( opt_parse_positional(ag) == 0) {
                    /* if not successful...then stop processing */
                    return ag;
                }
            } else if( ag_c(ag) == c ) {
                ag_c_advance(ag);          /* eat the '-' */
                if (ag_eow(ag)) {          /* if unadorned "--" */
                    ag_w_advance(ag);      /* eat the entire "--" */
                    return ag;             /* and stop processing */
                    
                    /* if "--something" */                    
                } else {
                    /* Special cases: --help, --optVersion, --menu, etc.. */
                    if (0==strcmp("help",ag_s(ag))) { 
                        long_usage();
                        if (!opt_menuflag) {
                            ag_free(ag);
                            optAbortRun();
                        }
                    } else if (0==strcmp("xml-opt-interface",ag_s(ag))) {
                        opt_write_xml_file(stdout);
                        if (!opt_menuflag)
                            optAbortRun();
                    } else if (0==strcmp("optVersion",ag_s(ag))) {
#ifdef OPT_VERSION                        
                        opt_mess_1("OPT Version %s\n",OPT_VERSION);
#endif
                        if (!opt_menuflag) {
                            ag_free(ag);
                            optAbortRun();
                        }
                    } else if (0==strcmp("menu",ag_s(ag))) { /* --menu */
			if (!opt_menu_enabled) {
			    opt_warning("menu not available for this application");
			    ag_w_advance(ag);
			}
		        else {
			    if ( !opt_menuflag )      /* If not called already */
				ag_w_advance(ag); /* Swallow the '--menu' */
			    opt_menu();       /* and call the menu */	
		        }
                    } else if (0==strcmp("version",ag_s(ag)) &&
                               !ISEMPTYSTRING(opt_pkg_version)) {
                      optWriteVersion();
                      if (!opt_menuflag) {
                          ag_free(ag);
                          optAbortRun();
                      }
                    } else {
                        /* if a generic "--longoption" */
                        opt_parse_longdelim(ag);
                    }
                }
                /* if "-something" */
            } else {
                iword = ag_w_number(ag);
                while( iword==ag_w_number(ag) && !ag_eow(ag) )
                    opt_parse_delim(ag);
            }
            break;

        case OPTFROMFILE:
            if (! opt_fromfname(argnext(ag))) {
                opt_warn_1("Cannot open options file [%s]\n",opt_filename);
            }
            break;
            
        case OPTTOFILE:
            opt_tofname(argnext(ag));
            if ( !opt_menuflag )
              opt_quit();
            break;

        case HELPCH:
	    if( !opt_fileflag)
		opt_help(argnext(ag));
            break;

        case INTERACT:
            if ( !opt_menuflag )		/* If not called already */
                opt_menu();				/* Call the menu */	
            else {						/* otherwise */ 
                opt_menuflag=OPT_FALSE;			/* turn off menu */
                ag_clear(ag);
            }
            break;
        case IGNOREEOL:
            ag_clear(ag);
            break;

        case RUN:
	    /* If we're on the command line then treat this as a positional option */
	    if (!(opt_menuflag || opt_fileflag))
		goto POSITIONAL;    /* I used goto! Hey, that felt ... good */

            /* Give precedence to optMain, then optRun */
            if ( opt_main_fcn != NULL ) {
                int retval;
                int ac;
                char **av;
                ac = ag_argc(ag);
                av = ag_argv(ag);
                ag_clear(ag);
                retval = optBeginMain(opt_main_fcn,ac,av);
                if (retval == OPT_ABORT)
                    optAbortMain();
                if (retval == OPT_QUIT)
                    opt_quit();
                if (retval == OPT_EXIT)
                    exit(opt_exit_number);
            }
            else if ( opt_run_fcn != NULL ) {
                int retval;
                retval = optBeginRun(opt_run_fcn);
                if (retval == OPT_ABORT)
                    optAbortRun();
                if (retval == OPT_QUIT)
                    opt_quit();
                if (retval == OPT_EXIT)
                    exit(opt_exit_number);
            } 
            else {
                /* if opt_run_fcn is not set, and 'RUN' is called
                 * from the menu, then exit the menu, so that
                 * the routine is run, and then program will exit
                 */
                if ( opt_menuflag ) {
                    opt_menuflag=OPT_FALSE;			/* turn off menu */
                    ag_clear(ag);
                } else {
                    /* RUN called from a file or the command line:
                     * There is no reason to be doing this.
                     */
                    opt_warning("No Run Function has been registered");
                }
                break;
				
	    /* Note that the QUITCH '.' is quite a common positional option! */
            /* So we only treat it as a quit signal if using the menu */
            case QUITCH:
		if (!opt_menuflag)
		    goto POSITIONAL;   /* There I go again... */
                opt_quit();
                break;
				
            case BANG:
                /* If there is a BANG as the first character of a
                 * menuline, then we shouldn't even get here.
                 */
                if (!opt_menuflag)
		    goto POSITIONAL;     /* Can't ... stop ... */
                else {
                    opt_warn_1("Shell to system invalid unless \'%c\' is the first character in the line\n",BANG);
                }
#if 0
                /* bang...shell to system */
                system(argnext(ag));
                break;
#endif

            case ' ':		/* blanks and empty characters, ignore */
            case '\t':
            case '\0':
                break;
		
	    POSITIONAL:
            default:
                /* in the default case, option may be positional */
                /* ---------------------------------------------- */
                if(DEBUG) ag_fprint(stderr,ag);
                ag_backspace(ag);
                if(DEBUG) {
                    opt_mess_1("opt_proc: in default section [%s]\n",ag_s(ag));
                    ag_fprint(stderr,ag);
                }
                /* Try to process it as an positional option */
                if( opt_parse_positional(ag) == 0) {
                    /* if not successful...then stop processing */
                    if(DEBUG) ag_fprint(stderr,ag);
                    return ag;
                }
                break;
                /* last case */
            }
        }
    }
    ag_free(ag);
    return (ARGVECTOR *)NULL;
}/*opt_process*/

/* opt_lineprocess()
 * process a string, by converting it first into an argument vector
 * then calling opt_process
 * return the number of processed arguments
 */

int
opt_lineprocess(char *line)
{
    int i,nargc;
    char **nargv;

    if( line==NULL || *line=='\0') return(0);
 
    line2argv(&nargc,&nargv,line);
    opt_process(nargc,nargv);
    /* Now free nargv */
    for (i=0; i<nargc; ++i) {
      free(nargv[i]);
    }
    free(nargv);    
    return(nargc);
}
/*	line2argv(): converts a one-line string into an argument vector */

static int
line2argv(int *pargc, char ***pargv, char *line)
{
    char *lincp;
    if(line==NULL) {	
	*pargc = 0;
	**pargv = NULL;
	return(0);
    }
	/*
	 * First thing to do is copy the line into
	 * a buffer ("lincp") so that input line will
 	 * not be corrupted -- also so that input line
 	 * doesn't have to be double null terminated
 	 */
    lincp = (char *)malloc((strlen(line)+2)*sizeof(char));
    strcpy(lincp,line);
	/*
	 * Next step is to double null terminate 
	 * the copied line
	 */
    lincp[strlen(line)+1]='\0';
    	/* count words in line */
    *pargc = l2a(lincp,NULL);
    	/* allocate array for arg vector */
    *pargv = (char **)malloc(((*pargc)+1)*sizeof(char *));
    	/* fill arg vector with words in line */
    l2a(lincp,*pargv);
        /* we are done with lincp, free it... */
    free(lincp);

    return(*pargc);
}

#define	QUOTE	'\"'
#define	BKSLASH	'\\'

/* l2a */
/*	l2a() is slave routine to line2argv()           */
/*	if argv==NULL then count the words in the line	*/
/*	if argv!=NULL then put those words into argv[]  */
/* 	WARNINGS:
 *	l2a assumes that input line is double-null terminated!!
 *	the line buffer is pointed to by argv[]'s so do not
 *	free the input line buffer
 */
	/* quoted material counts as a single word */
	/* so that -s"string with spaces"-g is parsed */
	/* into 
		-s
	        string with spaces
		-g
	*/
	/* Comment: should be able to escape with backslash */

static int
l2a(char *line, char **argv)
{
    int inquote=OPT_FALSE;	/* flag: inside quotation */
    int iarg;
    
    for(iarg=0; *line != '\0'; ++iarg) {
        if(!inquote) {
            while( isspace(*line) )
                ++line;	/* skip leading blank spaces */
            if( *line == QUOTE ) {
                inquote=OPT_TRUE;
                ++line;	/* skip past leading quote */
            }
        }
        if(argv!=NULL)			/* point to */
            argv[iarg] = line;	/* first character of line */
		
        while( !break_word(inquote,line) )
            ++line;
		
        if( *line==QUOTE ) {		/* toggle quote */ 
            inquote = (inquote==OPT_TRUE ? OPT_FALSE : OPT_TRUE );
            if(argv==NULL)	++line;	/* skip past quote */
        }
		
        if(argv!=NULL) {
            *line='\0';	/* Null terminate string */
            ++line;		/* and go to next character */
        }
        else 
            while( isspace(*line) )
                ++line;		/* skip trailing blanks */
		
    }
    /* We do the following to make the line free-able */
    if (argv != NULL) {
      int i;
      for (i=0; i<iarg; ++i) {
        argv[i] = strdup(argv[i]);
      }
    }
    /* Of course we will still need to free the argv's... */
    
    return(iarg);
}
static int
break_word(int inquote, char *line)
{
    if( *line == '\0' || *line==QUOTE )
        return(OPT_TRUE);
    if( !inquote && isspace(*line) )
        return(OPT_TRUE);
    
    return(OPT_FALSE);
}


	/* --------------------------- */
	/* write out a usage statement */
	/* --------------------------- */

/* A global function which can be used to programatically */
/* print the usage. */
void optPrintUsage() {
  long_usage();
}


void
short_usage(void)
{
  opt_message(opt_titlestring);
  opt_mess_1("Usage: %s ",opt_program_name);
  opt_mess_1("%s\n", opt_usagestring ? opt_usagestring : "[options] ");
  if (opt_menu_enabled)
      opt_mess_1("To invoke the menu, type:\n\t%s --menu\n", opt_program_name);
}	
void
long_usage(void)
{
    short_usage();
    opt_usage();
}

	/* --------------------------------- */
	/* put current options in named file */
	/* --------------------------------- */

int
opt_tofname(char *fname)
{
    FILE	*fp;
    
    if( *fname != OPTTOFILE ) {
	free((char *)opt_filename);
	opt_filename = strdup(fname);
    }
    backup_file(opt_filename);
    fp = fopen(opt_filename,"w");

    if(fp==NULL) {
	opt_mess_1("cant open opt file [%s] for writing\n",fname);
	return 0;
    }
    opt_tofile(fp);
    fclose(fp);
    return 1;
}
static int
backup_file(char *fname)
{
    int c;
    FILE *fp,*fp_bak;
    char *fname_bak;
    
    fp = fopen(fname,"r");
    if (fp == NULL)
        return 1;               /* file fname doesn't exist */
    fname_bak = append_string(fname,"~");
    fp_bak = fopen(fname_bak,"w");
    if (fp_bak == NULL) {
        fclose(fp);
        return 1;               /* cant open backup file */
    }
    /* copy file */
    while ((c = getc(fp)) != EOF) {
        putc(c,fp_bak);
    }
    opt_mess_1("Options file [%s] has been overwritten\n",fname_bak);
    fclose(fp_bak);
    fclose(fp);
    free((char *)fname_bak);    /* allocated in append_string */
    return 0;                   /* success */                   
}

	/* ----------------------------------------- */
	/* get and process options from a named file */
	/* ----------------------------------------- */

#define	MAXOPLINE	2560
int
opt_fromfname(char *fname)
{
    FILE	*fp;
    char	line[MAXOPLINE+2];	/* leave room for termination */
    int         tmp_fileflag;
    
    if (ISEMPTYSTRING(fname)) {
        opt_warn_1("%cfile.opt fails to specify file.opt\n",OPTFROMFILE);
        return 0;
    }
    if (*fname == OPTFROMFILE) {      /* if '@@' */
        fp = fopen(opt_filename,"r");
    } else {                          /* if '@file.opt' */
        fp = fopen(fname,"r");
    }
    if (fp == NULL) {
        return 0;
    }

    tmp_fileflag = opt_fileflag;
    opt_fileflag = OPT_TRUE;
    while( opt_getline(line,fp) > 0 )
    	opt_lineprocess(line);
    opt_fileflag = tmp_fileflag;
    
    fclose(fp);
    return 1;                   /* success */
}
int
opt_getline(char *line, FILE *fp)
{
    int	c;
    int	count=0;
	
    while( (c=getc(fp)) != '\n' && c!=EOF ) {
        *line = (char)c;
        ++line;
        if( ++count > MAXOPLINE-3 ) {
            opt_warning("Options line too long");
            break;
        }
    }
    *line++ = '\0';		/* double terminate line !! */
    *line++ = '\0';
    
    return(count);
}/*opt_getline*/

/* ------------------- */
/* basic help function */
/* ------------------- */

void
opt_help(char *s)
{
    if(s==NULL || *s=='\0') {
        opt_mess_2("\t%c %-20s\n",DELIM,"Options delimiter"); 
        opt_mess_2("\t%c %-20s\n",HELPCH,"Help");
        opt_mess_2("\t%c %-20s\n",RUN,"Run program and return to menu");
        opt_mess_2("\t%c %-20s\n",BANG,"Shell to Operating System"); 
        if (opt_menuflag)
            opt_mess_2("\t%c %-20s\n",INTERACT,"Exit menu"); 
        else
            opt_mess_2("\t%c %-20s\n",INTERACT,"Invoke menu"); 
        opt_mess_2("\t%c %-20s\n",ADDITIONAL_OPTS,"Additional options");
        opt_mess_3("\t%c<filename> %-20s [%s]\n",OPTFROMFILE,
                   "Get options from file",opt_filename);
        opt_mess_4("\t%c%c %-2s [%s]\n",OPTFROMFILE,OPTFROMFILE,
                   "Get options from file",opt_filename);
        opt_mess_2("\t%c<filename> %-20s\n",OPTTOFILE,
                   "Put options in file");
        opt_mess_4("\t%c%c %-2s [%s]\n",OPTTOFILE,OPTTOFILE,
                   "Put options in file",opt_filename);
        opt_mess_2("\t%c %-20s\n",QUITCH,"Quit");
    } else {
        if( s[1]=='\0')
            opt_get_help(s[0]); 
        else {
            opt_mess_1("Help unavailable for: [%s]\n",s);
        }
    }
}
/* ------ */
/* quit() */ 
/* ------ */
void
opt_quit(void)
{
    int retval;

    if( opt_quit_fcn != NULL )
        retval = (*opt_quit_fcn)();

    /* Let's free up memory */
    opt_free();
    /* Do we want to do anything with the retval? */
    /* eg, exit vs optAbortRun?? */
    exit(opt_exit_number);
}

/* OPT v3.19: options parsing tool */
/*
 * This SOFTWARE has been authored by an employee of the University of
 * California, operator of the Los Alamos National Laboratory under
 * Contract No. W-7405-ENG-36 with the U.S. Department of Energy. The
 * U.S. Government has rights to use, reproduce, and distribute this
 * SOFTWARE. Neither the Government nor the University makes any
 * warranty, express or implied, or assumes any liability or
 * responsibility for the use of this SOFTWARE. If SOFTWARE is modified
 * to produce derivative works, such modified SOFTWARE should be clearly
 * marked, so as not to confuse it with the version available from LANL.
 * 
 * Additionally, this program is free software; you can distribute it
 * and/or modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of the
 * License, or any later version.  Accordingly, this program is
 * distributed in the hope that it will be useful, but WITHOUT A
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 * for more details (http://www.gnu.org/copyleft/gpl.txt).
 * 
 * 
 * The OPT package is available from
 * http://nis-www.lanl.gov/~jt/Software
 */
/* $Id: opt.c,v 1.1.1.1 2005/10/29 19:20:46 ckemp Exp $ */
/* ag.c */

#ifndef	SINGLEFILE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "opt.h"                /* needed for __PROTO */
#include "ag.h"
#endif

/* primitive routines for manipulation argument vectors */

/*
  make-new-argvector              ag_new()
  reset-argvector                 ag_reset()
  flag-end-of-word                ag_eow()
  what-is-character               ag_c()
  what-is-string                  ag_s()

  read-character                  ag_c()
  read-next-character             ag_c_next()
  read-next-non-null-character    ag_cnn_next()
  read-character-advance          ag_c_advance()
  read-string                     ag_s()
  read-next-string                ag_s_next()
  read-string-advance             ag_s_advance()
  word-number                     ag_w_number()
  word-advance                    ag_w_advance()
  backspace                       ag_backspace()
  clear                           ag_clear()

  get-argc                        ag_argc()
  get-argv                        ag_argv()

  get-next-argument               argnext()
*/

ARGVECTOR *
ag_new(int argc, char **argv)
{
    ARGVECTOR *ag;
    ag = (ARGVECTOR *)malloc(sizeof(ARGVECTOR));
    /* Note...we do not actually make a COPY of argv; but this
     * shouldn't be necessary, since *none* of the following routines
     * should ever write to ag->c and *especially* not ag->v
     * ...um, actually ag_prepend_argv0 does
     */
    ag->v = argv;
    ag->c = argc;
    ag->iw = 0;
    ag->ic = 0;
    return(ag);
}
void
ag_free(ARGVECTOR *ag)
{
    free((char *)ag);
}
void
ag_fprint(FILE *fp, ARGVECTOR *ag)
{
        int i;
        if (ag == NULL) {
            fprintf(fp,"ag=<nul>\n");
            return;
        }
        for(i=0; i<ag->c; ++i)
            fprintf(fp,"[%s]",ag->v[i]);
        fprintf(fp,", <%d> <%d>\n",ag->iw,ag->ic);
}
void
ag_reset(ARGVECTOR *ag)
{
    ag->iw = 0;
    ag->ic = 0;
}
int
ag_w_number(ARGVECTOR *ag)
{
    return ag->iw;
}
void
ag_w_advance(ARGVECTOR *ag)
{
    ++(ag->iw);             /* advance to next word */
    ag->ic=0;               /* beginning of next word */
}
/*      ag_eow: end of word
        flag whether current position is at end of word
*/
int
ag_eow(ARGVECTOR *ag)
{
    if( ag->iw >= ag->c )
        return(TRUE);
    if( ag->ic >= (int)strlen(ag->v[ag->iw]) )
        return(TRUE);
    return(FALSE);
}
/*      ag_end: end of command line
                flag whether current position is at end of command line 
*/
int
ag_end(ARGVECTOR *ag)
{
        if( ag->iw >= ag->c )
                return(TRUE);
        if( ag_eow(ag) && ag->iw == (ag->c)-1 )
                return(TRUE);
        return(FALSE);
}
/*      ag_c:   return current character
                do not advance
*/
char
ag_c(ARGVECTOR *ag)
{
        return(ag->v[ag->iw][ag->ic]);
}
char
ag_c_next(ARGVECTOR *ag)
{
        return(ag->v[ag->iw][ag->ic+1]);
}
char
ag_cnn_next(ARGVECTOR *ag)
{
    if( ag_c_next(ag) == '\0' ) {
        if(ag->iw+1 >= ag->c)
            return('\0');
        else
            return(ag->v[ag->iw+1][0]);
    }
    else
        return( ag_c_next(ag) );
}
/*      ag_c_advance:   read current character, and advance to next
                        return '\0' if end of word
                        do not advance to next word
*/
char
ag_c_advance(ARGVECTOR *ag)
{
        char    c;              /* current character */

        if( ag_eow(ag) )
                return(c='\0'); /* return NULL to signal that current*/
                                /* character is past end of word     */
        c = ag->v[ag->iw][ag->ic];
        ++(ag->ic);             /* advance to next character */
        return(c);
}
void
ag_backword(ARGVECTOR *ag)
{
    ag->iw -= 1;
}

char
ag_backspace(ARGVECTOR *ag)
{
    if( --(ag->ic) < 0 ) {          /* if back past beginning of word */
        ag->ic=0;
        if(--(ag->iw) < 0) {        /* goto beginning of previous word */
            ag_reset(ag);           /* if no previous, reset */
        }
        else {
            while( !ag_eow(ag) )    /* goto end of prevous word */
                ag_c_advance(ag);
            ag_backspace(ag);       /* back to just before end */
        }
    }
    return(ag->v[ag->iw][ag->ic]);
}
/*      ag_s:   returns current string
                returns "" if current position is at end of word
                returns NULL if past end of argument vector
*/
char    *
ag_s(ARGVECTOR *ag)
{
        if( ag->iw < ag->c )
                return( ag->v[ag->iw]+ag->ic );
        else    return( NULL );
}
char    *
ag_s_next(ARGVECTOR *ag)
{
        if( ag->v[ag->iw][ag->ic+1] == '\0' )
        {
                if(ag->iw+1 >= ag->c)
                        return(NULL);
                else    return(ag->v[ag->iw+1]);
        }
        else    return( ag->v[ag->iw]+ag->ic+1 );
}
/*      ag_s_advance:   read current string and advance to next
                        returns NULL if current string is at end
                        does not check advanced string
*/
char    *
ag_s_advance(ARGVECTOR *ag)
{
    char    *s;             /* current string */
    
    if( ag_eow(ag) )        /* if end of word, go to next word */
        ag_w_advance(ag);
    if( ag_eow(ag) )        /* if still at end of word,  */
        s = NULL;       /* signify by returning NULL */
    else
        s = ag->v[ag->iw] + ag->ic;
    
    ag_w_advance(ag);       /* advance to next word */
    
    return(s);
}
int
ag_clear(ARGVECTOR *ag)
{
    while( !ag_end(ag) )
        argnext(ag);
    return 1;
}

int
ag_argc(ARGVECTOR *ag)
{
    return ag->c - ag->iw;
}
char **
ag_argv(ARGVECTOR *ag)
{
    return ag->v + ag->iw;
}
ARGVECTOR *
ag_copy(ARGVECTOR *ag)
{
    /* Here is how you would make a copy of an ARGVECTOR */
    /* To free this, do not use ag_free, but ag_freeall */
    int i;
    ARGVECTOR *agnew;
    agnew = (ARGVECTOR *)malloc(sizeof(ARGVECTOR));
    agnew->c = ag->c;
    agnew->v = (char **)malloc(sizeof(char *));
    for (i=0; i< agnew->c; ++i) {
        agnew->v[i] = strdup(ag->v[i]);
    }
    agnew->iw = 0;
    agnew->ic = 0;
    return agnew;
}
void
ag_freeall(ARGVECTOR *ag)
{
    int i;
    for (i=0; i<ag->c; ++i) {
        free(ag->v[i]);
    }
    free((char *)ag->v);
    free((char *)ag);
}
ARGVECTOR *
ag_prepend_argv0(ARGVECTOR *ag, char *argv0)
{
    int i;
    ARGVECTOR *agnew;
    agnew = (ARGVECTOR *)malloc(sizeof(ARGVECTOR));
    agnew->c = ag->c - ag->iw + 1;
    agnew->v = (char **)malloc(agnew->c * sizeof(char *));
    agnew->v[0] = argv0; /* NOTE: we are not strdup'ing */
    for (i=1; i< agnew->c; ++i) {
        agnew->v[i] = ag->v[i + ag->iw - 1];
    }
    agnew->iw = 0;
    agnew->ic = 0;
    return agnew;
}

    


	/* ------------------------ */
	/* return the next argument */
	/* ------------------------ */
char *
argnext(ARGVECTOR *ag)
{
    static char nullchar='\0';
    char *s;

    s = ag_s_advance(ag);
    if( s==NULL )	
    	s = &nullchar;
    return(s);
}




/* opt_regcc.cc */
/* Note: this file generated automatically from makeopth.pl */
/* Edits to this file may be lost! */
/* C++ bindings for opt */
#ifndef OPT_NOCPLUSINTERFACE
#ifdef __cplusplus
#include <cstring>
#include "opt.h"
using namespace std;
int OptRegister(char *v, char c, char *n, char *B)
{
    return optrega(v,OPT_CHAR,c,n,B);
}
int OptRegister(char *v, char *n, char *B)
{
    return optrega(v,OPT_CHAR,'\0',n,B);
}
int OptRegister(char *v, opt_TYPE o, char c, char *n, char *B)
{
    return optrega(v,o,c,n,B);
}
int OptRegister(char *v, opt_TYPE o, char *n, char *B)
{
    return optrega(v,o,'\0',n,B);
}
int OptRegister(char *v, opt_MODE m, char c, char *n, char *B)
{
    int i = optrega(v,OPT_CHAR,c,n,B);
    optmode_n(i, m);
    return i;
}
int OptRegister(char *v, opt_MODE m, char *n, char *B)
{
    int i = optrega(v,OPT_CHAR,'\0',n,B);
    optmode_n(i, m);
    return i;
}
int OptRegister(char *v, opt_TYPE o, opt_MODE m, char c, char *n, char *B)
{
    int i = optrega(v,o,c,n,B);
    optmode_n(i, m);
    return i;
}
int OptRegister(char *v, opt_TYPE o, opt_MODE m, char *n, char *B)
{
    int i = optrega(v,o,'\0',n,B);
    optmode_n(i, m);
    return i;
}
int OptRegister(short *v, char c, char *n, char *B)
{
    return optrega(v,OPT_SHORT,c,n,B);
}
int OptRegister(short *v, char *n, char *B)
{
    return optrega(v,OPT_SHORT,'\0',n,B);
}
int OptRegister(short *v, opt_TYPE o, char c, char *n, char *B)
{
    return optrega(v,o,c,n,B);
}
int OptRegister(short *v, opt_TYPE o, char *n, char *B)
{
    return optrega(v,o,'\0',n,B);
}
int OptRegister(short *v, opt_MODE m, char c, char *n, char *B)
{
    int i = optrega(v,OPT_SHORT,c,n,B);
    optmode_n(i, m);
    return i;
}
int OptRegister(short *v, opt_MODE m, char *n, char *B)
{
    int i = optrega(v,OPT_SHORT,'\0',n,B);
    optmode_n(i, m);
    return i;
}
int OptRegister(short *v, opt_TYPE o, opt_MODE m, char c, char *n, char *B)
{
    int i = optrega(v,o,c,n,B);
    optmode_n(i, m);
    return i;
}
int OptRegister(short *v, opt_TYPE o, opt_MODE m, char *n, char *B)
{
    int i = optrega(v,o,'\0',n,B);
    optmode_n(i, m);
    return i;
}
int OptRegister(unsigned short *v, char c, char *n, char *B)
{
    return optrega(v,OPT_USHORT,c,n,B);
}
int OptRegister(unsigned short *v, char *n, char *B)
{
    return optrega(v,OPT_USHORT,'\0',n,B);
}
int OptRegister(unsigned short *v, opt_TYPE o, char c, char *n, char *B)
{
    return optrega(v,o,c,n,B);
}
int OptRegister(unsigned short *v, opt_TYPE o, char *n, char *B)
{
    return optrega(v,o,'\0',n,B);
}
int OptRegister(unsigned short *v, opt_MODE m, char c, char *n, char *B)
{
    int i = optrega(v,OPT_USHORT,c,n,B);
    optmode_n(i, m);
    return i;
}
int OptRegister(unsigned short *v, opt_MODE m, char *n, char *B)
{
    int i = optrega(v,OPT_USHORT,'\0',n,B);
    optmode_n(i, m);
    return i;
}
int OptRegister(unsigned short *v, opt_TYPE o, opt_MODE m, char c, char *n, char *B)
{
    int i = optrega(v,o,c,n,B);
    optmode_n(i, m);
    return i;
}
int OptRegister(unsigned short *v, opt_TYPE o, opt_MODE m, char *n, char *B)
{
    int i = optrega(v,o,'\0',n,B);
    optmode_n(i, m);
    return i;
}
int OptRegister(int *v, char c, char *n, char *B)
{
    return optrega(v,OPT_INT,c,n,B);
}
int OptRegister(int *v, char *n, char *B)
{
    return optrega(v,OPT_INT,'\0',n,B);
}
int OptRegister(int *v, opt_TYPE o, char c, char *n, char *B)
{
    return optrega(v,o,c,n,B);
}
int OptRegister(int *v, opt_TYPE o, char *n, char *B)
{
    return optrega(v,o,'\0',n,B);
}
int OptRegister(int *v, opt_MODE m, char c, char *n, char *B)
{
    int i = optrega(v,OPT_INT,c,n,B);
    optmode_n(i, m);
    return i;
}
int OptRegister(int *v, opt_MODE m, char *n, char *B)
{
    int i = optrega(v,OPT_INT,'\0',n,B);
    optmode_n(i, m);
    return i;
}
int OptRegister(int *v, opt_TYPE o, opt_MODE m, char c, char *n, char *B)
{
    int i = optrega(v,o,c,n,B);
    optmode_n(i, m);
    return i;
}
int OptRegister(int *v, opt_TYPE o, opt_MODE m, char *n, char *B)
{
    int i = optrega(v,o,'\0',n,B);
    optmode_n(i, m);
    return i;
}
int OptRegister(unsigned long *v, char c, char *n, char *B)
{
    return optrega(v,OPT_ULONG,c,n,B);
}
int OptRegister(unsigned long *v, char *n, char *B)
{
    return optrega(v,OPT_ULONG,'\0',n,B);
}
int OptRegister(unsigned long *v, opt_TYPE o, char c, char *n, char *B)
{
    return optrega(v,o,c,n,B);
}
int OptRegister(unsigned long *v, opt_TYPE o, char *n, char *B)
{
    return optrega(v,o,'\0',n,B);
}
int OptRegister(unsigned long *v, opt_MODE m, char c, char *n, char *B)
{
    int i = optrega(v,OPT_ULONG,c,n,B);
    optmode_n(i, m);
    return i;
}
int OptRegister(unsigned long *v, opt_MODE m, char *n, char *B)
{
    int i = optrega(v,OPT_ULONG,'\0',n,B);
    optmode_n(i, m);
    return i;
}
int OptRegister(unsigned long *v, opt_TYPE o, opt_MODE m, char c, char *n, char *B)
{
    int i = optrega(v,o,c,n,B);
    optmode_n(i, m);
    return i;
}
int OptRegister(unsigned long *v, opt_TYPE o, opt_MODE m, char *n, char *B)
{
    int i = optrega(v,o,'\0',n,B);
    optmode_n(i, m);
    return i;
}
int OptRegister(double *v, char c, char *n, char *B)
{
    return optrega(v,OPT_DOUBLE,c,n,B);
}
int OptRegister(double *v, char *n, char *B)
{
    return optrega(v,OPT_DOUBLE,'\0',n,B);
}
int OptRegister(double *v, opt_TYPE o, char c, char *n, char *B)
{
    return optrega(v,o,c,n,B);
}
int OptRegister(double *v, opt_TYPE o, char *n, char *B)
{
    return optrega(v,o,'\0',n,B);
}
int OptRegister(double *v, opt_MODE m, char c, char *n, char *B)
{
    int i = optrega(v,OPT_DOUBLE,c,n,B);
    optmode_n(i, m);
    return i;
}
int OptRegister(double *v, opt_MODE m, char *n, char *B)
{
    int i = optrega(v,OPT_DOUBLE,'\0',n,B);
    optmode_n(i, m);
    return i;
}
int OptRegister(double *v, opt_TYPE o, opt_MODE m, char c, char *n, char *B)
{
    int i = optrega(v,o,c,n,B);
    optmode_n(i, m);
    return i;
}
int OptRegister(double *v, opt_TYPE o, opt_MODE m, char *n, char *B)
{
    int i = optrega(v,o,'\0',n,B);
    optmode_n(i, m);
    return i;
}
int OptRegister(long *v, char c, char *n, char *B)
{
    return optrega(v,OPT_LONG,c,n,B);
}
int OptRegister(long *v, char *n, char *B)
{
    return optrega(v,OPT_LONG,'\0',n,B);
}
int OptRegister(long *v, opt_TYPE o, char c, char *n, char *B)
{
    return optrega(v,o,c,n,B);
}
int OptRegister(long *v, opt_TYPE o, char *n, char *B)
{
    return optrega(v,o,'\0',n,B);
}
int OptRegister(long *v, opt_MODE m, char c, char *n, char *B)
{
    int i = optrega(v,OPT_LONG,c,n,B);
    optmode_n(i, m);
    return i;
}
int OptRegister(long *v, opt_MODE m, char *n, char *B)
{
    int i = optrega(v,OPT_LONG,'\0',n,B);
    optmode_n(i, m);
    return i;
}
int OptRegister(long *v, opt_TYPE o, opt_MODE m, char c, char *n, char *B)
{
    int i = optrega(v,o,c,n,B);
    optmode_n(i, m);
    return i;
}
int OptRegister(long *v, opt_TYPE o, opt_MODE m, char *n, char *B)
{
    int i = optrega(v,o,'\0',n,B);
    optmode_n(i, m);
    return i;
}
int OptRegister(unsigned int *v, char c, char *n, char *B)
{
    return optrega(v,OPT_UINT,c,n,B);
}
int OptRegister(unsigned int *v, char *n, char *B)
{
    return optrega(v,OPT_UINT,'\0',n,B);
}
int OptRegister(unsigned int *v, opt_TYPE o, char c, char *n, char *B)
{
    return optrega(v,o,c,n,B);
}
int OptRegister(unsigned int *v, opt_TYPE o, char *n, char *B)
{
    return optrega(v,o,'\0',n,B);
}
int OptRegister(unsigned int *v, opt_MODE m, char c, char *n, char *B)
{
    int i = optrega(v,OPT_UINT,c,n,B);
    optmode_n(i, m);
    return i;
}
int OptRegister(unsigned int *v, opt_MODE m, char *n, char *B)
{
    int i = optrega(v,OPT_UINT,'\0',n,B);
    optmode_n(i, m);
    return i;
}
int OptRegister(unsigned int *v, opt_TYPE o, opt_MODE m, char c, char *n, char *B)
{
    int i = optrega(v,o,c,n,B);
    optmode_n(i, m);
    return i;
}
int OptRegister(unsigned int *v, opt_TYPE o, opt_MODE m, char *n, char *B)
{
    int i = optrega(v,o,'\0',n,B);
    optmode_n(i, m);
    return i;
}
int OptRegister(char **v, char c, char *n, char *B)
{
    return optrega(v,OPT_STRING,c,n,B);
}
int OptRegister(char **v, char *n, char *B)
{
    return optrega(v,OPT_STRING,'\0',n,B);
}
int OptRegister(char **v, opt_TYPE o, char c, char *n, char *B)
{
    return optrega(v,o,c,n,B);
}
int OptRegister(char **v, opt_TYPE o, char *n, char *B)
{
    return optrega(v,o,'\0',n,B);
}
int OptRegister(char **v, opt_MODE m, char c, char *n, char *B)
{
    int i = optrega(v,OPT_STRING,c,n,B);
    optmode_n(i, m);
    return i;
}
int OptRegister(char **v, opt_MODE m, char *n, char *B)
{
    int i = optrega(v,OPT_STRING,'\0',n,B);
    optmode_n(i, m);
    return i;
}
int OptRegister(char **v, opt_TYPE o, opt_MODE m, char c, char *n, char *B)
{
    int i = optrega(v,o,c,n,B);
    optmode_n(i, m);
    return i;
}
int OptRegister(char **v, opt_TYPE o, opt_MODE m, char *n, char *B)
{
    int i = optrega(v,o,'\0',n,B);
    optmode_n(i, m);
    return i;
}
int OptRegister(float *v, char c, char *n, char *B)
{
    return optrega(v,OPT_FLOAT,c,n,B);
}
int OptRegister(float *v, char *n, char *B)
{
    return optrega(v,OPT_FLOAT,'\0',n,B);
}
int OptRegister(float *v, opt_TYPE o, char c, char *n, char *B)
{
    return optrega(v,o,c,n,B);
}
int OptRegister(float *v, opt_TYPE o, char *n, char *B)
{
    return optrega(v,o,'\0',n,B);
}
int OptRegister(float *v, opt_MODE m, char c, char *n, char *B)
{
    int i = optrega(v,OPT_FLOAT,c,n,B);
    optmode_n(i, m);
    return i;
}
int OptRegister(float *v, opt_MODE m, char *n, char *B)
{
    int i = optrega(v,OPT_FLOAT,'\0',n,B);
    optmode_n(i, m);
    return i;
}
int OptRegister(float *v, opt_TYPE o, opt_MODE m, char c, char *n, char *B)
{
    int i = optrega(v,o,c,n,B);
    optmode_n(i, m);
    return i;
}
int OptRegister(float *v, opt_TYPE o, opt_MODE m, char *n, char *B)
{
    int i = optrega(v,o,'\0',n,B);
    optmode_n(i, m);
    return i;
}
int OptRegister(unsigned char *v, char c, char *n, char *B)
{
    return optrega(v,OPT_UCHAR,c,n,B);
}
int OptRegister(unsigned char *v, char *n, char *B)
{
    return optrega(v,OPT_UCHAR,'\0',n,B);
}
int OptRegister(unsigned char *v, opt_TYPE o, char c, char *n, char *B)
{
    return optrega(v,o,c,n,B);
}
int OptRegister(unsigned char *v, opt_TYPE o, char *n, char *B)
{
    return optrega(v,o,'\0',n,B);
}
int OptRegister(unsigned char *v, opt_MODE m, char c, char *n, char *B)
{
    int i = optrega(v,OPT_UCHAR,c,n,B);
    optmode_n(i, m);
    return i;
}
int OptRegister(unsigned char *v, opt_MODE m, char *n, char *B)
{
    int i = optrega(v,OPT_UCHAR,'\0',n,B);
    optmode_n(i, m);
    return i;
}
int OptRegister(unsigned char *v, opt_TYPE o, opt_MODE m, char c, char *n, char *B)
{
    int i = optrega(v,o,c,n,B);
    optmode_n(i, m);
    return i;
}
int OptRegister(unsigned char *v, opt_TYPE o, opt_MODE m, char *n, char *B)
{
    int i = optrega(v,o,'\0',n,B);
    optmode_n(i, m);
    return i;
}
#endif /* __cpluscplus */
#endif /* OPT_NOCPLUSINTERFACE */
/* OPT v3.19: options parsing tool */
/*
 * This SOFTWARE has been authored by an employee of the University of
 * California, operator of the Los Alamos National Laboratory under
 * Contract No. W-7405-ENG-36 with the U.S. Department of Energy. The
 * U.S. Government has rights to use, reproduce, and distribute this
 * SOFTWARE. Neither the Government nor the University makes any
 * warranty, express or implied, or assumes any liability or
 * responsibility for the use of this SOFTWARE. If SOFTWARE is modified
 * to produce derivative works, such modified SOFTWARE should be clearly
 * marked, so as not to confuse it with the version available from LANL.
 * 
 * Additionally, this program is free software; you can distribute it
 * and/or modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of the
 * License, or any later version.  Accordingly, this program is
 * distributed in the hope that it will be useful, but WITHOUT A
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 * for more details (http://www.gnu.org/copyleft/gpl.txt).
 * 
 * 
 * The OPT package is available from
 * http://nis-www.lanl.gov/~jt/Software
 */
/* $Id: opt.c,v 1.1.1.1 2005/10/29 19:20:46 ckemp Exp $ */
/* isreal: returns whether a string corresponds to
 * a real number
 */

#ifndef SINGLEFILE
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <assert.h>
#include "opt.h"
#include "ag.h"
#include "opt_p.h"
#endif

#define IsSpace(t) isspace((int)(t))
#define IsDigit(t) isdigit((int)(t))
#define IsSign(t)  ((t)=='+' || (t)=='-')

int
opt_isvalidnumber(char *t)
{
  int status = OPT_NUM_VALID; /* innocent until proven guilty */
  
  if (!t) return FALSE;       /* null string is not a number */

  while( *t != '\0' && IsSpace(*t) ) ++t;   /* leading white space */
  
  if   ( *t != '\0' && IsSign(*t)  ) {      /* leading sign */
    if ( *t == '-' ) status |= OPT_NUM_NEG;
    ++t;  
  }

  if   (*t=='e' || *t=='E' ) return FALSE;  /* 'e5' not allowed! */
  
  while( *t != '\0' && IsDigit(*t) ) ++t;   /* many digits */
  if   ( *t != '\0' && ( *t=='.' ) ) {      /* at most one decimal point */
    status |= OPT_NUM_FLOAT;
    ++t;
  }
  while( *t != '\0' && IsDigit(*t) ) ++t;   /* more digits */

  /* check for exponent */
  if   (*t=='e' || *t=='E' ) {
    ++t;
    if ( *t == '\0' ) return FALSE;  /* '5e' isn't a number! */
    if ( *t != '\0' && IsSign(*t) ) {  
      if ( *t == '-' ) status |= OPT_NUM_FLOAT;  /* '1e-3' = 0.001 is FLOAT */
      ++t;
    }
    while( *t !=  '\0' && IsDigit(*t) ) ++t;
  }

  /* trailing space is okay */
  while( *t != '\0' && IsSpace(*t) ) ++t; 

  /* if made it to the end then its a real */
  return (*t=='\0') ? status : FALSE;

}

/* return numeric value of string: basically a wrapper on atof() */

double
opt_atof(char *s)
{
  int valid;
  double x;
  x = atof(s); /* Call atof() whether or not string is valid,
                * because some strings, eg '15x' can still be
                * interpreted as numeric.  But give a warning
                * unless the string really is valid
                */
  valid = opt_isvalidnumber(s);
  if (!valid) {
    opt_warn_2("String [%s] is not a valid number, will use [%g]",s,x);
  }
  return x;
}
long
opt_atoi(char *s)
{
  int valid;
  long x;
  x = (long)atof(s); /* Call atof() whether or not string is valid,
                      * because some strings, eg '15x' can still be
                      * interpreted as numeric.  But give a warning
                      * unless the string really is valid
                      */
  
  valid = opt_isvalidnumber(s);
  if (!valid || (valid & OPT_NUM_FLOAT)) {
    opt_warn_2("String [%s] is not a valid integer, will use [%ld]",s,x);
  }
  return x;
}

unsigned long
opt_atou(char *s)
{
  int valid;
  unsigned long x;
  x = (unsigned long)atof(s);
  valid = opt_isvalidnumber(s);
  if (!valid ||
      (valid & OPT_NUM_FLOAT) ||
      (valid & OPT_NUM_NEG)) {
    opt_warn_2("String [%s] is not a valid positive integer, will use [%lu]",s,x);
  }
  return x;
}

