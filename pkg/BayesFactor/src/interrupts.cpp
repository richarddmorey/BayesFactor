/*
 * Check for user interruption in C++ code interrupting execution of the current code
 *
 * This code has been written by Simon Urbanek
 * I took it from the R-devel mailing list
 * in the thread "[Rd] Interrupting C++ code execution"
 * The mail is dated: 25 april 2011
 *
 * It allows to check for user interruption without
 * leaving the c++ function that calls it.
 *
 * Potential drawbacks according to its author:
 * The problem with it is that it will eat all errors, even if they were not yours
 * (e.g. those resulting from events triggered the event loop), so I would not recommend it for general use.
 *
 * Taken from CRAN package RcppProgress 0.1
 */
#include "interrupts.h"
#include <R.h>
#include "bfcommon.h"

using namespace Rcpp;

static void chkIntFn(void *dummy) {
  R_CheckUserInterrupt();
}

// this will call the above in a top-level context so it won't longjmp-out of your context
bool checkInterrupt() {
	return (R_ToplevelExec(chkIntFn, NULL) == FALSE);
}

