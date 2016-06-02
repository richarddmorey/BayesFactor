/*
 * progress.cpp
 *
 * A Front-end class for InterruptableProgressMonitor.
 *
 * Author: karl.forner@gmail.com
 * Small edits: richarddmorey@gmail.com
 * Taken from CRAN package RcppProgress 0.1
 *
 */
#include "progress.h"


InterruptableProgressMonitor* Progress::_monitor_singleton = 0;
