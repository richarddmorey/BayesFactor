/*
 * interruptable_progress_monitor.hpp
 *
 * A class that monitor the progress of computations:
 *   - can display a progress bar
 *   - can handle user (R user level) or programmatic abortion
 *   - can be used in OpenMP loops
 *
 * Author: karl.forner@gmail.com
 * Small edits: richarddmorey@gmail.com
 * Taken from CRAN package RcppProgress 0.1
 */
#ifndef _RcppProgress_INTERRUPTABLE_PROGRESS_MONITOR_HPP
#define _RcppProgress_INTERRUPTABLE_PROGRESS_MONITOR_HPP

#include <Rcpp.h>
#include "interrupts.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

class InterruptableProgressMonitor {
public: // ====== LIFECYCLE =====

	/**
	 * Main constructor
	 *
	 * @param max the expected number of tasks to perform
	 * @param display_progress whether to display a progress bar in the console
	 */
	InterruptableProgressMonitor(unsigned long max = 1,  bool display_progress = true)  {
		reset(max, display_progress);
		display_progress_bar();
	}

	~InterruptableProgressMonitor() { }

public: // ===== ACCESSORS/SETTERS =====
	void set_display_status(bool on) { _display_progress = on; 	}
	bool is_display_on() const { return _display_progress; }
	unsigned long get_max() const { return _max; }
	bool is_aborted() const { return _abort; }



public: // ===== PBLIC MAIN INTERFACE =====
	/**
	 * increment the current progress.
	 *
	 * Iff called by the master thread, it will also update the display if needed
	 *
	 * @param amount the number of newly performed tasks to report
	 *
	 * @return false iff the computation is aborted
	 */
	bool increment(unsigned long amount=1) {
		if ( is_aborted() )
			return false;
		return is_master() ? update_master(_current + amount) : atomic_increment(amount);
	}

	/**
	 * set the current progress indicator
	 *
	 * Iff called by the master thread, it will also update the display if needed
	 *
	 * @param current the total number of performed tasks so far (by all threads)
	 *
	 * @return false iff the computation is aborted
	 */
	bool update(unsigned long current) {
		if ( is_aborted() )
			return false;
		return is_master() ? update_master(current) : atomic_update(current);
	}

	/**
	 * check that the no interruption has been requested and return the current status
	 *
	 * Iff called by the master thread, it will check for R-user level interruption.
	 *
	 * @return true iff the computation is aborted
	 */
	bool check_abort() {
		if ( is_aborted() )
			return true;

		if ( is_master() )  {
			check_user_interrupt_master();
			update_display();
		}
		return is_aborted();
	}

	/**
	 * request computation abortion
	 */
	void abort() {
#pragma omp critical
		_abort = true;

	}

	/**
	 * return true iff the thread is the master.
	 * In case of non-OpenMP loop, always return true
	 */
	bool is_master() const {
#ifdef _OPENMP
		return omp_get_thread_num() == 0;
#else
		return true;
#endif
	}

public: // ===== methods for MASTER thread =====

	/**
	 * set the current progress indicator and update the progress bar display if needed.
	 *
	 *
	 * @param current the total number of performed tasks
	 *
	 * @return false iff the computation is aborted
	 */
	bool update_master(unsigned long current) {
		// try to make it as fast as possible
		unsigned long last = _current;
		_current = current;
		if ( (current - last)*100 > _max )
			update_display();
		return ! is_aborted();
	}

	void check_user_interrupt_master() {
		if ( !is_aborted() && checkInterrupt() ) {
			//REprintf("detected User interruption...\n");
			abort();
		}
	}

public: // ===== methods for non-MASTER threads =====

	bool atomic_increment(unsigned long amount=1) {
#pragma omp atomic
		_current+=amount;
		return ! is_aborted();
	}

	bool atomic_update(unsigned long current) {
#pragma omp critical
		_current=current;
		return ! is_aborted();
	}

public: // ===== methods related to DISPLAY, should not be called directly =====

	void update_display() {
		if ( !is_display_on() )
			return;
		int nb_ticks = _compute_nb_ticks(_current) - _compute_nb_ticks(_last_displayed);
		if (nb_ticks > 0) {
			_last_displayed = _current;
			_display_ticks(nb_ticks);
		}

		if ( _current >= _max )
			end_display();
	}

	void end_display() {
		if ( !is_display_on() )
			return;
		if ( ! is_aborted() ) {
			// compute the remaining ticks and display them
			int remaining = 50 - _compute_nb_ticks(_last_displayed);
			_display_ticks(remaining);
		}
		REprintf("|\n");
	}

	void display_progress_bar() {
		if ( !is_display_on() )
			return;
		REprintf("0%   10   20   30   40   50   60   70   80   90   100%\n");
		REprintf("|----|----|----|----|----|----|----|----|----|----|\n");
	}

protected: // ==== other instance methods =====

	int _compute_nb_ticks(unsigned long current) {
		return current * 50 / _max;
	}

	void _display_ticks(int nb) {
		for (int i = 0; i < nb; ++i)
			REprintf("*");
	}

	/**
	 * reset the monitor.
	 *
	 * Currently not really useful
	 *
	 * @param max the expected number of tasks to perform
	 * @param display_progress whether to display a progress bar in the console
	 *
	 */
	void reset(unsigned long max = 1, bool display_progress = true) {
		_max = max;
		if ( _max <= 0 )
			_max = 1;
		_last_displayed = _current = 0;
		_display_progress = display_progress;
		_abort = false;
	}


private: // ===== INSTANCE VARIABLES ====
	unsigned long _max; 			// the nb of tasks to perform
	unsigned long _current; 		// the current nb of tasks performed
	unsigned long _last_displayed; 	// the nb of tasks last displayed
	bool _abort;					// whether the process should abort
	bool _display_progress;			// whether to display the progress bar

};

#endif
