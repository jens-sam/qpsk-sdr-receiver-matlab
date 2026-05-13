# QPSK SDR Receiver with Adaptive Equalization

This repository contains a portfolio snapshot of a MATLAB-based software-defined radio receiver project for QPSK digital communication signal recovery.

The project focused on building a complete receiver pipeline capable of recovering transmitted text from RF signal files under different channel impairment conditions, including intersymbol interference, carrier frequency offset, constellation rotation, and symbol-rate mismatch.

## Repository Note

This repository is intended for viewing and portfolio demonstration.

The full course framework and original signal data files may not be included in this public repository. The included files and report are intended to show the receiver design, signal processing workflow, MATLAB implementation approach, and project results.

## Project Overview

The receiver was developed in stages, beginning with RF signal analysis and baseband conversion, then extending into synchronization, equalization, carrier phase correction, timing recovery, and adaptive processing.

The project used MATLAB to implement and test a QPSK receiver across multiple received signal files with increasing levels of impairment.

## Receiver Pipeline

The receiver included the following major processing stages:

- RF spectrum analysis
- Carrier downconversion to baseband
- Matched filtering
- Symbol timing selection
- Preamble detection
- Payload extraction
- QPSK symbol decisions
- Bitstream reconstruction
- Text recovery
- Symbol-spaced LMS equalization
- Half-symbol-spaced equalization
- Carrier frequency offset estimation
- Decision-directed phase tracking
- Symbol-rate tracking
- Adaptive equalizer-based impairment compensation

## Technical Features

### Downconversion and Matched Filtering

The received RF signal was downconverted to complex baseband and filtered using a matched filter to improve symbol recovery.

### Timing Acquisition

The oversampled signal was evaluated across possible sampling phases, and the symbol timing phase was selected by maximizing decimated signal energy.

### Preamble Detection

A known repeated pilot sequence was used to detect the start of the packet. Preamble detection allowed the receiver to align the payload and train equalizers before message recovery.

### Adaptive Equalization

The project implemented LMS-based adaptive equalization to reduce intersymbol interference and improve constellation clustering.

Both symbol-spaced and half-symbol-spaced equalizer structures were explored.

### Carrier Phase Tracking

For longer received signals, residual carrier frequency offset caused constellation rotation. A decision-directed phase tracking loop was implemented to stabilize the constellation and improve decoding.

### Symbol-Rate Tracking

For signals with symbol-rate mismatch, timing drift accumulated over the packet. A timing recovery loop was used to track and compensate for symbol timing drift.

## Skills Demonstrated

- MATLAB
- Digital signal processing
- Software-defined radio
- QPSK modulation and demodulation
- Digital communication receiver design
- Matched filtering
- Synchronization
- Preamble detection
- LMS adaptive filtering
- Adaptive equalization
- Carrier frequency offset estimation
- Decision-directed phase tracking
- Symbol timing recovery
- Constellation analysis
- Eye diagram analysis
- RF signal processing
- Wireless communications

## Project Report

A technical project report is included to document the receiver design, processing stages, plots, results, and performance observations.

The report covers:

- demodulation and timing acquisition,
- preamble detection,
- symbol-spaced equalization,
- fractionally-spaced equalization,
- carrier phase tracking,
- symbol-rate tracking,
- adaptive equalizer performance,
- final receiver observations and conclusions.

## Portfolio Summary

This project demonstrates the design and implementation of a MATLAB-based digital communication receiver for QPSK RF signal recovery. The receiver combines synchronization, adaptive filtering, equalization, carrier tracking, and timing recovery to recover transmitted messages from impaired signal files.
