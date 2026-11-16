# STA496 Independent Reading Course

**Simultaneous Sieve Estimation and Inference for Time-Varying Nonlinear Regression with Locally Stationary Time Series**  
Fall 2025 · Department of Statistical Sciences, University of Toronto

## Overview

This repository documents my study of the simultaneous sieve inference method (SIMle) proposed in [Simultaneous Sieve Estimation and Inference for Time-Varying Nonlinear Time Series Regression](https://arxiv.org/pdf/2506.23069), as part of an independent reading course supervised by [Dr. Zhou Zhou](https://sites.google.com/view/zhouzhou-stat/home) at the University of Toronto.

It includes my implementation of SIMle, simulation scripts for reproducing the studies and results in the paper, and data processing code for a real-world application. It also contains my final report and the slides I presented at the department's [Reading Course Presentations](https://www.statistics.utoronto.ca/Undergraduate-Research#fall25-11).

## Presentation Abstract

Many interesting phenomena can naturally be viewed as time series, but they often evolve with time and in a way that is hard to parametrize. Such complex temporal dynamics make linear or parametric models less suitable. To this end, scientists have developed tools for nonlinear nonstationary time series. Yet, some gaps remain open-- uniform inference for the estimators and their practical implementation. In “Simultaneous Sieve Estimation and Inference for Time-Varying Nonlinear Regression with Locally Stationary Time Series” (Ding & Zhou, 2025), the authors propose a sieve-based estimator and present both theoretical and empirical results for uniform inference of this estimator. In this talk, I will review the theoretical framework in this proposal as well as its empirical performance across a range of scenarios.

- [Presentation slides](Katrina%20Sha%20STA496_Final_Presentation.pdf)
- [Final report](STA496_Final_Report.pdf)

## Repository Structure

| Directory                                  | Contents                                                                                         |
| ------------------------------------------ | ------------------------------------------------------------------------------------------------ |
| [`data/`](data/)                           | Hong Kong respiratory disease dataset and data processing code                                   |
| [`method/`](method/)                       | My implementation of SIMle                                                                       |
| [`sim_res/`](sim_res/)                     | Selected results from simulations run locally                                                    |
| [`simulation_script/`](simulation_script/) | Scripts for running simulations locally and on a remote high-performance computing (HPC) cluster |
