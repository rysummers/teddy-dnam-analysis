# TEDDY Longitudinal DNA Methylation Analysis

## Overview

This project investigates developmental DNA methylation (DNAm) trajectories in 
the TEDDY (The Environmental Determinants of Diabetes in the Young) cohort.

Longitudinal mixed-effects spline models are used to estimate CpG-specific
methylation trajectories across early childhood. Predicted trajectories are 
subsequently evaluated for evidence of latent methylation strata, filtered to
remove non-informative patterns, and clustered to identify common developmental
methylation profiles.

The overall goal is to characterize large-scale patterns of DNAm change over 
time and identify distinct classes of developmental methylation trajectories.

## Objectives

* Model longitudinal methylation trajectories across early childhood.
* Identify CpGs exhibiting sub-trajectory patterns (strata).
* Remove non-informative ("flat") trajectories.
* Cluster CpGs according to developmental trajectory shape.
* Characterize the dynamics and biological features of resulting trajectory classes.

## Analysis Workflow

```text
M-value Matrix
        ↓
Preprocessing
        ↓
CpG-wise Mixed-Effects Spline Models
        ↓
Predicted Methylation Trajectories
        ↓
Trajectory Filtering
    • Flat CpGs
    • Sub-trajectory CpGs (GMM)
        ↓
Trajectory Clustering
        ↓
Trajectory Annotation (characterization)

```

### Preprocessing

Current preprocessing steps include:

* Removal of sex chromosome CpGs
* Selection of control subjects (current direction - may change in downstream analysis)
* Evaluation of outlier detection strategies

### Longitudinal Modeling

A mixed-effects spline model is fit separately for each CpG to estimate
age-associated methylation trajectories.

Predicted methylation values are generated across a predefined age grid and 
used as trajectory features for downstream analyses.

### Gaussian Mixture Modeling (GMM)

Gaussian mixture models are applied to identify CpGs exhibiting sub-trajectory
patterns (multiple methylation strata).

These CpGs are separated from the primary trajectory clustering workflow and 
may be investigated independently.

### Trajectory Clustering

Remaining CpGs are clustered according to their predicted developmental
methylation profiles.

The objective is to identify common classes of methylation change, including:

* increasing trajectories
* decreasing trajectories
* complex non-linear patterns (early change then plateau)

### Trajectory Annotation

Trajectory classes will be summarized using:

* overall trajectory shape
* rates of methylation change
* inflection points
* genomic annotation
* functional annotation


## Downstream analyses

- exposure association analyses
- hypothesis generation related to T1D risk 

## Repository Structure

| Directory  | Description                                     |
| ---------- | ----------------------------------------------- |
| `R/`       | Analysis functions and modeling code            |
| `scripts/` | Standalone analysis scripts and HPC batch jobs  |
| `results/` | Intermediate and final analysis outputs         |
| `reports/` | Analysis reports and summaries                  |
| `figures/` | Figures and visualizations                      |
| `docs/`    | Project documentation and supporting literature |


## Project Status

Current development efforts focus on:

* refining trajectory representations
* identifying multi-modal CpGs using GMMs
* defining robust criteria for flat trajectories
* evaluating trajectory clustering approaches
* characterizing developmental methylation patterns
