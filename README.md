# The Potential of Controversy Affecting Viewership Numbers in Modern Media: The Try Guys Case Study

A numerical analysis of the impact of public controversy on YouTube viewership, 
using The Try Guys as a case study following the Ned Fulmer scandal in October 2022.

## Overview

In October 2022, The Try Guys, one of YouTube's most prominent creator groups, 
experienced a high-profile public scandal following the removal of founding member 
Ned Fulmer. This project investigates whether and how that event measurably impacted 
their viewership, using a before/after comparative analysis centered on the scandal date.

## Methodology

- **Data collection:** Manually documented view counts for all Try Guys videos 
  published one year before and one year after the October 2022 scandal.
- **Feature engineered:** Calculated average monthly views per video for each month 
  in the two-year window to normalize for upload frequency.
- **Analysis:** Applied regression analysis, logistic differential equation modeling, 
  and interpolation techniques to model audience behavioral dynamics over time.
- **Event study design:** Treated the scandal date as an intervention point and 
  compared pre- and post-period trends to quantify the viewership impact.

## Tools & Libraries

- Python
- NumPy — numerical computation and data manipulation
- Matplotlib — data visualization and trend plotting
- SymPy — symbolic mathematics for differential equation modeling

## Key Findings

Viewership dropped significantly following the scandal. Average monthly views per 
video declined measurably in the post-scandal period compared to the pre-scandal 
baseline, suggesting that public controversy had a lasting negative impact on 
audience engagement beyond the immediate news cycle.

## Files

- `FINAL_CODE.py` — full analysis code
- `501 Final Project.pdf` — full technical report detailing methodology and findings

## Context

This project was completed as part of graduate coursework in Computational Applied 
Mathematics at California State University, Fullerton (Spring and Summer 2024).
