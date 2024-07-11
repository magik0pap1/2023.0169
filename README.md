# 2023.0169
files for Lancia Vidoni paper on sub-quadratic 2-OPT moves for TSP

The repository contains all files necessary to replicate computational experiments
as those reported in the paper "Average case sub-quadratic  exact and heuristic procedures
for the traveling salesman 2-OPT neighborhood", by Giuseppe Lancia and Paolo Vidoni,
INFORMS Journal on Computing, 2024.

The programs are in c. They can be compiled by running

$ makefile KoptLS

Two pdf documents are included that describe:

(i) the suite KoptLS for TSP local search via k-OPT moves for (k=2,3,4).
Only k=2 is relevant to the JOC paper.

(ii) a detailed step-by-step guide on how to use
KoptLS in order to simulate the experiment relative
to all tables and figures reported in the paper.
