# 2023.0169

This archive is distributed in association with the INFORMS Journal on Computing under the MIT License.

The software and data in this repository are a snapshot of the software and data that were used in the research reported on in the paper  "Average case sub-quadratic  exact and heuristic procedures for the traveling salesman 2-OPT neighborhood" by G. Lancia and P. Vidoni. 

Cite

To cite the contents of this repository, please cite both the paper and this repo, using their respective DOIs.

https://doi.org/10.1287/ijoc.2019.0169

https://doi.org/10.1287/ijoc.2019.0169.cd

Description

The repository contains all files necessary to replicate computational experiments
as those reported in the paper "Average case sub-quadratic  exact and heuristic procedures
for the traveling salesman 2-OPT neighborhood", by Giuseppe Lancia and Paolo Vidoni,
INFORMS Journal on Computing, 2024.

Building

The programs are in c. They can be compiled by running
$ makefile KoptLS
This creates the suite KoptLS, i.e., a program for TSP local search via k-OPT moves for (k=2,3,4). Only k=2 is relevant to the paper.
Please, find further instructions on how to use the suite KoptLS on the file KoptLSmanual.pdf.

Results

Please, find all instructions on how to replicate the paper's computational experiments on the file KoptLS_JOC.pdf.
This is a detailed step-by-step guide on how to use KoptLS in order to simulate the experiment relative to all tables and figures reported in the paper.

Tables and figures can be found on tablesAndFigures.pdf
