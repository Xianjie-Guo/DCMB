DCMB - Dual Correction strategy based MB learning software

This software implements the dual-correction-strategy-based MB learning (causal feature selection) algorithm, DCMB, as described in [1].


Compilation
===========
./make


Usage
=====
./executable_file_name data_file_path label_index

Example
========
find causal features of label attribute (or find Markov blanket for variable 37):
./main data/Alarm_s500.txt 37


The input data file should contain one data sample per row, each sample consisting of one integer value for each variable. Values on a row should be separated by whitespace (tabs or spaces). For an example data file with 37 binary variables (columns) and 500 samples (rows) see "Alarm_s500.txt". label_index is used to specify which column in the data matrix is the label (target) attribute. "label_index=37" denotes that the last column in the data matrix ("Alarm_s500.txt") is label attribute.


The default output files are "mb.out" (running time) and "time.out" (Markov blanket (causal features) of label).

References
==========
[1] Xianjie Guo, Kui Yu, Lin Liu, Fuyuan Cao, and Jiuyong Li. Causal Feature Selection With Dual Correction[J]. IEEE Transactions on Neural Networks and Learning Systems, 2022.

