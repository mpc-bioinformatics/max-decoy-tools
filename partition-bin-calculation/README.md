# Calculate partition/bin boundaries for the database

This script calculates the (database) bin/partition borders (in Da*1000000) for peptides

As input a TSV file with columns for "peptides", "bin_start" and "bin_stop" is needed, in this order.
This file can be generated from the database entries of an earlier digest. The smaller the bins in
the given TSV file, the better.
