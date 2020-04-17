# this script calculates the (database) bin/partition borders (in Da*1000000) for peptides

# As input a TSV file with columns for "peptides", "bin_start" and "bin_stop" is needed, in this order.
# This file can be generated from the database entries of an earlier digest. The smaller the bins in
# the given TSV file, the better.

# adjust nr_bins and the bin_file
nr_bins <- 100
bin_file <- "/mnt/data/uni/MaxDecoy/20200417-bin_calculation/peptides_per_1Da_bin.tsv"
output_file <- "/mnt/data/uni/MaxDecoy/20200417-bin_calculation/peptides_binned.tsv"

bin_data <- read.table(bin_file, header = TRUE, sep = "\t")
bin_data <- as.matrix(bin_data)

nr_all_peptides <- sum(bin_data[,1])
peptides_per_bin <- nr_all_peptides / nr_bins

bin_matrix <- matrix(0, nrow=nr_bins, ncol=3)
binned_peptides <- 0
bin_start <- 0
bin_peptides <- 0
peptides_fill_bin <- peptides_per_bin
c_bin <- 1
for (i in 1:nrow(bin_data)) {
    bin_peptides <- bin_peptides + bin_data[i,1]
    binned_peptides <- binned_peptides + bin_data[i,1]

    if (binned_peptides >= peptides_fill_bin) {
        peptides_fill_bin <- peptides_fill_bin + peptides_per_bin
        
        bin_end <- bin_data[i-1,3]
        bin_matrix[c_bin,] <- c(bin_start, bin_data[i-1,3], bin_peptides)
        c_bin <- c_bin + 1
        bin_peptides <- 0
        bin_start <- bin_data[i,2]
    }
}

bin_matrix[c_bin,] <- c(bin_start, 999999999999, bin_peptides)
colnames(bin_matrix) <- c("bin_start", "bin_end", "peptides_in_bin")


if (sum(bin_data[,1]) == sum(bin_matrix[,3])) {
    bin_matrix
    write.table(format(bin_matrix, scientific = FALSE), file = output_file, append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
    print("ok")
} else {
    print("wrong sums")
}
