# kmer
The k-mer features of a set of DNA/ Genomic sequences are of great significance on revealing the hidden patterns in that sequence population. Here, *k* is a integer constant and can range from 2 to several dozens, depending on the real application requirements. The k-mer features therefore are widely used in many applications of Bioinformatics, such as for building new prediction methods, etc. 

This piece of Python code is for generating k-mer features for a list of DNA/ genomic sequences. Given a list of *m* DNA sequences, it returns a 2-d array with shape (*m*, *4<sup>k</sup>*) for the 1-hot representation of the kmer features. For a specific *k*, the total number of k-mer features is *4<sup>k</sup>*. For a DNA sequence, the value of each k-mer feature could be the number of occurrences of this k-mer, or its percentage of occurrences compared to all the other different k-mer features.

This code is carefully written to ensure efficiency. It runs sufficiently fast for a set of a very large number of sequences.

How to use this code:

Example:

```python
seq_list = ['ATCGA', 'TCGAC']  # a list of DNA sequences

k = 6  # choose the value for k
obj = kmer_featurization()  # initialize a kmer_featurization object
kmer_features = obj.obtain_kmer_feature_for_a_list_of_sequences(seq_list, write_number_of_occurrences=False)
# If you would like the k-mer features to be the percentage of occurrences (ranging from 0 to 1) as stated above, then leave write_number_of_occurrences as False (the default). If you prefer the features to be the counts for each k-mer occurrence, then set it to True.
```
