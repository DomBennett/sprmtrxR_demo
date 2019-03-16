
# TODO:
# - fix error with different numbers of bps
# - fix error with other groups
# - ensure valid tip names at start of process

# outsider::module_install('dombennett/om..raxml')
raxml <- outsider::module_import('raxml', 'dombennett/om..raxml')

input_file <- 'n304.fasta'
partition_file <- 'partition.txt'
raxml('-f', 'a', '-m', 'GTRGAMMA', '-T', '1', '-#', '10', '-x', sample(1:10000000, 1),
      '-n', 'primates', '-s', input_file, ' -q ', partition_file, '-p', '100')
