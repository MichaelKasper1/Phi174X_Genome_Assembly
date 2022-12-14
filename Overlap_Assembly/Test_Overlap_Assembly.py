#python3

# Test will output TRUE if the test genome assembler is able to assemble this test genome.
def run_test():
    
    #small 'genome'
    test_genome = "gagttttatcgcttcca"
    
    #the MER_LEN variable needs to be changed to 5 from 100 before implementing this test.
    #data is random 5mers from within the test string
    reads = [
        "cagag", "tatcg", "tttat", "gtttt", "tcgct",
        "gcttc", "ttcca", "agttt", "gagtt",
    ]
    # pass to class
    t = Overlap_Assembly_Error_Free(reads)
    genome = t.overlap_assembly()
    print(sorted(genome) == sorted(test_genome)
