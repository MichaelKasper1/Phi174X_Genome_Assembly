Python 3.8.10 (v3.8.10:3d8993a744, May  3 2021, 08:55:58) 
[Clang 6.0 (clang-600.0.57)] on darwin
Type "help", "copyright", "credits" or "license()" for more information.
>>> # python 3


class AssembleSimpleGenome:
    """
    Genome Assembler for error free, constant length reads
    """
    # Initializing data
    def __init__(self, data):
        self.data = data
        self.n = len(self.data)

    # Create assembly using simple overlaps
    def overlap_assembly(self):
        next_s = [(0, None) for _ in range(self.n)]

        for i in range(self.n - 1):
            s1 = self.data[i]
            
            for j in range(i + 1, self.n):
                s2 = self.data[j]
                prev_overlap_size_i, _ = next_s[i]
                prev_overlap_size_j, _ = next_s[j]
                overlap_size_i = self.find_bigger_overlap_size(s1, s2, prev_overlap_size_i)
                overlap_size_j = self.find_bigger_overlap_size(s2, s1, prev_overlap_size_j)
                if overlap_size_i is not None:
                    next_s[i] = (overlap_size_i, j)
                if overlap_size_j is not None:
                    next_s[j] = (overlap_size_j, i)
                    
        str_parts = [self.data[0]]
        overlap_size, cur = next_s[0]
        while cur != 0:
            str_parts.append(self.data[cur][overlap_size:])
            overlap_size, cur = next_s[cur]
        s = "".join(str_parts)

        # Slice when overlaps are identical or circularly identical
        s = s[:-overlap_size]
        return s
    
    #length of reads for this assignment was fixed at 100
    READ_LEN = 5
    
    #staticmethod is used 
    @staticmethod
    def find_bigger_overlap_size(s1, s2, prev_overlap_size):
        res = None
        for overlap_size in range(AssembleSimpleGenome.READ_LEN - 1, prev_overlap_size, -1):
            if s2.startswith(s1[(AssembleSimpleGenome.READ_LEN - overlap_size):]):
                res = overlap_size
                break
        return res

def assemble():
    #number of rows from specific assignment
    n_rows = 1618
    data = [input().strip()]
    for _ in range(n_rows - 1):
        s = input().strip()
        if s != data[-1]:
            data.append(s)

    t = AssembleSimpleGenome(data)
    ans = t.solve()
    print(ans)

# Automatically run
if __name__ == "__main__":
    assemble()
    
# Test will output TRUE if the test genome assembler is able to assemble this test genome.
def run_test():
    
    #small 'genome'
    s = "gagttttatcgcttcca"
    
    #the MER_LEN variable needs to be changed to 5 from 100 before implementing this test.
    #data is random 5mers from within the test string
    data = [
        "cagag",  # 0
        "tatcg",  # 1
        "tttat",  # 2
        "gtttt",  # 3
        "tcgct",  # 4
        "gcttc",  # 5
        "ttcca",  # 6
        "agttt",  # 7
        "gagtt",  # 8
    ]
    t = AssembleSimpleGenome(data)
    genome = t.overlap_assembly()
    print(sorted(genome) == sorted(s)

    
#In further projects I will make genome assemblers using:
# De Bruijn graphs
# Bubble detection algorithms to detect sequencing errors
# Data with non uniform read lengths