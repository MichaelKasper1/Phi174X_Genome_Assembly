#python3

class Overlap_Assembly_Error_Free:
    """
    Genome Assembler for error free, constant length reads
    by overlaps and a hamilton
    
    This genome assembler has problems when: 
    - the reads are variable lengths
    - there are possible errors in reads
    - the genomes are really large
    
    In further projects I will make genome assemblers using:
    De Bruijn graphs
    Bubble detection algorithms to detect sequencing errors
    Data with non uniform read lengths
    """
    # Initializing data
    def __init__(self, reads):
        self.reads = reads
        self.n = len(self.reads)

    # Create assembly using simple overlaps
    def overlap_assembly(self):
        next_s = [(0, None) for _ in range(self.n)]
        
        #iterate through mer length
        for i in range(self.n - 1):
            mer_a = self.reads[i]
            
            #os means overlap size
            for j in range(i + 1, self.n):
                mer_b = self.reads[j]
                previous_os_i, _ = next_s[i]
                previous_os_j, _ = next_s[j]
                os_i = self.find_bigger_os(mer_a, mer_b, prev_os_i)
                os_j = self.find_bigger_os(mer_b, mer_a, prev_os_j)
                if overlap_size_i is not None:
                    next_s[i] = (os_i, j)
                if overlap_size_j is not None:
                    next_s[j] = (os_j, i)
                    
        str_parts = [self.data[0]]
        os, current = next_s[0]
        
        #if current overlap is 0, append them together
        while current != 0:
            str_parts.append(self.reads[cur][os:])
            os, cur = next_s[cur]
        genome = "".join(str_parts)

        # Slice when overlaps are identical or circularly identical
        genome = genome[:-os]
        return genome
    
    #length of reads for this assignment was fixed at 100
    Read_len = 5
    
    #find which two strings have the bigger overlap
    @staticmethod
    def find_bigger_overlap_size(mer_a, mer_b, prev_os):
        res = ""
        for overlap_size in range(Overlap_Assembly_Error_Free.Read_len - 1, prev_os, -1):
            if mer_b.startswith(mer_a[(Overlap_Assembly_Error_Free.Read_len - os):]):
                res = overlap_size
                break
        return res

def assemble():
    #number of rows from specific assignment
    n_rows = 1618
    data = [input().strip()]
    #data preprocessing
    for _ in range(n_rows - 1):
        s = input().strip()
        if s != reads[-1]:
            reads.append(s)

    tmp = AssembleSimpleGenome(reads)
    ans = tmp.solve()
    print(ans)
