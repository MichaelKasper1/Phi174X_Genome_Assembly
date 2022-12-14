# Phix174_Genome_Assembly

I am currently working locally on this project. I have uploaded a simple genome assembler but will be adding to it soon.

At Washington University school of medicine, I was working on a project using data from Illumina DNA sequencers. Since read allignment is quite complicated, we chose to use external softwares to allign our reads. We used Bowtie. I became interested in how genome alignment worked and found a data structures and algorithms course on coursera whose capstone was a genome assembly project!

The aim project is to recreate the genome of a specific E Coli outbreak from Germany. My goals from the project were to implement knowledge of algorithms and to continue building my understanding of working with large genetic datasets. I follow the guidelines for the project from the course.

In the work, I select an optimal kmer size, write code to find sequencing error via bubble detection, and assemble the Phi174X genome.

In April 2011, hundreds of people in Germany were hospitalized with a deadly disease that often started as food poisoning with bloody diarrhea. It was the beginning of the deadliest outbreak in recent history, caused by a mysterious bacterial strain that we will refer to as E. coli X. Within a few months, the outbreak had infected thousands and killed 53 people. To prevent the further spread of the outbreak, computational biologists all over the world had to answer the question “What is the genome sequence of E. coli X?” in order to figure out what new genes it acquired to become pathogenic. The 2011 German outbreak represented an early example of epidemiologists collaborating with computational biologists to stop an outbreak. In this Genome Assembly Programming Challenge, you will follow in the footsteps of the bioinformaticians investigating the outbreak by developing a program to assemble the genome of the deadly E. coli X strain. However, before you embark on building a program for assembling the E. coli X strain, we have to explain some genomic concepts and warm you up by having you solve a simpler problem of assembling a small virus.
