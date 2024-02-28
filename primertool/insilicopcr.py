import re
from io import StringIO
from Bio import SeqIO
from urllib.request import urlopen


class InSilicoPCR:
    def __init__(self, forward_primer, reverse_primer, chromosome):
        self.forward_primer = forward_primer
        self.reverse_primer = reverse_primer
        if not chromosome.startswith("chr"):
            raise ValueError("Chromosome must be of the form 'chr1'")
        self.chromosome = chromosome

        # Get the results from UCSC In-Silico PCR
        url = (f"https://genome.ucsc.edu/cgi-bin/hgPcr?"
               f"hgsid=1700687194_KVgKjWQCEx7I4aHWdxjPXjja0ZGc"
               f"&org=Human&db=hg38"
               f"&wp_target=genome"
               f"&wp_f={forward_primer}"
               f"&wp_r={reverse_primer}"
               f"&Submit=submit"
               f"&wp_size=4000"
               f"&wp_perfect=15"
               f"&wp_good=15"
               f"&boolshad.wp_flipReverse=0"
               f"&boolshad.wp_append=0")
        page = urlopen(url)
        html_bytes = page.read()
        html = html_bytes.decode("utf-8")

        # Get the html output from the pcr
        start_idx = html.find("<PRE>") + len("<PRE>")
        end_idx = html.find("</PRE>")
        output_html = html[start_idx:end_idx]

        # remove <A> tags
        pattern = re.compile(r'<A.*?>')
        output_html = pattern.sub('', output_html)
        pattern = re.compile(r'</A>')
        output_html = pattern.sub('', output_html)

        # translate output_html (fasta string) to biopython fasta
        fasta_pcr = SeqIO.parse(StringIO(output_html), "fasta")

        # remove all entries not corresponding to the given chromosome
        self.fasta_pcr = [entry for entry in fasta_pcr if entry.id.startswith(chromosome)]

    def has_unique_primers(self) -> bool:
        return len(self.fasta_pcr) == 1

# test
# in_silico_pcr = InSilicoPCR('CCTGGGCAACAAAGCAAGAC', 'TGCGCTTGTAATGTCAATAGCT', 'chr7')
# print(in_silico_pcr.has_unique_primers())
