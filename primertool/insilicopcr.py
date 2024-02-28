import re
from io import StringIO
from Bio import SeqIO
from urllib.request import urlopen


class InSilicoPCR:
    """
    This class provides a python interface to the UCSC In-Silico PCR tool. It takes a forward and a reverse primer
    and returns the PCR product.
    """

    def __init__(self, forward_primer: str, reverse_primer: str, max_product_size: int = 4000,
                 min_perfect_match: int = 15, min_good_match: int = 15, flip_reverse_primer: bool = False):
        """
        :param forward_primer: forward primer sequence (required)
        :param reverse_primer: reverse primer sequence (required)
        :param max_product_size: maximum product size (optional, default=4000)
        :param min_perfect_match: minimum perfect match (optional, default=15)
        :param min_good_match: minimum good match (optional, default=15)
        :param flip_reverse_primer: flip reverse primer (optional, default=False)
        """
        self.forward_primer = forward_primer
        self.reverse_primer = reverse_primer

        # Get the results from UCSC In-Silico PCR (https://genome.ucsc.edu/cgi-bin/hgPcr)
        url = (f"https://genome.ucsc.edu/cgi-bin/hgPcr?"
               f"hgsid=1700687194_KVgKjWQCEx7I4aHWdxjPXjja0ZGc"
               f"&org=Human"
               f"&db=hg38"
               f"&wp_target=genome"
               f"&wp_f={forward_primer}"
               f"&wp_r={reverse_primer}"
               f"&Submit=submit"
               f"&wp_size={max_product_size}"
               f"&wp_perfect={min_perfect_match}"
               f"&wp_good={min_good_match}"
               f"&boolshad.wp_flipReverse={int(flip_reverse_primer)}"
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
        self.fasta_pcr = list(SeqIO.parse(StringIO(output_html), "fasta"))

        # remove all entries not corresponding to the given chromosome
        # self.fasta_pcr = [entry for entry in fasta_pcr if entry.id.startswith(chromosome)]

    def is_uniquely_binding(self) -> bool:
        """
        Check if the PCR product is uniquely binding to just one site.
        """
        return len(self.fasta_pcr) == 1

# test
# in_silico_pcr = InSilicoPCR('CCTGGGCAACAAAGCAAGAC', 'TGCGCTTGTAATGTCAATAGCT')
# print(in_silico_pcr.is_uniquely_binding())
