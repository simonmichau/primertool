import streamlit as st
from primertool import primertool as pt

st.title("Primertool")

st.sidebar.markdown("# Primertool")

st.markdown('---')

kuerzel = st.text_input(label='', placeholder='Initials', help='Enter your initials (e.g. "XY")')
genome_assembly = st.selectbox('Genome assembly', ('hg38', 'hg19'), index=0, help='Select genome assembly')
generator = st.selectbox('Generate from:', ('Variant', 'Exon', 'Gene', 'Genomic Position'), index=0,
                         help='Select what input you want to use to generate primers')

st.markdown('---')

if generator == 'Variant':
    st.subheader("VariantPrimerGenerator")
    hgvs_variant = st.text_input(label='', placeholder='HGVS variant',
                                 help='Enter variant in HGVS notation (e.g. "NM_000410.3:c.845G>A")')

    if genome_assembly and hgvs_variant:
        if not kuerzel:
            st.warning("Please enter your initials")
            kuerzel = None
        df_ordertable = pt.VariantPrimerGenerator(hgvs_variant, genome_assembly, kuerzel=kuerzel).ordertable
        st.markdown(f"Primers for variant {hgvs_variant} in genome assembly {genome_assembly} are:")
        st.dataframe(df_ordertable.style.hide_index())
elif generator == 'Exon':
    st.subheader("ExonPrimerGenerator")
    nm_number = st.text_input(label='', placeholder='NM number', help='Enter NM number (e.g. "NM_000451")')
    exon_number = st.number_input(label='', min_value=0, step=1, help='Enter exon number')

    if genome_assembly and nm_number and exon_number:
        if not kuerzel:
            st.warning("Please enter your initials")
            kuerzel = None
        df_ordertable = pt.ExonPrimerGenerator(nm_number, exon_number, genome_assembly, kuerzel=kuerzel).ordertable
        st.markdown(f"Primers for exon {exon_number} in gene {nm_number} in genome assembly {genome_assembly} are:")
        st.dataframe(df_ordertable.style.hide_index())
elif generator == 'Gene':
    st.subheader("GenePrimerGenerator")
    nm_number = st.text_input(label='', placeholder='NM number', help='Enter NM number (e.g. "NM_000451")')

    if genome_assembly and nm_number:
        if not kuerzel:
            st.warning("Please enter your initials")
            kuerzel = None
        df_ordertable = pt.GenePrimerGenerator(nm_number, genome_assembly, kuerzel=kuerzel).ordertable
        st.markdown(f"Primers for gene {nm_number} in genome assembly {genome_assembly} are:")
        st.dataframe(df_ordertable.style.hide_index())
elif generator == 'Genomic Position':
    st.subheader("GenomicPositionPrimerGenerator")
    chromosome = st.text_input(label='', placeholder='Chromosome', help='Enter chromosome (e.g. "chr1")')
    start_position = st.number_input(label='', min_value=0, step=1, help='Enter start position')
    end_position = st.number_input(label='', min_value=0, step=1, help='Enter end position')

    if genome_assembly and chromosome and start_position and end_position:
        if not kuerzel:
            st.warning("Please enter your initials")
            kuerzel = None
        df_ordertable = pt.GenomicPositionPrimerGenerator(chromosome, start_position, end_position, genome_assembly, kuerzel=kuerzel).ordertable
        st.markdown(f"Primers for genomic position {chromosome}:{start_position}-{end_position} in genome assembly {genome_assembly} are:")
        st.dataframe(df_ordertable.style.hide_index())

#if button_clicked:
#    if not kuerzel:
#        st.warning("Please enter your initials")
#    if kuerzel and genome_assembly and hgvs_variant:
#        st.markdown(f"Primers for variant {hgvs_variant} in genome assembly {genome_assembly} are:")
#        st.dataframe(df_ordertable.style.hide_index())
