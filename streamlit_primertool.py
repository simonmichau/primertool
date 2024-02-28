import streamlit as st
from primertool import primertool as pt

st.set_page_config(page_title="Primertool", page_icon="logo.png")  # Dna icons created by Freepik - Flaticon

st.title("Primertool")
st.markdown(f"This tool is built to facilitate the process of finding primers for non-technical users. "
            f"Primers can be found for a given variant, exon, entire gene or genomic position.")
st.markdown("To get started, please")
st.markdown("- enter your initials (Kürzel); Default is None")
st.markdown(f"- select the genome assembly you are working with (default is hg38)")
st.markdown(f"- select the primer generation method that corresponds to the type of primer you are looking for and "
            f"enter the required information")
st.markdown(f'- hit the enter button and wait until the order table is generated (depending on your input '
            f'this may take some time)')
st.markdown(f"Note: When running the tool for the first time, it may take several minutes to download the selected "
            f"reference genome. Subsequent runs will be faster.")

st.sidebar.markdown("# Primertool")

st.markdown('---')

kuerzel = st.text_input(label='', placeholder='Initials', help='Enter your initials (e.g. "XY")')
genome_assembly = st.selectbox('Genome assembly', ('hg38', 'hg19'), index=0, help='Select genome assembly')
generator = st.selectbox('Generate from:', ('Variant', 'Exon', 'Gene', 'Genomic Position'),
                         help='Select what input you want to use to generate primers')

st.markdown('---')

if generator == 'Variant':
    st.subheader("VariantPrimerGenerator")
    hgvs_variant = st.text_input(label='', placeholder='HGVS variant',
                                 help='Enter variant in HGVS notation (e.g. "NM_000410.3:c.845G>A")')

    if genome_assembly and hgvs_variant:
        if not kuerzel:
            st.warning("Don't forget to enter your initials")
            kuerzel = None
        with st.spinner('Generating primers...'):
            df_ordertable = pt.VariantPrimerGenerator(hgvs_variant, genome_assembly, kuerzel=kuerzel).ordertable
        st.success(f"Primers for variant {hgvs_variant} in genome assembly {genome_assembly} are:")
        st.dataframe(df_ordertable.style.hide_index())
elif generator == 'Exon':
    st.subheader("ExonPrimerGenerator")
    nm_number = st.text_input(label='', placeholder='NM number', help='Enter NM number (e.g. "NM_000451")')
    exon_number = st.number_input(label='', min_value=0, step=1, help='Enter exon number')

    if genome_assembly and nm_number and exon_number:
        if not kuerzel:
            st.warning("Don't forget to enter your initials")
            kuerzel = None
        with st.spinner('Generating primers...'):
            df_ordertable = pt.ExonPrimerGenerator(nm_number, exon_number, genome_assembly, kuerzel=kuerzel).ordertable
        st.success(f"Primers for exon {exon_number} in gene {nm_number} in genome assembly {genome_assembly} are:")
        st.dataframe(df_ordertable.style.hide_index())
elif generator == 'Gene':
    st.subheader("GenePrimerGenerator")
    nm_number = st.text_input(label='', placeholder='NM number', help='Enter NM number (e.g. "NM_000451")')

    if genome_assembly and nm_number:
        if not kuerzel:
            st.warning("Don't forget to enter your initials")
            kuerzel = None
        with st.spinner('Generating primers...'):
            df_ordertable = pt.GenePrimerGenerator(nm_number, genome_assembly, kuerzel=kuerzel).ordertable
        st.success(f"Primers for gene {nm_number} in genome assembly {genome_assembly} are:")
        st.dataframe(df_ordertable.style.hide_index())
elif generator == 'Genomic Position':
    st.subheader("GenomicPositionPrimerGenerator")
    chromosome = st.text_input(label='', placeholder='Chromosome', help='Enter chromosome (e.g. "chr1")')
    start_position = st.number_input(label='', min_value=0, step=1, help='Enter start position')
    end_position = st.number_input(label='', min_value=0, step=1, help='Enter end position')

    if genome_assembly and chromosome and start_position and end_position:
        if not kuerzel:
            st.warning("Don't forget to enter your initials")
            kuerzel = None
        with st.spinner('Generating primers...'):
            df_ordertable = pt.GenomicPositionPrimerGenerator(chromosome, start_position, end_position, genome_assembly, kuerzel=kuerzel).ordertable
        st.success(f"Primers for genomic position {chromosome}:{start_position}-{end_position} in genome assembly {genome_assembly} are:")
        st.dataframe(df_ordertable.style.hide_index())


# st.balloons()
# st.snow()

# if button_clicked:
#    if not kuerzel:
#        st.warning("Please enter your initials")
#    if kuerzel and genome_assembly and hgvs_variant:
#        st.markdown(f"Primers for variant {hgvs_variant} in genome assembly {genome_assembly} are:")
#        st.dataframe(df_ordertable.style.hide_index())
