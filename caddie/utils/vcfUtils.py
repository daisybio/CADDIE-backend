# import vcf   # VCF 4.1, 4-2
# import vcfpy  # VCF 4.3
# import vap
from sqlalchemy import create_engine


class VCFUtils:

    @staticmethod
    def __read_pyvcf(file):
        """
        :param file:
        :return: vcf reader object
        """
        # return vcf.Reader(fsock=file)

    @staticmethod
    def __annotate(vcf_reader):
        """
        Needs polyphen-2.2.2-whess-2011_12.sqlite
        (ftp://genetics.bwh.harvard.edu/pph2/whess/polyphen-2.2.2-whess-2011_12.sqlite.bz2)

        :param vcf_reader: vcf reader object
        :return annotations: list of tuples with (vcf reader record, vap annotation)
        vcf record has attributes "CHROM", "POS", "REF", "ALT"
        vap annotation has attributes "gene", "protein", "aa_change", "hvar_pred", "hvar_prob", "hdiv_pred", "hdiv_prob"
        """

        engine = create_engine('sqlite:///data/polyphen2/polyphen-2.2.2-whess-2011_12.sqlite', pool_pre_ping=True)
        conn = engine.connect()

        annotations = []
        for record in vcf_reader:
            for alt in record.ALT:
                # annotation = vap.annotate_variant(conn, f'chr{record.CHROM}', record.POS, record.REF, alt)
                annotation = None
                if annotation is not None:
                    annotations.append((record, annotation))
        return annotations

    @staticmethod
    def __filter(annotations, threshold=0.85):
        """
        "While PolyPhen-2 HDIV uses alleles encoding human proteins and their closely related mammalian homologs as
        TN observations, PolyPhen-2 HVAR applies common human nsSNVs as TN observations"
        https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4375422/

        Filter by probability value threshold
        https://ionreporter.thermofisher.com/ionreporter/help/GUID-57A60D00-0654-4F80-A8F9-F6B6A48D0278.html

        :param annotations:
        :param threshold:
        :return annotations_filtered:list of tuples with (vcf reader record, vap annotation)
        vcf record has attributes "CHROM", "POS", "REF", "ALT"
        vap annotation has attributes "gene", "protein", "aa_change", "hvar_pred", "hvar_prob", "hdiv_pred", "hdiv_prob"
        """
        annotations_filtered = []
        for record, annotation in annotations:
            if annotation.hvar_prob > threshold:
                annotations_filtered.append(annotation.gene)
        return annotations_filtered

    @staticmethod
    def filter_file(file, threshold):
        """
        Apply pipeline of VCF utils to load and extract data from vcf file

        :param file:
        :param threshold:
        :return: list of significant gene names in HUGO format (https://www.genenames.org/)
        """
        try:
            file = file.replace(' ', '')
            vcf_reader = VCFUtils.__read_pyvcf(file.splitlines())
            annotations = VCFUtils.__annotate(vcf_reader)
            seeds = VCFUtils.__filter(annotations, threshold)
        except Exception as e:
            print(e)

        # # Open file, this will read in the header
        # print("here")
        #
        # try:
        #     stream = io.StringIO(file)
        #     reader = vcfpy.Reader(stream=stream)
        #     print(reader.header.samples.names)
        #     for record in reader:
        #         if not record.is_snv():
        #             continue
        #         print(record)
        # except Exception as e:
        #     print(e)

        return seeds
