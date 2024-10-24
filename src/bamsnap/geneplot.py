from PIL import ImageFont, Image, ImageDraw, ImageOps
from .conf import COLOR, GENE_ANNOT_FILE, REFER_SEQ_VERSION
from .util import getTemplatePath, getDataPath, comma, gzopen, decodeb, convert_int_list, getrgb, get_scale
import tabix

LINETYPE = [
    'exon',
    'CDS',
    'start_codon',
    'stop_codon',
    'five_prime_utr',
    'three_prime_utr',
    'Selenocysteine'
]

class GeneAnnot():
    def __init__(self, rec, header):
        self.data = {}
        for i in range(len(header)):
            self.data[header[i]] = rec[i] if i < len(rec) else ''

        self.chrom = self.data.get('CHROM')
        self.spos = int(self.data.get('SPOS')) + 1  # Convert 0-based to 1-based
        self.epos = int(self.data.get('EPOS'))
        self.strand = self.data.get('STRAND', '+')  # Default to '+'
        self.is_negative = self.strand == '-'
        self.gene_id = self.data.get('NAME', 'Gene')
        self.gene_name = self.gene_id
        self.gene_biotype = 'gene'
        self.transcripts = []
        self.visible_transcripts = []
        self.load_transcript()
        self.level = None

    def load_transcript(self):
        # Treat the gene as a single transcript for BED files
        t1 = TranscriptAnnot(self.gene_id, self.gene_biotype, self.spos, self.epos)
        t1.subregion['exon_spos'] = [self.spos]
        t1.subregion['exon_epos'] = [self.epos]
        self.transcripts.append(t1)

    def set_visible_transcript(self, region_spos, region_epos):
        # Determine if the transcript overlaps with the region of interest
        self.visible_transcripts = []
        for transcript in self.transcripts:
            if transcript.epos >= region_spos and transcript.spos <= region_epos:
                self.visible_transcripts.append(transcript)


class TranscriptAnnot():
    def __init__(self, tid, biotype, spos, epos):
        self.transcript_id = tid
        self.biotype = biotype
        self.subregion = {}
        self.spos = int(spos)
        self.epos = int(epos)

    def set_subregion(self):
        # For BED data, subregions are already set
        pass


class GenePlot():
    def __init__(self, chrom, spos, epos, xscale, w, refversion="hg38", show_transcript=True, gene_annot_file=None, opt={}):
        self.chrom = chrom
        self.nchrom = chrom.replace('chr', '')
        self.spos = max(spos, 1)  # Ensure spos is at least 1
        self.epos = max(epos, self.spos)  # Ensure epos >= spos
        self.g_len = self.epos - self.spos + 1
        self.fontsize = opt.get('gene_fontsize', 10)
        self.font = ImageFont.truetype(getTemplatePath('VeraMono.ttf'), self.fontsize)
        self.show_transcript = show_transcript
        self.w = w
        self.h = 0
        self.im = None
        self.bgcolor = "FFFFFF"
        self.noline = 0
        self.margin = opt.get('gene_margin', 5)
        self.lineheight = opt.get('gene_lineheight', 2)  # Adjust as needed
        self.gene_pos_color = opt.get('gene_pos_color', "ffac9c")
        self.gene_neg_color = opt.get('gene_neg_color', "A19Cff")
        self.xscale = xscale
        self.opt = opt  # Store options

        # Use the provided gene annotation file if specified
        if gene_annot_file:
            self.gene_annot_file = gene_annot_file
        else:
            self.gene_annot_file = getDataPath(GENE_ANNOT_FILE.replace("#REFSEQVERSION#", REFER_SEQ_VERSION[refversion]))

        self.gene_annot_tb = tabix.open(self.gene_annot_file)
        self.gene_annot_header = []
        self.gene_annot = []

        self.load_gene_structure()

    def assign_levels(self):
        levels = []  # List of lists, where each sublist represents a level with genes

        # Sort genes by start position
        sorted_genes = sorted(self.gene_annot, key=lambda ga: ga.spos)

        for ga in sorted_genes:
            placed = False
            for level_index, level in enumerate(levels):
                # Check if the gene overlaps with any gene in the current level
                overlap = False
                for other_ga in level:
                    if self.genes_overlap(ga, other_ga):
                        overlap = True
                        break
                if not overlap:
                    # Place gene in this level
                    level.append(ga)
                    ga.level = level_index
                    placed = True
                    break
            if not placed:
                # Create a new level
                levels.append([ga])
                ga.level = len(levels) - 1

        self.total_levels = len(levels)
        self.noline = self.total_levels

    def genes_overlap(self, ga1, ga2):
        # Returns True if ga1 and ga2 overlap in genomic positions
        return not (ga1.epos <= ga2.spos or ga2.epos <= ga1.spos)

    def load_gene_structure(self):
        if not self.gene_annot_header:
            for line in gzopen(self.gene_annot_file):
                line = decodeb(line)
                if line.startswith("#"):
                    self.gene_annot_header = line[1:].strip().split('\t')
                    break
            else:
                # No header found; set default BED header
                self.gene_annot_header = ['CHROM', 'SPOS', 'EPOS', 'NAME']

        # Adjust positions to be within valid ranges
        spos = max(self.spos, 1)
        epos = max(self.epos, spos)  # Ensure epos >= spos

        pos_str = f"{self.chrom}:{spos}-{epos}"
        self.gene_annot = []
        try:
            genes_found = 0
            for rec in self.gene_annot_tb.querys(pos_str):
                genes_found += 1
                ga = GeneAnnot(rec, self.gene_annot_header)
                ga.set_visible_transcript(self.spos, self.epos)
                self.gene_annot.append(ga)
            print(f"Number of genes found in region {pos_str}: {genes_found}")
        except tabix.TabixError:
            # Handle cases where no annotations are found for the region
            print(f"No genes found in region {pos_str}")

        # Assign levels after loading all genes
        self.assign_levels()

    def draw(self, dr):
        # Calculate per-level height
        per_level_height = self.margin + (self.font.getbbox('C')[3] - self.font.getbbox('C')[1]) + self.margin + self.lineheight + self.margin

        for ga in self.gene_annot:
            level = ga.level
            yi = level * per_level_height

            yi += self.margin  # Add top margin within the level

            for t1 in ga.visible_transcripts:
                # Construct the text label conditionally
                if t1.transcript_id != ga.gene_name:
                    txt = f"{ga.gene_name} ({t1.transcript_id})"
                else:
                    txt = ga.gene_name

                # Draw the text label
                x_center = self.xscale.get_x((t1.spos + t1.epos) / 2)['cpos']
                text_width = self.font.getlength(txt)
                x_text = min(max(x_center - text_width / 2, 0), self.w - text_width)
                dr.text((x_text, yi), txt, font=self.font, fill=COLOR['COORDINATE'])

                yi += (self.font.getbbox('C')[3] - self.font.getbbox('C')[1]) + self.margin  # Move yi below the text label

                # Draw the gene line below the text
                x1_line = max(self.xscale.get_x(t1.spos)['spos'], 0)
                x2_line = min(self.xscale.get_x(t1.epos)['epos'], self.w)
                col1 = getrgb(self.gene_neg_color) if ga.is_negative else getrgb(self.gene_pos_color)
                dr.line([(x1_line, yi), (x2_line, yi)], fill=col1, width=self.lineheight)

                # Optionally, draw arrows to indicate strand direction
                if self.show_transcript:
                    arrow_size = 5
                    if ga.is_negative:
                        dr.polygon([
                            (x2_line, yi),
                            (x2_line - arrow_size, yi - arrow_size),
                            (x2_line - arrow_size, yi + arrow_size)
                        ], fill=col1)
                    else:
                        dr.polygon([
                            (x1_line, yi),
                            (x1_line + arrow_size, yi - arrow_size),
                            (x1_line + arrow_size, yi + arrow_size)
                        ], fill=col1)


    def get_image(self):
        if self.im is None:
            # Calculate the image height based on the number of levels
            text_bbox = self.font.getbbox('C')
            text_height = text_bbox[3] - text_bbox[1]
            per_level_height = self.margin + text_height + self.margin + self.lineheight + self.margin
            self.h = int(self.total_levels * per_level_height)

            self.im = Image.new('RGBA', (self.w, self.h), getrgb(self.bgcolor))
            dr = ImageDraw.Draw(self.im)
            self.draw(dr)
        return self.im
