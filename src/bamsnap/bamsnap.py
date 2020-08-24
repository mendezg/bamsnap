import os
import time
import multiprocessing as mp
from PIL import Image, ImageDraw, ImageColor, ImageOps, ImageFont
from pyfaidx import Fasta
from zipfile import ZipFile
from .drawreadset import DrawReadSet, CoveragePlot, CoverageHeatmap
from .coordinates import COORDINATES
from .geneplot import GenePlot
from .basetrack import BaseTrack
from .scale import Xscale
from .bam import BAM
from .util import get_url, getTemplatePath, fileOpen, fileSave, check_dir, getrgb, get_scale, is_exist
from .conf import COLOR, IMAGE_MARGIN_BOTTOM


class BamSnap():

    def __init__(self, opt):
        self.opt = opt
        self.font = {}
        self.bamlist = []
        self.outfnamelist = []
        self.has_opt_error = False
        self.is_single_image_out = False
        self.drawplot = self.opt['draw']
        self.bamplot = self.opt['bamplot']
        self.process = {}
        self.split_poslist = self.get_split_poslist(self.opt['poslist'], self.opt['process'])
        self.set_ref_fasta()
        self.set_is_single_image_out()

    def get_split_poslist(self, poslist, no_process):
        split_poslist = {}
        icore = 0
        for pos1 in poslist:
            try:
                split_poslist[icore].append(pos1)
            except KeyError:
                split_poslist[icore] = [pos1]
            icore += 1
            if icore >= no_process:
                icore = 0
        return split_poslist

    def set_ref_fasta(self):
        self.fasta = None
        if is_exist(self.opt['ref']):
            self.fasta = Fasta(self.opt['ref'], rebuild=False)

    def get_font(self, font_size, font_type='regular'):
        try:
            font = self.font[font_size]
        except KeyError:
            if font_type == 'bold':
                font = ImageFont.truetype(getTemplatePath('VeraMono-Bold.ttf'), font_size)
            else:
                font = ImageFont.truetype(getTemplatePath('VeraMono.ttf'), font_size)
            self.font[font_size] = font
        return font

    def set_is_single_image_out(self):
        flag = False
        if len(self.opt['poslist']) == 1:
            if not self.opt['separated_bam']:
                flag = True
            else:
                fnameuppper = self.opt['out'].upper()
                if fnameuppper.endswith('.PNG') or fnameuppper.endswith('.JPG'):
                    flag = True
        self.is_single_image_out = flag

    def load_bamlist(self):
        ks = self.opt.keys()
        if 'bamlist' in ks and self.opt['bamlist'] is not None:
            for line in open(self.opt['bamlist']):
                line = line.strip()
                if len(line) > 0 and line[0] != '#':
                    arr = line.split('\t')
                    if len(arr) == 2:
                        self.bamlist.append(BAM(arr[0], arr[1]))
                    else:
                        self.bamlist.append(BAM(arr[0]))
        elif 'bam' in ks and len(self.opt['bam']) > 0:
            for idx, bamfile in enumerate(self.opt['bam']):
                try:
                    title = self.opt['title'][idx]
                except IndexError:
                    title = ""
                self.bamlist.append(BAM(bamfile, title))

    def get_refseq(self, pos1):
        refseq = {}
        if self.opt['ref'] == "":
            refseq = self.get_refseq_from_ucsc(pos1)
        else:
            refseq = self.get_refseq_from_localfasta(pos1)
        return refseq

    # def set_refseq_from_ncbiapi(self):
    #     url = "http://www.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide"
    #     url += "&id=NC_000001&seq_start=97533596&seq_stop=97533606&rettype=fasta&retmode=text"
    #     cont = get_url(url)
    #     seq = ""
    #     for line in cont.strip().split('\n'):
    #         if line[0] != '>':
    #             seq += line.strip()
    #     print ('seq:',seq)
    #     i = 0
    #     for gpos in range(self.opt['g_spos']-1000,self.opt['g_epos']+1000+1):
    #         self.refseq[gpos] = seq[i]
    #         i += 1

    def get_refseq_from_ucsc(self, pos1):
        spos = pos1['g_spos']-self.opt['margin'] - 500
        epos = pos1['g_epos']+self.opt['margin'] + 1 + 500
        # seqver = "hg38"
        seqver = self.opt['refversion']
        if not pos1['chrom'].startswith('chr'):
            chrom = 'chr' + pos1['chrom']
        else:
            chrom = pos1['chrom']
        url = "http://genome.ucsc.edu/cgi-bin/das/" + seqver + \
            "/dna?segment=" + chrom + ":" + str(spos) + "," + str(epos)
        # url = "http://genome.ucsc.edu/cgi-bin/das/" + seqver + \
        #     "/dna?segment=chr"+pos1['chrom']+":"+str(spos+1)+","+str(epos+1)
        cont = get_url(url)
        seq = ""
        for line in cont.strip().split('\n'):
            if line[0] != '<':
                seq += line.strip().upper()
        i = 0
        refseq = {}
        for gpos in range(spos, epos):
            refseq[gpos] = seq[i]
            i += 1
        return refseq

    def get_refseq_from_localfasta(self, pos1):
        spos = pos1['g_spos']-self.opt['margin'] - 500
        epos = pos1['g_epos']+self.opt['margin']+1 + 500
        seq = self.get_refseq_from_fasta(pos1['chrom'], spos, epos, self.opt['ref_index_rebuild'])
        i = 0
        refseq = {}
        for gpos in range(spos, epos):
            refseq[gpos+1] = seq[i]
            i += 1
        return refseq

    def get_refseq_from_fasta(self, chrom, spos, epos, rebuild_index=False):
        f = self.fasta
        fastachrommap = {}
        for c1 in list(f.keys()):
            arr = c1.split(' ')
            tchrom = arr[0]
            fastachrommap[tchrom] = c1
        refseq = f[fastachrommap[chrom]][spos:epos+1]
        return str(refseq)

    def add_margin_to_image(self, ia, margin_left=0, margin_top=0, margin_right=0, margin_bottom=0):
        if max(margin_left, margin_top, margin_right, margin_bottom) > 0:
            ia = ImageOps.expand(ia, (margin_left, margin_top, margin_right, margin_bottom), fill=getrgb(self.opt['bgcolor']))
        return ia

    def save_image(self, ia, bam, pos1):
        imgtitle = ""
        if self.is_single_image_out:
            outfname = self.opt['out']
            imgtitle = outfname
            u = outfname.upper()
            if not (u.endswith('.JPG') or u.endswith('.JPEG') or u.endswith('.PNG')):
                outfname += '.' + self.opt['imagetype']
        else:
            if self.opt['save_image_only']:
                path = os.path.join(self.opt['out'])
            else:
                path = os.path.join(self.opt['out'], self.opt['image_dir_name'])

            bamtitle = bam.title.replace(' ', '_').replace('#', '_')
            
            if self.opt['separated_bam']:
                imgtitle = bamtitle + '_'

            if 't_pos' in pos1:
                imgtitle += pos1['chrom'] + ':' + str(pos1['t_pos'])
            else:
                imgtitle += pos1['chrom'] + ':' + str(pos1['t_spos']) + '-' + str(pos1['t_epos'])
            
            outfname = os.path.join(path, imgtitle) + '.' + self.opt['imagetype']

        check_dir(outfname)
        if self.opt['imagetype'] == "jpg":
            ia = ia.convert("RGB")
            ia.save(outfname, "JPEG", quality=100, optimize=True)
        else:
            ia.save(outfname, "PNG", quality=1, optimize=True)

        # self.opt['log'].info("(" + mp.current_process().name + ") Saved " + outfname)
        self.outfnamelist.append({'t': imgtitle, 'img': outfname.split('/')[-1]})
        return outfname

    def save_html(self):
        cont = fileOpen(getTemplatePath('template.html'))
        out = os.path.join(self.opt['out'], 'bamsnap_index.html')
        cont_bamsnapimages = ""
        for imgfile in self.outfnamelist:
            imgpath = os.path.join(self.opt['image_dir_name'], imgfile['img'])
            # print (imgpath)
            cont_bamsnapimages += imgfile['t'] + "<br>"
            cont_bamsnapimages += "<img src='"+imgpath+"'><br>"
        cont = cont.replace('##BAMSNAPIMAGES##', cont_bamsnapimages)
        fileSave(out, cont, 'w')
        self.opt['log'].info('Saved '+out)

    def get_title_image(self, title, w, fontsize):
        margin = 10
        font = self.get_font(fontsize, 'bold')
        fontsize = font.getsize(title)
        h = fontsize[1] + 3 + margin * 2
        im = Image.new('RGB', (w, h), getrgb(self.opt['bgcolor']))
        dr = ImageDraw.Draw(im)

        if self.opt['border']:
            h1 = int(h/2)
            dr.line([(0, h1), (w, h1)], fill=getrgb('000000'), width=1)

        # dr.rectangle([(margin, margin - 1 ),(margin * 3 + fontsize[0], margin + fontsize[1] + 1)], fill=(255,255,255,255), outline=(120,120,120,255))
        dr.rectangle([(margin, margin - 1 ),(margin * 3 + fontsize[0], margin + fontsize[1] + 1)], fill=(255,255,255,255), outline=(255,255,255,255))
        dr.text((margin * 2, margin), title, font=font, fill=COLOR['LABEL'])
        return im, h

    def append_image(self, ia, ia_sub):
        offset_top = ia.height
        ia = ImageOps.expand(ia, (0, 0, 0, ia_sub.height), fill=getrgb(self.opt['bgcolor']))
        ia.paste(ia_sub, (0, offset_top))
        return ia

    def get_image_seperator(self, w, h):
        im = Image.new('RGB', (w, h), getrgb(self.opt['bgcolor']))
        dr = ImageDraw.Draw(im)
        h1 = int(h/2)
        dr.line([(0, h1), (w, h1)], fill=getrgb('000000'), width=1)
        return im

    def init_image(self, image_w, bgcolor="FFFFFF"):
        ia = Image.new('RGB', (image_w, 0), getrgb(bgcolor))
        return ia

    def get_bamplot_image(self, bam, pos1, image_w, xscale, refseq):
        rset = DrawReadSet(bam, pos1['chrom'], pos1['g_spos'], pos1['g_epos'], xscale, refseq)
        rset.read_gap_w = self.opt['read_gap_width']
        rset.read_gap_h = self.opt['read_gap_height']
        rset.read_thickness = self.opt['read_thickness']
        rset.coverage_vaf = self.opt['coverage_vaf']
        rset.opt = self.opt
        rset.xscale = xscale
        rset.calculate_readmap(is_strand_group=True)

        im = self.init_image(image_w, self.opt['bgcolor'])
        border_top = 0
        if not self.opt['no_title']:
            ia_sub, title_height = self.get_title_image(bam.title, image_w, self.opt['title_fontsize'])
            border_top = im.height + int(title_height/2)
            im = self.append_image(im, ia_sub)

        for pidx, plot1 in enumerate(self.bamplot):
            if plot1 == "heatmap":
                covhmplot = CoverageHeatmap(rset, xscale)
                covhmplot.font = self.get_font(self.opt['coverage_fontsize'])
                ia_sub = covhmplot.get_image(image_w, self.opt['heatmap_height'], self.opt['heatmap_bgcolor'])
                im = self.append_image(im, ia_sub)

            if plot1 == "coverage":
                covplot = CoveragePlot(rset, xscale, self.opt['coverage_vaf'])
                covplot.coverage_color = self.opt['coverage_color']
                covplot.font = self.get_font(self.opt['coverage_fontsize'])
                ia_sub = covplot.get_image(image_w, self.opt['coverage_height'], self.opt['coverage_bgcolor'])
                im = self.append_image(im, ia_sub)

            if plot1 == "read":
                if self.opt['read_group'] == "":
                    h_all = rset.get_estimated_height('all')
                    ia_sub = rset.get_image(
                        image_w, h_all, 'all', self.opt['read_color'], self.opt['read_bgcolor'], self.opt['read_color_by'])
                    im = self.append_image(im, ia_sub)

                elif self.opt['read_group'] == "strand":
                    # self.opt['read_bgcolor'] = "F0F000"
                    h_pos = rset.get_estimated_height('pos_strand') 
                    ia_sub = rset.get_image(image_w, h_pos, 'pos_strand',
                                            self.opt['read_pos_color'], self.opt['read_bgcolor'], self.opt['read_color_by'])
                    im = self.append_image(im, ia_sub)

                    # self.opt['read_bgcolor'] = "00F0F0"
                    h_neg = rset.get_estimated_height('neg_strand') 
                    ia_sub = rset.get_image(image_w, h_neg, 'neg_strand',
                                            self.opt['read_neg_color'], self.opt['read_bgcolor'], self.opt['read_color_by'])
                    im = self.append_image(im, ia_sub)

            if plot1 == "coordinates":
                im = self.append_coordinates_image(im, pos1, image_w, xscale)

            if plot1 == "gene":
                im = self.append_geneplot_image(im, pos1, image_w, xscale)

            if plot1 == "base":
                im = self.append_baseplot_image(im, pos1, image_w, xscale, refseq)

        if self.opt['border']:
            padding_bottom = 15
            margin_top = 25
            margin_bottom = 10
            im = self.add_margin_to_image(im, 0,margin_top,0,margin_bottom+padding_bottom)

            dr = ImageDraw.Draw(im)
            dr.line([(0, border_top + margin_top), (0, im.height-1-margin_bottom)], fill=getrgb('000000'), width=1)
            dr.line([(im.width-1, border_top + margin_top), (im.width-1, im.height-1-margin_bottom)], fill=getrgb('000000'), width=1)
            dr.line([(0, im.height-1-margin_bottom), (im.width-1, im.height-1-margin_bottom)], fill=getrgb('000000'), width=1)

        return im

    def append_coordinates_image(self, ia, pos1, image_w, xscale):
        coord = COORDINATES(pos1['chrom'], pos1['g_spos'], pos1['g_epos'], xscale, image_w, self.opt['coordinates_height'], self.opt['debug'])
        coord.font = self.get_font(self.opt['coordinates_fontsize'])
        coord.axisloc = self.opt['coordinates_axisloc']
        coord.bgcolor = self.opt['coordinates_bgcolor']
        coord.labelcolor = self.opt['coordinates_labelcolor']
        ia_sub = coord.get_image()
        ia = self.append_image(ia, ia_sub)
        return ia

    def append_geneplot_image(self, ia, pos1, image_w, xscale):
        geneplot = GenePlot(pos1['chrom'], pos1['g_spos'], pos1['g_epos'], xscale, image_w, self.opt['refversion'], show_transcript=True)
        geneplot.font = self.get_font(self.opt['gene_fontsize'])
        geneplot.gene_pos_color = self.opt['gene_pos_color']
        geneplot.gene_neg_color = self.opt['gene_neg_color']
        ia_sub = geneplot.get_image()
        ia = self.append_image(ia, ia_sub)
        return ia

    def append_baseplot_image(self, ia, pos1, image_w, xscale, refseq):
        baseplot = BaseTrack(pos1['chrom'], pos1['g_spos'], pos1['g_epos'], refseq, xscale, image_w)
        baseplot.font = self.get_font(self.opt['base_fontsize'])
        ia_sub = baseplot.get_image(self.opt['base_margin_top'], self.opt['base_margin_bottom'])
        ia = self.append_image(ia, ia_sub)
        return ia

    def drawplot_bamlist(self, pos1, image_w, bamlist, xscale, refseq):
        ia = self.init_image(image_w, self.opt['bgcolor'])
        drawA = None

        for pidx, plot1 in enumerate(self.drawplot):
            if plot1 == "coordinates":
                ia = self.append_coordinates_image(ia, pos1, image_w, xscale)

            if plot1 == "bamplot":
                for bidx, bam in enumerate(bamlist):
                    if (bidx + 1) % 10 == 0:
                        self.opt['log'].info("..processing " + bam.filename + " (" + str(bidx + 1) + ")")
                    ia_sub = self.get_bamplot_image(bam, pos1, image_w, xscale, refseq)
                    ia = self.append_image(ia, ia_sub)

                    if not self.opt['border'] and self.opt['separator_height'] > 0:
                        ia_sub = self.get_image_seperator(image_w, self.opt['separator_height'])
                        ia = self.append_image(ia, ia_sub)

            if plot1 == "gene":
                ia = self.append_geneplot_image(ia, pos1, image_w, xscale)

            if plot1 == "base":
                ia = self.append_baseplot_image(ia, pos1, image_w, xscale, refseq)

        if self.opt['grid'] > 0:
            drawA = ImageDraw.Draw(ia)
            for posi in range(pos1['g_spos'], pos1['g_epos'] + 1, self.opt['grid']):
                x1 = xscale.xmap[posi]['spos'] - 1
                drawA.line([(x1, 0), (x1, ia.height)], fill=COLOR['GRID_COLOR'], width=1)

        if self.opt['center_line']:
            drawA = ImageDraw.Draw(ia)
            x1 = ((self.opt['g_epos']-self.opt['g_spos'])/2)*self.scale_x
            for h1 in range(0, ia.height, 7):
                drawA.line([(x1, h1), (x1, h1+2)], fill=COLOR['CENTER_LINE'], width=1)

        if not self.opt['no_target_line']:
            drawA = ImageDraw.Draw(ia)
            x1 = xscale.xmap[pos1['t_spos']]['spos'] - 1
            x2 = xscale.xmap[pos1['t_spos']]['epos'] + 1
            for h1 in range(0, ia.height, 7):
                drawA.line([(x1, h1), (x1, h1+2)], fill=COLOR['CENTER_LINE'], width=1)
                drawA.line([(x2, h1), (x2, h1+2)], fill=COLOR['CENTER_LINE'], width=1)

            if self.opt['debug']:
                for xi in range(pos1['g_spos'], pos1['g_epos']+1):
                    x1 = xscale.xmap[xi]['spos']
                    drawA.line([(x1, 20), (x1, ia.height)], fill=(0,0,0,255), width=1)
                drawA.line([(int(image_w/2), 0), (int(image_w/2), ia.height)], fill=(255,0,0,255), width=1)

        if 'highlight' in self.opt.keys():
            drawA = ImageDraw.Draw(ia)
            for h1 in self.opt['highlight'].split(','):
                arr = h1.split(':')
                arr2 = arr[1].split('-')
                t_spos = int(arr2[0])
                t_epos = int(arr2[1])
                x1 = (t_spos-self.opt['g_spos'])*self.scale_x - 4
                x2 = (t_epos-self.opt['g_spos'])*self.scale_x - 2
                for h1 in range(0, ia.height, 7):
                    drawA.line([(x1, h1), (x1, h1+2)], fill=COLOR['CENTER_LINE'], width=1)
                    drawA.line([(x2, h1), (x2, h1+2)], fill=COLOR['CENTER_LINE'], width=1)

        ia = self.add_margin_to_image(ia, self.opt['plot_margin_left'], self.opt['plot_margin_top'], self.opt['plot_margin_right'], self.opt['plot_margin_bottom'])
        imagefname = self.save_image(ia, bamlist[0], pos1)
        return imagefname
    
    def generate_zipfile(self):
        outzip = self.opt['out'] + '.zip'
        zo = ZipFile(outzip, 'w')
        for folderName, subfolders, filenames in os.walk(self.opt['out']):
            for filename in filenames:
                filePath = os.path.join(folderName, filename)
                zo.write(filePath, filePath)
        zo.close()
        self.opt['log'].info("(" + mp.current_process().name + ") Saved " + outzip)
    
    def start_process_drawplot(self, image_w, bamlist):
        for tno in range(self.opt['process']):
            self.process[tno] = mp.Process(target=self.start_process_drawplot_bamlist, args=(image_w, bamlist, self.split_poslist[tno]), name='proc ' + str(tno+1))
            self.process[tno].start()

    def start_process_drawplot_bamlist(self, image_w, bamlist, poslist):
        for pos1 in poslist:
            t11 = time.time()
            refseq = self.get_refseq(pos1)
            xscale = Xscale(pos1['g_spos'], pos1['g_epos'], image_w)
            imagefname = self.drawplot_bamlist(pos1, image_w, bamlist, xscale, refseq)
            t12 = time.time()
            self.opt['log'].info("(" + mp.current_process().name + ") Saved " + imagefname + " : " + str(round(t12-t11, 5)) + " sec")

    def run(self):
        t0 = time.time()
        timemap = {'set_refseq': 0}
        self.load_bamlist()

        if not self.has_opt_error and len(self.opt['poslist']) > 0:
            image_w = self.opt['width'] - self.opt['plot_margin_left'] - self.opt['plot_margin_right']

            if self.opt['separated_bam']:
                for bidx, bam in enumerate(self.bamlist):
                    self.start_process_drawplot(image_w, [bam])
            else:
                self.start_process_drawplot(image_w, self.bamlist)

            for tno in range(self.opt['process']):
                self.process[tno].join()

        t2 = time.time()
        if not self.is_single_image_out:
            if not self.opt['save_image_only']:
                self.save_html()
            if self.opt['zipout']:
                self.generate_zipfile()

        self.opt['log'].debug('Total running time for getting reference sequence (set_refseq): ' +
                              str(round(timemap['set_refseq'], 3))+' sec')
        self.opt['log'].info('Total running time: ' + str(round(t2-t0, 1))+' sec')
