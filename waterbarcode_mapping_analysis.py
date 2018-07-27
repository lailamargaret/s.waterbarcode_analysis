import subprocess as sp
#from Bio import SeqIO

def process(in_fastq1):

	_DB_hg19 = '/data/applications/stamp/pipelines/production/resources/genome/hg19.fa'
	_DB_phix = '/StampFileShare/test2/pipeline_data/heme-stamp/production/runs/HEME0035-water_barcode/phix/phix_genome.fa' 

	list = in_fastq1.rsplit('/', 1)
	directory = list[0]
	list2 = in_fastq1.split('/')
	run_name = list2[-2]
	in_fastq2 = in_fastq1.replace('R1', 'R2')

	hg19_aligned_sam = '%s/hg19_aligned.sam' % directory
	with open (hg19_aligned_sam, 'w') as f:
		sp.call(['bwa' , 'mem', _DB_hg19, in_fastq1, in_fastq2], stdout = f)

	hg19_aligned_bam = '%s/hg19_aligned.bam' % directory
	with open (hg19_aligned_bam, 'w') as f:
		sp.call(['samtools', 'view', '-S', '-b', hg19_aligned_sam], stdout = f)

	hg19_unmapped_bam = '%s/hg19_unmapped.bam' % directory
	with open (hg19_unmapped_bam, 'w') as f:
		sp.call(['samtools', 'view', '-S', '-u', '-f', '12', '-F', '256', hg19_aligned_sam], stdout = f)

	#create 2 fastqs
	fastq1 = '%s/hg19_unmapped_R1.fastq' % directory
	fastq2 = fastq1.replace('R1', 'R2')

	sp.call(['bamToFastq', '-i', hg19_unmapped_bam, '-fq', fastq1, '-fq2', fastq2])

	hg19_unmapped_phix_sam = '%s/hg19_unmapped_phix.sam' % directory
	with open(hg19_unmapped_phix_sam, 'w') as f:
		sp.call(['bwa', 'mem', _DB_phix, fastq1, fastq2], stdout = f)

	unmapped_bam = '%s/unmapped.bam' % directory
	with open (unmapped_bam, 'w') as f:
	        sp.call(['samtools', 'view', '-S', '-u', '-f', '12', '-F', '256', hg19_unmapped_phix_sam], stdout = f)


	fastq1 = '%s/unmapped_R1.fastq' % directory
	fastq2 = fastq1.replace('R1', 'R2')
	sp.call(['bamToFastq', '-i', unmapped_bam, '-fq', fastq1, '-fq2', fastq2])

	fasta1 = '%s/unmapped_R1.fasta' % directory
	fasta2 = fasta1.replace('R1', 'R2')
	temp = '%s/temp.txt' % directory			

	with open(fastq1, 'r') as f1, open(temp, 'w+') as t:
			sp.call(['paste', '-', '-', '-', '-'], stdin = f1, stdout = t)
	with open(fasta1, 'w') as f2, open(temp, 'r') as t:
			lines = t.readlines()
			for line in lines:
				content = line.split('\t')
				angle_bracket  = content[0].replace('@', '>')
				f2.write(angle_bracket + '\n' + content[1] + '\n')
	with open(fastq1, 'r') as f1, open(temp, 'w+') as t:
                        sp.call(['paste', '-', '-', '-', '-'], stdin = f1, stdout = t)
        with open(fasta2, 'w') as f2, open(temp, 'r') as t:
                        lines = t.readlines()
                        for line in lines:
                                content = line.split('\t')
                                angle_bracket  = content[0].replace('@', '>')
                                f2.write(angle_bracket + '\n' + content[1] + '\n')					
	sp.call(["rm", temp])		


	#Each item to show on the txt file 
	r_original_bam = sp.check_output(['samtools', 'flagstat', hg19_aligned_bam]).split(" ")[0]
	r_nm_hg19 = sp.check_output(['samtools', 'flagstat', hg19_unmapped_bam]).split(" ")[0]
	r_m_hg19 = int(r_original_bam) - int(r_nm_hg19)
	r_nm_phix_hg19 = sp.check_output(['samtools', 'flagstat', unmapped_bam]).split(" ")[0]
	r_m_phix = int(r_nm_hg19) - int(r_nm_phix_hg19)
	
	#Write to unmapped_reads.txt which contains information about what mapped
        chart = '%s/unmapped_reads.txt' % directory
	with open (chart, 'w') as ch:
		ch.write("Run: %s\n" % run_name)
		ch.write("Reads in original file: %s\n" % r_original_bam)
		ch.write("Reads mapped to hg19: %s\n" % r_m_hg19)
		ch.write("Remaining unmapped reads: %s\n" % r_nm_hg19)
		ch.write("Reads mapped to phix: %s\n" % r_m_phix)
		ch.write("Remaining unmapped reads: %s (after removing reads that map to either genome)\n" % r_nm_phix_hg19)
	
	#write to a master database file
	waterbarcode_analysis_master = "/home/upload/waterbarcodeanalysis.txt"
	with open (waterbarcode_analysis_master, 'a') as f:
		f.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (run_name, r_original_bam, r_m_hg19, r_nm_hg19, r_m_phix, r_nm_phix_hg19))


def main():
	input = raw_input("Enter the full filepath of the FIRST fastq file: " )
	process(input)      

main()
