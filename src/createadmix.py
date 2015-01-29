"""

Create synthetic bulk RNA-seq samples from single-cell data
with an idealised mixing (heterogeneity)

Currently supports .fastq.gz

kieran.campbell@dpag.ox.ac.uk

"""

import os
import numpy as np

class SyntheticAdmixCreator:
	"""Take single-cell data and create synthetic bulk data """

	def __init__(self, cell_directories, mixing_coefficients, 
		total_reads, output_fastq, paired_end = True):
		"""
		Create new SyntheticAdmixCreator

		@param cell_directories: List of directories holding fastq files, where each
		directory corresponds to one 'type' of file
		@param mixing_coefficients: A numpy array where each entry is the proportion of the
		bulk data made up of cells of the corresponding type given by the same index
		in cell_directories
		@total_reads: Total number of reads to dump to output fastq file
		@output_fastq: The output fastq file
		@paired_end: Are the input files paired end? If so, matching files xxx_1.fastq.gz and
		xxx_2.fastq.gz will be searched for in each 'cell_directories'
		"""
		self._check_inputs(cell_directories, mixing_coefficients)

		# declare all class variables
		self.cell_directories = cell_directories
		self.mixing_coefficients = mixing_coefficients
		self.total_reads = total_reads
		self.output_fastq = output_fastq
		self.paired_end = paired_end

		self.assigned_reads = mixing_coefficients * total_reads # number of reads for each cell type
		self.assigned_reads_per_cell = None

		""" Dictionary holding list of fastq files for each cell type """
		self.cell_fastq_dict = {}
		self.cell_dict = {}
		self.cells_per_type = None


		self._get_ncells()
		if self.paired_end:
			self._find_paired_ends()
		else:
			self._trim_dict()

	def _check_inputs(self, cell_directories, mixing_coefficients):

		for d in cell_directories:
			assert os.path.exists(d), "Directory %s does not exit" %d

		assert mixing_coefficients.sum() == 1, "Mixing coefficients must sum to 1"

		return		

	def _get_ncells(self):
		""" Crawls each directory counting the number of fastq files to work 
		out how many reads per cell we should take (at random) """

		for directory in self.cell_directories:
			self.cell_fastq_dict[directory] = [f for f in os.listdir(directory) if f.endswith(".fastq.gz")]

		self.cells_per_type = np.array([len(self.cell_fastq_dict[d]) for d in self.cell_directories])
		if self.paired_end:
			assert np.any(self.cells_per_type % 2 == 0), "Paired end reads must have even number of fastq.gz files"
			self.cells_per_type = self.cells_per_type / 2			

		self.assigned_reads_per_cell = self.assigned_reads / self.cells_per_type

	def  _find_paired_ends(self):
		""" Iterates over self.cell_fastq_dict and sorts files into pairs based on _1 and _2 """

		for celltype, fastq_list in self.cell_fastq_dict.iteritems():
			forward_strand = [x for x in fastq_list if "_1" in x]
			reverse_strand = [x for x in fastq_list if "_2" in x]

			forward_strand_names = sorted([x.replace("_1.fastq.gz","") for x in forward_strand])
			reverse_strand_names = sorted([x.replace("_2.fastq.gz","") for x in reverse_strand])

			assert(forward_strand_names == reverse_strand_names), "Paired end reads must have matching names"
			self.cell_dict[celltype] = forward_strand_names			


	def _trim_dict(self):
		""" Creates self.cell_dict from self.cell_fastq_dict by trimming .fastq.gz off end """
		for celltype, fastq_list in self.cell_fastq_dict.iteritems():
			self.cell_dict[celltype] = [x.replace("fastq.gz","") for x in fastq_list]


