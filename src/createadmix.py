"""

Create synthetic bulk RNA-seq samples from single-cell data
with an idealised mixing (heterogeneity)

Currently supports .fastq.gz

kieran.campbell@dpag.ox.ac.uk

"""

import os, sys
import numpy as np
import gzip

class SyntheticAdmixCreator:
	"""Take single-cell data and create synthetic bulk data 

	NB: have to be very careful that we always access the dictionary in the order
	of cell directories, as order not preserved in python dictionaries """

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
		self.output_file = None

		self._get_ncells()
		if self.paired_end:
			self._find_paired_ends()
			self.output_file = [self.output_fastq.replace(".fastq.gz","") + x + ".fastq.gz" for x in ('_1','_2')]
		else:
			self._trim_dict()
			self.output_file = [self.output_fastq]

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

	def cellIO(self):

		outfilestreams = [gzip.open(f,'wb') for f in self.output_file]

		for i in range(len(self.cell_directories)):
			""" iterate over each type of cell """
			cell_type = self.cell_directories[i]
			reads_per_cell = self.assigned_reads_per_cell[i] # select reads_per_cell from 
			for cell in self.cell_dict[cell_type]:
				self._write_one_cell(cell_type, cell, reads_per_cell, outfilestreams)


		[f.close() for f in outfilestreams]

	def _write_one_cell(self, directory, cell_file, reads_per_cell, outfilestream):
		""" select reads_per_cell reads at random from cell file in directory, and output
		to outfilestream """
		if self.paired_end:
			print "Reading %s" % cell_file
			fs = [gzip.open(os.path.join(directory, f),'rb') for f in [cell_file + x + ".fastq.gz" for x in ("_1","_2")]]
			
			lines = [f.readlines() for f in fs]
			n_lines = [len(l) for l in lines]

			assert n_lines[0] == n_lines[1], "Paired end reads must have equal number of lines in _1 and _2 files"
			n_lines = n_lines[0]
			assert n_lines % 4 == 0, "Number of lines in %s  and %s is not a multiple of 4" % cell_file
			

			""" now select reads_per_cell from n_lines """
			chosen_lines = np.random.choice(n_lines / 4,  reads_per_cell) * 4

			for l in chosen_lines:
				[outfilestream[i].writelines(lines[i][l:(l+4)]) for i in (0,1)]
			lines = None
			[f.close() for f in fs]

		else:
			print  "Non-paired-end currently not supported"
			sys.exit(1)













