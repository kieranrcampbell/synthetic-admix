"""

Create synthetic bulk RNA-seq samples from single-cell data
with an idealised mixing (heterogeneity)

kieran.campbell@dpag.ox.ac.uk

"""

import os
import numpy as np

class SyntheticAdmixCreator:
	"""Take single-cell data and create synthetic bulk data """

	def __init__(self, cell_directories, mixing_coefficients, total_reads, output_fastq):
		"""
		Create new SyntheticAdmixCreator

		@param cell_directories: List of directories holding fastq files, where each
		directory corresponds to one 'type' of file
		@param mixing_coefficients: A numpy array where each entry is the proportion of the
		bulk data made up of cells of the corresponding type given by the same index
		in cell_directories
		@total_reads: Total number of reads to dump to output fastq file
		@output_fastq: The output fastq file
		"""
		_check_inputs(cell_directories, mixing_coefficients)

		# declare all class variables
		self.cell_directories = cell_directories
		self.mixing_coefficients = mixing_coefficients
		self.total_reads = total_reads
		self.output_fastq = output_fastq

		self.assigned_reads = mixing_coefficients * total_reads # number of reads for each cell type
		self.assigned_reads_per_cell = None

		""" Dictionary holding list of fastq files for each cell type """
		self.cell_fastq_dict = {}
		self.cells_per_type = None


		self._get_ncells()

	def _check_inputs(self, cell_directories, mixing_coefficients):

		for d in cell_directories:
			assert os.path.exists(d), "Directory %s does not exit" %d

		assert mixing_coefficients.sum() == 1, "Mixing coefficients must sum to 1"

		return		

	def _get_ncells(self):
		""" Crawls each directory counting the number of fastq files to work 
		out how many reads per cell we should take (at random) """

		for directory in self.cell_directories:
			self.cell_fastq_dict[directory] = [f for f in os.listdir(directory) if f.endswith(".fastq")]

		self.cells_per_type = np.array([len(self.cell_fastq_dict[d]) for d in self.cell_directories])
		self.assigned_reads_per_cell = self.assigned_reads / self.cells_per_type

	



