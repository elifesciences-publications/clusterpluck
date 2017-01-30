from setuptools import setup, find_packages
from glob import glob
import os

__author__ = "Knights Lab"
__copyright__ = "Copyright (c) 2016--, %s" % __author__
__credits__ = ["Robin Shields-Cutler", "Benjamin Hillmann", "Dan Knights"]
__email__ = "cutler.robin@gmail.com"
__license__ = "MIT"
__maintainer__ = "Robin Shields-Cutler"
__version__ = "0.0.1-dev"

long_description = ''

setup(
	name='clusterpluck',
	version=__version__,
	packages=find_packages(),
	url='',
	license=__license__,
	author=__author__,
	author_email=__email__,
	description='',
	long_description=long_description,
	# scripts=glob(os.path.join('clusterpluck', 'scripts', '*py')),
	keywords='',
	install_requires=['numpy', 'scipy'],
	entry_points=dict(console_scripts=[
		'orf_matrix_collapse = clusterpluck.scripts.orf_matrix_collapse:main',
		'extract_cluster_aaseq = clusterpluck.parsers.extract_cluster_aaseq:main',
		'list_cluster_types = clusterpluck.parsers.list_cluster_types:main',
		'blastp_to_matrix = clusterpluck.scripts.blastp_to_matrix:main',
		'combine_metrics = clusterpluck.scripts.combine_metrics:main',
		'parallel_in_common = clusterpluck.scripts.parallel_in_common:main',
		'parallel_orf_collapse = clusterpluck.scripts.parallel_orf_collapse:main',
		'matrix_dicer = clusterpluck.tools.matrix_dicer:main',
		'matrix_serial_dicer = clusterpluck.tools.matrix_serial_dicer:main',
		'matrix_original_dicer = clusterpluck.tools.matrix_original_dicer:main',
		'orfs_in_common = clusterpluck.scripts.orfs_in_common:main',
		'compile_mpfa = clusterpluck.tools.compile_mpfa:main',
		'ofu_matrix = clusterpluck.scripts.ofu_matrix:main',
		'otu_x_ofu = clusterpluck.scripts.otu_x_ofu:main',
		'cluster_lookup = clusterpluck.tools.cluster_lookup:main',
		'shOFUn = clusterpluck.wrappers.shOFUn:main',
	]),
)
