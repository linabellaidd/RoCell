from setuptools import setup, find_packages

setup(
    name = 'RoCell',
    packages = find_packages(exclude=[]),
    include_package_data = True,
    version = '1.0',
    description = 'RoCell: an end-to-end workflow for scRNAseq analysis',
    url = 'https://github.com/linabellaidd/RoCell',
    packages = ['distutils', 'distutils.command'],
    install_requires = [
        'anndata',
        'numpy',
        'pandas',
        'scanpy',
        'scipy',
        'sklearn',
        'torch']
)
