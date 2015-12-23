from setuptools import setup, find_packages

setup(
    name = 'CorrelatedVariants',
    version = '0.1.0',
    author = 'Pacific Biosciences',
    author_email = 'devnet@pacificbiosciences.com',
    license = open('LICENSE.txt').read(),
    packages = find_packages('.'),
    package_dir = {'':'.'},
    zip_safe = False,
    scripts=[
        'bin/correlatedVariants',
        'bin/rareCaller'
    ],
    install_requires = [
        'pbcore >= 0.2',
        'GenomicConsensus',
        'numpy >= 1.6.0',
        'scipy >= 0.9.0',
        'h5py >= 1.3.0'
    ]
)
