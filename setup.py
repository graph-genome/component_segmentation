from setuptools import setup, find_packages

setup(
    name='component_segmentation',
    version='1.0.0',
    description='Segmenting syntenic blocks of shared sequence.',
    author='Josiah Seaman, Dmytro Trybushnyi, Artem Tarasov, Christian Kubica, Simon Heumos, Thomas Townsley, Pantograph Team',
    author_email='josiah@newline.us',
    packages=find_packages(exclude=('data',)),
    include_package_data=True,
    scripts=['matrixcomponent/segmentation.py'],
    install_requires=[
        'nose>=1.3.7',
        'rdflib>=4.2.2',
        'DNASkittleUtils>=1.0.13',
        'numpy>=1.18.2',
        'pytest>=5.4.1',
        'sortedcontainers>=2.1.0',
        'joblib>=0.14.1',
        'numba>=0.48.0',
        'psutil>=5.7.0',
        'recordclass>=0.13.2'
    ],
    url='https://graphgenome.org',
    download_url='https://github.com/graph-genome/component_segmentation',
    keywords=['bioinformatics', 'dna', 'fasta', 'graph genome', 'GFA', 'MSA'],
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: Apache Software License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
)