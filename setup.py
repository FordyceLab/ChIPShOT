#!/usr/bin/env python3

from distutils.core import setup

setup(name="ChIPShOT",
      version="0.1",
      description="Chromatin Immunoprecipitation Shape/Occupancy Toolset: \
                   command line tools and python modules to assist in the \
                   analysis of ChIP-seq data for transcription factor binding",
      author="Tyler Shimko",
      author_email="tshimko@stanford.edu",
      url="https://github.com/FordyceLab/ChIPShOT",
      packages=["chipshot"],
      license='MIT',
      include_package_data=True,
      install_requires=[
        "biopython",
        "scipy",
        "tqdm",
      ],
      entry_points={
        "console_scripts": [
            "download_sra=download_sra:main",
            "align=align:main",
            "call_peaks=call_peaks:main"
        ],
    },
      zip_safe=False)
