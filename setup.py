import os
from setuptools import setup

# version_py = os.path.join(os.path.dirname(__file__), 'monas', 'version.py')
# version = open(version_py).read().strip().split(
#     '=')[-1].replace('"', '').strip()
long_description = """
"""

version=0.5

HERE = os.path.dirname(__file__)

install_requires = ['biopython']

# with open(os.path.join(HERE, "requirements.txt"), "r") as f:
#     install_requires = [x.strip() for x in f.readlines() if x.strip() not in require_exemption]

def package_files(directory):
    paths = []
    for (path, directories, filenames) in os.walk(directory):
        for filename in filenames:
            paths.append(os.path.join('..', path, filename))
    return paths

extra_files = package_files('monas/misc')
extra_files.extend(package_files('monas/references'))
extra_files.append('bin_path.json')
extra_files.append('kdr_list.json')
print(extra_files)

setup(
    name="monas",
    version=version,
    install_requires=install_requires,
    requires=['python (>=3.10)'],
    packages=['monas'],
    author="Kentaro Itokawa",
    description='Genotyping and annotate VGSC gene mutations in insects',
    long_description=long_description,
    url="https://github.com/ItokawaK/MoNaS",
    package_dir={'monas': "monas"},
    package_data={
        'monas': extra_files,
        },
    zip_safe=False,
    include_package_data=True,
    entry_points={
        'console_scripts': [
            'monas=monas.genotype:main',
            'monas_sanger=monas.genotype_sanger:main',
            'bed2gff3=monas.bed2gff3:main',
            'extract_exons=monas.extract_exons:main',
            'finalize_table=monas.finalize_table:main',
            'make_AA_alignment=monas.make_AA_alignment:main',
            'make_jbrowse_config=monas.make_jbrowse_config:main'


        ],
    },
    author_email="itokawa@niid.go.jp",
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ]
)