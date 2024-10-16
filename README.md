# BamSnap (Fork)

**BamSnap**: A lightweight command-line visualization tool for sequencing reads in BAM files.

**This is a fork of the original [BamSnap](https://github.com/parklab/bamsnap) tool developed by the Park Lab.**

---

## About This Fork

This fork was created to address specific issues with the original BamSnap tool when working with short reference sequences and to update the codebase for compatibility with newer versions of the Pillow library (version 10 and above).

### Why Fork?

- **Handling Short References**: The original BamSnap encountered errors when processing BAM files aligned to small reference sequences (e.g., viral genomes, small plasmids). This fork implements fixes to allow BamSnap to handle short reference sequences gracefully without errors.

- **Compatibility with Pillow 10+**: The original code used methods from the Pillow library that have been deprecated or removed in version 10 and above. This fork updates the code to be compatible with Pillow 10+, ensuring continued functionality with the latest versions of the library.

### Changes Made

- **Bounds Checking for Short References**: Implemented bounds checking and position adjustments throughout the code to ensure genomic positions stay within valid ranges of the reference sequence.

- **Updated Pillow Method Calls**: Replaced deprecated methods such as `getsize()` with their modern equivalents like `getbbox()` to maintain compatibility with Pillow 10+.

- **Code Refactoring**: Minor refactoring and cleanup to improve code readability and maintainability.

---

## Installation

### Prerequisites

- **Python 3.4+**
- [**Pillow (Python Imaging Library)**](https://pypi.org/project/Pillow/) (version 10.0.0 or higher)
- [**pysam**](https://pypi.org/project/pysam/)
- [**pyfaidx**](https://pypi.org/project/pyfaidx/)
- [**pytabix**](https://pypi.org/project/pytabix/)

### Install from GitHub

Clone this forked repository and install BamSnap:

```
*bash*
git clone https://github.com/mendezg/bamsnap.git
cd bamsnap
pip install -e .
*/bash*
```

### Note on Installation

This fork is not currently available via PyPI. Please install it directly from this repository using the instructions above.

---

## Usage

### Simple Usage

```
*bash*
bamsnap -bam test.bam -pos 1:7364529 -out test.png
*/bash*
```

For more details, refer to the BamSnap [**Documentation**](http://bamsnap.readthedocs.io/en/latest).

*Note*: The documentation refers to the original BamSnap. While most functionalities remain the same, some details may differ due to the changes in this fork.

---

## Example Use Case

- [**1000 Genome Data**](https://bamsnap-1kg.s3.amazonaws.com/index.html): 1000 genomic loci in 2504 individuals
- [**BamSnap Plot Gallery**](https://bamsnap.readthedocs.io/en/latest/gallery.html)

*Note*: These resources are from the original BamSnap and are provided for reference.

---

## License

BamSnap is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

---

## Acknowledgments

This fork is based on the original BamSnap tool developed by the Park Lab. We thank the original authors for their work and contributions to the bioinformatics community.

---

## Contact

For questions or comments about this fork, please contact [mendezg](mailto:mendezg@umd.edu).

---

## Contributing

Contributions to improve this fork are welcome. Please submit pull requests or open issues on GitHub.

---

## Changelog

### Version 0.3.0 (Fork)

- **Fixed**: Issues with processing BAM files aligned to short reference sequences.
- **Updated**: Codebase for compatibility with Pillow 10+.
- **Improved**: Code readability and maintainability.
