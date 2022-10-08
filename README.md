# SpectralBP

[![license: MIT](https://img.shields.io/github/license/slashdotfield/SpectralBP)](LICENSE)
[![latest release](https://img.shields.io/github/v/release/slashdotfield/SpectralBP)](https://github.com/slashdotfield/SpectralBP/releases)
![compatibility](https://img.shields.io/badge/Mathematica-_11.x_12.x-brightgreen.svg)
[![DOI](https://zenodo.org/badge/186220358.svg)](https://zenodo.org/badge/latestdoi/186220358)



A Mathematica package for the numerical solution of ODE eigenvalue problems via a pseudospectral method using the Bernstein basis.

## Installation

SpectralBP was developed using Mathematica 11.0.

To install, [download the latest release](https://github.com/slashdotfield/SpectralBP/releases), which is distributed as a paclet file.

Open a Mathematica Notebook, and evaluate

			PacletInstall["dir"]

where dir is the corresponding directory of the downloaded paclet. For example, "C:\\Users\\user\\Downloads\\SpectralBP-1.0.1.paclet".

Alternatively, one may save the paclet to Mathematica's default notebook directory, which may be found out by running the command \$UserDocumentsDirectory

One may then evaluate

			PacletInstall["SpectralBP-1.0.1.paclet"]

If installation is successful, one would be prompted by an output such as Paclet\[SpectralBP-,1.0.1,<>\].

## Use

To use SpectralBP, one may simply evaluate

			Needs["SpectralBP`"]

To find documentation on each command, simply prefix a command with a question mark, such as

			?GetModes

For more details, click the '>>' hyperlink that appears. The documentation window that pops up also allows you to find the tutorial notebooks, on [quantum mechanics](https://github.com/slashdotfield/SpectralBP/blob/master/SpectralBP/Documentation/English/Tutorials/QuantumMechanicsTutorial.nb) and [quasinormal mode](https://github.com/slashdotfield/SpectralBP/blob/master/SpectralBP/Documentation/English/Tutorials/QuasinormalModesTutorial.nb) calculations.

## Uninstallation or Upgrading

To uninstall, simply run the command

			PacletUninstall["SpectralBP"]

To upgrade, Mathematica always uses the latest version of a package, so it is always safe to install a newer version over an older one.

To list all installed versions, simply evaluate

			PacletFind["SpectralBP"]

## License

Sean Julian Fortuna

SpectralBP is covered by the MIT License. For details, see [LICENSE.md](LICENSE).
