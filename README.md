# Mark0 in Python

## About The Project

This small project presents a python interface for the simple Mark-0 Agent-based Model. The model is based on the work of Stanislao Gualdi, Jean-Philippe Bouchaud, Marco Tarzia, Francesco Zamponi, and Dhruv Sharma. The relevant references are given below. 

The model is written using C++, and the code was initially developed by Stanislao Gualdi, and extended by Dhruv Sharma. This version of the code has made some additional modifications relating to the reading of parameters (from Python) as well as the exporting of the resulting outputs. Furthermore, the Python interface allows for easy parallelization using the `multiprocessing` module. Further expansions in this direction are planned for a future project. 

**References**<br>
- Gualdi et al. (2015) [_Tipping Points in Macroeconomic Agent-based Models](https://doi.org/10.1016/j.jedc.2014.08.003)
- Gualdi et al. (2017) [_Monetary policy and dark corners in a stylized agent-based model_](http://link.springer.com/10.1007/s11403-016-0174-z)
- Bouchaud et al. (2018) [_Optimal inflation target: insights from an agent-based model_](http://www.economics-ejournal.org/economics/journalarticles/2018-15)
- Sharma et al. (2020) [_V-, U-, L-, or W-shaped recovery after COVID? Insights from an Agent Based Model_](https://doi.org/10.1371/journal.pone.0247823)

<p align="right">(<a href="#top">back to top</a>)</p>

## Getting Started

1. **Compilation**: To be able to use the C++ code underlying Mark-0, it will need to be compiled. The C++ code does make use of the `libgsl` library. A sample compile statement on the linux machine would involve entering the directory with the mark0 cpp file and executing: `g++ mark0_covid.cpp -l:libgsl.so.23.1.0 -o mark0_covid`. 
2. **Using the model**: The model, once compiled, can be accessed using `from mark0_covid import Mark0_COVID`, where the Mark0_COVID object has functions to simulate and classify phases. An example of this is provided in the jupyter noteboo, where the Mark0_COVID object has functions to simulate and classify phases. An example of this is provided in the jupyter notebook.
3. **Requirements**: the python environment requires `pandas`, and `yaml`


<p align="right">(<a href="#top">back to top</a>)</p>

## License
Distributed under the MIT License. See `LICENSE.txt` for more information.

<p align="right">(<a href="#top">back to top</a>)</p>

## Contact

Karl Naumann-Woleske - [@KNaumannWoleske](https://twitter.com/KNaumannWoleske) - karl [at] karlnaumann.com

Project Link: [https://github.com/KarlNaumann/Mark0](https://github.com/KarlNaumann/Mark0)

<p align="right">(<a href="#top">back to top</a>)</p>

<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/github/contributors/KarlNaumann/Mark0.svg?style=for-the-badge
[contributors-url]: https://github.com/KarlNaumann/Mark0/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/KarlNaumann/Mark0.svg?style=for-the-badge
[forks-url]: https://github.com/KarlNaumann/Mark0/network/members
[stars-shield]: https://img.shields.io/github/stars/KarlNaumann/Mark0.svg?style=for-the-badge
[stars-url]: https://github.com/KarlNaumann/Mark0/stargazers
[issues-shield]: https://img.shields.io/github/issues/KarlNaumann/Mark0.svg?style=for-the-badge
[issues-url]: https://github.com/KarlNaumann/Mark0/issues
[license-shield]: https://img.shields.io/github/license/KarlNaumann/Mark0.svg?style=for-the-badge
[license-url]: https://github.com/KarlNaumann/Mark0/LICENSE.txt
