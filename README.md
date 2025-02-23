## Table of Contents

1. [Introduction](#introduction)
   1. [What is DynPy?](#1-what-is-dynpy)
   2. [Key Features](#2-key-features)
   3. [Getting Started on CoCalc](#3-getting-started-on-cocalc)
2. [How to Start / Basic Usage](#how-to-start--basic-usage)
   1. [Example Scripts](#1-example-scripts)
   2. [Creating First Document / Report](#2-creating-first-document--report)
3. [Reporting Module](#reporting-module)
   1. [Overview of the Reporting Module](#1-overview-of-the-reporting-module)
   2. [Creating Reports](#2-creating-reports)
      - [Setting Up the Reporting Environment](#setting-up-the-reporting-environment)
      - [Defining Report Content](#defining-report-content)
      - [Creating Sections and Subsections](#creating-sections-and-subsections)
      - [Adding Content (Text, Images, Equations)](#adding-content-text-images-equations)
      - [Incorporating Simulation Results](#incorporating-simulation-results)
      - [Adding Visualizations and Data Tables](#adding-visualizations-and-data-tables)
   3. [Exporting Reports](#3-exporting-reports)
      - [Supported Formats](#supported-formats)
      - [Exporting Procedures](#exporting-procedures)
   4. [Practical Examples](#4-practical-examples)
      - [Generating a Simple Report](#generating-a-simple-report)
      - [Advanced Reporting Features](#advanced-reporting-features)
   5. [Customization Options (Advanced)](#5-customization-options-advanced)
      - [Formatting Text and Equations](#formatting-text-and-equations)
      - [Customizing Layout and Styles](#customizing-layout-and-styles)
      - [Utilizing Templates for Consistency](#utilizing-templates-for-consistency)
4. [Simulation Engine](#simulation-engine)
5. [Data Handling](#data-handling)
6. [Dynamic Modeling](#dynamic-modeling)
7. [Visualization Tools](#visualization-tools)
8. [Installation & Setup (Optional, for Local Development)](#installation--setup-optional-for-local-development)
   1. [Requirements](#1-requirements)
   2. [Manual Installation](#2-manual-installation)
9. [Usage Examples](#usage-examples)
   1. [Simulating a Dynamic System](#1-simulating-a-dynamic-system)
   2. [Data Import & Export](#2-data-import--export)
   3. [Running a Custom Model](#3-running-a-custom-model)
   4. [Generating Reports](#4-generating-reports)
10. [Licensing Information](#licensing-information)

# Introduction

## 1. What is DynPy?

DynPy is a Python module designed for engineering calculations on dynamical systems. It provides a comprehensive framework for modeling, simulating, and analyzing dynamic mechanical systems.

## 2. Key Features

- **Dynamics Module:** Tools for modeling mechanical systems and their dynamics.
- **Mechanical Models:** A collection of predefined mechanical models developed by experts.
- **Symbolic and Numeric Solvers:** Tools for solving Ordinary Differential Equations (ODEs) using symbolic and numerical methods.
- **Reporting Module:** A structured reporting system for generating and exporting reports.

## 3. Getting Started on CoCalc

To begin working with DynPy, you need an account on [CoCalc](https://cocalc.com/).

1. Create an account on CoCalc.
2. Accept the project invitation using this [link](https://cocalc.com/app?project-invite=hXnPFLqokQsoK6TG).
3. Open the [README FIRST](https://cocalc.com/projects/b51ce971-5b39-4911-ad97-ef59f15f0039/files/README%20FIRST.ipynb) file.
4. Follow the instructions in the introductory guide.

---

# How to Start / Basic Usage

## 1. Example Scripts

To view exemplary capabilities of dynpy, run the following example script:

```python
from dynpy.models.mechanics.pendulum import Pendulum

Pendulum().interactive_preview()
```

## 2. Creating First Document / Report
```python
from dynpy.utilities.report import *
from dynpy.utilities.documents.guides import Guide

doc = Guide('./output/sample_report', title="Sample Report")

section = Section('Exemplary section name')
CurrentContainer(section)

display(Markdown(''' Exemplary Markdown text in the section '''))
display(ReportText(' Exemplary text appended into section '))

doc.append(section)

doc.generate_pdf(clean_tex=True)
```

## 3. Looking for some help
```python
from dynpy.utilities.creators import list_of_guides
list_of_guides()
```
### Defining Report Content

#### Creating a section of document
```pyton
section = Section('Sample section name')
```

#### Creating a subsection of document
```pyton
subsection = Subsection('Sample subsection name');
```

#### Selecting a section or subsection to add content to
```pyton
section = Section('Sample section name')
CurrentContainer(section);
```

#### Adding text to section via ReportText
```pyton
display(ReportText('Sample text'));
```

#### Adding text to section via Markdown
```pyton
display(Markdown(
'''
Sample text 
'''
))
```

#### Adding an image into the section
```pyton
Picture('/route/to/image', caption = 'Sample caption')
```

#### Appending sections and subsections into the document
```python
doc.append(sample_section_name)
```

### Incorporating Simulation Results
```python
# Add simulated data here
```

### Adding Visualizations, Formulas and Data Tables
Creating plot and adding it to document
```python
import matplotlib.pyplot as plt

def create_plot():
    plt.plot([0, 1, 2], [0, 1, 4])
    plt.savefig("./plot.png")

Picture('./plot.png', caption='Sample plot')
```

Adding formula to the document
```python
d, r, fib, fia, thetaa, thetab = symbols('d r varphi_B varphi_A phi_A phi_B');

harvestine_formula = Eq(d, 2 * r * sp.asin(sp.sqrt(sp.sin((fib - fia) / 2)**2 + (cos(fia) * cos(fib) * sp.sin((thetab - thetaa) / 2)**2))))

display(SympyFormula(harvestine_formula))
```

Creating table and adding it to document
```python
from dynpy.utilities.adaptable import *

predicted_travel_time = Subsection('Predicted Travel Time');
CurrentContainer(predicted_travel_time);

time_s = Symbol('time_s', positive=True)
time_h = Symbol('time_h')
length = Symbol('length', positive=True)
velocity = Symbol('velocity', positive=True)

dane = {
    'Start': ['NYC', 'Albany', 'Syracuse', 'Buffalo', 'Cleveland', 'Toledo', 'Total line'],
    'Stop': ['Albany', 'Syracuse', 'Buffalo', 'Cleveland', 'Toledo', 'Chicago', ''],
    time_s: [3348, 3386, 2782, 4362, 2606, 4824, 21308],
    time_h: ['00:55:48', '00:56:26', '00:46:22', '01:12:42', '00:43:26', '01:20:24', '05:55:08'],
    length: [215.981, 219.844, 225.822, 295.54, 172.905, 369.093, 1499.185],
    velocity: [232.24, 233.74, 292.22, 243.91, 238.86, 275.44, 253.29]
}

unit_dict = {
    time_s: ureg.second,
    time_h: ureg.hour,
    length: ureg.kilometer,
    velocity: ureg.kilometer / ureg.hour
}

LatexDataFrame.set_default_units(unit_dict)

def format_cell(x):
    if isinstance(x, str):
        return x
    else:
        return f'${latex(x)}$'

tabelka = LatexDataFrame.formatted(
    data = dane,
).map(format_cell)

tabelka.columns = ['Start', 'Stop', 'Time [s]', 'Time [h]', 'Length [km]', 'Velocity [km/h]']

predicted_travel_time.append(tabelka.reported(caption="Travel Time Data Table"))
```

## 3. Exporting Reports
### Supported Formats
- **PDF**
- **LaTeX**
- **Markdown**
### Exporting Procedures
```python
doc.generate_pdf(clean_tex=True)
```

## 4. Practical Examples

### Generating a Simple Report
```python
from dynpy.utilities.report import *
from dynpy.utilities.templates.document import Guide

doc = Guide('./reports/sample-report')

sample_section = Section('Sample Section')
CurrentContainer(sample_section)

display(ReportText('Sample text'));

second_sample_section = Section('Second Sample Section')
CurrentContainer(second_sample_section)

display(ReportText('Sample text'));

sample_subsection = Subsection('Sample Subsection');
CurrentContainer(sample_subsection);

display(ReportText('Sample text'));

second_sample_subsection = Subsection('Second sample Subsection');
CurrentContainer(second_sample_subsection);

display(ReportText('Sample text'));

doc.append(sample_section)
doc.append(second_sample_section)
doc.append(sample_subsection)
doc.append(second_sample_subsection)

doc.generate_pdf(clean_tex=True)
```

## 5. Customization Options (Advanced)

### Formatting Text and Equations
```python
///
```

### Customizing Layout and Styles
```python
///
```

# Custom styles

### Utilizing Templates for Consistency
```python
///
```

### Use predefined templates
```python
///
```

# Simulation Engine
```python
///
```

# Data Handling
```python
///
```

# Dynamic Modeling
```python
///
```
# Installation & Setup (Optional, for Local Development)
## Requirements
Python Version: **Python 3.8+**. Required Libraries:
- **numpy**
- **pylatex**
- **sympy**
- **pandas**
- **matplotlib**
- **scipy**
- **pint**
- **pypandoc**
- **wand**
## Manual Installation
```bash
pip install numpy pylatex sympy pandas matplotlib scipy pint pypandoc wand

git clone https://github.com/bogumilchilinski/dynpy.git
```

# Licensing Information
DynPy is distributed under an open-source license. Refer to the LICENSE file for details.
