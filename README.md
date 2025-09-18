## Table of Contents

- [Introduction](#introduction)
  - [1. What is DynPi?](#1-what-is-dynpi)
  - [2. Key Features](#2-key-features)
  - [3. Getting Started on CoCalc](#3-getting-started-on-cocalc)
- [How to Start / Basic Usage](#how-to-start--basic-usage)
  - [1. Example Scripts](#1-example-scripts)
  - [2. Creating First Document / Report](#2-creating-first-document--report)
  - [3. Looking for some help](#3-looking-for-some-help)
    - [Defining Report Content](#defining-report-content)
      - [Creating a section of document](#creating-a-section-of-document)
      - [Creating a subsection of document](#creating-a-subsection-of-document)
      - [Selecting a section or subsection to add content to](#selecting-a-section-or-subsection-to-add-content-to)
      - [Adding text to section via ReportText](#adding-text-to-section-via-reporttext)
      - [Adding text to section via Markdown](#adding-text-to-section-via-markdown)
      - [Adding an image into the section](#adding-an-image-into-the-section)
      - [Appending sections and subsections into the document](#appending-sections-and-subsections-into-the-document)
    - [Incorporating Simulation Results](#incorporating-simulation-results)
    - [Adding Formulas, Data Tables and Visualizations](#adding-formulas-data-tables-and-visualizations)
  - [3. Exporting Reports](#3-exporting-reports)
    - [Supported Formats](#supported-formats)
    - [Exporting Procedures](#exporting-procedures)
  - [4. Practical Examples](#4-practical-examples)
    - [Generating a Simple Report](#generating-a-simple-report)
- [ODESystem Class Usage](#odesystem-class-usage)
    - [ODESystem class usage preview based on a damped mehanical oscillator](#odesystem-class-usage-preview-based-on-a-damped-mehanical-oscillator)
    - [ODESystem operations on dynamical systems](#odesystem-operations-on-dynamical-systems)
    - [More information about ODESystem class:](#more-information-about-odesystem-class)
  - [5. Customization Options (Advanced)](#5-customization-options-advanced)
    - [Formatting Text and Equations](#formatting-text-and-equations)
    - [Customizing Layout and Styles](#customizing-layout-and-styles)
- [Custom styles](#custom-styles)
    - [Utilizing Templates for Consistency](#utilizing-templates-for-consistency)
    - [Use predefined templates](#use-predefined-templates)
- [Simulation Engine](#simulation-engine)
- [Data Handling](#data-handling)
- [Dynamic Modeling](#dynamic-modeling)
- [Installation \& Setup (Optional, for Local Development)](#installation--setup-optional-for-local-development)
  - [Requirements](#requirements)
  - [Manual Installation](#manual-installation)
- [Licensing Information](#licensing-information)

# Introduction

## 1. What is DynPi?

DynPi is a Python module designed for engineering calculations on dynamical systems. It provides a comprehensive framework for modeling, simulating, and analyzing dynamic mechanical systems. DynPy is the first tool of its kind to simplify performing scientific activities and enginnering calculations! It’s ready-to-use code that generates reports in engineering mechanics. You can create your own mechanical simulations, use tools for modeling mechanical systems and their dynamics, and take advantage of a built-in solver for ordinary differential equations (ODEs). All of this is available in just 10 code cells – with every step and hint clearly described below and within the code itself (marked with #). We’re also including an instructional video showing how to create a sample report:
Writing your own diploma still looking for tools to use? Word, DynPy and .... ChatGPt. These tools have its pros and cons - all depends on how do you use it and how do u like to work. Word it's classic - everyone knows and can use it. Starter for wrtining, formatting, adding adnotations and commentaries. Good for theory, literature reviews, humanities papers. However, when it comes to simulations and calculations, then sorry you're obliged to open other programs.
DynPy is an engineering combine. Simulations, calculations and modelling it's hard and DynPy is superior in this especailly if you are into mechanics and automation. DynPy allows you to create professional looking raports. Moreover they are created in PDF format. Cons? You must pick up a little coding. Once you learn it, it works wonders.
ChatGPT? It's your go-to buddy for everything – it'll sometimes help write an introduction, check your spelling, suggest how to calculate something, or explain what a given piece of code does. However, it won't always be accurate, it won't substitute DynPy for simulations and reporting.

## 2. Key Features

- **Dynamics Module:** Tools for modeling mechanical systems and their dynamics.
- **Mechanical Models:** A collection of predefined mechanical models developed by experts.
- **Symbolic and Numeric Solvers:** Tools for solving Ordinary Differential Equations (ODEs) using symbolic and numerical methods.
- **Reporting Module:** A structured reporting system for generating and exporting reports.

## 3. Getting Started on CoCalc

To begin working with DynPi, you need an account on [CoCalc](https://cocalc.com/).

1. Create an account on CoCalc.
2. Accept the project invitation using this [link](https://cocalc.com/app?project-invite=hXnPFLqokQsoK6TG).
3. Open the [README](https://cocalc.com/projects/b51ce971-5b39-4911-ad97-ef59f15f0039/files/READme.ipynb) file.
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
from dynpy.utilities.documents.document import Report
SympyFormula._break_mode = 'eq'

doc = Report('./output/sample_report', title="Sample Report")

section = Section('Exemplary section name')
CurrentContainer(section)

display(Markdown(''' Exemplary Markdown text in the section '''))
display(ReportText(' Exemplary text appended into section '))

doc.append(section)

doc.generate_pdf(clean_tex=False)
```

## 3. Looking for some help

Documentclasses with generic content of exemplary document:

```python
from dynpy.utilities.documents.document import *

Report.base_setup()
BeamerPresentation.base_setup()
Guide.base_setup()
WutThesis.base_setup()
```

Guides that provides step by step instructions:

```python
from dynpy.utilities.creators import list_of_guides
list_of_guides()
```

### Defining Report Content

#### Creating a section of document

```python
section = Section('Sample section name')
```

#### Creating a subsection of document

```python
subsection = Subsection('Sample subsection name');
```

```result
```

#### Selecting a section or subsection to add content to

```python
section = Section('Sample section name')
CurrentContainer(section);
```

#### Adding text to section via ReportText

```python
from dynpy.utilities.report import *

display(ReportText('Sample text'));
```

#### Adding text to section via Markdown

```python
from dynpy.utilities.report import *

display(Markdown(
'''
Sample text
'''
))
```

#### Adding an image into the section

```python
Picture('./images_folder/image_name.PNG', caption = 'Sample caption') # './images_folder/image_name.PNG' is certain '/route/to/file' path
```

#### Appending sections and subsections into the document

```python
doc.append(sample_section_name)
```

### Incorporating Simulation Results

```python
# Add simulated data here
```

### Adding Formulas, Data Tables and Visualizations

Adding formula to the document

```python
from dynpy.utilities.report import *
from sympy import *


d, r, fib, fia,  = symbols('d r varphi_B varphi_A') #many symbols at once

thetaa = Symbol('thetaa') #separate definition of the symbol
thetab = Symbol('thetab') #separate definition of the symbol


harvestine_formula = Eq(d, 2 * r * asin(sqrt(sin((fib - fia) / 2)**2 + (cos(fia) * cos(fib) * sin((thetab - thetaa) / 2)**2))))
display(SympyFormula(harvestine_formula))
```

Creating table and adding it to document

```python
from dynpy.utilities.report import *
from dynpy.utilities.adaptable import *

predicted_travel_time = Subsection('Predicted Travel Time');
CurrentContainer(predicted_travel_time)

time_s = Symbol('t_s', positive=True)
time_h = Symbol('t_h')
length = Symbol('l', positive=True)
velocity = Symbol('v', positive=True)

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



report_table = LatexDataFrame.formatted(
    data = dane)


display(report_table.reported(caption="Travel Time Data Table"))
```



Creating plot and adding it to document

```python
import matplotlib.pyplot as plt
from dynpy.utilities.report import PltPlot


plt.plot([0, 1, 2], [0, 1, 4],label='linear')
plt.xlabel('x-axis')
plt.ylabel('y-axis')
plt.title('Sample Plot')
plt.legend()

data_plot=PltPlot(caption='Sample plot')
display(data_plot)
```


```python
from dynpy.utilities.report import *
from dynpy.utilities.adaptable import *

LatexDataFrame.set_default_units(unit_dict)


results_df = LatexDataFrame.formatted(
    data = dane)


#METHOD #1 (optional way) - via matplotlib
results_df.to_latex_dataframe.plot()
dataplot_plt=PltPlot(caption='Sample plot')
display(dataplot_plt)

#METHOD #2 (preferred way) - via matplotlib
dataplot_tikz=results_df.to_pylatex_tikz(height=8).in_figure(caption='Sample plot 2')#comment for preview
display(dataplot_tikz)#comment for preview
```

There is a possibility of providing analitycal solution

```python
from dynpy.solvers.linear import AnalyticalSolution
from sympy import *
import numpy as np

t =Symbol('t')
x = Function('x')(t)



results_df=AnalyticalSolution.from_vars_and_rhs([x],[cos(t)+sin(t)**2]).compute_solution(np.linspace(0,100,10001))#.plot()

#METHOD #1 (optional way) - via matplotlib
results_df.to_latex_dataframe.plot()
dataplot_plt=PltPlot(caption='Sample plot')
display(dataplot_plt)

#METHOD #1 (preferred way) - via matplotlib
dataplot_tikz=results_df.to_pylatex_tikz(height=8).in_figure(caption='Sample plot 2')#comment for preview
display(dataplot_tikz)#comment for preview
```

There is a possibility of providing analitycal solution

```python
from dynpy.solvers.linear import AnalyticalSolution
from sympy import *
import numpy as np

t =Symbol('t')
x = Function('x')(t)
y = Function('y')(t)

results_signals_df=AnalyticalSolution.from_vars_and_rhs([x,y],[cos(t)+sin(t)**2,1+exp(cos(t/2)+sin(t/2)**2)]).compute_solution(np.linspace(0,100,10001))#.plot()


### PREVIEW
#results_signals_df.plot(title='Sample plot', xlabel='Time [s]', ylabel='Value', figsize=(14, 3))#uncomment for preview, !!!!DON"T USE IN RAPORT!!!!

data_plot=results_signals_df.to_pylatex_tikz(height=8).in_figure(caption='Time history of the process')#comment for preview
display(data_plot)

```

## 3. Exporting Reports

### Supported Formats

- **PDF**
- **LaTeX**
- **Markdown**

### Exporting Procedures

```python

#for LaTeX report (only LaTeX distribution is needed)
from dynpy.utilities.creators import PdfLatexGenerator

PdfLatexGenerator(doc).generate_file()

#for LaTeX report (LaTeX distribution and perl are needed)
doc.generate_pdf(clean_tex=False)

#for `.docx` file
import pypandoc
doc.generate_tex('./output/sample_report')
pypandoc.convert_file('./output/sample_report.tex',to='docx',format='tex',outputfile="./output/sample_report.docx")

```

## 4. Practical Examples

### Generating a Simple Report

```python
from dynpy.utilities.report import *
from dynpy.utilities.documents.document import Report

doc = Report('./output/sample-report')

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

doc.generate_pdf(clean_tex=False)
```

# ODESystem Class Usage

### ODESystem class usage preview based on a damped mehanical oscillator

```python
from sympy import *
from dynpy.solvers.linear import ODESystem
from sympy.physics.mechanics import dynamicsymbols

m = Symbol('m',positive = True)
c = Symbol('c',positive = True)
k = Symbol('k',positive = True)
t=Symbol('t')
x=dynamicsymbols('x')

 
eq1 = Eq(m*x.diff(t,2)+c*x.diff(t)+k*x,0)


odesys = ODESystem(odes = Matrix([eq1.lhs-eq1.rhs]),dvars =Matrix([x],ode_order=2))

odesys.solution
```

### ODESystem operations on dynamical systems

```python
from dynpy.models.mechanics.trolley import ForcedSpringMassSystem
import numpy as np
from sympy import *
dsys =ForcedSpringMassSystem()
t =dsys.ivar
z = dsys.z
k=dsys.k
F=dsys.F
g=dsys.g

m=dsys.m

dane={
  k:1,
  m:1,
  F:0,
  g:9.81,

}
t_span=np.linspace(0.0,10,1000)

ode=dsys.eoms


ode.subs(dane).numerized().compute_solution(t_span, [0.1, 0.0]).plot()

ode.solution.with_ics([0.1,0]).subs(dane).numerized().compute_solution(t_span).plot()
```

### More information about ODESystem class:

```python
from dynpy.utilities.documents.guides import BasicsOfODESystemGuide
BasicsOfODESystemGuide()


```

```result
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
- **pygithub**
- **wand**
- **pymupdf**

## Manual Installation

```bash
pip install numpy pylatex sympy pandas matplotlib scipy pint pypandoc wand pymupdf

pip install dynpi
```

Installing the Development Environment for Engineering Analysis

1. Install Python (Microsoft Store)
2. Install Visual Studio Code (Microsoft Store)
2.5 Install Git
3. Install the DynPy library (https://github.com/bogumilchilinski/dynpy)
Git clone command https://github.com/bogumilchilinski/dynpy
4. Install the library - instructions on GitHub - bogumilchilinski/dynpy (https://github.com/bogumilchilinski/dynpy?tab=readme-ov-file#manual-installation)
4.5 Installing the plugin in VS Code (git extension package + latex workshop)
5. Creating a virtual environment in VSCode
Working folder on the main drive + subfolders (output, images)
Set the kernel and Jupyter Notebook environment
Set Git Autofetch: True in VSCode settings

Support LaTeX:
6. Install a separate Latex distribution - e.g., MikaTex (https://miktex.org/download#dok)
7. (optional) Strawberry Pearl environment
Installation: https://strawberryperl.com/ + remaining steps to obtain the PDF (admin required)

# Licensing Information

DynPy is distributed under an open-source license. Refer to the LICENSE file for details.

