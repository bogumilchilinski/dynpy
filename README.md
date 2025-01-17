# Introduction

It's a module that allows to engineering calculations on dynamical systems. 

There are four main parts of the entire project:

- dynamics module

- mechanical models - lead of development and maintenance: Amadeusz Radomski (@amvdek); Grzegorz Długopolski (@grzegorzdl);

- symbolic and numeric solvers for ODE systems;

- reporting module.

First step for starting a project is to create an account in [COCALC](https://cocalc.com/). 

Then, using the following [LINK](https://cocalc.com/app?project-invite=hXnPFLqokQsoK6TG), accept the invitation.

Afterwards, you will be directed to the page, where you should click the [README FIRST](https://cocalc.com/projects/b51ce971-5b39-4911-ad97-ef59f15f0039/files/README%20FIRST.ipynb) file (you can click this link if you have trouble seeing the page). There, you have access to the introductory code, which is prepared for you.

In this file, you will find the essential information on how to create a blank Jupiter (where you will run the codes), use Cocalc, access usefull commands and more...

# Help and guides for DynPy

You can access the introductory guide with the following code:

```python {kernel="python3"}
from dynpy.utilities.documents.guides import IntroToCocalcGuide, UsageOfDynamicSystemsGuide

IntroToCocalcGuide();
```

You can list all of the available guides with the following call:

```python {kernel="python3"}
from dynpy.utilities.creators import list_of_guides
list_of_guides()
```

If you are looking for information on reporting and creating a PDF file, we can use the command below to view the tutorial:

```python {kernel="python3"}
from dynpy.utilities.documents.guides import BasicsOfReportingGuide
BasicsOfReportingGuide();
```

# Dynamic systems

Next for an example, run the codes below and you will see how it works:

You can preview the pendulum using such a function.

```python {kernel="python3"}
import sympy 
from sympy import Symbol

from dynpy.models.mechanics.pendulum import Pendulum

Pendulum().interactive_preview()
```

# Reporting Module Usage

This reporting module is designed to create, compile and export a report into a PDF document. Below is a breakdown of the main tools and functions used to make exporting the report possible, along with instructions on how to install the required dependencies.

## Dependencies

The module relies on the following Python libraries:

1. **numpy** - Used for numerical computations.
2. **pylatex** - Used to create and generate LaTeX documents in Python, enabling PDF export.
3. **sympy** - For symbolic mathematics and algebraic calculations.
4. **pandas** - For data manipulation and analysis.
5. **matplotlib** - For creating visualizations and plots.
6. **scipy** - For scientific computing and technical computations.
7. **pint** - A library for defining, operating, and converting physical quantities and units.
8. **pypandoc** - A simple interface for Pandoc, a universal document converter. It allows you to convert files between different markup formats (e.g., Markdown to LaTeX, HTML to DOCX) directly from Python.
9. **wand** - A Python binding for ImageMagick, a powerful image manipulation tool. It is used to process and manipulate images, such as resizing, cropping, and format conversion, directly from Python scripts.

### Installation Instructions

To install all the dependencies, you can either install them individually or use a `requirements.txt` file. Below are both methods.

#### Option 1: Installing Individually

Run the following commands in your terminal or command prompt to install each dependency:

```bash
pip install numpy
pip install pylatex
pip install sympy
pip install pandas
pip install matplotlib
pip install scipy
pip install pint
pip install pypandoc
pip install wand
```

### Option 2: Using requirements.txt

Alternatively, You can install all the dependencies at once by running the following command:

```bash
pip install -r requirements.txt
```

### Verify Installation

After installation, you can verify that the libraries have been successfully installed by running the following Python script:

```python
import numpy
import pylatex
import sympy
import pandas
import matplotlib
import scipy
import pint
import pypandoc
import wand

print("All libraries are installed successfully!")
```

# Creating a Simple Document

To create a document, follow the steps below.

## Import the Required Modules

Begin by importing the necessary libraries from dynpy and sympy:

```python
from dynpy.utilities.report import *
from dynpy.utilities.templates.document import Guide
```

## Initialize the Document

To start creating a document, use the Guide class to specify the output file location and filename:

```python
doc = Guide('./path/to/save/the/document/exemplary-file-name')
```

This initializes the document and sets the location where the final PDF file will be saved.

## Creating Sections and Subsections

Documents are structured into sections and subsections. Use the Section and Subsection classes to create these.

### Creating a Section:

```python
exemplary_section = Section('Exemplary section name')
```

This creates a section with the title "Exemplary section name."

To make the section active (i.e., the current container for adding content), use:

```python
CurrentContainer(exemplary_section)
```

### Creating a Subsection:

Similarly, subsections can be created within sections:

```python
exemplary_subsection = Subsection('Exemplary subsection name')
CurrentContainer(exemplary_subsection)
```

This creates a subsection and makes it the current container.

## Adding Content to Sections

Various types of content, such as text, Markdown, images, and symbolic formulas, can be added to the sections and subsections.

### Adding Plain Text

To add plain text to a section or subsection, use the ReportText class with display:

```python
display(ReportText('Exemplary text'))
```

### Adding Markdown Text

You can add Markdown-formatted text using the Markdown class:

```python
display(Markdown('''Exemplary Markdown text'''))
```

### Adding Images

To add an image, use the Picture class. You can also add optional captions and set custom widths for the image:

```python
Picture('path/to/file/', caption='Exemplary caption.', width=NoEscape('0.3\\textwidth'))
```

### Adding Sympy Formulas

You can insert symbolic mathematical formulas using SympyFormula. For example, to display a symbolic formula defined in sympy, use:

```python
import sympy as sp

a, b, c = sp.symbols('a b c')

exemplary_formula = sp.Eq(
    a, 
    sp.asin(sp.sqrt(
        sp.sin((b - c) / 2)**2 )
    ))
)

display(SympyFormula(exemplary_formula))
```

This will render the symbolic formula in LaTeX format.

## Appending Sections to the Document

Once you've created sections and added content, append the sections to the document:

```python
doc.append(exemplary_section)
doc.append(exemplary_subsection)
```

## Generating the PDF

After adding all the content, you can generate the PDF using the generate_pdf() function. The clean_tex=True option automatically removes the intermediate LaTeX files:

```python
doc.generate_pdf(clean_tex=True)
```

This will create a PDF file at the location specified when initializing the document.

