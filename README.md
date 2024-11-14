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

