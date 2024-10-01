# Introduction

It's a module that allows to engineering calculations on dynamical systems. 

There are four main parts of the entire project:

- dynamics module

- mechanical models - lead of development and maintenance: Amadeusz Radomski (@amvdek); Grzegorz Długopolski (@grzegorzdl);

- symbolic and numeric solvers for ODE systems;

- reporting module.

Using the code below in Jupyter enviroment on [Free Access Project](https://cocalc.com/app?project-invite=hXnPFLqokQsoK6TG) <- ([CLICK LINK](https://cocalc.com/app?project-invite=hXnPFLqokQsoK6TG)) we can learn more about how to and what to use Python in engineering calculations:

```python {kernel="python3"}
from dynpy.utilities.documents.guides import IntroToCocalcGuide, UsageOfDynamicSystemsGuide

IntroToCocalcGuide();
```

Run this code in the blank Jupyter you have created.

You will see the guide in Output after running it, i.e. a CELL\-by\-CELL \(step\-by\-step\) procedure.

If you are looking for information on reporting and creating a PDF file, we can use the command below to view the tutorial:

```python {kernel="python3"}
from dynpy.utilities.documents.guides import BasicsOfReportingGuide
BasicsOfReportingGuide();
```

# Help and guides for DynPy

You can list all of the available guides with the following call:

```python {kernel="python3"}
from dynpy.utilities.creators import list_of_guides
list_of_guides()
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

