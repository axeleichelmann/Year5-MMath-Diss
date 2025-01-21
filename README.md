# University of Edinburgh Master's Dissertation - 'Modelling Interseasonal Energy Storage'

![UoE banner](assets/school-of-maths.png)

This repository contains the files and scripts used to mathematically model the implementation of various long term Thermal Energy Storage (TES) technologies into the UK's energy grid.

The main contents of this directory are as follows:

* [sPLAN guide.pdf](https://github.com/axeleichelmann/Year5-MMath-Diss/blob/main/sPLAN%20guide.pdf) - Contains user guide for the original integrated system environment of energy investment planning problems (Planner) written by [Dr. Rodrigo Garcia Nava](https://www.linkedin.com/in/rodrigogarciana/), upon which the models of the UK energy grid used in this dissertaion were built. The formulation in Julia of the optimisation model being solved by the Planner can be found in the **
* [Modelling_ISES_Final_Report.pdf](https://github.com/axeleichelmann/Year5-MMath-Diss/blob/main/Modelling_ISES_Final_Report.pdf) - The final 52 page dissertation report detailing the research problem, goals, literature review, method, results, and conclusions. This report achieved a 1st - the highest possible grade.
* [Final Models/BTES_Model](https://github.com/axeleichelmann/Year5-MMath-Diss/tree/main/Final%20Models/BTES_Model) - Contents for the Planner that models the implementation of Borehole Thermal Energy Storage (BTES) into the UK energy grid.
* [Final Models/PTES_Model](https://github.com/axeleichelmann/Year5-MMath-Diss/tree/main/Final%20Models/PTES_Model) - Contents for the Planner that models the implementation of Pit Thermal Energy Storage (PTES) into the UK energy grid.
* [Final Models/TTES_Model](https://github.com/axeleichelmann/Year5-MMath-Diss/tree/main/Final%20Models/TTES_Model) - Directory for the Planner that models the implementation of Tank Thermal Energy Storage (TTES) into the UK energy grid.



## Basic Description of the Optimization Problem
The optimisation problem consists of two parts: an investment master problem that decides how much new capacity to install for each type of generation/storage technology and an operational subproblem that finds the optimal schedule for the usage of these technologies given their capacities. These two models are linked with the overall objective of minimizing the combined sum of the investment and operational costs of meeting the UK’s energy demands for one year. This two-part problem is solved using a standard Benders decomposition algorithm.

**A more detailed explanation of the optimisation problem can be found in Chapter 3 of the dissertation report.**



## Basic Description of the Planner Structure
The planner directories contain various files which work in conjunction to produce a mathematical optimisation problem in Julia. As far as the contents of the Planner directories go, the key components to pay attention to are:
* **OREIA3_2017/load_models.jl** - Contains the formulation of the investment master problem and operational subproblem in Julia.
* **OREIA3_2017/data/investment.xlsx** - Contains data regarding the existing capacity and investment costs/limits of technologies included in the optimisation model.
* **OREIA3_2017/data/operation.xlsx** - Contains data regarding how the UK energy grid & demands was modelled including half-hourly demand data, energy nodes, lines, renewable energy availability, etc.

A more detailed explanation of the contents of the Planner directory can be found in the [**Planner User Guide**](https://github.com/axeleichelmann/Year5-MMath-Diss/blob/main/sPLAN%20guide.pdf)**
