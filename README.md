Elizabeth Line × Retail RV — Reproducible Workflow

This repository reproduces the dissertation results on pre-opening capitalisation of the Elizabeth line for retail properties. The pipeline is Python + SQL for data preparation and R (Quarto) for modelling and analysis.
1) Repository structure

data/                

upload final dataset (cleaned); Please download the raw dataset follow the URL in the paper's appendix.And follow the script(data_cleaning.ipynb and twfe_data_cleaning.ipynb) to get the final dataset if need repeat the process.

notebook/            

analysis notebooks and Quarto file

01_data_clean.ipynb（rv,rings,ptal,imd,caz,tc）

02_descriptive.ipynb/qmd

03_main_model_twfe.qmd

outputs/             

key outputs (tables, figures)

README.md

Run order by filename:

01_data_clean.ipynb → 02_descriptive_statistics.ipynb/qmd → 03_main_model_twfe.qmd

2) Software and environments

Python: JupyterLab, pandas, geopandas, numpy, matplotlib, sqlalchemy, a database driver (e.g., psycopg2 for Postgres or sqlite), and any geospatial libs you use.
SQL: a running RDBMS or a local SQLite file. The notebooks read from and write to it via SQLAlchemy.
R: R ≥ 4.2 with packages for modelling and reporting, for example:
fixest (TWFE and robust vcov incl. Conley), sandwich, lmtest
sf, dplyr, data.table, readr, ggplot2

3) Data and credentials

Place required inputs in data/.Before use please unzip the dataset.And set up file path.

4) How to run

Step 1 — Data cleaning (Python + SQL)

Open Jupyter and run the notebook top to bottom:

jupyter lab


In Jupyter, open notebook/01_data_clean.ipynb.
This step:

reads raw inputs from data/

queries and writes to the SQL database defined by SQLALCHEMY_URL

produces analysis-ready tables for retail properties and stations

saves cleaned exports in data/ or outputs/ as specified in the notebook


Step 2 — Descriptive statistics (Python & R)

Run notebook/02_descriptive_statistics.ipynb/qmd.

This step:

reads cleaned data

generates descriptive tables and figures used in Chapter 3 and 4.1

writes results to outputs/tables/ and outputs/figures/


Step 3 — Modelling and analysis (R + Quarto)

Render the main analysis:

quarto render notebook/main_model_twfe.qmd

performs diagnostics including spatial autocorrelation checks and Conley HAC

runs robustness and heterogeneity analyses

exports publication-ready tables and figures


5) Reproducibility conventions

Paths: notebooks read DATA_DIR and OUTPUT_DIR . Keep relative paths stable.

Randomness: set seeds where applicable and record them in the notebooks or qmd.

Rounding: display rounding is applied in exported tables; full-precision results are available in intermediate objects.

Session info: each notebook prints package versions at the end. The Quarto doc includes an R session summary.

Make it portable: avoid absolute paths, keep file names lower_snake_case, and write outputs only under outputs/.

6) Mapping outputs to the dissertation

Chapter 3: variables dictionary and descriptive tables from 02_descriptive_statistics.ipynb/qmd in outputs/tables/; maps and histograms in outputs/figures/.

Chapter 4: main regression table and coefficient plots from main_model_twfe.qmd in outputs/tables/ and outputs/figures/. Diagnostics and robustness summaries are exported with matching filenames.

Each table and figure used in the dissertation should cite the exact export filename and run date in the caption notes.


7) License and data access

Code: choose a license (e.g., MIT) and add LICENSE.

Data: All dataset are public, please use them follow the related regulation.



modelsummary or broom for tables

Quarto: for rendering main_model_twfe.qmd.
