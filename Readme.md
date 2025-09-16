# Introduction
In this exercise, we are going to use Pandas to parse the results from a miRAW run. The run was to predict the targets of the genes in the Reactome pathway `R-HSA-111931`.

There is a Python file `target_prediction_parsing.py` which loads and processes the data. You can get an idea about what the program does by looking at the `main` function at the bottom of the file. (we use a `main` function because it's harder to work out what the program does if the function calls are scattered throughout the file.

First of all, we can do some preliminary investigation of the data.

**Question 1**: How many predictions are there in the file?

**Question 2**: What are the different columns presented in the predictions?

**Question 3**: What pathway is `R-HSA-111931` and how many genes are there in the pathway?

**Question 4**: If you look at the code, you will see there is an `parse_args` function. Which parameters can be specified by the program? and which columns do they correspond to in the data?

Hint, as the file is large you can view the first few lines of the file using the `head` command from a terminal window

in the function `summariseMiRAWResults()` it reports the number of targets for each `miRNA`, and summarises the prediction probability.

**Task 1**: 
Add code to report the number of targets for each `gene`

**Task 2**:
Add code to report the number of `canonical` vs `non-canonical` predictions

**Task 3**:
Add code to report the length distribution of the `miRNA` targets. Hint, you need to add a new column in the same way you did for the code to read the `rno.gff3` file.

**Question 6**:
Do longer `miRNA` transcripts have more targets? Can you generate a plot to show `miRNA length` versus `no of targets`?