
  DTU 02225 DRTS, Mini Project 1: DM vs EDF Scheduling Analysis
  README
================================================================================

1. CONTENTS

  mini_project1.py            Main analysis script (Python 3)
  console_output.txt          Full terminal output from the final run
  custom_testcases/           15 custom-generated task sets (12 implicit-deadline + 3 constrained-deadline)
  results/                    All generated plots

  NOTE: The course-provided test cases (test_examples/) and the batch
  datasets (automotive, uunifast) are not included in this archive.
  The test_examples/ folder is available in the course-provided
  Taskset-Generator-Exercise repository. The batch datasets can be
  generated using the tool described in Section 4. All results from
  these datasets are pre-computed in results/ and console_output.txt.

2. PREREQUISITES

  Python 3.8 or later
  matplotlib (pip install matplotlib)

3. HOW TO RUN

  To run the full analysis:

  1. Clone the course-provided Taskset-Generator-Exercise repository and
     place the test_examples/ folder next to the script.

  2. Generate or place the batch datasets next to the script:
       automotive-utilDist/automotive-perDist/1-core/25-task/0-jitter/
       uunifast-utilDist/uniform-discrete-perDist/1-core/25-task/0-jitter/


  3. Run:
       python3 mini_project1.py

  4. To save the console output:
       python3 mini_project1.py 2>&1 | tee console_output.txt

================================================================================
