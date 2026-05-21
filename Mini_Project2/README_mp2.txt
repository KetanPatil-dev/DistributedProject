
  DTU 02225 DRTS, Mini Project 2: CBS WCD Analysis and Simulation for TSN
  README
================================================================================

1. CONTENTS

  mini_project2.py            Main analysis script (Python 3)
  console_output.txt          Full terminal output from the final run
  plots/                      All generated plots:
    wcd_cbs_vs_sp.png             Analytical WCD: CBS vs Strict Priority (5 scenarios)
    cbs_analytical_vs_sim.png     Analytical WCD vs simulation max/avg (5 scenarios)
    cdf_scenario2.png             Response time CDFs for Scenario 2
    credit_trace_scenario1.png    CBS credit evolution trace for Scenario 1
    starvation_demo_scenario5.png BE starvation demo: CBS vs SP CDFs (Scenario 5)
    idleslope_sensitivity_scenario6.png  WCD vs idleSlope allocation (Scenario 6)

  NOTE: The professor's JSON test cases (tsn-test-cases/) are not included
  in this archive. Clone them separately (see Section 4). All results from
  the JSON test case are pre-computed in console_output.txt.

2. PREREQUISITES

  Python 3.10 or later
  matplotlib  (pip install matplotlib)
  numpy       (pip install numpy)

3. HOW TO RUN

  To run the full analysis (Scenarios 1-6 + JSON test case if available):

    python3 mini_project2.py

  To save the console output:

    python3 mini_project2.py 2>&1 | tee console_output.txt

  Plots are saved automatically to a plots/ subfolder next to the script.

4. JSON TEST CASE (OPTIONAL)

  The script automatically analyses the professor's JSON Test Case 1 if the
  tsn-test-cases repository is found. To enable this:

    git clone https://github.com/paulpop/tsn-test-cases ~/Desktop/tsn-test-cases

  The script searches for the test case in the following locations:
    - Same folder as the script
    - Parent folder of the script
    - ~/Desktop/tsn-test-cases/examples/test_case_1

  Note: Per professor clarification (29 Apr 2026), all links use
  default_bandwidth_mbps = 100 Mbps. Per-link overrides are ignored.

5. SCENARIOS

  Scenario 1  1 stream per class (1 hop)        — Baseline / sanity check
  Scenario 2  3A, 3B, 2BE (1 hop)               — Realistic mixed load
  Scenario 3  High-load stress test (1 hop)     — Near-capacity, varied frame sizes
  Scenario 4  Multi-hop (2 hops)                — WCD scaling with number of hops
  Scenario 5  BE Starvation Demo (1 hop)        — CBS vs SP starvation protection
  Scenario 6  idleSlope Sensitivity (1 hop)     — Bandwidth allocation trade-off

6. NETWORK MODEL

  Line topology: ES -> SW1 -> ... -> SWn -> ES (all streams same direction)
  Link rate:     100 Mb/s
  CBS queues:    AVB Class A (highest), AVB Class B, Best Effort (lowest)
  idleSlope:     alpha+ = alpha- = 0.5 * BW = 50 Mb/s per CBS class

================================================================================
